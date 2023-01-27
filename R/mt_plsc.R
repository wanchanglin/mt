#' =========================================================================
#' wll-02-10-2007: tune the best number of components
tune.plsc <- function(x, y, pls = "simpls", ncomp = 10, tune.pars, ...) {
  if (missing(tune.pars)) {
    tune.pars <- valipars(sampling = "rand", niter = 1, nreps = 10)
  }

  cat("ncomp tune (", ncomp, "):", sep = "")
  res <- sapply(1:ncomp, function(i) {
    cat(" ", i, sep = "")
    flush.console()
    accest(x, y,
      pars = tune.pars, method = "plsc", pls = pls, ncomp = i,
      tune = FALSE, ...
    )$acc
  })
  cat("\n")
  list(ncomp = which.max(res), acc.tune = res)
}

#' ========================================================================
#' PLS for classification
#' History:
#'   wll-01-10-2007: commence
#'   lwc-21-05-2012: use wrapper function of "mvr".
#'   lwc-21-05-2012: It should not be difficult to get R2 for "accest" with
#'                   methods of plsc and plslda.
plsc.default <- function(x, y, pls = "simpls", ncomp = 10, tune = FALSE,
                         ...) {
  #' arguments validity checking
  if (missing(x) || missing(y)) {
    stop("data set or class are missing")
  }
  if (nrow(x) != length(y)) stop("x and y don't match.")
  if (length(unique(y)) < 2) {
    stop("Classification needs at least two classes.")
  }
  if (any(is.na(x)) || any(is.na(y))) {
    stop("NA is not permitted in data set or class labels.")
  }

  pls <- match.arg(pls, c("kernelpls", "simpls", "oscorespls"))
  pls.fit <- paste(pls, ".fit", sep = "")

  #' initialisation
  x <- as.matrix(x)
  y <- as.factor(y)
  n <- nrow(x)
  p <- ncol(x)

  if (ncomp < 1 || ncomp > min(n - 1, p)) {
    ncomp <- min(n - 1, p)
    #' stop("Invalid number of components, ncomp")
  }

  #' find the best number of components
  if (tune) {
    val <- tune.plsc(x, y, pls = "simpls", ncomp = ncomp, ...)
    ncomp <- val$ncomp
  }

  #' Convert response to be suitable format for PLS
  y.1 <- class.ind(y)

  #' Call PLS for regression
  #' pls.out <- do.call(pls.fit, c(list(X=x, Y=y.1, ncomp=ncomp), list(...)))

  #' lwc-21-05-2012: use wrapper function of "mvr".
  pls.out <- plsr(y.1 ~ x, method = pls, ncomp = ncomp, ...)
  names(pls.out$Ymeans) <- colnames(y.1)
  #' NOTE: New version of pls strip off the names of Ymeans

  #' Predict with models containing ncomp components.
  B <- pls.out$coefficients[, , ncomp, drop = T]
  if (length(pls.out$Xmeans) > 1) { #' wll: 07-01-2008
    B0 <- pls.out$Ymeans - pls.out$Xmeans %*% B
  } else {
    B0 <- pls.out$Ymeans - pls.out$Xmeans * B
  }
  B0 <- rep(B0, each = nrow(x))
  B1 <- x %*% B
  pred <- B1 + B0
  #' softmax for convert values to probabilities
  poste <- exp(pred)
  poste <- poste / rowSums(poste)
  #' generates the predict class
  nm <- names(pls.out$Ymeans)
  pred.cl <- factor(nm[max.col(poste)], levels = levels(y))
  dimnames(poste) <- list(rownames(x), nm)

  conf <- table(y, pred.cl)
  acc <- round(sum(diag(conf)) * 100 / n, 2)

  lc <- unclass(pls.out$scores) #' latent components
  colnames(lc) <- paste("LC", 1:ncol(lc), sep = "")

  res <- list(
    x = lc, cl = y, pred = pred.cl, posterior = poste, conf = conf, acc = acc,
    ncomp = ncomp, pls.method = pls, pls.out = pls.out
  )

  if (tune) res$acc.tune <- val$acc.tune
  res$call <- match.call()
  res$call[[1]] <- as.name("plsc")
  class(res) <- "plsc"
  return(res)
}

#' =========================================================================
#' wll-01-10-2007: predict method. See predict.mvr and coef.mvr in package
#'                 pls for details.
#' NOTE: More comments are in function pred.pls which can be call directly
#'       by simpls.fit etc.
predict.plsc <- function(object, newdata, ...) {
  if (!inherits(object, "plsc")) stop("object not of class \"plsc\"")
  if (missing(newdata)) {
    return(list(
      class = object$pred, posterior = object$posterior,
      x = unclass(object$pls.out$scores)
    ))
  }
  if (is.null(dim(newdata))) {
    dim(newdata) <- c(1, length(newdata))
  } #' a row vector

  newdata <- as.matrix(newdata)
  if (ncol(newdata) != length(object$pls.out$Xmeans)) {
    stop("wrong number of variables")
  }

  #' rotated data (projection)
  x <-
    scale(newdata, center = object$pls.out$Xmeans, scale = FALSE) %*%
    object$pls.out$projection

  #' Predict with models containing ncomp components.
  ncomp <- object$ncomp
  B <- object$pls.out$coefficients[, , ncomp, drop = T]
  if (length(object$pls.out$Xmeans) > 1) { #' wll: 07-01-2008
    B0 <- object$pls.out$Ymeans - object$pls.out$Xmeans %*% B
  } else {
    B0 <- object$pls.out$Ymeans - object$pls.out$Xmeans * B
  }
  B0 <- rep(B0, each = nrow(newdata))
  B1 <- newdata %*% B
  pred <- B1 + B0
  poste <- exp(pred)
  poste <- poste / rowSums(poste)

  #' predict classification
  nm <- names(object$pls.out$Ymeans)
  cl <- factor(nm[max.col(poste)], levels = levels(object$cl))
  #' NOTE: The levels of factor must be consistent with train data set.
  dimnames(poste) <- list(rownames(x), nm)

  list(class = cl, posterior = poste, x = x)
}

#' ========================================================================
#' wll-23-05-2007: print method for plsc
print.plsc <- function(x, ...) {
  alg <- switch(x$pls.method,
    kernelpls  = "kernel",
    simpls     = "simpls",
    oscorespls = "orthogonal scores",
    stop("Unknown fit method.")
  )
  cat("Partial least squares classification, fitted with the", alg, "algorithm.")
  #' cat("\nCall:\n", deparse(x$call), "\n")
  cat("\nCall:\n")
  dput(x$call)

  cat("\nConfusion matrix of training data:\n")
  print(x$conf)
  cat("\nAccuracy rate of training data:\n")
  print(x$acc)
  invisible(x)
}

#' ========================================================================
#' wll-23-05-2007: summary method for plsc
summary.plsc <- function(object, ...) {
  structure(object, class = "summary.plsc")
}

#' ========================================================================
#' lwc-23-05-2007: summary method for plsc
print.summary.plsc <- function(x, ...) {
  print.plsc(x)

  cat("\nNumber of components considered:", x$ncomp)
  cat("\nProportion of components:\n")
  evar <- 100 * x$pls.out$Xvar / x$pls.out$Xtotvar
  names(evar) <- paste("LC", 1:length(evar), sep = "")
  print(evar)

  lev <- levels(x$cl)
  cat("\nNumber of Classes: ", length(lev), "\n\n")
  cat("Levels:", if (is.numeric(lev)) "(as integer)", "\n", lev)
  cat("\n\n")
}

#' =========================================================================
plsc <- function(x, ...) UseMethod("plsc")

#' =========================================================================
plsc.formula <- function(formula, data = NULL, ..., subset,
                         na.action = na.omit) {
  call <- match.call()
  if (!inherits(formula, "formula")) {
    stop("method is only for formula objects")
  }
  m <- match.call(expand.dots = FALSE)
  if (identical(class(eval.parent(m$data)), "matrix")) {
    m$data <- as.data.frame(eval.parent(m$data))
  }
  m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m$na.action <- na.action
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  x <- model.matrix(Terms, m)
  y <- model.extract(m, "response")
  attr(x, "na.action") <- attr(y, "na.action") <- attr(m, "na.action")

  ret <- plsc.default(x, y, ..., na.action = na.action)

  ret$call <- call
  ret$call[[1]] <- as.name("plsc")
  ret$terms <- Terms
  if (!is.null(attr(m, "na.action"))) {
    ret$na.action <- attr(m, "na.action")
  }
  class(ret) <- c("plsc.formula", class(ret))
  return(ret)
}

#' ========================================================================
#' wll-13-12-2007: plot method for plsc using lattice.
plot.plsc <- function(x, dimen, ...) {
  lc.names <- function(object, comps) {
    labs <- paste("LC", 1:length(object$Xvar), sep = "")
    if (missing(comps)) {
      comps <- seq(along = labs)
    } else {
      labs <- labs[comps]
    }
    evar <- 100 * object$Xvar / object$Xtotvar
    evar <- evar[comps]

    labs <- paste(labs, " (", format(evar, digits = 2, trim = TRUE),
      " %)",
      sep = ""
    )
    return(labs)
  }

  if (missing(dimen)) {
    dimen <- seq(along = colnames(x$x))
  } else {
    #' check validity
    if (!all(dimen %in% c(1:ncol(x$x)))) {
      stop("dimen is not valid")
    }
  }
  dfn <- lc.names(x$pls.out, dimen)

  y <- x$cl
  x <- data.frame(x$x[, dimen, drop = FALSE])
  names(x) <- dfn

  #' call group plot
  p <- grpplot(x, y, plot = "pairs", ...)
  p
}

#' =========================================================================
#' wll-22-05-2007: plot method for plsc. It plot PLS latent components.
plot.plsc.1 <- function(x, panel = panel.plsc, cex = 0.7, dimen,
                        abbrev = FALSE, ...) {
  panel.plsc <- function(x, y, ...) {
    text(x, y, as.character(g.nlda), cex = tcex, col = unclass(g), ...)
  }

  lc.names <- function(object, comps) {
    labs <- paste("LC", 1:length(object$Xvar), sep = "")
    if (missing(comps)) {
      comps <- seq(along = labs)
    } else {
      labs <- labs[comps]
    }
    evar <- 100 * object$Xvar / object$Xtotvar
    evar <- evar[comps]

    labs <- paste(labs, " (", format(evar, digits = 2, trim = TRUE),
      " %)",
      sep = ""
    )
    return(labs)
  }

  xval <- x$x
  g <- x$cl

  if (abbrev) levels(g) <- abbreviate(levels(g), abbrev)
  assign("g.nlda", g)
  assign("tcex", cex)

  if (missing(dimen)) {
    dimen <- seq(along = colnames(xval))
  } else {
    #' check validity
    if (!all(dimen %in% c(1:ncol(xval)))) {
      stop("dimen is not valid")
    }
  }

  xval <- xval[, dimen, drop = FALSE]
  varlab <- lc.names(x$pls.out, dimen)
  nDimen <- length(dimen)

  if (nDimen <= 2) {
    if (nDimen == 1) { #' One component
      ldahist(xval[, 1], g, ...)
      #' ldahist(xval, g, xlab=varlab,...)
    } else { #' Second component versus first
      xlab <- varlab[1]
      ylab <- varlab[2]
      eqscplot(xval, xlab = xlab, ylab = ylab, type = "n", ...)
      panel(xval[, 1], xval[, 2], ...)
    }
  } else { #' Pairwise scatterplots of several components
    pairs(xval, labels = varlab, panel = panel, ...)
  }
  invisible(NULL)
}

#'  1) tune.plsc
#'  2) plsc.default
#'  3) predict.plsc
#'  4) print.plsc
#'  5) summary.plsc
#'  6) print.summary.plsc
#'  7) plsc
#'  8) plsc.formula
#'  9) plot.plsc
#' 10) plot.plsc.1
