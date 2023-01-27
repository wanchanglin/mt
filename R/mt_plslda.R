#' ========================================================================
#' wll-02-10-2007: tune the best number of components
tune.plslda <- function(x, y, pls = "simpls", ncomp = 10, tune.pars, ...) {
  if (missing(tune.pars)) {
    tune.pars <- valipars(sampling = "rand", niter = 1, nreps = 10)
  }

  cat("ncomp tune (", ncomp, "):", sep = "")
  res <- sapply(1:ncomp, function(i) {
    cat(" ", i, sep = "")
    flush.console()
    accest(x, y,
      pars = tune.pars, method = "plslda", pls = pls, ncomp = i,
      tune = FALSE, ...
    )$acc
  })
  cat("\n")
  list(ncomp = which.max(res), acc.tune = res)
}

#' ========================================================================
#' PLS+LDA for classification
#' History:
#'   wll-21-05-2007: commence
#'   lwc-21-05-2012: use wrapper function of "mvr".
#'   lwc-21-05-2012: It should not be difficult to get R2 for "accest" with
#'                   methods of plsc and plslda.
plslda.default <- function(x, y, pls = "simpls", ncomp = 10, tune = FALSE,
                           ...) {

  #' Generates Class Indicator Matrix from a Factor.
  #' A matrix which is zero except for the column corresponding to the class.
  #' note:from package NNET
  class.ind <- function(cl) {
    n <- length(cl)
    cl <- as.factor(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (unclass(cl) - 1)] <- 1
    dimnames(x) <- list(names(cl), levels(cl))
    x
  }

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
    val <- tune.plslda(x, y, pls = "simpls", ncomp, ...)
    ncomp <- val$ncomp
  }

  #' Use PLS for dimension reduction
  #' pls.out <- do.call(pls.fit, c(list(X=x, Y=class.ind(y), ncomp=ncomp),
  #'                    list(...)))

  #' lwc-21-05-2012: use wrapper function of "mvr".
  pls.out <- plsr(class.ind(y) ~ x, method = pls, ncomp = ncomp, ...)

  #' Use latent variables as input data for LDA.
  x.lv <- unclass(pls.out$scores)

  #' Transform test data using weight matrix (projection)(Xt = X*W)
  #' Ztest <- scale(Xtest,center=pls.out$Xmeans,scale=FALSE) %*%
  #'          pls.out$projection

  lda.out <- lda(x.lv, y)
  pred <- predict(lda.out, x.lv)
  conf <- table(y, pred$class)
  acc <- round(sum(diag(conf)) * 100 / n, 2)

  lc <- unclass(pls.out$scores) #' latent components
  colnames(lc) <- paste("LC", 1:ncol(lc), sep = "")

  res <- list(
    x = lc, cl = y, pred = pred, posterior = pred$posterior,
    conf = conf, acc = acc, ncomp = ncomp, pls.method = pls,
    pls.out = pls.out, lda.out = lda.out
  )

  if (tune) res$acc.tune <- val$acc.tune
  res$call <- match.call()
  res$call[[1]] <- as.name("plslda")
  class(res) <- "plslda"
  return(res)
}

#' ========================================================================
predict.plslda <- function(object, newdata, ...) {
  if (!inherits(object, "plslda")) stop("object not of class \"plslda\"")
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
  x <- scale(newdata, center = object$pls.out$Xmeans, scale = FALSE) %*%
    object$pls.out$projection

  pred <- predict(object$lda.out, x)

  list(class = pred$class, posterior = pred$posterior, x = x)
}

#' ========================================================================
#' lwc-23-05-2007: print method for plslda
print.plslda <- function(x, ...) {
  alg <- switch(x$pls.method,
    kernelpls  = "kernel",
    simpls     = "simpls",
    oscorespls = "orthogonal scores",
    stop("Unknown fit method.")
  )
  cat(
    "Partial least squares classification, fitted with the", alg,
    "algorithm."
  )
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
#' lwc-23-05-2007: summary method for plslda
summary.plslda <- function(object, ...) {
  structure(object, class = "summary.plslda")
}

#' ========================================================================
#' lwc-23-05-2007: summary method for plslda
print.summary.plslda <- function(x, ...) {
  print.plslda(x)

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
plslda <- function(x, ...) UseMethod("plslda")

#' =========================================================================
plslda.formula <- function(formula, data = NULL, ..., subset,
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

  ret <- plslda.default(x, y, ..., na.action = na.action)

  ret$call <- call
  ret$call[[1]] <- as.name("plslda")
  ret$terms <- Terms
  if (!is.null(attr(m, "na.action"))) {
    ret$na.action <- attr(m, "na.action")
  }
  class(ret) <- c("plslda.formula", class(ret))
  return(ret)
}

#' ========================================================================
#' wll-13-12-2007: plot method for plsc using lattice.
plot.plslda <- function(x, dimen, ...) {
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
#' lwc-22-05-2007: plot method for plslda. It plot PLS latent components.
plot.plslda.1 <- function(x, panel = panel.plslda, cex = 0.7, dimen,
                          abbrev = FALSE, ...) {
  panel.plslda <- function(x, y, ...) {
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

#'  1) tune.plslda
#'  2) plslda.default
#'  3) predict.plslda
#'  4) print.plslda
#'  5) summary.plslda
#'  6) print.summary.plslda
#'  7) plslda
#'  8) plslda.formula
#'  9) plot.plslda
#' 10) plot.plslda.1
