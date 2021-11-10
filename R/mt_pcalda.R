#' =========================================================================
#' wll-02-10-2007: tune the best number of components
tune.pcalda <- function(x, y, ncomp = NULL, tune.pars, ...) {
  n <- nrow(x)
  g <- length(levels(y))

  if (is.null(ncomp)) {
    ncomp <- max(g, n - g)
  } else {
    if (ncomp < 1 || ncomp > max(g, n - g)) {
      ncomp <- max(g, n - g)
    }
  }

  if (missing(tune.pars)) {
    tune.pars <- valipars(sampling = "rand", niter = 1, nreps = 4)
  }

  cat("ncomp tune (", ncomp, "):", sep = "")
  res <- sapply(1:ncomp, function(i) {
    cat(" ", i, sep = "")
    flush.console()
    accest(x, y,
      pars = tune.pars, method = "pcalda", ncomp = i,
      tune = F, ...
    )$acc
  })
  cat("\n")

  list(ncomp = which.max(res), acc.tune = res)
}

#' =========================================================================
#' wll-02-10-2007: Tune the best number of components
#' NOTE: This is a debug version of tune.pcalda, which can take user defined
#'       range of ncomp.
tune.pcalda.1 <- function(x, y, ncomp = NULL, tune.pars, ...) {
  n <- nrow(x)
  g <- length(levels(y))

  if (is.null(ncomp)) {
    ncomp <- 1:max(g, n - g)
  } else {
    if (max(ncomp) < 1 || max(ncomp) > max(g, n - g)) {
      ncomp <- 1:max(g, n - g)
    }
  }

  if (missing(tune.pars)) {
    tune.pars <- valipars(sampling = "rand", niter = 1, nreps = 4)
  }

  #'  res <- sapply(1:ncomp, function(i) {
  #'    accest(x, y, pars = tune.pars, method = "pcalda", ncomp=i, tune=F,...)$acc
  #'  })
  cat("ncomp tune (", max(ncomp), "):", sep = "")
  func <- function(i) {
    cat(" ", i, sep = "")
    flush.console()
    accest(
      dat = x, cl = y, pars = tune.pars, method = "pcalda", ncomp = i,
      tune = F, ...
    )$acc
  }
  res <- sapply(ncomp, FUN = func)
  names(res) <- ncomp
  cat("\n")

  list(ncomp = which.max(res), acc.tune = res)
}

#' =========================================================================
#' wll-02-10-2007: tune the best number of components
#'  NOTE: This is a test version of tune ncomp, which uses LDA directly and
#'        runs fast than tune.pcalda.
tune.pcalda.2 <- function(x, y, ncomp = NULL, center = TRUE, scale. = FALSE,
                          tune.pars, ...) {
  n <- nrow(x)
  g <- length(levels(y))

  if (is.null(ncomp)) {
    ncomp <- max(g, n - g)
  } else {
    if (ncomp < 1 || ncomp > max(g, n - g)) {
      ncomp <- max(g, n - g)
    }
  }

  pca.out <- prcomp(x, center = center, scale. = scale., ...)
  ncomp <- min(ncomp, length(pca.out$center))

  if (missing(tune.pars)) {
    tune.pars <- valipars(sampling = "rand", niter = 1, nreps = 10)
  }

  cat("ncomp tune (", ncomp, "):", sep = "")
  res <- sapply(1:ncomp, function(i) {
    cat(" ", i, sep = "")
    flush.console()
    acc <- accest(pca.out$x[, 1:i, drop = F], y,
      pars = tune.pars,
      method = "lda"
    )$acc
  })
  cat("\n")

  list(ncomp = which.max(res), acc.tune = res)
}

#' =========================================================================
#' PCA+LDA for classification
#' History:
#'   wll-22-06-2007: commence
#'   wll-01-07-2007: try number of PCs as n - g. Over-fitting.
#'   wll-24-01-2008: stip off constant PCs within group
pcalda.default <- function(x, y, center = TRUE, scale. = FALSE, ncomp = NULL,
                           tune = FALSE, ...) {

  #' arguments validity checking
  if (missing(x) || missing(y)) {
    stop("data set or class are missing")
  }
  x <- as.matrix(x)
  if (nrow(x) != length(y)) stop("x and y don't match.")
  y <- as.factor(y)
  if (any(table(y) == 0)) stop("Can't have empty classes in y.")
  #' lwc-NOTE: Or simple apply: y <- factor(y), which will drop factor levels

  if (length(unique(y)) < 2) {
    stop("Classification needs at least two classes.")
  }
  if (any(is.na(x)) || any(is.na(y))) {
    stop("NA is not permitted in data set or class labels.")
  }

  n <- nrow(x)
  g <- length(levels(y))

  #' The singularity problem of the within-class scatter matrix is overcome
  #' if number of retained PCs varies at least g to at most n-g. Here g and
  #' n is the number of classes and training data, respectively.
  #' if(is.null(ncomp)) ncomp <- max(g,n - g)

  #' NOTE-wll: Too good for training, which means overfits the training data.

  #' if(is.null(ncomp)) ncomp <- max(g,round(n/2))

  #' Check the number of components
  if (is.null(ncomp)) {
    ncomp <- if (tune) max(g, n - g) else max(g, round(n / 2))
  } else {
    if (ncomp < 1 || ncomp > max(g, n - g)) {
      ncomp <- max(g, n - g)
    }
  }

  #' find the best number of components
  if (tune) {
    val <- tune.pcalda(x, y, ncomp, ...)
    ncomp <- val$ncomp
  }

  #' dimension reduction by PCA
  pca.out <- prcomp(x, center = center, scale. = scale., ...)
  ncomp <- min(ncomp, length(pca.out$center))
  x.tmp <- pca.out$x[, 1:ncomp, drop = F]

  #' stip off PCs constant within groups
  x.tmp <- preproc.const(x.tmp, y)
  ncomp <- ncol(x.tmp)
  #' NOTE-28-01-2008: If the variables being strippt off is not in the end
  #'   of columns, they positions should be sotred somewhere. But this
  #'   situation is rare in PCs. Refer to predict.pcalda where the ncomp is
  #'   used.

  lda.out <- lda(x.tmp, y)
  pred <- predict(lda.out, x.tmp)
  conf <- table(y, pred$class)
  acc <- round(sum(diag(conf)) * 100 / n, 2)

  x <- scale(x.tmp, center = colMeans(lda.out$means), scale = FALSE) %*%
    lda.out$scaling

  res <- list(
    x = x, cl = y, pred = pred$class, posterior = pred$posterior,
    conf = conf, acc = acc, ncomp = ncomp, pca.out = pca.out,
    lda.out = lda.out
  )

  if (tune) res$acc.tune <- val$acc.tune
  res$call <- match.call()
  res$call[[1]] <- as.name("pcalda")
  class(res) <- "pcalda"
  return(res)
}

#' =========================================================================
predict.pcalda <- function(object, newdata, ...) {
  if (!inherits(object, "pcalda")) stop("object not of class \"pcalda\"")
  if (missing(newdata)) {
    res <- list(class = object$pred, posterior = object$posterior, x = object$x)
    return(res)
  }

  if (is.null(dim(newdata))) {
    dim(newdata) <- c(1, length(newdata))
  } #' a row vector

  newdata <- as.matrix(newdata)
  if (ncol(newdata) != length(object$pca.out$center)) {
    stop("wrong number of variables")
  }

  #' rotated data (projection) by PCA
  x <- predict(object$pca.out, newdata)
  x <- x[, 1:object$ncomp, drop = F]

  #' predict using LDA
  pred <- predict(object$lda.out, x)
  #' list(class=pred$class, posterior = pred$posterior, x = pred$x)

  return(pred)
}

#' ========================================================================
#' wll-22-06-2007: print method for pcalda
#' wll-15-01-2008: add ratio(svd values)
print.pcalda <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n")
  cat("\nNumber of components considered:", x$ncomp)

  cat("\nConfusion matrix of training data:\n")
  print(x$conf)
  cat("\nRatio of between- and within-group s.d.:\n")
  df <- colnames(x$x)
  ratio <- x$lda.out$svd
  names(ratio) <- df
  print(ratio)
  cat("\nAccurary rate of training data:\n")
  print(x$acc)
  invisible(x)
}

#' ========================================================================
#' lwc-22-06-2007: summary method for pcalda
summary.pcalda <- function(object, ...) {
  structure(object, class = "summary.pcalda")
}

#' ========================================================================
#' lwc-22-06-2007: summary method for pcalda
print.summary.pcalda <- function(x, ...) {
  print.pcalda(x)

  lev <- levels(x$cl)
  cat("\nNumber of Classes: ", length(lev), "\n\n")
  cat("Levels:", if (is.numeric(lev)) "(as integer)", "\n", lev)
  cat("\n\n")
}

#' ========================================================================
#' wll-13-12-2007: plot method for pcalda using lattice.
#' wll-11-10-2008: Definition of svd in lda documents:
#'   svd is the singular values, which give the ratio of the between- and
#'   within-group standard deviations on the linear discriminant variables.
#'   Their squares are the canonical F-statistics.
#' wll-15-01-2008: add svd listed in the plot.
plot.pcalda <- function(x, dimen, ...) {
  ld.names <- function(object, comps) {
    labs <- paste("LD", 1:length(object$means), sep = "")
    if (missing(comps)) {
      comps <- seq(along = labs)
    } else {
      labs <- labs[comps]
    }
    svd <- object$svd
    svd.p <- 100 * svd / sum(svd)
    #' svd.p <- 100 * svd^2/sum(svd^2)
    #' wll-11-01-2008: check lda's print method: Proportion of trace
    svd <- svd[comps]
    svd.p <- svd.p[comps]

    #' labs <- paste(labs, " (", format(svd.p, digits = 2, trim = TRUE),
    #'               " %)", sep = "")
    labs <- paste(labs, " (", format(svd, digits = 2, trim = TRUE), ", ",
      format(svd.p, digits = 2, trim = TRUE), "%)",
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
  dfn <- ld.names(x$lda.out, dimen)

  y <- x$cl
  x <- data.frame(x$x[, dimen, drop = FALSE])
  names(x) <- dfn

  #' call group plot
  p <- grpplot(x, y, plot = "pairs", ...)
  p
}

#' =========================================================================
#' lwc-22-06-2007: plot method for pcalda. It plot LDA scores.
plot.pcalda.1 <- function(x, panel = panel.pcalda, cex = 0.7, dimen,
                          abbrev = FALSE, ...) {
  panel.pcalda <- function(x, y, ...) {
    text(x, y, as.character(g.nlda), cex = tcex, col = unclass(g), ...)
  }

  ld.names <- function(object, comps) {
    labs <- paste("LD", 1:length(object$means), sep = "")
    if (missing(comps)) {
      comps <- seq(along = labs)
    } else {
      labs <- labs[comps]
    }
    svd <- object$svd
    svd <- 100 * svd^2 / sum(svd^2)
    evar <- svd[comps]

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
  varlab <- ld.names(x$lda.out, dimen)
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

#' =========================================================================
pcalda <- function(x, ...) UseMethod("pcalda")

#' =========================================================================
pcalda.formula <- function(formula, data = NULL, ..., subset,
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

  ret <- pcalda.default(x, y, ..., na.action = na.action)

  ret$call <- call
  ret$call[[1]] <- as.name("pcalda")
  ret$terms <- Terms
  if (!is.null(attr(m, "na.action"))) {
    ret$na.action <- attr(m, "na.action")
  }
  class(ret) <- c("pcalda.formula", class(ret))
  return(ret)
}

#'  1) tune.pcalda
#'  2) tune.pcalda.1
#'  3) tune.pcalda.2
#'  4) pcalda.default
#'  5) predict.pcalda
#'  6) print.pcalda
#'  7) summary.pcalda
#'  8) print.summary.pcalda
#'  9) plot.pcalda
#' 10) plot.pcalda.1
#' 11) pcalda
#' 12) pcalda.formula
