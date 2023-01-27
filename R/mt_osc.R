#' =======================================================================
#' Orthogonal signal correction (OSC) for data pre-processing
#' History:
#'   wll-04-06-2007: commence
#'   wll-12-06-2007: minor changes
osc.default <- function(x, y, method = "wold", center = TRUE, osc.ncomp = 4,
                        pls.ncomp = 10, tol = 1e-3, iter = 20, ...) {
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

  method <- match.arg(method, c("wold", "sjoblom", "wise"))

  #' initialisation
  x <- as.matrix(x)
  y <- as.factor(y)
  n <- nrow(x)
  p <- ncol(x)

  if (pls.ncomp < 1 || pls.ncomp > min(n - 1, p)) {
    pls.ncomp <- min(n - 1, p)
    #' stop("Invalid number of components, ncomp")
  }

  #' Select OSC algorithm:
  oscFunc <- switch(method,
    wold = osc_wold,
    sjoblom = osc_sjoblom,
    wise = osc_wise
  )

  #' call OSC algorithm
  res <- oscFunc(x, y,
    center = center, osc.ncomp = osc.ncomp, pls.ncomp = pls.ncomp,
    tol = tol, iter = iter, ...
  )

  #' Build and return the object:
  res$call <- match.call()
  res$call[[1]] <- as.name("osc")
  res$center <- center
  res$osc.ncomp <- osc.ncomp
  res$pls.ncomp <- pls.ncomp
  res$method <- method

  class(res) <- "osc"
  return(res)
}

#' ========================================================================
#' wll-04-06-2007: predict method for OSC
predict.osc <- function(object, newdata, ...) {
  #' if(!inherits(object, "osc")) stop("object not of class \"osc\"")
  if (missing(newdata)) {
    return(object$x)
  }
  if (is.null(dim(newdata))) {
    dim(newdata) <- c(1, length(newdata))
  } #' a row vector

  newdata <- as.matrix(newdata)
  if (ncol(newdata) != ncol(object$x)) stop("wrong number of variables")

  if (object$center) {
    newdata <- sweep(newdata, 2, colMeans(newdata), "-")
  } #' column-wise center

  if (F) { #' Tom Fearn, On OSC, Chmom. Intell. Lab. Syst. 50(2000):47-52
    t <- newdata %*% object$w
    p <- t(newdata) %*% t %*% ginv(t(t) %*% t)
    x <- newdata - t %*% t(p)
  } else {
    x <- newdata - newdata %*% object$w %*% t(object$p)
  }

  #' calculate the removed variance of X
  Q2 <- sum(x^2) / sum(newdata^2) * 100

  return(list(x = x, Q2 = Q2))
}

#' ========================================================================
#' wll-04-06-2007: print method for osc
print.osc <- function(x, ...) {
  alg <- switch(x$method,
    wold     = "Wold et al approach",
    sjoblom  = "Sjoblom et al approach",
    wise     = "Wise and Gallagher approach",
    stop("Unknown approach.")
  )
  cat("Orthogonal signal correction (OSC), fitted with the", alg, ".")
  cat("\nCall:\n", deparse(x$call), "\n")

  cat("\nR2 (percentage):", x$R2)
  cat("\nAngle (degree):\t", x$angle)
  cat("\n")

  invisible(x)
}

#' ========================================================================
#' wll-04-06-2007: summary method for osc
summary.osc <- function(object, ...) {
  structure(object, class = "summary.osc")
}

#' =======================================================================
#' wll-04-06-2007: summary method for osc
print.summary.osc <- function(x, ...) {
  print.osc(x)

  cat("\nNumber of OSC components:\t", x$osc.ncomp)
  cat("\nNumber of PLS components:\t", x$pls.ncomp)

  cat("\nData dimension:\t\t\t", dim(x$x))
  cat("\nWeight dimension:\t\t", dim(x$w))
  cat("\nLoading dimension:\t\t", dim(x$p))
  cat("\nScore dimension:\t\t", dim(x$t))
  cat("\n")
}

#' =======================================================================
osc <- function(x, ...) UseMethod("osc")

#' =======================================================================
osc.formula <- function(formula, data = NULL, ..., subset,
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

  ret <- osc.default(x, y, ..., na.action = na.action)

  ret$call <- call
  ret$call[[1]] <- as.name("osc")
  ret$terms <- Terms
  if (!is.null(attr(m, "na.action"))) {
    ret$na.action <- attr(m, "na.action")
  }
  class(ret) <- c("osc.formula", class(ret))
  return(ret)
}

#' ========================================================================
#' wll-04-06-2007: Wold algorithm for OSC
osc_wold <- function(x, y, center = TRUE, osc.ncomp = 4, pls.ncomp = 10,
                     tol = 1e-3, iter = 20, ...) {
  if (center) x <- sweep(x, 2, colMeans(x), "-") #' column-wise centre
  x.ori <- x
  y <- class.ind(y) #' convert class labels to numeric values

  np <- nw <- nt <- list()
  for (i in 1:osc.ncomp) {
    pc <- prcomp(x, scale = F)
    t <- pc$x[, 1] #' PC1 as initial score
    dif <- 1
    k <- 0
    while (dif > tol & k < iter) {
      k <- k + 1
      #' Orthogonalize t to y
      tnew <- t - (y %*% ginv(t(y) %*% y) %*% t(y) %*% t)
      #' calculate weight vector using PLS
      pls <- simpls.fit(x, tnew, ncomp = pls.ncomp, ...)
      w <- pls$coefficients[, , ncomp = pls.ncomp, drop = FALSE]
      w <- w / sqrt(sum(w^2))

      tnew <- x %*% w
      #' Check for convergence
      dif <- sqrt(sum((tnew - t)^2) / sum(tnew^2))
      t <- tnew
    }

    p <- t(x) %*% t %*% ginv(t(t) %*% t)
    x <- x - t %*% t(p)

    nw[[i]] <- w
    np[[i]] <- p
    nt[[i]] <- t
  }
  nw <- do.call(cbind, nw)
  np <- do.call(cbind, np)
  nt <- do.call(cbind, nt)
  #' OSC-correct the original data set
  x <- x.ori - x.ori %*% nw %*% t(np)

  #' Calculate the fraction of the variation in X. (also called the removed
  #' /remained variance of X)
  R2 <- sum(x^2) / sum(x.ori^2) * 100
  #' R2 <- var(as.vector(x^2))/var(as.vector(x.ori^2))

  #' Calculate angle which assesses that t vector is orthogonal to y.
  angle <- t(nt) %*% y
  norm <- ginv(sqrt(apply(nt^2, 2, sum) * sum(y^2)))
  angle <- t(angle) %*% t(norm)
  angle <- mean(acos(angle) * 180 / pi)

  res <- list(x = x, R2 = R2, angle = angle, w = nw, p = np, t = nt,
              center = center)
  return(res)
}

#' ========================================================================
#' wll-03-06-2007: Sjoblom algorithm
#' wll-03-06-2007: Fix a algorithm mis-understanding in updating weights.
#' wll-03-06-2007: Orthogonalize t to y in the last step. This improve the
#'                 performance of OSC.
osc_sjoblom <- function(x, y, center = TRUE, osc.ncomp = 4, pls.ncomp = 10,
                        tol = 1e-3, iter = 20, ...) {
  if (center) x <- sweep(x, 2, colMeans(x), "-") #' column-wise centre
  x.ori <- x
  y <- class.ind(y) #' convert class labels to numeric values

  np <- nw <- nt <- list()
  for (i in 1:osc.ncomp) {
    pc <- prcomp(x, scale = F)
    t <- pc$x[, 1] #' PC1 as initial score
    dif <- 1
    k <- 0
    while (dif > tol & k < iter) {
      k <- k + 1
      #' Orthogonalize t to y (by tnew = t - y*inv(y'*y)*y'*t).
      tnew <- t - (y %*% ginv(t(y) %*% y) %*% t(y) %*% t)
      #' Update weights and scores
      w <- t(x) %*% tnew %*% ginv(t(tnew) %*% tnew)

      #' w    <- t(x) %*% tnew                  #' NOTE: not this one!

      w <- w / sqrt(sum(w^2))
      tnew <- x %*% w
      #' Check for convergence
      dif <- sqrt(sum((tnew - t)^2) / sum(tnew^2))
      t <- tnew
    }

    #' fit PLS model
    pls <- simpls.fit(x, tnew, ncomp = pls.ncomp, ...)

    #' extract the coefficients as weights in OSC
    w <- pls$coefficients[, , ncomp = pls.ncomp, drop = FALSE]

    #' Update scores, loads and corrected data
    t <- x %*% w
    #' Orthogonalize t to y
    t <- t - y %*% ginv(t(y) %*% y) %*% t(y) %*% t

    p <- t(x) %*% t %*% ginv(t(t) %*% t)
    x <- x - t %*% t(p)

    nw[[i]] <- w
    np[[i]] <- p
    nt[[i]] <- t
  }
  nw <- do.call(cbind, nw)
  np <- do.call(cbind, np)
  nt <- do.call(cbind, nt)
  #' OSC-correct the original data set
  x <- x.ori - x.ori %*% nw %*% t(np)

  #' Calculate the fraction of the variation in X. (also called the removed
  #' /remained variance of X)
  R2 <- sum(x^2) / sum(x.ori^2) * 100
  #' R2 <- var(as.vector(x^2))/var(as.vector(x.ori^2))

  #' Calculate angle which assesses that t vector is orthogonal to y.
  angle <- t(nt) %*% y
  norm <- ginv(sqrt(apply(nt^2, 2, sum) * sum(y^2)))
  angle <- t(angle) %*% t(norm)
  angle <- mean(acos(angle) * 180 / pi)

  res <- list(
    x = x, R2 = R2, angle = angle, w = nw, p = np, t = nt,
    center = center
  )

  return(res)
}

#' ========================================================================
#' wll-03-06-2007: Wise algorithm
osc_wise <- function(x, y, center = TRUE, osc.ncomp = 4, pls.ncomp = 10,
                     tol = 1e-3, iter = 20, ...) {
  if (center) x <- sweep(x, 2, colMeans(x), "-")
  x.ori <- x
  y <- class.ind(y)

  np <- nw <- nt <- list()
  for (i in 1:osc.ncomp) {
    pc <- prcomp(x, scale = F)
    told <- pc$x[, 1] #' initial score
    p <- pc$rotation[, 1] #' initial loadings

    dif <- 1
    k <- 0
    while (dif > 1e-5 & k < iter) {
      k <- k + 1
      #' Calculate scores from loads (by t = x*p/(p'*p) ).
      t <- (x %*% p) %*% solve(t(p) %*% p)
      #' Othogonalize t to y  (by tnew = t - y*inv(y'*y)*y'*t).
      tnew <- t - (y %*% solve(t(y) %*% y) %*% t(y) %*% t)
      #' Compute a new loading  (by pnew = x'*tnew/(tnew'*tnew) ).
      pnew <- (t(x) %*% tnew) %*% solve(t(tnew) %*% tnew)
      #' Check for convergence
      dif <- sqrt(sum((tnew - told)^2) / sum(tnew^2))
      told <- tnew
      p <- pnew
    }

    #' fit PLS model
    nc <- min(pls.ncomp, qr(x, tol = 1e-9)$rank)
    pls <- simpls.fit(x, tnew, ncomp = nc, ...)
    w <- pls$coefficients[, , ncomp = nc, drop = FALSE]
    w <- w / sqrt(sum(w^2))
    #' pls <- plsr(tnew~x,ncomp=ncomp,method="simpls")
    #' w   <- coef.mvr(pls,ncomp=nc)

    #' Calculate new scores vector
    t <- x %*% w
    #' Othogonalize t to y
    t <- t - y %*% ginv(t(y) %*% y) %*% t(y) %*% t
    #' Compute new p
    p <- t(x) %*% t %*% ginv(t(t) %*% t)
    #' Remove orthogonal signal from x
    x <- x - t %*% t(p)

    nw[[i]] <- w
    np[[i]] <- p
    nt[[i]] <- t
  }
  nw <- do.call(cbind, nw)
  np <- do.call(cbind, np)
  nt <- do.call(cbind, nt)
  #' OSC-correct the original data set
  x <- x.ori - x.ori %*% nw %*% t(np)

  #' Calculate the fraction of the variation in X. (also called the removed
  #' /remained variance of X)
  R2 <- sum(x^2) / sum(x.ori^2) * 100
  #' R2 <- var(as.vector(x^2))/var(as.vector(x.ori^2))

  #' Calculate angle which assesses that t vector is orthogonal to y.
  angle <- t(nt) %*% y
  norm <- ginv(sqrt(apply(nt^2, 2, sum) * sum(y^2)))
  angle <- t(angle) %*% t(norm)
  angle <- mean(acos(angle) * 180 / pi)

  res <- list(
    x = x, R2 = R2, angle = angle, w = nw, p = np, t = nt,
    center = center
  )

  return(res)
}

#'  1) osc.default
#'  2) predict.osc
#'  3) print.osc
#'  4) summary.osc
#'  5) print.summary.osc
#'  6) osc
#'  7) osc.formula
#'  8) osc_wold
#'  9) osc_sjoblom
#' 10) osc_wise
