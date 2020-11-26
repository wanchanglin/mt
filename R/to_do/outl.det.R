#' lwc-07-01-2014: outlier detection by univariate and multivariate.
#' lwc-08-01-2014: Need more efforts. Refer to pca.outlier, pca.outlier.1,
#'  pcaplot and outl.det in my another package FIEmspro.
#' wl-14-06-2019, Fri: have a check and will finish soon. should use some
#' functions such as such as `influencePlot` and `outlierTest` in 'car'.
#' also check Cook's distance.

#' ========================================================================
#' lwc-30-01-2013: Multivariate outlier detection
#' lwc-04-02-2013: Major changes.
#' Note: "mve" and "mcd" are based on approximate search. User need to set up
#'       random seed by set.seed. For details, see ?cov.rob
#' Usages:
if (F) {
  set.seed(134)
  x <- cbind(rnorm(80), rnorm(80), rnorm(80))
  y <- cbind(rnorm(10, 5, 1), rnorm(10, 5, 1), rnorm(10, 5, 1))
  x <- rbind(x, y)
  outl <- pca.outlier(iris[, 1:4], adj = -0.5)
  outl$outlier
  outl.1 <- outl.det.m(x, method = "mcd", conf.level = 0.95)

  #' wll-01-12-2015: the results bewtween pca.outlier and outl.det.m are
  #' different. Need to be careful.
}
outl.det.m <- function(x, method = "mcd", conf.level = 0.95) {
  require(MASS)
  method <- match.arg(method, c("mve", "mcd", "classical"))
  x <- as.matrix(x)
  covr <- cov.rob(x, method = method) #' MASS. NAs are not allowed.
  dis <- sqrt(mahalanobis(x, center = covr$center, cov = covr$cov))
  cutoff <- sqrt(qchisq(conf.level, ncol(x)))
  outlier <- which(dis > cutoff)

  if (!is.null(names(outlier))) outlier <- names(outlier)
  return(outlier)
}

#' ========================================================================
#' lwc-30-01-2013: Multivariate outlier detection
#' Note: 1.) check out.det in package FIEmspro
#'       2.) check out pca.outlier in package mt
#'       3.) check aq.plot in package mvutlier.
outl.det.m.1 <- function(x, conf.level = 0.95) {
  require(MASS)
  require(robustbase)

  #' covr    <- cov.rob(x, method="mcd")   #' MASS. NAs are not allowed.
  #' covr <- cov.mcd(x)   #' MASS. NAs are not allowed.
  covr <- covMcd(x) #' robustbase
  dis <- sqrt(mahalanobis(x, center = covr$center, cov = covr$cov))
  cutoff <- sqrt(qchisq(conf.level, ncol(x)))
  outlier <- which(dis > cutoff)

  if (!is.null(names(outlier))) outlier <- names(outlier)
  return(outlier)
}

#' ========================================================================
#' lwc-20-01-2013: Univariate outlier detection
#' Note: Specific. Use outl.det.u.1.R instead.
outl.det.u <- function(df) {
  out <- boxplot.stats(df$Value)$out
  sel <- df$Value %in% out
  df <- df[!sel, ]
  return(df)
}

#' ========================================================================
#' lwc-23-01-2013: Outlier detection based on boxplot. Univariate outlier
#' removing based on one variable of a data frame (df).
#' lwc-17-06-2013: add value for general purpose
#' Usage:
#'   outl.det.u.1(iris, value="Sepal.Width")
outl.det.u.1 <- function(df, value, range = 1.5) {
  out <- boxplot.stats(df[[value]], coef = range)$out
  sel <- df[[value]] %in% out
  idx <- which(sel) #' Do we need to output outlier index?

  df <- df[!sel, ]
  return(df)
}

#' ========================================================================
#' lwc-23-01-2013: outlier detection based on boxplot. Package car
#' lwc-07-01-2013: minor changes.
#' Usage:
#'   outl.det.u.2(iris, value="Sepal.Width")
outl.det.u.2 <- function(df, value) {
  require(car)
  rownames(df) <- NULL
  idx <- Boxplot(~ get(value), data = df, id.n = 10) #' Do we need to output outlier index?

  if (!is.null(idx)) df <- df[-as.numeric(idx), ]
  return(df)
}

#' =======================================================================
#' lwc-06-02-2013: Modification of aq.plot in package mvoutlier.
#' Modified:
#'  1.) change delta behivour. It is easily adjusted by the confidence level.
#'  2.) use prcomp instead of princomp
#'  3.) Remove arw(adaptive reweighted estimator for multivariate location and
#'      scatter)
#'  4.) remain two subplots.
#'  5.) Fix a bug of plotting pchisq line
#'  6.) use MASS::cov.mcd instead of robustbase:::covMcd
#' Usages:
if (F) {
  set.seed(134)
  x <- cbind(rnorm(80), rnorm(80), rnorm(80))
  y <- cbind(rnorm(10, 5, 1), rnorm(10, 5, 1), rnorm(10, 5, 1))
  x <- rbind(x, y)
  outl <- aq.plot.1(x, conf.level = 0.95)
}
aq.plot.1 <- function(x, conf.level = 0.975, plotting = TRUE, ...) {
  if (is.vector(x) == TRUE || ncol(x) == 1) {
    stop("x must be at least two-dimensional")
  }

  delta <- qchisq(conf.level, df = ncol(x))

  #' Note that NAs are not allowed.
  covr <- MASS:::cov.mcd(x, ...)
  #' covr <- robustbase:::covMcd(x,...)
  #' lwc: covMcd is not robust.

  dist <- mahalanobis(x, center = covr$center, cov = covr$cov)
  s <- sort(dist, index = TRUE)

  #' --------------
  if (plotting) {
    z <- x
    if (ncol(x) > 2) {
      p <- prcomp(x, center = TRUE, scale. = TRUE)
      z <- data.frame(p$x[, 1:2])
      #' p <- princomp(x,covmat=covr)
      #' z <- p$scores[,1:2]
    }

    X11(width = 8, height = 5)
    par(mfrow = c(1, 2))
    #' -----------
    plot(s$x, (1:length(dist)) / length(dist),
      col = 3,
      xlab = "Ordered squared robust distance",
      ylab = "Cumulative probability", type = "n"
    )
    text(s$x, (1:length(dist)) / length(dist), as.character(s$ix), col = 3, cex = 0.8)

    t <- seq(0, max(dist), length.out = 100)
    ## t <- seq(0,max(dist), by=0.01)    #' lwc: too many points
    lines(t, pchisq(t, df = ncol(x)), col = 6)

    abline(v = delta, col = 5)
    xrange <- par("usr")
    xrange <- xrange[2] - xrange[1]
    text(
      x = delta + xrange / 20, y = 0.4, paste(100 * (pchisq(delta, df = ncol(x))), "% Quantile", sep = ""),
      col = 5, pos = 2, srt = 90, cex = 0.8
    )
    #' -----------
    plot(z,
      col = 3, type = "n",
      main = paste("Outliers based on ", 100 * (pchisq(delta, df = ncol(x))),
        "% quantile",
        sep = ""
      ),
      xlab = "", ylab = ""
    )
    if (any(dist > delta)) {
      text(z[dist > delta, 1], z[dist > delta, 2],
        dimnames(as.data.frame(x)[dist > delta, ])[[1]],
        col = 2, cex = 0.8
      )
      #' lwc: "red", potential outliers
    }
    if (any(dist <= delta)) {
      text(z[dist <= delta, 1], z[dist <= delta, 2],
        dimnames(as.data.frame(x)[dist <= delta, ])[[1]],
        col = 3, cex = 0.8
      )
    }
  }
  outliers <- which(sqrt(dist) > sqrt(delta))
  return(outliers)
}

#' =======================================================================
#' lwc-06-02-2013: wrapper function for outlier detection. Very specific.
#' df <- subset(data.w, Disease=="Sham" & Treatment=="TETA")
#' tmp <- outl.func(df)
outl.func <- function(df, conf.level = 0.975, plotting = TRUE, ...) {
  tmp <- df[, -c(1:3)]
  #' filling missing values
  tmp <- as.data.frame(sapply(tmp, function(x) {
    m <- mean(x, na.rm = TRUE)
    x[is.na(x)] <- m
    x
  }))
  #' Remove outliers
  set.seed(123)
  out <- aq.plot.1(tmp, conf.level = conf.level, plotting = plotting, ...)
  #' require(mvoutlier)
  #' out <- aq.plot(tmp)
  #' out <- which(out[[1]])
  tmp <- df[-out, ]
  res <- list(outliers = out, data = tmp)
}
