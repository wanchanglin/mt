#' lwc-07-09-2010: Functions not documented
#' wl-20-11-2020, Fri: tidy up

#' ========================================================================
#' wll-23-06-2015: Get correlation coefficient and p-values
#' Note:
#'   This file is modified from 'cor.table' of package 'picante'
#' and 'corr.test' of package 'psych'.
#'   The original implementation is from Bill Venables, the author of R great
#' book MASS. For details, see
#' https://stat.ethz.ch/pipermail/r-help/2001-November/016201.html
#' Arguments:
#'  x - a data frame or matrix for correlation analysis column-wise
#'  cor.method - method for correlation
#'  adj.method - p-value correction method
#'  ... - other parameter for correlation.
#' Values:
#'  r - correlation coefficient
#'  p - statistics matrix, in which the lower trianular is p-values and the
#'      upper triangular is adjusted p-values
cor.tab <- function(x, cor.method = c("pearson", "kendall", "spearman"),
                    adj.method = c(
                      "holm", "hochberg", "hommel",
                      "bonferroni", "BH", "BY", "fdr", "none"
                    ),
                    ...) {
  R <- cor(x, method = cor.method, ...)
  df <- dim(x)[1] - 2

  if (T) {
    t <- R * sqrt(df / (1 - R^2))
    P <- 2 * (1 - pt(abs(t), df)) #' from corr.test: two-tailed
    #' P <- 2*pt(t, df)           #' from cor.table: right-tailed (greater)
    #' P[R>0] <- 2*pt(t[R>0], df,lower.tail=FALSE)
  } else { #' from Bill Venables
    F <- R^2 * df / (1 - R^2)
    P <- 1 - pf(F, 1, df)
  }
  diag(P) <- NA

  #' get adjusted p-values
  idx <- upper.tri(P, diag = FALSE)
  pval <- P[idx]
  padj <- p.adjust(pval, method = adj.method)
  P[upper.tri(P, diag = FALSE)] <- padj

  list(r = R, p = P)
}

#' ========================================================================
#' wll-24-06-2015: Convert a symmetric table(short format) to long format
#' Arguments:
#'   x     - A symmetric matrix-like data set
#'   tri   - Triangular being used
#' Returns:
#'   A data frame of paire-wise comparision
#' Usages
if (F) {
  #' library(pysch)
  #' co <- corr.test(mtcars, method="spearman",adjust="BH")
  #' From pysch: For symmetric matrices, p values adjusted for multiple tests
  #' are reported above the diagonal.
  co <- cor.tab(mtcars, cor.method = "spearman", adj.method = "BH")
  names(co)

  corr <- sym2long(co$r, tri = "upper")
  pval <- sym2long(t(co$p), tri = "upper")
  padj <- sym2long(co$p, tri = "upper")

  tmp.1 <- data.frame(corr, pval, padj)
}
sym2long <- function(x, tri = c("upper", "lower")) {
  tri <- match.arg(tri)

  if (tri == "upper") {
    ind <- lower.tri(x)
    x[ind] <- NA
    diag(x) <- NA
    x <- x[-nrow(x), -1, drop = F]
  } else if (tri == "lower") {
    ind <- upper.tri(x)
    x[ind] <- NA
    diag(x) <- NA
    x <- x[-1, -ncol(x), drop = F]
  } else { #' It never reaches here because of match.arg.
    stop("Invalid method")
  }

  if (F) {
    #' require(reshape)
    res <- melt(x)
    res <- res[complete.cases(res), ]
    colnames(res) <- c("com1", "com2", "var")
  } else {
    idx <- which(!is.na(x), arr.ind = T)
    fs1 <- rownames(x)[idx[, 1]]
    fs2 <- colnames(x)[idx[, 2]]
    res <- data.frame(
      com1 = fs1, com2 = fs2,
      var = x[idx], stringsAsFactors = FALSE
    )
  }
  return(res)
}

#' ========================================================================
#' wll-23-06-2015: get long format of correlation analysis
#' wll-24-11-2015: function 'sym2long' is more general.
#' Arguments:
#'   x     - A matrix-like data set
#' Returns:
#'   A data frame of paire-wise correleation coeficients
cor.long <- function(x, use = "pairwise.complete.obs", method = "pearson") {
  co <- cor(x, use = use, method = method)
  co[upper.tri(co)] <- NA
  diag(co) <- NA
  co <- co[-1, -ncol(co), drop = F]

  if (T) {
    #' remove NAs
    idx <- which(!is.na(co), arr.ind = T)
    fs1 <- rownames(co)[idx[, 1]]
    fs2 <- colnames(co)[idx[, 2]]
    res <- data.frame(
      com1 = fs1, com2 = fs2, cor = co[idx],
      stringsAsFactors = FALSE
    )
  } else {
    require(reshape)
    res <- melt(co)
    res <- res[complete.cases(res), ]
    colnames(res) <- c("com1", "com2", "cor")
  }
  return(res)
}

#' ========================================================================
#' lwc-07-07-2011: batch shifting: remove mean withing each batch/block
#' Arguments:
#'  x - data matrix
#'  y - categorical data for batch/block information
#' References:
#' Silvia Wagner, et.al, Tools in Metabonomics: An Integrated Validation
#' Approach for LC-MS Metabolic Profiling of Mercapturic Acids in Human
#' Urine Anal. Chem., 2007, 79 (7), pp 2918-2926, DOI: 10.1021/ac062153w
#' Publication Date (Web): February 23, 2007
batch.shift <- function(x, y, type = "mean") {
  x <- as.data.frame(x)

  g.mean <- sapply(x, function(x) tapply(x, y, type, na.rm = T))
  g.mean <- sapply(1:ncol(x), function(i) g.mean[, i][y])
  x <- x - g.mean

  return(x)
}

#' ========================================================================
#' lwc-02-06-2011: Relative standard deviation of matrix/data frame in
#' column
rsd <- function(x) {
  mn <- colMeans(x, na.rm = TRUE)
  std <- apply(x, 2, sd, na.rm = TRUE)
  #' std <- sd(x,na.rm=TRUE)           #' sd(<data.frame>) is deprecated.
  res <- 100 * std / mn
  return(res)
}

#' ========================================================================
#' lwc-02-06-2011: Relative standard deviation.
#' lwc-16-06-2011: Re-write according to function sd.
#' Note: The code style steal from function 'sd', i.e.,  provide a vector 
#' function and use apply/sapply for matrix/data frame.
rsd.1 <- function(x, na.rm = TRUE) {
  if (is.matrix(x)) {
    apply(x, 2, rsd.1, na.rm = na.rm)
  } else if (is.vector(x)) {
    100 * (sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm))
  } else if (is.data.frame(x)) {
    sapply(x, rsd.1, na.rm = na.rm)
  } else {
    100 * sd(as.vector(x), na.rm = na.rm) / mean(as.vector(x), na.rm = na.rm)
  }
}

#' ========================================================================
#' lwc-07-12-2010: Order a data frame by its colums. Only returns ord.
#' Note: Hacked from arrange function from package plyr
arrange_ord <- function(df, ...) {
  ord <- eval(substitute(order(...)), df, parent.frame())
  #' unrowname(df[ord, ])
  return(ord)
}

#' ========================================================================
#' lwc-01-12-2011: General matrix plot by lattice
#' 1.) Modified from function plot.ranef.mer in package lme4.
#' 2.) The combination of eval and substitute for generalisation is useful;
#' 3.) Beware the use of as.name.
#' Usages:
#' plot.mat(iris)
#' plot.mat(iris[,1:2])
#' plot.mat(iris[,1, drop=F])
plot.mat <- function(mat, ...) {
  mat <- as.data.frame(mat, stringsAsFactors = FALSE)
  #' to-do: any converting prevention?  (type.convert)
  mat <- Filter(is.numeric, mat)

  cn <- lapply(colnames(mat), as.name)
  switch(min(ncol(mat), 3),
    qqmath(eval(substitute(~x, list(x = cn[[1]]))), mat, ...),
    xyplot(
      eval(substitute(y ~ x, list(y = cn[[1]], x = cn[[2]]))),
      mat, ...
    ),
    splom(~mat, ...)
  )
}

#' ========================================================================
#' lwc-22-06-2011: check whether or not matrix has negative value and return
#'   index of negative columns.
check_neg <- function(x) {
  #' tmp <- unlist(x)
  #' idx <- which(tmp < 0)
  #' tmp[idx]

  tmp <- x < 0
  #' tmp  <- is.nan(y)
  tmp[is.na(tmp)] <- FALSE
  ind <- apply(tmp, 2, function(x) any(x))
  names(ind) <- NULL
  ind <- which(ind)
  if (length(ind) == 0) ind <- NULL
  return(ind)
}

#' ========================================================================
#' lwc-23-01-2013: panel function for qqmathline with confidence interval.
#'   The implementation is done by adding some confidence codes into
#'   panel.qqmathline. The cofidence interval code segment is modified from
#'   qqPlot.default in package car.  Also refer to some QQ functions such as
#'   qqnorm.default(base), qqline(base) and panel.qqmath(lattice)
#' usages:
if (F) {
  qqmath(~ height | voice.part,
    aspect = "xy", data = singer,
    prepanel = prepanel.qqmathline,
    panel = function(x, ...) {
      panel.qqconf(x, ...)
      #' panel.qqmathline(x, ...)
      panel.qqmath(x, ...)
    }
  )

  vp.comb <-
    factor(sapply(strsplit(as.character(singer$voice.part), split = " "), "[", 1),
      levels = c("Bass", "Tenor", "Alto", "Soprano")
    )
  vp.group <-
    factor(sapply(strsplit(as.character(singer$voice.part), split = " "), "[", 2))
  qqmath(~ height | vp.comb,
    data = singer,
    groups = vp.group, auto.key = list(space = "right"),
    aspect = "xy",
    prepanel = prepanel.qqmathline,
    panel = function(x, ...) {
      #' panel.qqmathline(x, ...)
      panel.qqconf(x, ...)
      panel.qqmath(x, ...)
    }
  )
}

panel.qqconf <- function(x, y = x, distribution = qnorm,
                         probs = c(0.25, 0.75),
                         qtype = 7, groups = NULL, conf = 0.95, ...,
                         identifier = "qqmathline") {
  y <- as.numeric(y)
  stopifnot(length(probs) == 2)
  distribution <- lattice:::getFunctionOrName(distribution)
  nobs <- sum(!is.na(y))
  if (!is.null(groups)) {
    panel.superpose(
      x = y, y = NULL, distribution = distribution,
      probs = probs, qtype = qtype, groups = groups, conf = conf,
      panel.groups = panel.qqmathline, ...
    )
  } #' panel.groups = panel.qqconf,...) #' lwc: have problems.
  else if (nobs > 0) {
    yy <- quantile(y, probs, names = FALSE, type = qtype, na.rm = TRUE)
    xx <- distribution(probs)
    r <- diff(yy) / diff(xx)
    panel.abline(c(yy[1] - xx[1] * r, r), ..., identifier = identifier)

    #' get/draw confidence interval.
    b <- r
    a <- yy[1] - xx[1] * b
    P <- ppoints(nobs)
    z <- distribution(P)
    zz <- qnorm(1 - (1 - conf) / 2)
    SE <- (b / dnorm(z)) * sqrt(P * (1 - P) / nobs)
    #' Note: density function should not be fixed. It should consistent
    #' with 'distribution'. For details, see qqPlot.default of car
    #' package.
    fit.value <- a + b * z
    upper <- fit.value + zz * SE
    lower <- fit.value - zz * SE

    panel.lines(z, upper, lty = 2, col = "red", lwd = 2)
    panel.lines(z, lower, lty = 2, col = "red", lwd = 2)
    #' panel.lines(z, upper, lty=2,col="red",lwd=2,...)
    #' panel.lines(z, lower, lty=2,col="red",lwd=2,...)
  }
}

#' ======================================================================
#' lwc-17-06-2010: Normality check by Shapiro test.
#'                 Note that H0 (not H1) is normal distribution.
#'  Arguments:
#'  dat   - data matrix
#'  alpha - confidence level. Default is 0 which means returning all
#'          variables
normality.test <- function(dat, alpha = 0) {
  tmp <- apply(dat, 2, function(y) shapiro.test(y)$p.value)
  tmp <- sort(tmp)
  idx <- tmp >= alpha #' H0 is normal distribution
  idx <- idx[!is.na(idx)] #' lwc-05-02-2010: in case p-val is NaN
  p <- format(tmp[idx], digits = 3) #' pval <- round(tmp[idx], 4)
  res <- list(vars = names(tmp)[idx], pval = p)
  #' res <- do.call(cbind, res)
  #' Note:  should use the above line.
}

#' ========================================================================
#' wll-05-12-2007: Calculate the pattern of missing values.
#' Details:
#' This function is useful for investigating any structure of missing
#' observation in the data.
#' Value:
#'   A matrix with (nrow(x)+1, ncol(x)+1) dimension. Except the last row and
#'   column, each row corresponds to a missing data pattern
#'   (1=observed, 0=missing). The row names shows the number of pattern.
#'   The last row contains the number of mising values
#'   with respect to each column and the last column represent the counts of
#'   each row.
#' See Also:
#'   md.pattern in package mice and prelim.norm in package norm.
#' NOTE: 1.The motivaton of the function is that Ted Harding mentioned that
#'       that prelim.norm can only encode NA-patterns in an R integer for up
#'       to 31 columns. More than that, and it will not work properly or at
#'       all. (http://article.gmane.org/gmane.comp.lang.r.general/55185).
#'       Function md.pattern has also this problem since it modified from
#'       prelim.norm. 2. The function is not sorted at current stage.
#' Usage:
#'   library(mice)
#'   data(nhanes)
#'   md.pattern(nhanes)     #' from mice
#'   mv.pattern(nhanes)
mv.pattern <- function(x) {
  #' The fowwling function is taken from
  #' http://article.gmane.org/gmane.comp.lang.r.general/16575
  "%all.==%" <- function(a, b) apply(b, 2, function(x) apply(t(a) == x, 2, all))

  if (!(is.matrix(x) | is.data.frame(x))) {
    stop("Data should be a matrix or dataframe")
  }

  #' get the pattern of missing values
  mat <- 1 * !is.na(x)
  pattern <- unique(mat)
  counts <- colSums(mat %all.==% t(unique(mat)))
  rownames(pattern) <- counts

  #' -- add some statistics -----
  #' number of missing values with respect to column (variable)
  nmis <- apply(1 * is.na(x), 2, sum)
  #' number of missing values in the pattern
  pmis <- ncol(pattern) - apply(pattern, 1, sum)

  pattern <- rbind(pattern, c(nmis)) #' a trick to take off the row name
  #' using c
  pattern <- cbind(pattern, c(pmis, sum(nmis)))
  pattern
}

#' =========================================================================
#' lwc-22-05-2012: Wrapper function for logistic regression
LogReg <- function(x, y, ...) {
  y <- factor(y)
  if (nlevels(y) != 2) {
    stop("Logistic Regression is only for two class problems")
  }

  xy <- data.frame(x, y)
  res <- glm(y ~ ., data = xy, family = binomial("logit"), ...)

  class(res) <- c("LogReg", class(res))
  return(res)
}

#' =========================================================================
#' lwc-22-05-2012: Predict method for LogReg.
#' Note: 1.) The response variable is a factor with two levels, where the
#'           last level is considered the event (success).
#'       2.) glm models the second factor level. See Details in ?glm.
#'       3.) newdata can be missing and NULL.
predict.LogReg <- function(object, newdata = NULL, ...) {
  #' get the original levels of categorical data
  lev <- levels(factor(object$model[["y"]]))

  #' predict with response
  res <- predict.glm(object, newdata, type = "response", ...)
  cls <- factor(ifelse(res < .5, lev[1], lev[2]))
  prob <- cbind(1 - res, res)
  colnames(prob) <- lev

  res <- list(class = cls, posterior = prob)
  return(res)
}

#' ========================================================================
#' lwc-17-05-2012: Wrapper function for Logistic/Multinomial Regression.
LogRegAnn <- function(x, y, ...) {
  require(nnet, quietly = TRUE)
  xy <- data.frame(x, y)
  res <- multinom(y ~ ., data = xy, ...)
  class(res) <- c("LogRegAnn", class(res))
  return(res)
}

#' ========================================================================
#' lwc-22-05-2012: Predict method for LogRegAnn.
#' Note: 1.) newdata can be missing, but not allowed as NULL.
predict.LogRegAnn <- function(object, newdata, ...) {
  cls <- nnet:::predict.multinom(object, newdata, type = "class")
  prob <- nnet:::predict.multinom(object, newdata, type = "probs")

  #' reshape probs to be consistent with other classifier such as lda and
  #' qda.
  if (nlevels(cls) == 2) {
    prob <- cbind(prob, 1 - prob)
    colnames(prob) <- levels(cls)
  }

  res <- list(class = cls, posterior = prob)
  return(res)
}

#' ========================================================================
#' Predict method for direct use of PLS methods: kernelpls, simpls and
#'   oscorespls, for classification.
#' History:
#'   wll-01-10-2007: commence
#' Note: 1). The coefficients in PLS package are the cumulative
#'           coefficients. It should be taken off Y-means and X-means before
#'           calculating the regression output. The method of 'coef' in PLS
#'           also takes off the two means. For details, see predict.mvr and
#'           coef.mvr.
#'       2.) The block comments code lines are used to predict with models
#'       containing ncomp[1] components, ncomp[2] components, etc. if ncomp
#'       takes form like: ncomp=1:10.
#' Test code for pred.pls
if (F) {
  data(iris)
  x <- as.matrix(subset(iris, select = -Species))
  y <- iris$Species
  y.1 <- mt:::class.ind(y)

  #' ncomp <- min(n - 1, p)

  #' Direct call simpls.fit
  res <- simpls.fit(X = x, Y = y.1, ncomp = 4)
  names(res)
  val <- pred.pls(res, x, ncomp = 4)
  tmp <- exp(val)
  tmp <- tmp / rowSums(tmp)
  pr <- mda:::softmax(tmp)
  #' table(y, pr)
  class.rate(y, pr)

  #' Call wrapper function
  res <- mvr(y.1 ~ x, ncomp = 4)
  val <- predict(res, x, ncomp = 4)
  val <- as.data.frame(val)
  tmp <- exp(val)
  tmp <- tmp / rowSums(tmp)
  colnames(tmp) <- levels(y)
  pr <- mda:::softmax(tmp)
  class.rate(y, pr)
  (z <- plsc(x, y, ncomp = 4))
}

pred.pls <- function(object, newdata, ncomp) {
  nobs <- nrow(newdata)

  if (F) {
    #' Predict with models containing ncomp[1] components,
    #' ncomp[2] components, etc.
    B <- res$coefficients
    dPred <- dim(B)
    dPred[1] <- dim(newdata)[1]
    dnPred <- dimnames(B)
    dnPred[1] <- dimnames(newdata)[1]
    pred <- array(dim = dPred, dimnames = dnPred)
    for (i in seq(along = 1:ncomp)) {
      B0 <- object$Ymeans - object$Xmeans %*% B[, , i]
      B0 <- rep(B0, each = nobs)
      B1 <- newdata %*% B[, , i]
      pred[, , i] <- B1 + B0
    }
  }

  #' Predict with models containing ncomp components.
  B <- res$coefficients[, , ncomp, drop = T]
  B0 <- object$Ymeans - object$Xmeans %*% B
  B0 <- rep(B0, each = nobs)
  B1 <- newdata %*% B
  pred <- B1 + B0

  return(pred)
}

#' ========================================================================
#' lwc-17-04-2007:Computes the p-value for the two sample t-test using a
#'  permutation test. The permutation density can also be plotted.
#' NOTE: 1.) This function is modified from twotPermutation in package DAAG.
#'       2.) There are several packages in R to do permutation test:
#'           exactRankTests, coin, DAAG, BHH2, multtest
#'       3.) This function has been modified to fit S3 method.
test.perm <- function(x, y, alternative = c("two.sided", "less", "greater"),
                      n.perm = 999, plotting = TRUE) {
  alternative <- match.arg(alternative)

  if (missing(x) || missing(y)) {
    stop("x or y missing!")
  }
  if (!is.vector(x) || !is.vector(y)) {
    stop(" x or y must be vector")
  }

  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2
  z <- c(x, y)

  meano <- mean(x) - mean(y) #' original mean difference
  meanp <- rep(0, n.perm) #' permutation mean difference
  for (i in 1:n.perm) {
    idx <- sample(n, n1, replace = FALSE)
    meanp[i] <- mean(z[idx]) - mean(z[-idx])
  }
  #' pval <- (sum(meanp >= abs(meano)) + sum(meanp <= -abs(meano)))/n.perm
  pval <- switch(alternative,
    "two.sided" = {
      (sum(meanp >= abs(meano)) + sum(meanp <= -abs(meano)) + 1) / (n.perm + 1)
    }, "greater" = {
      (sum(meanp >= meano) + 1) / (n.perm + 1)
    }, "less" = {
      (sum(meanp <= meano) + 1) / (n.perm + 1)
    }
  )

  #' FIX-ME: I am not sure that the p-values for the three situations are
  #' correct.

  if (plotting) {
    plot(density(meanp), xlab = "", main = "", yaxs = "i", cex.axis = 0.8)
    abline(v = meano)
    abline(v = -meano, lty = 2)
    mtext(side = 3, line = 0.5, text = expression(bar(x) - bar(y)), at = meano)
    mtext(side = 3, line = 0.5, text = expression(- (bar(x) - bar(y))), at = -meano)
  }
  pval
}

#' ===================================================================
#' lwc-02-02-2007: randomly select k samples from each factor/class
#' Arguments:
#'   y - a factor indicating class info
#'   k - number of samples selected from each class
#'   n - number of replicates
#' Test:
#'   data(abr1)
#'   cl   <- factor(abr1$fact$class)
#'   (tmp <- sub.samp(cl,6))
#'   cl[tmp[[1]]]
#'   table(cl[tmp[[1]]])
#' NOTE: It should extend to chose different number of samples for different
#'       classes.
sub.samp <- function(y, k, n = 10) {
  if (!is.factor(y)) stop("y is not of class factor")
  idx <- 1:length(y)
  g <- levels(y)
  ng <- length(g)

  nidx <- list()
  for (j in 1:ng) {
    nidx <- c(nidx, list(idx[which(y == g[j])]))
  }

  out.idx <- list()
  for (i in 1:n) {
    kidx <- c()
    for (j in 1:ng) {
      kidx <- c(kidx, sample(nidx[[j]], k))
    }
    kidx <- sample(kidx) #' shuffling
    out.idx <- c(out.idx, list(kidx))
  }
  return(out.idx)
}

#' =========================================================================
#' wll-07-08-2014: Mahalanobis distances between groups.
#'  Modified from Chapter 4 Multidimensional Scaling, Applied Multivariate
#'  Analysis by Brian Everitt and Torsten Hothorn, Springer 2011.
#' Usages
if (F) {
  data("skulls", package = "HSAUR2")
  skulls

  mah <- mah.dis(skulls[, -1], skulls[, 1])
  #' mah <- mah.dis.1(skulls[,-1], skulls[,1])
  mds <- cmdscale(as.dist(mah))
  lim <- range(mds) * 1.2
  plot(mds,
    xlab = "Coordinate 1", ylab = "Coordinate 2",
    xlim = lim, ylim = lim, type = "n"
  )
  text(mds, labels = levels(skulls$epoch), cex = 0.7)
}

mah.dis <- function(x, y) {
  #' Get covariance matrices of the data in each group (Not use here)
  grp_var <- tapply(1:nrow(x), y, function(i) var(x[i, ]))

  #' Get common covariance matrix
  com_var <- var(x)

  #' Get group mean
  grp_cen <- sapply(x, function(i) tapply(i, y, mean))

  #' Get Mahalanobis distances between groups
  grp_mah <-
    apply(grp_cen, 1, function(cen) mahalanobis(grp_cen, cen, com_var))

  return(grp_mah)
}

#' =========================================================================
#' wll-28-08-2008: Mahalanobis distances between groups by Gobor
#' wll-06-08-2014: fix a bug
mah.dis.1 <- function(x, y) {
  stopifnot(is.data.frame(x), !missing(y))
  stopifnot(dim(x)[1] != dim(y)[1])
  y <- as.factor(y)
  man <- manova(as.matrix(x) ~ y)
  E <- summary(man)$SS[2] #' Matrix E
  S <- as.matrix(E$Residuals) / man$df.residual
  InvS <- solve(S)

  #' Get group mean
  #' wll-06-08-2014: by(x,y,mean) does not work in the new version of R.
  #' mds  = matrix(unlist(by(x, y, mean)), byrow=T, ncol=ncol(x))
  mds <- sapply(x, function(x) tapply(x, y, mean))
  #' or
  #' tmp <- aggregate(x, list(y), mean)
  #' mds <- tmp[,-1]
  #' rownames(mds) <- tmp[,1]

  f <- function(a, b) {
    mapply(function(a, b) mahalanobis(mds[a, ], mds[b, ], InvS, TRUE), a, b)
  }
  seq. <- seq(length = nrow(mds))
  names(seq.) <- levels(y)
  res <- outer(seq., seq., f)

  #' From Chapter 4 Multidimensional Scaling, Applied Multivariate Analysis
  #' res <- apply(mds, 1, function(cen) mahalanobis(mds, cen, InvS))
}

#' ========================================================================
#' wll-07-08-2014: Mahalanobis distance among units in a dataset or among
#'  observations in two distinct datasets. From R package StatMatch.
#' Usages:
if (F) {
  md1 <- mah.dist.2(iris[1:6, 1:4])
  md2 <- mah.dist.2(data.x = iris[1:6, 1:4], data.y = iris[51:60, 1:4])
  vv <- var(iris[, 1:4])
  md1a <- mah.dist.2(data.x = iris[1:6, 1:4], vc = vv)
  md2a <- mah.dist.2(data.x = iris[1:6, 1:4], data.y = iris[51:60, 1:4], vc = vv)
}

mah.dist.2 <- function(data.x, data.y = NULL, vc = NULL) {
  xx <- as.matrix(data.x)
  if (is.null(data.y)) {
    yy <- as.matrix(data.x)
  } else {
    yy <- as.matrix(data.y)
  }

  if (is.null(vc)) {
    if (is.null(data.y)) {
      vc <- var(xx)
    } else {
      vc <- var(rbind(xx, yy))
    }
  }

  ny <- nrow(yy)
  md <- matrix(0, nrow(xx), ny)
  for (i in 1:ny) {
    md[, i] <- mahalanobis(xx, yy[i, ], cov = vc)
  }
  if (is.null(data.y)) {
    dimnames(md) <- list(rownames(data.x), rownames(data.x))
  } else {
    dimnames(md) <- list(rownames(data.x), rownames(data.y))
  }
  sqrt(md)
}

#' ======================================================================
#' wll-04-02-2008: Venn diagram functions
#' NOTE: venn and incidence.table are from Duncan Murdoch's venn package
#' ======================================================================
venn <- function(id, category, cutoff = 1, duplicates = FALSE, tab, main) {
  if (missing(tab)) {
    #' Create incidence table from id and category
    tab <- incidence.table(as.character(id), category,
      cutoff = cutoff,
      duplicates = duplicates
    )

    if (missing(main)) {
      main <- paste(
        "Count of", deparse(substitute(id)),
        "by", deparse(substitute(category))
      )
    }
  } else if (missing(main)) {
    main <- paste("Venn diagram of", deparse(substitute(tab)))
  }

  #' Convert rows to binary numbers and count them
  index <- tab %*% 2^ (1:ncol(tab) - 1)
  itab <- table(index)

  save <- par(
    pty = "s", # Must be square to get labels right
    mar = c(1, 0, 1, 0) * par("mar")
  ) # Don"t need side margins

  on.exit(par(save))

  if (ncol(tab) == 2) {
    #' Set up coordinates.  xlim and ylim must be same length or coordinate
    #' system will not be isometric
    plot(1, 1,
      xlim = c(-1.3, 2.3), ylim = c(-1.8, 1.8), bty = "n", axes = FALSE,
      type = "n", xlab = "", ylab = "", main = main
    )

    if (!is.na(zero <- itab[as.character(0)])) {
      title(sub = paste(zero, "not shown"))
    }

    #' Plot 2 circles
    cx <- c(0, 1.1)
    cy <- c(0, 0)
    mx <- mean(cx)
    my <- mean(cy)
    symbols(cx, cy, circles = rep(1, 2), inches = FALSE, add = TRUE)

    #' Put counts in the regions
    text(
      c(mx + 2 * (cx[1] - mx), mx + 2 * (cx[2] - mx)),
      c(my + 2 * (cy[1] - my), my + 2 * (cy[2] - my)),
      itab[as.character(c(1, 2))]
    )
    text(mx, my, itab["3"])

    #' Label the circles

    #' text(c(cx[1] -1, cx[2] + 1), c(0, 0),
    #'     pos=c(2, 4), colnames(tab))

    text(c((mx + 3 * (cx[1] - mx) + cx[1]) / 2, (mx + 3 * (cx[2] - mx) + cx[2]) / 2),
      c((my + 3 * (cy[1] - my) + cy[1] - 1.8) / 2, (my + 3 * (cy[2] - my) + cy[2] - 1.8) / 2),
      pos = c(2, 4), colnames(tab)
    )
  } else if (ncol(tab) == 3) {
    #' Set up coordinates.  xlim and ylim must be same length or coordinate
    #' system will not be isometric
    plot(1, 1,
      xlim = c(-1.5, 2.6), ylim = c(-1.5, 2.6), bty = "n", axes = FALSE,
      type = "n", xlab = "", ylab = "", main = main
    )

    if (!is.na(zero <- itab[as.character(0)])) {
      mtext(paste(zero, "not shown"), side = 1)
    }

    #' Plot 3 circles
    cx <- c(0, 1.1, 0.55)
    cy <- c(0, 0, 1.1 * sqrt(3) / 2)
    mx <- mean(cx)
    my <- mean(cy)
    symbols(cx, cy, circles = rep(1, 3), inches = FALSE, add = TRUE)

    #' Put counts in the regions
    text(
      c(mx + 2 * (cx[3] - mx), mx + 2 * (cx[1] - mx), mx + 2 * (cx[2] - mx)),
      c(my + 2 * (cy[3] - my), my + 2 * (cy[1] - my), my + 2 * (cy[2] - my)),
      itab[as.character(c(1, 2, 4))]
    ) #' wll-04-02-2008: fix a bug
    text(
      c(
        mx + (cx[1] + cx[3] - 2 * mx), mx + (cx[2] + cx[3] - 2 * mx),
        mx + (cx[2] + cx[1] - 2 * mx)
      ),
      c(
        my + (cy[1] + cy[3] - 2 * my), my + (cy[2] + cy[3] - 2 * my),
        my + (cy[2] + cy[1] - 2 * my)
      ),
      itab[as.character(c(3, 5, 6))]
    ) #' wll-04-02-2008: fix a bug
    text(mx, my, itab["7"])

    #' Label the circles
    text(c(
      mx + 2.6 * (cx[3] - mx), (mx + 3 * (cx[1] - mx) + cx[1]) / 2,
      (mx + 3 * (cx[2] - mx) + cx[2]) / 2
    ),
    c(
      my + 2.6 * (cy[3] - my), (my + 3 * (cy[1] - my) + cy[1] - 1.2) / 2,
      (my + 3 * (cy[2] - my) + cy[2] - 1.2) / 2
    ),
    pos = c(3, 2, 4), colnames(tab)
    )
  } else {
    stop("Can only Venn 2 or 3 categories")
  }
}

#' ======================================================================
incidence.table <- function(id, category, names = NULL, cutoff = 1,
                            duplicates = FALSE) {
  if (!duplicates) {
    #' Count combinations and convert to TRUE/FALSE
    tab <- table(as.character(id), category)
    tab >= cutoff
  } else {
    #' Count combinations
    tab <- table(as.character(id), category)

    #' Set up matrix with one row per id
    result <- matrix(FALSE, length(id), ncol(tab))

    #' Set appropriate entries TRUE
    for (i in 1:ncol(tab)) {
      result[, i] <- tab[as.character(id), i] >= cutoff
    }

    #' Return nice looking matrix
    rownames(result) <- as.character(names)
    colnames(result) <- colnames(tab)
    result
  }
}

#' ========================================================================
#' lwc-17-05-2011: Write list of data frame or matrix to Excel's XLS.
list2xls <- function(x, filename = "tmp.xls", FreezeRow = 1,
                     row.names = TRUE) {
  tmp <- names(x)
  for (i in tmp) assign(i, x[[i]])
  mt:::WriteXLS(tmp,
    ExcelFileName = filename, FreezeRow = FreezeRow,
    row.names = row.names
  )
  rm(tmp, list = tmp)
  #' Remove tmp and variables whose names are included in tmp
  invisible()
}

#' ========================================================================
#' lwc-06-04-2011: Slight modification of WriteXLS.R (2.1.0,2010-09-18) 
#' Internal function to write data frames to an Excel file using Perl
WriteXLS <- function(x, ExcelFileName = "R.xls", SheetNames = NULL,
                     perl = "perl", verbose = FALSE,
                     Encoding = c("UTF-8", "latin1"), row.names = FALSE,
                     AdjWidth = FALSE, AutoFilter = FALSE,
                     BoldHeaderRow = FALSE, FreezeRow = 0, FreezeCol = 0,
                     envir = parent.frame()) {
  require("WriteXLS", quietly = TRUE) #' added by LWC

  #' If 'x' is a single name, it is either a single data frame or a list of
  #' data frames. If 'x' is >1 names in a character vector, it is presumed to
  #' be a vector of data frame names. If not a list name, create a list of
  #' data frames from the vector, for consistency in subsequent processing.
  if (length(x) == 1) {
    TMP <- get(as.character(x), envir = envir)

    #' is TMP a list and not single data frame
    if ((is.list(TMP)) & (!is.data.frame(TMP))) {
      DF.LIST <- TMP
    } else {
      DF.LIST <- list(TMP)
      names(DF.LIST) <- x
    }
  } else {
    DF.LIST <- sapply(as.character(x), function(x) get(x, envir = envir),
      simplify = FALSE
    )
    names(DF.LIST) <- x
  }

  #' Check to be sure that each element of DF.LIST is a data frame

  #' lwc-06-04-2010: comment the following line for allowing matrix format.
  #' if (!all(sapply(DF.LIST, is.data.frame))) stop("One or more of the
  #' objects named in 'x' is not a data frame or does not exist")

  #' Additional checks for Excel 2003 limitations
  #' 256 columns, including rownames, if included
  #' 65,536 rows (including header row)
  if (!all(sapply(DF.LIST, function(x) (nrow(x) <= 65535) & (ncol(x) <= 256)))) {
    stop("One or more of the data frames named in 'x' exceeds 65535 rows or 256 columns")
  }

  Encoding <- match.arg(Encoding)

  #' -------------------------------------------------------------------------
  #'  Check to see if SheetNames is specified and if so: check for
  #'  duplications they are same length as the number of dataframes CHECK
  #'  TO see if any SheetNames are >31 chars, which is the Excel Limit
  #'  check for invalid characters: []:*?/\ ELSE check to see if first 31
  #'  characters of data frame names are unique
  if (!is.null(SheetNames)) {
    if (any(duplicated(SheetNames))) {
      message("At least one entry in 'SheetNames' is duplicated. Excel worksheets must have unique names.")
      return(invisible(FALSE))
    }

    if (length(DF.LIST) != length(SheetNames)) {
      message("The number of 'SheetNames' specified does not equal the number of data frames in 'x'")
      return(invisible(FALSE))
    }

    if (any(nchar(SheetNames) > 31)) {
      message("At least one of 'SheetNames' is > 31 characters, which is the Excel limit")
      return(invisible(FALSE))
    }

    if (any(grep("\\[|\\]|\\*|\\?|:|/|\\\\", SheetNames))) {
      message("Invalid characters found in at least one entry in 'SheetNames'. Invalid characters are: []:*?/\\")
      return(invisible(FALSE))
    }
    names(DF.LIST) <- SheetNames
  } else {
    if (any(duplicated(substr(names(DF.LIST), 1, 31)))) {
      message("At least one data frame name in 'x' is duplicated up to the first 31 characters. Excel worksheets must have unique names.")
      return(invisible(FALSE))
    }

    if (any(grep("\\[|\\]|\\*|\\?|:|/|\\\\", names(DF.LIST)))) {
      message("Invalid characters found in at least one data frame name in 'x'. Invalid characters are: []:*?/\\")
      return(invisible(FALSE))
    }
  }

  #' Get path to WriteXLS.pl
  Perl.Path <- file.path(path.package("WriteXLS"), "Perl")
  Fn.Path <- file.path(Perl.Path, "WriteXLS.pl")

  #' Get path for Tmp.Dir for CSV files
  Tmp.Dir <- file.path(tempdir(), "WriteXLS")

  #' Remove Tmp.Dir and Files
  clean.up <- function() {
    if (verbose) {
      cat("Cleaning Up Temporary Files and Directory\n\n")
    }
    unlink(Tmp.Dir, recursive = TRUE)
  }

  #' Clean up on function exit
  on.exit(clean.up())

  #' Cleanup now, in case Tmp.Dir still exists from a prior run
  if (file.exists(Tmp.Dir)) {
    if (verbose) {
      cat("Cleaning Up Temporary Files and Directory From Prior Run\n\n")
    }

    unlink(Tmp.Dir, recursive = TRUE)
  }

  #' Create Tmp.Dir for new run
  if (verbose) {
    cat("Creating Temporary Directory for CSV Files: ", Tmp.Dir, "\n\n")
  }

  dir.create(Tmp.Dir, recursive = TRUE)

  #' ---------------------------------------------------------------------
  #'  Write Comma Delimited CSV files
  for (i in seq(along = DF.LIST)) {
    if (verbose) {
      cat("Creating CSV File: ", i, ".csv", "\n", sep = "")
    }

    if (row.names) {
      write.table(DF.LIST[[i]],
        file = paste(Tmp.Dir, "/", i, ".csv", sep = ""),
        sep = ",", quote = TRUE, na = "",
        row.names = TRUE, col.names = NA
      )
    } else {
      write.table(DF.LIST[[i]],
        file = paste(Tmp.Dir, "/", i, ".csv", sep = ""),
        sep = ",", quote = TRUE, na = "", row.names = FALSE
      )
    }
  }

  #' Write 'x' (character vector of data frame names) to file
  #' appending Tmp.Dir and ".csv" to each x
  x <- paste(Tmp.Dir, "/", seq(length(DF.LIST)), ".csv", sep = "")
  write(as.matrix(x), file = paste(Tmp.Dir, "/FileNames.txt", sep = ""))

  #' ---------------------------------------------------------------------
  if (verbose) {
    cat("Creating SheetNames.txt\n")
  }
  write(as.matrix(names(DF.LIST)),
    file = paste(Tmp.Dir, "/SheetNames.txt", sep = "")
  )

  if (verbose) {
    cat("\n")
  }

  #' Call Perl script
  cmd <- paste(perl,
    " -I", shQuote(Perl.Path),
    " ", shQuote(Fn.Path),
    " --CSVPath ", shQuote(Tmp.Dir),
    " --verbose ", verbose,
    " --AdjWidth ", AdjWidth,
    " --AutoFilter ", AutoFilter,
    " --BoldHeaderRow ", BoldHeaderRow,
    " --FreezeRow ", FreezeRow,
    " --FreezeCol ", FreezeCol,
    " --Encoding ", Encoding,
    " ", shQuote(ExcelFileName),
    sep = ""
  )

  #' Call the external Perl script and get the result of the call
  Result <- system(cmd)

  #' Check to see if Result != 0 in the case of the failure of the Perl
  #' script This should also raise an error for R CMD check for package
  #' testing on R-Forge and CRAN
  if (Result != 0) {
    message("The Perl script 'WriteXLS.pl' failed to run successfully.")
    return(invisible(FALSE))
  } else {
    return(invisible(TRUE))
  }
}
#' ==================================================
#' TOC on 25-11-2015
#' ==================================================
#' (1). cor.tab
#' (2). sym2long
#' (3). cor.long
#' (13). batch.shift
#' (14). rsd
#' (15). rsd.1
#' (16). arrange_ord
#' (17). plot.mat
#' (18). check_neg
#' (19). panel.qqconf
#' (20). normality.test
#' (21). mv.pattern
#' (23). LogReg
#' (24). predict.LogReg
#' (25). LogRegAnn
#' (26). predict.LogRegAnn
#' (27). pred.pls
#' (28). test.perm
#' (29). sub.samp
#' (30). mah.dis
#' (31). mah.dis.1
#' (33). mah.dist.2
#' (34). venn
#' (35). incidence.table
