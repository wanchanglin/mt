#' =========================================================================
#' Mult-Classifier accuracy estimation with results of accuracy and
#' significant test using ANOVA plus TukeyHSD test.
#' History:
#' 05-12-06: commence
#' 10-01-07: minor changes.
#' 31-01-07: prepare RD file and put it into package 'mt'
#' 24-03-07: Re-calculate improved bootstrap error.
#' 28-03-07: Add overall confusion matrix
#' 01-07-07: Change largely. Produced directly from accest.
maccest.default <- function(dat, cl, method = "svm", pars = valipars(),
                            tr.idx = NULL, comp = "anova", ...) {
  if (missing(dat) || missing(cl)) {
    stop("data set or class are missing")
  }
  if (missing(method)) {
    stop("'method' is missing")
  }

  #' if (!is.factor (cl)) stop("cl must be a factor.")
  cl <- as.factor(cl) #' some classifier need it as factor, such as SVM.

  if (nrow(dat) != length(cl)) stop("mat and cl don't match.")
  if (length(unique(cl)) < 2) {
    stop("Classification needs at least two classes.")
  }
  if (any(is.na(dat)) || any(is.na(cl))) {
    stop("NA is not permitted in data set or class labels.")
  }

  dat <- as.matrix(dat)
  rownames(dat) <- NULL
  n <- nrow(dat)

  #' construct index of train data
  if (is.null(tr.idx)) {
    if (pars$sampling == "cv" && pars$nreps > n) {
      pars$sampling <- "loocv"
    }
    tr.idx <- trainind(cl, pars = pars)
  }
  pars$niter <- length(tr.idx)
  pars$nreps <- length(tr.idx[[1]])

  #' apply single accest for maccest
  res <- lapply(method, function(m) {
    cat("\n--Classifier = :", m)
    flush.console()
    accest(dat, cl, method = m, pars = pars, tr.idx = tr.idx, ...)
  })
  names(res) <- method

  acc <- sapply(res, function(x) x$acc)
  acc.iter <- sapply(res, function(x) x$acc.iter)
  #' acc.std  <- sapply(res, function(x) x$acc.std)
  acc.std <-
    sapply(res, function(x) ifelse(!is.null(x$acc.std), x$acc.std, NA))
  conf <- lapply(res, function(x) x$conf)

  mar <- sapply(res, function(x) ifelse(!is.null(x$mar), x$mar, NA))
  mar.iter <- sapply(res, function(x) x$mar.iter)
  auc <- sapply(res, function(x) ifelse(!is.null(x$auc), x$auc, NA))
  auc.iter <- sapply(res, function(x) x$auc.iter)

  #' significant test
  if (length(method) < 2 || pars$niter < 2) {
    h.test <- NULL
    gl.pval <- NULL
    mc.pval <- NULL
  } else {
    comp <- match.arg(comp, c("anova", "fried"))
    h.test <- switch(comp,
      "anova"  = mc.anova(acc.iter),
      "fried"  = mc.fried(acc.iter)
    )
    gl.pval <- h.test$gl.pval
    mc.pval <- h.test$mc.pval
    comp <- paste(names(h.test)[1], "+", names(h.test)[2], sep = " ")
  }

  #' prepare the returned values
  ret <- list(
    method = method,
    acc = acc,
    acc.std = acc.std,
    acc.iter = acc.iter,
    mar = mar,
    mar.iter = mar.iter,
    auc = auc,
    auc.iter = auc.iter,
    comp = comp,
    h.test = h.test,
    gl.pval = gl.pval,
    mc.pval = mc.pval,
    sampling = switch(pars$sampling,
      "loocv"  = "leave-one-out cross-validation",
      "cv"     = "cross validation",
      "boot"   = "bootstrap",
      "rand"   = "randomised validation (holdout)"
    ),
    niter = pars$niter,
    nreps = pars$nreps,
    tr.idx = tr.idx,
    conf = conf
  )
  if (pars$sampling == "boot") {
    ret$acc.boot <- lapply(res, function(x) x$acc.boot)
  }

  class(ret) <- "maccest"
  return(ret)
}

#' ========================================================================
print.maccest <- function(x, digits = 3, ...) {
  cat("\nMethod:\t\t\t", x$method)
  cat("\nAccuracy:\t\t", round(x$acc, digits))
  #' cat("\nSTD of Accuracy:\t",round(x$acc.std,digits))
  cat("\nAUC:\t\t\t", round(x$auc, digits))
  cat("\nMargin:\t\t\t", round(x$mar, digits))

  cat("\n\nNo. of iteration:\t", x$niter)
  cat("\nSampling:\t\t", x$sampling)
  cat("\nNo. of replications:\t", x$nreps)

  cat("\n")
  if (!is.null(x$h.test)) {
    cat("\nComparison:\t\t", x$comp)
    cat("\nGlobal test p-value:\t", round(x$gl.pval, digits))
    cat("\n\nMultiple comparison p-values:\n")
    print(round(x$mc.pval, digits))
  }

  cat("\n")
  invisible(x)
}

#' ========================================================================
summary.maccest <- function(object, ...) {
  structure(object, class = "summary.maccest")
}

#' ========================================================================
print.summary.maccest <- function(x, digits = 3, ...) {
  print.maccest(x)
  cat("\n\nAccuracy on each iteration:\n")
  print(round(x$acc.iter, digits))
  cat("\nSummary of accuracy on each iteration:\n")
  #' print(summary(x$acc.iter))
  print(apply(x$acc.iter, 2, summary)) #' nicely formatted summary

  invisible(x)
}

#' ========================================================================
boxplot.maccest <- function(x, ...) {
  if (x$niter == 1) {
    stop("Number of iteration (niter) must be greater than 1")
  }

  col <- "lightgray"
  xlab <- "Classifier"
  ylab <- "Accuracy Rate"

  ylim <- c(0, 1.0)

  main <- "Classifier Accuracy"

  boxplot(data.frame(x$acc.iter),
    main = main, col = col, xlab = xlab,
    ylab = ylab, ylim = ylim
  )
}

#' =======================================================================
#' wll-02-12-2006: user defined x-ticks
#' wll-04-12-2006: std bar
#' wll-03-07-2007: Check validity of acc.std
plot.maccest <- function(x, main = NULL, xlab = NULL, ylab = NULL, ...) {
  dots <- list(...)
  #' ylim <- if("ylim" %in% names(dots)) dots$ylim else
  #'        c(min(x$acc - x$acc.std) - 0.1, max(x$acc + x$acc.std) + 0.1)
  ylim <- if ("ylim" %in% names(dots)) {
    dots$ylim
  } else if (!any(is.na(x$acc.std))) {
    c(min(x$acc - x$acc.std) - 0.1, max(x$acc + x$acc.std) + 0.1)
  } else {
    c(min(x$acc) - 0.1, max(x$acc) + 0.1)
  }

  if (is.null(main)) {
    main <- paste("Performance of classifier (Sampling: ", x$sampling, " )",
      sep = ""
    )
  }

  if (is.null(xlab)) xlab <- "Classifier"
  if (is.null(ylab)) ylab <- "Accuracy"

  plot(x$acc,
    type = "o", main = main, xlab = xlab, ylab = ylab, col = "blue",
    ylim = ylim, xaxt = "n", ...
  )

  ns <- 1:length(x$method)
  if (!any(is.na(x$acc.std))) {
    segments(ns, x$acc - x$acc.std, ns, x$acc + x$acc.std)
  }

  # Now draw the x axis with text labels
  axis(1, at = seq(1, length(x$method), by = 1), labels = x$method)
}

#' ========================================================================
maccest <- function(dat, ...) UseMethod("maccest")

maccest.formula <- function(formula, data = NULL, ..., subset,
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
  m$scale <- NULL
  m[[1]] <- as.name("model.frame")
  m$na.action <- na.action
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  x <- model.matrix(Terms, m)
  y <- model.extract(m, "response")
  attr(x, "na.action") <- attr(y, "na.action") <- attr(m, "na.action")

  ret <- maccest.default(x, y, ..., na.action = na.action)

  ret$call <- call
  ret$call[[1]] <- as.name("maccest")
  ret$terms <- Terms
  if (!is.null(attr(m, "na.action"))) {
    ret$na.action <- attr(m, "na.action")
  }
  class(ret) <- c("maccest.formula", class(ret))
  return(ret)
}

#' =====================================================================
#' Estimates the accuracy of pairwise comparison using multi-classifiers.
#' History:
#'   03-12-06: Create
#'   18-12-06: keep all the results
#'   03-07-07: add auc and margin
#' NOTE: It is difficult to provide user-defined data partition before
#'       data extracting for pairwise comparison.
mbinest <- function(dat, cl, choices = NULL, method,
                    pars = valipars(), ...) {
  dat.bin <- .dat.sel(dat, cl, choices = choices)

  #' get length of each pairwise comparison
  len <- sapply(dat.bin$cl, length)

  if (pars$sampling == "cv" && pars$nreps > min(len)) {
    pars$sampling <- "loocv"
  }
  if (pars$sampling == "loocv") pars$niter <- 1

  res <- lapply(names(dat.bin$cl), function(x) {
    cat("\nRun = :", x, "\n")
    flush.console() #' for Windows
    maccest(dat.bin$dat[[x]], dat.bin$cl[[x]],
      method = method,
      pars = pars, ...
    )
  })

  acc <- t(sapply(res, function(x) x$acc))
  auc <- t(sapply(res, function(x) x$auc))
  mar <- t(sapply(res, function(x) x$mar))

  #' get comparison names
  com <- apply(dat.bin$com, 1, paste, collapse = "-")
  names(res) <- rownames(acc) <- rownames(auc) <- rownames(mar) <- com

  ret <- list(
    all = res, com = dat.bin$com, acc = acc, auc = auc,
    mar = mar, method = method, niter = pars$niter,
    sampling = pars$sampling
  )
  if (pars$sampling != "loocv") ret$nreps <- pars$nreps
  return(ret)
}

#' ====================================================================
#' Friedman test + Wilcoxon test for Multi-classifier significant test
#' lwc-14-12-2006: commence
mc.fried <- function(x, p.adjust.method = p.adjust.methods, ...) {
  p.adjust.method <- match.arg(p.adjust.method)

  #' significant test using Friedman test
  f.htest <- friedman.test(x, ...)
  #' global null hypothesis test p value
  gl.pval <- f.htest$p.value

  #' post-hoc test by Wilcoxon test
  dname <- colnames(x)
  n <- nrow(x)
  acc <- as.vector(x)
  algo <- factor(rep(dname, each = n))

  #' t.htest <- pairwise.t.test(acc, algo, p.adj = "bonf",pool.sd = T,...)
  w.htest <- pairwise.wilcox.test(acc, algo, 
	                                p.adjust.method = p.adjust.method, ...)
  lo <- lower.tri(w.htest$p.value, T)
  mc.pval <- w.htest$p.value[lo]
  #' pairwise comparison names
  dname <- dimnames(w.htest$p.value)
  tmp <- outer(dname[[1]], dname[[2]], paste, sep = "-")
  names(mc.pval) <- tmp[lower.tri(tmp, T)]

  ret <- list(
    fried = f.htest, wilcox = w.htest, gl.pval = gl.pval,
    mc.pval = mc.pval
  )
  ret
}

#' ==============================================================
#' ANOVA + TukeyHSD for Multi-classifier significant test.
#' lwc-14-12-2006: commence
mc.anova <- function(x, ...) {
  #' prepare for ANOVA
  dname <- colnames(x)
  n <- nrow(x)
  acc <- as.vector(x)
  algo <- factor(rep(dname, each = n))

  #' ANOVA for the global null hypothesis test
  aov.tab <- summary(fm1 <- aov(acc ~ algo, ...))
  gl.pval <- aov.tab[[1]][1, 5]

  #' post-hoc test using Tukey HSD
  t.htest <- TukeyHSD(fm1, "algo", ...) #' fm1 must be an output of 'aov'
  mc.pval <- t.htest$algo[, 4]
  #' plot(t.htest)

  ret <- list(
    anova = aov.tab, tukey = t.htest, gl.pval = gl.pval,
    mc.pval = mc.pval
  )
  ret
}

#' ========================================================================
#' lwc-19-12-2006: normality test using shpiro.test, and plot boxplot and
#' histogram
#' lwc-26-02-2010: Using lattice. 
#' wl-12-11-2021, Fri: return two lattice objects.
mc.norm <- function(x, ...) {
  x <- data.frame(x)
  #' normality test
  s.test <- lapply(x, function(x) shapiro.test(x))

  rownames(x) <- NULL
  x <- stack(x)
  p.bw <- bwplot(~ values | ind,  data = x, as.table = T, xlab = "", 
                 pch = "|", scales = list(cex = .75, relation = "free"), ...)
  p.hist <- 
  histogram(~ values | ind,
    data = x, as.table = T,
    scales = list(cex = .75, relation = "free"),
    panel = function(x, ...) {
      panel.histogram(x, ...)
      panel.mathdensity(
        dmath = dnorm, col = "black",
        args = list(mean = mean(x), sd = sd(x))
      )
    }, ...
  )

  #' densityplot(~ values | ind, data=x, as.table=T,
  #'             scales=list(cex =.75,relation="free"), plot.points = F)
  #' qqmath(~ values | ind, data=x, as.table=T, #' f.value = ppoints(100),
  #'        scales=list(cex =.75,relation="free"),
  #'        xlab = "Standard Normal Quantiles")

  return(list(s.test = s.test, bwplot = p.bw, histogram = p.hist))
}

#'  1) maccest.default
#'  2) print.maccest
#'  3) summary.maccest
#'  4) print.summary.maccest
#'  5) boxplot.maccest
#'  6) plot.maccest
#'  7) maccest
#'  8) maccest.formula
#'  9) mbinest
#' 10) mc.fried
#' 11) mc.anova
#' 12) mc.norm
