#' =======================================================================
#' Feature frequency and stability of feature ranking.
#' History:
#'   01-02-2007: commence
#'   04-02-2007: feature ranking frequency and overlap rate
#'   10-02-2007: add the stability
#'   09-07-2007: output all frequency
#' Arguments:
#'   x           - a matrix or data frame of feature order
#'   rank.cutoff - top feature order cut-off
#'   freq.cutoff - feature frequency cut-off
#' References:
#' 1.) Chad A Davis, et al., (2006), Reliable gene signatures for microarray
#' classification: assessment of stability and performance. Bioinformatics,
#' vol.22, no.19, pages 2356 - 2363.
#' 2.) Michiels, S.; Koscielny, S. & Hill, C. Prediction of cancer outcome
#' with microarrays: a multiple random validation strategy. Lancet, 365,
#' 488-492
feat.freq <- function(x, rank.cutoff = 50, freq.cutoff = 0.5) {
  x <- as.matrix(x)
  x <- x[1:rank.cutoff, , drop = F]
  tmp <- table(x)
  tmp <- sort(tmp, decreasing = T) / ncol(x)
  fre <- as.vector(tmp) #' strip off a character 'x'
  names(fre) <- names(tmp)
  #' overlap.rate <- sum(fre==1)/rank.cutoff

  res <- list(
    freq.all = fre,
    freq = fre[fre >= freq.cutoff],
    rank.cutoff = rank.cutoff,
    freq.cutoff = freq.cutoff,
    stability = sum(fre) / length(fre)
    #' overlap.fs   = if (overlap.rate) names(fre[fre==1]) else NULL,
  )
  return(res)
}

#' ========================================================================
#' lwc-16-02-2007: Calculate the consensus of feature selection by different
#' methods
#' Internal function.
#' Arguments:
#'  freq - a list consisting of feature frequency more than a threshold
feat.cons <- function(freq, disp = TRUE) {
  #' lwc-18-02-2007: If only retrieve the names, the following code line
  #'                 is enough.
  #' fs.names <- unique(unlist(lapply(freq, names)))

  fs.names <- lapply(freq, function(x) names(x))
  fs.names <- do.call("c", fs.names)
  fs.names <- sort(table(fs.names), decreasing = T)

  mat <- matrix(
    data = NA, nrow = length(fs.names), ncol = length(freq),
    dimnames = list(names(fs.names), names(freq))
  )

  fs.tab <- sapply(names(freq), function(x) {
    dname <- names(freq[[x]])
    if (disp) mat[dname, x] <- freq[[x]] else mat[dname, x] <- 1
    #' mat[dname,x] <- ifelse(disp, freq[[x]],1) #' do not work
    return(mat[, x])
  })

  #' print(fs.tab, na.print="")
  return(fs.tab)
}

#' ========================================================================
#' wll-04-07-2007: Multiple feature selectors with resampling procedures.
#'                 Will give rank table based on the statistics and feature
#'                 consensus table.
#' wll-12-12-2007: add output of fs.freq and R doc
#' wll-29-04-2008: Wrapper function for multiple feature selection with or
#'                 without re-sampling.
#' wll-29-05-2008: Add no-resampling code;
#' wll-21-10-2009: Integrate re-sampling and no-resampling
#' lwc-15-02-2010: drop fs.tab.
#' lwc-25-02-2010: minor changes
#' Arguments:
#'   x         - data matirx or data fram
#'   y         - class labels (factor)
#'   method    - a set of feature selections methods
#'   pars      - validation control parameters.
#'   is.resam  - Boolean indicator of re-sampling.
feat.mfs <- function(x, y, method, pars = valipars(), is.resam = TRUE,
                     ...) {
  res <- lapply(method, function(m) {
    cat("\nFeature Selector = :", m, "\n")
    flush.console()
    if (is.resam) {
      feat.rank.re(x, y, method = m, pars = pars, ...)
    } else {
      #'    model <- if (is.function(m)) m
      #'             else if (is.character(m)) get(m)
      #'             else eval(m)
      model <- get(m)
      model(x, y, ...)
    }
  })
  names(res) <- method

  fs.rank <- sapply(res, function(x) x$fs.rank)
  fs.order <- sapply(res, function(x) x$fs.order)
  fs.order <- as.data.frame(fs.order, stringsAsFactors = F)

  fs.stats <- if (is.resam) {
    sapply(res, function(x) x$fs.stats)
  } else {
    sapply(res, function(x) x$stats)
  }

  feat.res <- list(
    fs.order = fs.order, fs.rank = fs.rank, fs.stats = fs.stats,
    all = res
  )

  return(feat.res)
}

#' ======================================================================
#' lwc-25-02-2010: Calculate frequency and consensus from results of 
#'  feat.mfs' with 'is.resam=TRUE'. 
#' wll-05-12-2015: Should give an error checking here.
feat.mfs.stab <- function(fs.res, rank.cutoff = 20, freq.cutoff = 0.5) {
  order.list <- lapply(fs.res$all, function(x) x$order.list)

  freq.all <- lapply(order.list, function(x) {
    feat.freq(x, rank.cutoff = rank.cutoff, freq.cutoff = freq.cutoff)
  })

  fs.freq <- lapply(freq.all, function(x) x$freq)
  fs.subs <- lapply(freq.all, function(x) names(x$freq))
  fs.stab <- lapply(freq.all, function(x) x$stability)
  fs.cons <- feat.cons(fs.freq)
  #' print(fs.cons,digits=2,na.print="")

  fs <- list(
    fs.freq = fs.freq, fs.subs = fs.subs,
    fs.stab = fs.stab, fs.cons = fs.cons
  )
  return(fs)
}

#' =======================================================================
#' lwc-03-09-2010: Plot the stats values of multiple feature selection.
#' lwc-06-09-2010: add fs.tab
#' wll-05-12-2015: should use this function in data analysis
#' Note:
#' Arguments:
#'  fs.stats  - Stats value of features
#'  cumu.plot - A logical value indicating the cumulative scores should be
#'              plotted.
feat.mfs.stats <- function(fs.stats, cumu.plot = FALSE, main = "Stats Plot",
                           ylab = "Values", xlab = "Index of variable", ...) {

  fs.stats <- as.matrix(fs.stats)
  nam <- colnames(fs.stats)
  fs.stats <- lapply(nam, function(x) { #'  x = nam[1]
    val <- fs.stats[, x] #' if stats is data frame, no names for val.
    val <- sort(val, decreasing = T, na.last = T)
  })
  names(fs.stats) <- nam

  #' Get feautes tab based on stats
  fs.tab <- lapply(fs.stats, function(x) {
    list(fs = names(x), val = x)
  })
  #' fs.tab <- list2df(un.list(fs.tab))

  #' Note-09-03-2010: If you use cumulative scores, you can easily calculate
  #' the numbers, such as fix at 80%.
  if (cumu.plot) {
    st <- lapply(fs.stats, function(x) cumsum(x / sum(x, na.rm = T)))
  } else {
    st <- fs.stats
  }

  #' reshape data for plotting
  st <- do.call(cbind, st)
  st <- cbind(st, idx = 1:nrow(st))
  st <- data.frame(st, stringsAsFactors = F)

  #' wl-06-11-2021, Sat: use base R function to get long format
  st_long <- with(st, data.frame(idx = st$idx, stack(st, select = -idx)))
  # st_long <- data.frame(idx = st$idx,stack(st,select = -idx))
  st_long <- st_long[c("idx", "ind", "values")]
  names(st_long) <- c("idx", "variable", "value")

  #' st_l <- reshape::melt(st, id = "idx")
  st.p <- xyplot(value ~ idx | variable,
    data = st_long,
    as.table = T, type = c("g", "l"),
    scales = list(cex = .75, relation = "free"),
    par.strip.text = list(cex = 0.65),
    ylab = ylab, xlab = xlab, main = main, ...
  )
  st.p

  res <- list(stats.tab = fs.tab, stats.long = st_long, stats.p = st.p)
  return(res)
}

#' =======================================================================
#' wll-06-11-2008: Use Borda count to get the final feature order
#' Note: Previous name is fs.agg
feat.agg <- function(fs.rank.list) {
  fs.score <- apply(fs.rank.list, 1, sum)
  fs.order <- order(fs.score, decreasing = F) #' order from best to worst.
  fs.rank <- order(fs.order, decreasing = F) #' feature rank score.
  names(fs.rank) <- rownames(fs.rank.list)
  temp <- names(fs.rank[fs.order])
  if (!is.null(temp)) fs.order <- temp
  return(list(fs.order = fs.order, fs.rank = fs.rank))
}

#' ======================================================================
#' wll-20-03-2007: resampling-based on feature ranking/selection
#' wll-23-07-2008: some change in loops handling
feat.rank.re <- function(x, y, method, pars = valipars(), tr.idx = NULL,
                         ...) {
  #' validity checking
  if (missing(x) || missing(y)) {
    stop("data set or class are missing")
  }
  if (length(dim(x)) != 2) {
    stop("'x' must be a matrix or data frame")
  }
  y <- as.factor(y) #' some classifier need it as factor, such as SVM.

  if (nrow(x) != length(y)) stop("x and y don't match.")
  if (length(unique(y)) < 2) {
    stop("Classification needs at least two classes.")
  }
  if (any(is.na(x)) || any(is.na(y))) {
    stop("NA is not permitted in data set or class labels.")
  }

  n <- nrow(x) #' number of samples
  p <- ncol(x) #' size of feature

  #' construct index of train data
  if (is.null(tr.idx)) {
    if (pars$sampling == "cv" && pars$nreps > n) {
      pars$sampling <- "loocv"
      pars$niter <- 1
    }
    if (pars$sampling == "cv" && pars$nreps < 2) {
      stop("Number of fold (nreps) for cv must greater than 1")
    }
    tr.idx <- trainind(y, pars = pars)
  } else {
    pars$sampling <- c("user")
  }
  pars$niter <- length(tr.idx)
  pars$nreps <- length(tr.idx[[1]])

  #' feature selection with re-sampling
  cat("Iter (", pars$niter, "):", sep = "")
  res.all <- lapply(1:pars$niter, function(i) {
    cat(" ", i, sep = "")
    flush.console()
    train.ind <- tr.idx[[i]]
    res <- lapply(1:pars$nreps, function(j) {
      x.tr <- x[train.ind[[j]], , drop = F]
      y.tr <- y[train.ind[[j]]]
      do.call(method, c(list(x = x.tr, y = y.tr), list(...)))
    })
    names(res) <- paste("Reps", 1:pars$nreps, sep = "_")
    res
  })
  cat("\n")
  names(res.all) <- paste("Iter", 1:pars$niter, sep = "_")

  rank.list <-
    lapply(res.all, function(x) as.data.frame(sapply(x, function(y) y$fs.rank)))

  order.list <-
    lapply(res.all, function(x) as.data.frame(sapply(x, function(y) y$fs.order)))
  stats.list <-
    lapply(res.all, function(x) as.data.frame(sapply(x, function(y) y$stats)))
  rank.list <- do.call("cbind", rank.list)
  order.list <- do.call("cbind", order.list)
  stats.list <- do.call("cbind", stats.list)

  fs.stats <- apply(stats.list, 1, mean)

  #' Use Borda count to get the final feature order
  fs.score <- apply(rank.list, 1, sum)
  fs.order <- order(fs.score, decreasing = F) #' feature order from best to worst.
  fs.rank <- order(fs.order, decreasing = F) #' feature rank score.
  names(fs.rank) <- rownames(rank.list)
  temp <- names(fs.rank[fs.order])
  if (!is.null(temp)) {
    fs.order <- noquote(temp)
  } #' lwc-16-02-2010: Should we remove noquote?

  res <- list(
    method = method,
    fs.order = fs.order, #' feature order
    fs.rank = fs.rank, #' feature rank
    fs.stats = fs.stats, #' means of stats
    rank.list = rank.list, #' full feature rank list
    order.list = order.list, #' full feature order list
    pars = pars, #' resampling parameters
    tr.idx = tr.idx, #' index of training samples.
    all = res.all
  ) #' all results of re-sampling
  return(res)
}

#' =======================================================================
#' wll-31-10-2007: feature selection using VIP of PLS.
#' NOTE: This function supports multiple response, Y (i.e. dummy matrix for
#' discriminat). The Mahalanobis distance of VIP is computed as the values
#' for selection of feature/variables.
#' References:
#' 1. Svante Wold et. al., PLS-regression: a basic tool of chemometrics.
#'    Chemometrics and Intelligent Laboratory Systems 58(2001), 109-130.
#' 2. Sarah J. Dixon, et.al., Pattern recognition of gas chromatography mass
#'    spectrometry of human volatiles in sweat to distinguish the sex of
#'    subjects and determine potential discriminatory marker peaks.
#'    Chemometrics and Intelligent Laboratory Systems 87(2007), 161-172.
#' 3. Chong, Il-Gyo and Jun, Chi-Hyuck, Performance of some variable
#'    selection methods when multicollinearity is present, Chemometrics and
#'    Intelligent Laboratory Systems 78(2005), 103-112.
fs.plsvip <- function(x, y, ncomp = 10, ...) {
  if (!is.data.frame(x)) x <- as.data.frame(x)
  if (length(y) != nrow(x)) {
    stop("x and y is not consistent.")
  }

  val <- plsc(x, y, pls = "oscorespls", ncomp = ncomp)
  #' NOTE: Only NIPLS supports VIP values.
  pls <- val$pls.out

  #' calculate SS
  y.lo <- (unclass(pls$Yloadings))^2
  x.sc <- colSums(pls$scores^2)
  x.sc <- matrix(x.sc, nrow = nrow(y.lo), ncol = ncol(y.lo), byrow = T)
  SS <- y.lo * x.sc #' not matrix product %*%.

  #' calculate normalised squared weight
  W <- (unclass(pls$loading.weights))^2
  sumW <- colSums(W)
  sumW <- matrix(sumW, nrow = nrow(W), ncol = ncol(W), byrow = T)
  W <- W / sumW

  SSW <- W %*% t(SS)

  sumSS <- apply(SS, 1, sum)
  sumSS <- matrix(sumSS, nrow = nrow(SSW), ncol = ncol(SSW), byrow = T)

  vip <- sqrt(nrow(SSW) * (SSW / sumSS))
  #' vip <- rowSums(abs(vip))

  #' Mahalanobis distances
  val <- sqrt(mahalanobis(vip, colMeans(vip), cov(vip), inverted = T))

  #' feature rank and feature order
  fs.order <- order(val, decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(val)
  nam <- names(val[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(fs.order = fs.order, fs.rank = fs.rank, stats = val)
  return(res)
}

#' =======================================================================
#' wll-31-10-2007: feature selection using VIP of PLS.
#' NOTE: This function supports multiple response, Y (i.e. dummy matrix for
#'       discriminat). The final VIP is the means of absolute value of VIP.
fs.plsvip.1 <- function(x, y, ncomp = 10, ...) {
  if (!is.data.frame(x)) x <- as.data.frame(x)
  if (length(y) != nrow(x)) {
    stop("x and y is not consistent.")
  }

  val <- plsc(x, y, pls = "oscorespls", ncomp = ncomp)
  #' NOTE: Only NIPLS supports VIP values.
  pls <- val$pls.out

  #' calculate SS
  y.lo <- (unclass(pls$Yloadings))^2
  x.sc <- colSums(pls$scores^2)
  x.sc <- matrix(x.sc, nrow = nrow(y.lo), ncol = ncol(y.lo), byrow = T)
  SS <- y.lo * x.sc #' not matrix product %*%.

  #' calculate normalised squared weight
  W <- (unclass(pls$loading.weights))^2
  sumW <- colSums(W)
  sumW <- matrix(sumW, nrow = nrow(W), ncol = ncol(W), byrow = T)
  W <- W / sumW

  SSW <- W %*% t(SS)

  sumSS <- apply(SS, 1, sum)
  sumSS <- matrix(sumSS, nrow = nrow(SSW), ncol = ncol(SSW), byrow = T)

  vip <- sqrt(nrow(SSW) * (SSW / sumSS))

  val <- rowMeans(abs(vip))

  #' Mahalanobis distances
  #' val <- sqrt(mahalanobis(vip, colMeans(vip), cov(vip), inverted=T))

  #' feature rank and feature order
  fs.order <- order(val, decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(val)
  nam <- names(val[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(fs.order = fs.order, fs.rank = fs.rank, stats = val)
  return(res)
}

#' =======================================================================
#' wll-29-10-2007: feature selection using VIP of PLS.
#' NOTE: To calculate VIP, two conditions needs to satisfy:
#'       1.) PLS algorithm is NIPLS;
#'       2.) Y must be single vector, not multiple vector, i.e matrix. Hence
#'       for classification, the coding of Y as a single vector is not
#'       efficient, especially for the multi-class problem. For two-class
#'       problem, coding of 0 and 1 or 1 and -1 may be OK for a single y
#'       vector.
fs.plsvip.2 <- function(x, y, ncomp = 10, ...) {
  if (!is.data.frame(x)) x <- as.data.frame(x)
  if (length(y) != nrow(x)) {
    stop("x and y is not consistent.")
  }

  #' convert to numeric, especially for factor.
  y <- as.numeric(y)
  #' NOTE: need to consider a nice way to convert

  pls <- oscorespls.fit(as.matrix(x), y, ncomp = ncomp)
  #' NOTE: Only NIPLS supprots VIP values.

  #' VIP values (taken from http://mevik.net/work/software/pls.html)
  SS <- c(pls$Yloadings)^2 * colSums(pls$scores^2)
  Wnorm2 <- colSums(pls$loading.weights^2)
  SSW <- sweep(pls$loading.weights^2, 2, SS / Wnorm2, "*")

  val <- sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
  val <- val[ncomp, ] #' extract VIP values for ncomp components

  #' feature rank and feature order
  fs.order <- order(val, decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(val)
  nam <- names(val[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(fs.order = fs.order, fs.rank = fs.rank, stats = val)
  return(res)
}

#' =========================================================================
#' wll-29-10-2007: feature selection using regression coefficient of PLS.
#' wll-19-11-2015: add 'pls:::' in the front of 'coef.mvr' since it 'hides'
#'   in the new version of 'pls'
#' NOTE: I try to use robust estimation of center and covarian matrix by
#'       cov.rob in package MASS. But it seems the collinearity problem to
#'       appear. Therefore, the simple Mahalanobis distance is used.
#' NOTES: 1.) Mahalanobis distance and leverage are often used to detect
#'       outliers especially in the development of linear regression models.
#'       A point that has a greater Mahalanobis distance from the rest of
#'       the sample population of points is said to have higher leverage
#'       since it has a greater influence on the slope or coefficients of
#'       the regression equation. (From
#'       http://en.wikipedia.org/wiki/Mahalanobis_distance)
fs.pls <- function(x, y, pls = "simpls", ncomp = 10, ...) {
  if (!is.data.frame(x)) x <- as.data.frame(x)
  if (length(y) != nrow(x)) {
    stop("x and y is not consistent.")
  }

  val <- plsc(x, y, pls = pls, ncomp = ncomp, ...)
  #' wl-08-11-2021, Mon: Use this one
  coe <- drop(coef(val$pls.out, ncomp = val$ncomp))
  # coe <- drop(pls:::coef.mvr(val$pls.out, ncomp = val$ncomp))


  #' lwc-14-06-2010: After runing plsc, ncomp may change; hence ncomp here
  #'                 use val$ncomp

  #' Mahalanobis distances
  val <- sqrt(mahalanobis(coe, colMeans(coe), cov(coe), inverted = T))
  #' val <- sapply(as.data.frame(t(coe)), function(x) sqrt(sum(x^2)))
  #' val <- sapply(as.data.frame(t(coe)), function(x) sum(abs(x)))

  #' feature rank and feature order
  fs.order <- order(val, decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(val)
  nam <- names(val[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(fs.order = fs.order, fs.rank = fs.rank, stats = val)
  return(res)
}

#' ======================================================================
#' wll-30-10-2007: feature selection using loadings of PCA.
#' lwc-22-09-2011: a bug fixed.
#' NOTE: 1. To check the eignenvalue, use screeplot().
#'       2. If it combines with other supervised methods such as fs.rf and
#'       fs.anove as the input methods for feat.mfs and feat.mfs.1, the
#'       'thres' should be given explicitly in case of conflicting with 'y'.
#'       e.g., fs.m <- c("fs.anova", "fs.rf", "fs.pca") feat.mfs.1(dat, cls,
#'       method=fs.m, is.resam=F, thres=0.8)
fs.pca <- function(x, thres = 0.8, ...) {
  x <- as.matrix(x)

  obj <- prcomp(x, ...)
  vars <- obj$sdev^2
  vars <- vars / sum(vars) #' Proportion of Variance
  cumvars <- cumsum(vars) #' Cumulative Proportion
  names(cumvars) <- colnames(obj$rotation)

  id <- which(cumvars >= thres)[1]
  if (id == 1) id <- 2 #' lwc-22-09-2011:
  lo <- obj$rotation[, 1:id] #' selected loadings

  #' Mahalanobis distances
  #' rob <- cov.rob(lo, method="mve")    #' c("mve", "mcd", "classical")
  #' val <- sqrt(mahalanobis(lo, rob$center, rob$cov,tol = 1e-7))

  val <- sqrt(mahalanobis(lo, colMeans(lo), cov(lo), inverted = T))
  #' val <- sapply(as.data.frame(t(lo)), function(x) sqrt(sum(x^2)))
  #' val <- sapply(as.data.frame(t(lo)), function(x) sum(abs(x)))

  #' feature rank and feature order
  fs.order <- order(val, decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(val)
  nam <- names(val[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(fs.order = fs.order, fs.rank = fs.rank, stats = val)
  return(res)
}

#' =======================================================================
#' lwc-12-04-2007: feature selection using ratio of between-group
#'                 to within-group sums of squres (BW).
#' References: Dudoit S, Fridlyand J, Speed T.P. Comparison of
#'   discrimination methods for the classification of tumors using gene
#'   expression data. J Amer Statist Assoc 2002, 97:7.
#' NOTE: Someone claims that BW ratio for multiclass classification is a
#'       modification of the F-ratio statistics for one-way ANOVA.
fs.bw <- function(x, y, ...) {
  if (!is.data.frame(x)) x <- as.data.frame(x)
  if (length(y) != nrow(x)) {
    stop("x and y is not consistent.")
  }

  bw <- sapply(x, function(z) {
    #' z <- x[,1]      #' for debug
    mn.all <- mean(z)
    mn.grp <- tapply(z, y, mean)
    tmp.1 <- tmp.2 <- 0
    for (i in 1:length(z)) {
      cls <- y[i] #' which class
      tmp.1 <- tmp.1 + (mn.grp[[cls]] - mn.all)^2
      tmp.2 <- tmp.2 + (z[i] - mn.grp[[cls]])^2
    }
    tmp.1 / tmp.2
  })

  fs.order <- order(bw, decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(bw)
  nam <- names(bw[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(fs.order = fs.order, fs.rank = fs.rank, stats = bw)
  return(res)
}

#' =========================================================================
#' lwc-06-04-07: Feature selection using RELIEF
#' wll-06-04-07: According to the original algorithm, a random sample should
#' be drawn in each computation. But result of each run will be different.
#' So I change the idea and each sample will be used to update the weight.
#' wll-15-10-07: 1. Extend to ReliefF, in which main ideas are that there
#'               are k (k>=1, default as 10) nearest hits/misses and all the
#'               hits and misses are averaged. 2. Add the user defined
#'               number of instances to sample. Default is all instances will
#'               be used.
#' References:
#' 1.) KIRA, K. and RENDEL, L. (1992). The Feature Selection Problem :
#' Traditional Methods and a new algorithm. Proc. Tenth National Conference
#' on Artificial Intelligence, MIT Press, 129-134.
#' 2.) KONONENKO, I., SIMEC, E., and ROBNIK-SIKONJA, M. (1997). Overcoming
#' the myopia of induction learning algorithms with RELIEFF. Applied
#' Intelligence Vol7, 1, 39-55.
#' 3.) Igor Kononenko, Estimating Attributes: Analysis and Extensions of
#' RELIEF, European Conference on Machine Learning, Ed. Francesco Bergadano
#' and Luc De Raedt, 1994, 171-182, Springer
#' 4.) MARKO ROBNIK-SIKONJA and IGOR KONONENKO, Theoretical and Empirical
#' Analysis of ReliefF and RReliefF, Machine Learning, 53, 23<U+FFFD>C69, 2003
fs.relief <- function(x, y, m = NULL, k = 10, ...) {

  #' Find the nearest neighbors from a matrix
  nearest <- function(x, mat, k = 10) {
    #' Euclidean distance
    dis <- sapply(as.data.frame(t(mat)), function(y) sqrt(sum((x - y)^2)))
    k <- min(k, length(dis)) #' wll-21-03-2008: fix a bug spotted by Ian Scott.
    #' ind  <- which.min(dis)
    ind <- sort.list(dis)[1:k]

    return(mat[ind, , drop = F])
  }

  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.factor(y)) y <- as.factor(y)
  if (length(y) != nrow(x)) {
    stop("x and y is not consistent.")
  }

  n <- nrow(x)
  p <- ncol(x)
  gp <- levels(y)
  prio <- table(y) / n #' Computing the prior

  #' Calculating the range of each feature. range = Max - Min
  rng <- sapply(as.data.frame(x), function(x) diff(range(x)))

  if (is.null(m)) {
    m <- n
  } else {
    m <- min(m, n)
  }
  idx <- sample(1:n, m, replace = F)

  weight <- rep(0, p)
  for (i in idx) {
    #' split x by group
    dat <- split.data.frame(x[-i, , drop = F], y[-i])

    #' find nearest neighbours
    near <- lapply(dat, function(z) nearest(x[i, ], z, k = k))

    hit <- gp[gp == y[i]]
    miss <- gp[gp != y[i]]

    delta <- rep(0, p)
    for (j in 1:p) {
      diff.hit <- -mean(abs(x[i, ][j] - near[[hit]][, j, drop = T]))
      diff.miss <- lapply(miss, function(z) {
        prio[z] * mean(abs(x[i, ][j] - near[[z]][, j, drop = T]))
      })
      diff.miss <- do.call("sum", diff.miss)
      diff.miss <- diff.miss / (1 - prio[hit])
      delta[j] <- (1 / m) * ((diff.hit + diff.miss) / rng[j])
    }
    #' updat weight
    weight <- weight + delta
  }

  names(weight) <- colnames(x)
  fs.order <- order(weight, decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(weight)
  nam <- names(weight[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(fs.order = fs.order, fs.rank = fs.rank, stats = weight)
  return(res)
}

#' ========================================================================
#' RFE-SVM feature selection
#' History:
#'   lwc-02-11-2006:
#'   lwc-15-11-2006: fix fs.len as "power2"
#'   lwc-12-01-2007: re-write
#'   lwc-17-01-2007: sequence of number of features in decreasing order
#'   lwc-18-01-2008: dots argument has a problem. Bu I think it is problem of
#'                   svm, not dots' problem.  So I have to strip off ... here.
#'                   Does svm treat ... and list(...) differently?
fs.rfe <- function(x, y, fs.len = "power2", ...) {
  #' lwc-16-01-2007: avoid multiple actual arguments in calling svm.
  #' dots <- list(...)
  #' if(hasArg(kernel)) dots$kernel <- NULL

  #' LWC-24-10-2006: Calculates the primal variables w which are stored in
  #' warray
  wts.svm <- function(x) {
    #' warray[k,l,] is the weight vector for the binary pb class k against
    #' class l
    ncl <- length(x$labels)
    classk <- rep(1:ncl, x$nSV)
    p <- dim(x$SV)[2]

    #' array of the weight vectors
    warray <- array(0, dim <- c(ncl, ncl, p))
    #' loop to use the coefs
    for (s in 1:dim(x$SV)[1]) {
      for (co in 1:(ncl - 1)) {
        #' find the two class problem
        k <- classk[s]
        l <- ((1:ncl)[-k])[co]
        warray[k, l, ] <- warray[k, l, ] + x$coefs[s, co] * x$SV[s, ]
        warray[l, k, ] <- warray[l, k, ] + x$coefs[s, co] * x$SV[s, ]
      }
    }
    #' return twice the sum of the absolute value of primal variables w
    wts <- apply(abs(warray), 3, sum)
    return(wts)
  }

  y <- as.factor(y) #' 31-03-2007: must be factor if for classification.
  p <- ncol(x)
  fs.order <- seq(1, p)

  #' get feature lengths for SVM-RFE computation.
  len <- get.fs.len(p, fs.len = fs.len)
  len <- sort(len, decreasing = T) #' must be decreasing order for SVM-RFE

  nlen <- length(len)
  for (i in 1:nlen) {
    #' extract index of feature for this length.
    sel <- fs.order[1:len[i]]
    #' call SVM with linear kernel.
    model <- svm(x[, sel, drop = F], y, kernel = "linear")
    #' model <- svm(x[,sel,drop=F], y, kernel = "linear",dots)
    #' calculate the weights
    wts <- wts.svm(model)
    #' sort the feature based on the weights
    ord <- order(wts, decreasing = TRUE)
    #' update the feature set
    fs.order[1:len[i]] <- sel[ord]
  }

  fs.rank <- order(fs.order)

  names(fs.rank) <- colnames(x)
  nam <- colnames(x)[fs.order]
  if (!is.null(nam)) fs.order <- nam

  #'  wll-05-07-2007: add stats for consistent with other methods.
  #'                  No other purpose.
  stats <- length(fs.rank) - fs.rank

  res <- list(fs.order = fs.order, fs.rank = fs.rank, stats = stats)
  return(res)
}

#' ========================================================================
#' SNR for feature selection
#' wll-20-03-2007: SNR is only for two-class classification.
#' wll-29-10-2008: is similar to Fisher criterion rate: (m1-m2)^2/(v1+v2)
fs.snr <- function(x, y, ...) {
  y <- as.factor(y)
  if (length(levels(y)) != 2) {
    stop("'y' must have two classes")
  }

  g.mn <- sapply(data.frame(x), function(x) tapply(x, y, mean))
  g.sd <- sapply(data.frame(x), function(x) tapply(x, y, sd))
  snr <- abs(g.mn[1, ] - g.mn[2, ]) / (g.sd[1, ] + g.sd[2, ])

  fs.order <- order(abs(snr), decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(snr)
  nam <- names(snr[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(fs.order = fs.order, fs.rank = fs.rank, stats = abs(snr))
  return(res)
}

#' ========================================================================
#' AUC for feature selection
#' wll-20-03-2007:  AUC is only for two-class classification.
fs.auc <- function(x, y, ...) {
  y <- as.factor(y)
  if (length(levels(y)) != 2) {
    stop("'y' must have two classes")
  }

  levels(y) <- c(0, 1) #' change levels as 1,0
  y <- as.numeric(levels(y))[as.integer(y)]

  auc <- sapply(as.data.frame(x), function(x) {
    y <- y[order(x, decreasing = TRUE)]
    tmp <- cumsum(y) / sum(y)
    mean(tmp[y == 0])
  })

  #' library(limma)
  #' auc  <- sapply(as.data.frame(x),function(z) auROC(y,z))
  #' library(verification)
  #' auc  <- sapply(as.data.frame(x),function(z) roc.area(y,z)$A)
  #' library(ROCR)
  #' auc  <- sapply(as.data.frame(x),function(z)
  #'           as.numeric(performance(prediction(z,y),measure="auc")@y.values))

  auc[auc < 0.5] <- 1 - auc[auc < 0.5]

  fs.order <- order(abs(auc), decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(auc)
  nam <- names(auc[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(fs.order = fs.order, fs.rank = fs.rank, stats = abs(auc))
  return(res)
}

#' ========================================================================
#' randomForest for feature selection
#' lwc-20-03-2007:Commence
#' wll-08-12-2015:Use scaled important scores.
#' Note-08-12-2015:
#'  1.) If there are large quantity of zeros in the important score which in
#'      turn lead to ties of rank list, the results of reampling in which
#'      rank aggregation is used are not reasonable.
#'  2.) Random Forest is random method, which leads to different results for
#'      different runs even if the random seed has been set by set.seed().
#'  3.) The application of 'fs.rf' should be limited.
fs.rf <- function(x, y, ...) {
  tmp <- randomForest(x, y, importance = T, ...)
  meas <- tmp$importance[, ncol(tmp$importance) - 1]
  meas[meas <= 0] <- 0
  #' Or use the following two lines
  if (F) {
    meas <- importance(tmp, type = 1, scale = TRUE)
    meas <- meas[, 1]
  }

  fs.order <- order(meas, decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(meas)
  nam <- names(meas[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(fs.order = fs.order, fs.rank = fs.rank, stats = meas)
  return(res)
}

#' ========================================================================
#' wll-29-10-2008: Another version of RF for feature selection based on
#'                 successively eliminating the least important variables.
fs.rf.1 <- function(x, y, fs.len = "power2", ...) {
  y <- as.factor(y)
  p <- ncol(x)
  fs.order <- seq(1, p) #' initialisation

  #' get feature lengths
  len <- get.fs.len(p, fs.len = fs.len)
  len <- sort(len, decreasing = T) #' must be decreasing order

  nlen <- length(len)
  for (i in 1:nlen) {
    #' extract index of feature for this length.
    sel <- fs.order[1:len[i]]
    #' call randomForest
    rf <- randomForest(x[, sel, drop = F], y,
      importance = T,
      keep.forest = FALSE, ...
    )
    imp <- importance(rf, type = 1, scale = TRUE)
    imp <- imp[, 1]
    #' sort the feature based on the scaled important scores
    ord <- order(imp, decreasing = T, na.last = T)
    #' update the feature set
    fs.order[1:len[i]] <- sel[ord]
  }

  fs.rank <- order(fs.order)
  names(fs.rank) <- colnames(x)
  nam <- colnames(x)[fs.order]
  if (!is.null(nam)) fs.order <- nam

  #'  Add stats for consistent with other methods.  No other purpose.
  stats <- length(fs.rank) - fs.rank
  res <- list(fs.order = fs.order, fs.rank = fs.rank, stats = stats)
  return(res)
}

#' =========================================================================
#' Welch test for feature selection
#' lwc-05-01-2007: commence
#' lwc-16-03-2007: minor change. Note that do not use sort with na.last=T
#' lwc-19-05-2008: Change oneway.test as t.test with dot arguments. So this
#'                 method supports paired t.test (Welch). And also strip off
#'                 noquote.
#' wll-17-06-2008: change data.frame as as.data.frame in case of data frame
#'                 names changed. (e.g., names of a data frame is
#'                 numeric-like characteristic)
fs.welch <- function(x, y, ...) {
  tmp <- sapply(as.data.frame(x), function(x) {
    tmp <- t.test(x ~ y, var.equal = F, ...)
    #' tmp <- oneway.test(x ~ y,var.equal=F)
    c(tmp$statistic, tmp$p.value)
  })

  stats <- tmp[1, ]
  pval <- tmp[2, ]

  fs.order <- order(abs(stats), decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(stats)
  nam <- names(stats[fs.order])
  if (!is.null(nam)) fs.order <- nam #'  fs.order <- noquote(nam)

  res <- list(
    fs.order = fs.order, fs.rank = fs.rank, stats = abs(stats),
    pval = pval
  )
  return(res)
}

#' =========================================================================
#' wll-19-05-2008: Welch test for feature selection
#'   NOTE: This function selects features based on the p-values rather than
#'   absolute values of statistics. And also supports additional arguments
#'   passing, such as paired test or not, and the alternative hypothesis.
fs.welch.1 <- function(x, y, ...) {
  tmp <- sapply(as.data.frame(x), function(x) {
    tmp <- t.test(x ~ y, var.equal = F, ...)
    c(tmp$statistic, tmp$p.value)
  })

  stats <- tmp[1, ]
  pval <- tmp[2, ]

  fs.order <- order(pval, decreasing = F, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(pval)
  nam <- names(pval[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(
    fs.order = fs.order, fs.rank = fs.rank, pval = pval,
    stats = stats
  )
  return(res)
}

#' =========================================================================
#' Wilcoxon test for feature selection
#' lwc-21-06-2010: commence
#' Note: No doc and not export
fs.wilcox <- function(x, y, ...) {
  tmp <- sapply(as.data.frame(x), function(x) {
    tmp <- wilcox.test(x ~ y, ...)
    c(tmp$statistic, tmp$p.value)
  })

  stats <- tmp[1, ]
  pval <- tmp[2, ]

  fs.order <- order(abs(stats), decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(stats)
  nam <- names(stats[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(fs.order = fs.order, fs.rank = fs.rank, stats = abs(stats),
              pval = pval)
  return(res)
}

#' =========================================================================
#' ANOVA for feature selection
#' lwc-05-01-2007: commence
#' lwc-16-03-2007: minor change
fs.anova <- function(x, y, ...) {
  tmp <- sapply(as.data.frame(x), function(x) {
    tmp <- oneway.test(x ~ y, var.equal = T)
    c(tmp$statistic, tmp$p.value)
  })

  stats <- tmp[1, ]
  pval <- tmp[2, ]

  fs.order <- order(abs(stats), decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(stats)
  nam <- names(stats[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(
    fs.order = fs.order, fs.rank = fs.rank, stats = abs(stats),
    pval = pval
  )
  return(res)
}

#' =========================================================================
#' Kruskal-Wallis test for feature selection
#' lwc-05-01-2007: commence
#' lwc-16-03-2007: minor change
fs.kruskal <- function(x, y, ...) {
  tmp <- sapply(as.data.frame(x), function(x) {
    tmp <- kruskal.test(x ~ y)
    c(tmp$statistic, tmp$p.value)
  })

  stats <- tmp[1, ]
  pval <- tmp[2, ]

  fs.order <- order(abs(stats), decreasing = T, na.last = T)
  fs.rank <- order(fs.order)

  names(fs.rank) <- names(stats)
  nam <- names(stats[fs.order])
  if (!is.null(nam)) fs.order <- nam

  res <- list(
    fs.order = fs.order, fs.rank = fs.rank, stats = abs(stats),
    pval = pval
  )
  return(res)
}

#' =========================================================================
#' Feature ranking validation by error estimation.
#' History:
#'   08-11-2006: commence
#'   09-11-2006: output different errors.
#'   26-11-2006: Borda count for selecting final feature order
#'   10-01-2007: Add user-defined data partitioning

#'   12-10-2007: Add fs.order as a argument. This allows the user to
#'               estimate error using a feature order calculated somewhere
#'               else. Actually this function replaces rankvali (deprecated).
frankvali.default <- function(dat, cl, cl.method = "svm",
                              fs.method = "fs.auc", fs.order = NULL,
                              fs.len = "power2", pars = valipars(),
                              tr.idx = NULL, ...) {
  #' validity checking
  if (missing(dat) || missing(cl)) {
    stop("data set or class are missing")
  }
  if (length(dim(dat)) != 2) {
    stop("'dat' must be a matrix or data frame")
  }
  if (!is.factor(cl)) {
    stop("cl must be a factor.")
  }
  if (nrow(dat) != length(cl)) stop("dat and cl don't match.")
  if (length(unique(cl)) < 2) {
    stop("Classification needs at least two classes.")
  }
  if (any(is.na(dat)) || any(is.na(cl))) {
    stop("NA is not permitted in data set or class labels.")
  }

  dat <- as.matrix(dat)
  rownames(dat) <- NULL # strip off the row names

  n <- nrow(dat) #' number of samples
  p <- ncol(dat) #' size of feature

  len <- get.fs.len(p, fs.len = fs.len)
  nlen <- length(len)

  #' construct index of train data
  if (is.null(tr.idx)) {
    if (pars$sampling == "cv" && pars$nreps > n) {
      pars$sampling <- "loocv"
    }
    tr.idx <- trainind(cl, pars = pars)
  }
  pars$niter <- length(tr.idx)
  pars$nreps <- length(tr.idx[[1]])

  err.all <- list()
  fs.list <- list()

  for (i in 1:pars$niter) {
    train.ind <- tr.idx[[i]]
    res <- list()
    #' generic loop for loocv, cv, scv,random and bootstrap.
    for (j in 1:length(train.ind)) {
      dat.tr <- dat[train.ind[[j]], , drop = F]
      cl.tr <- cl[train.ind[[j]]]
      dat.te <- dat[-train.ind[[j]], , drop = F]
      cl.te <- cl[-train.ind[[j]]]
      #' Error estimation of feature selection with fs.order or with rank
      #' method.
      res[[j]] <- frank.err(dat.tr, cl.tr, dat.te, cl.te,
        cl.method = cl.method,
        fs.method = fs.method, fs.order = fs.order,
        fs.len = fs.len, ...
      )
    } #' end of j
    #' feature ranking list 
    fs.list[[i]] <- sapply(res, function(x) cbind(x$fs.rank))
    rownames(fs.list[[i]]) <- colnames(dat)

    #' error estimation
    err.all[[i]] <- t(sapply(res, function(x) cbind(x$error)))
    colnames(err.all[[i]]) <- len
    #' colnames(err.all[[i]]) <- paste("Len_", len, sep="")
  } #' End of i

  names(err.all) <- paste("Iter_", seq(1, pars$niter), sep = "")
  names(fs.list) <- paste("Iter_", seq(1, pars$niter), sep = "")

  err.iter <- t(sapply(err.all, function(x) apply(x, 2, mean)))
  err.avg <- apply(err.iter, 2, mean)
  #' or try
  #' err.mat <- do.call("rbind",err.all)
  #  err.avg <- apply(err.mat,2,mean)

  if (is.null(fs.order)) {
    #' final feature ranking
    #' Use Borda count to get the final feature order
    fs.mat <- do.call("cbind", fs.list)
    fs.score <- apply(fs.mat, 1, sum)
    fs.order <- order(fs.score, decreasing = F) #' fs order from best to worst.
    fs.rank <- order(fs.order, decreasing = F) #' fs rank score.
    names(fs.rank) <- rownames(fs.mat)
    temp <- names(fs.rank[fs.order])
    if (!is.null(temp)) {
      fs.order <- noquote(temp)
    }
  } else {
    fs.rank <- order2rank(fs.order)
    fs.method <- "User defined"
  }

  ret <- list(
    fs.method = fs.method,
    cl.method = cl.method,
    fs.len = len, #' computational levels
    err.avg = err.avg, #' average error
    err.iter = err.iter, #' error matrix on each iteration
    err.all = err.all, #' all error matrix
    fs.order = fs.order, #' final feature order
    fs.rank = fs.rank, #' final feature rank
    sampling = switch(pars$sampling,
      "cv"     = "cross validation",
      "loocv"  = "leave-one-out cross-validation",
      "boot"   = "bootstrap",
      "random" = "randomised validation (holdout)"
    ),
    niter = pars$niter, #' number of iteration
    nreps = pars$nreps
  )

  if (is.null(fs.order)) {
    ret$fs.list <- fs.list #' feature list of all computation
  }
  class(ret) <- "frankvali"
  return(ret)
}

#' ========================================================================
frankvali <- function(dat, ...) UseMethod("frankvali")

#' ========================================================================
#' lwc-12-11-2006:
frankvali.formula <- function(formula, data = NULL, ..., subset,
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

  ret <- frankvali.default(x, y, ..., na.action = na.action)

  ret$call <- call
  ret$call[[1]] <- as.name("frankvali")
  ret$terms <- Terms
  if (!is.null(attr(m, "na.action"))) {
    ret$na.action <- attr(m, "na.action")
  }
  class(ret) <- c("frankvali.formula", class(ret))
  return(ret)
}

#' =======================================================================
#' lwc-13-11-2006: boxplot the error rate on each iteration or computation.
boxplot.frankvali <- function(x, ...) {
  col <- "lightgray"
  xlab <- "Numbers of Feature"
  ylab <- "Error Rate"
  ylim <- c(0, 1.0)

  if (x$niter > 1) {
    main <- "Error rate on each iteration"
    #' tmp <- data.frame(x$err.iter)     #' why does data.frame change names?
    tmp <- as.data.frame(x$err.iter) #' wll-1706-2008: Is it OK?
    colnames(tmp) <- colnames(x$err.iter)
  } else {
    main <- "Error rate on each computation"
    tmp <- as.data.frame(x$err.all) #' tmp <- data.frame(x$err.all)
    colnames(tmp) <- colnames(x$err.all)
  }

  boxplot(tmp, main = main, col = col, xlab = xlab, ylab = ylab, ylim = ylim)
}

#' ========================================================================
print.frankvali <- function(x, digits = 3, ...) {
  cat("\nFeature selection method:\t\t", x$fs.method)
  cat("\nClassification method:\t\t", x$cl.method)
  cat("\nSampling:\t\t", x$sampling)
  cat("\nNo. of iteration.:\t", x$niter)
  cat("\nNo. of replications:\t", x$nreps)

  cat("\nFeature numbers:\t", x$fs.len)
  cat("\nAverage error:\t\t", round(x$err.avg, digits))

  cat(
    "\nFeature order (top 10):\t",
    x$fs.order[1:min(10, length(x$fs.order))]
  ) #' best to worst
  cat("\n")
  invisible(x)
}

#' ========================================================================
summary.frankvali <- function(object, ...) {
  structure(object, class = "summary.frankvali")
}

#' ========================================================================
print.summary.frankvali <- function(x, ...) {
  print.frankvali(x)
  cat("\nError of each iteration:\n")
  print(x$err.iter)
  #' cat("\nError of all computation:\n")
  #' print(x$err.all)
  if (!is.null(x$fs.list)) {
    cat("\nFeature ranking List:\n")
    print(x$fs.list)
  }
}

#' =========================================================================
#' Error estimation of feature ranking on a single data set.
#' History:
#'   08-11-2006: commence
#'   10-11-2006: generalise cl.method.
#'   10-01-2007: minor modification
#'   24-09-2007: Generalise fs.method
#'   12-10-2007: Add fs.order as a argument.
frank.err <- function(dat.tr, cl.tr, dat.te, cl.te, cl.method = "svm",
                      fs.method = "fs.auc", fs.order = NULL,
                      fs.len = "power2", ...) {
  if (missing(dat.tr) || missing(cl.tr) || missing(dat.te) || missing(cl.te)) {
    stop(" training and test data missing")
  }

  #' feature ranking
  if (is.null(fs.order)) {
    tmp <- do.call(fs.method, c(list(x = dat.tr, y = cl.tr), list(...)))
    fs.order <- tmp$fs.order
    fs.rank <- tmp$fs.rank
  } else {
    fs.rank <- order2rank(fs.order)
    fs.method <- "User defined"
  }

  #' generate feature length for error estimation
  p <- ncol(dat.tr)
  len <- get.fs.len(p, fs.len = fs.len)
  nlen <- length(len)

  #' error estimation
  error <- numeric(length = nlen)
  names(error) <- len
  #' names(error) <- paste("Len_", len, sep="")
  for (i in 1:nlen) {
    #' feature selection
    sel <- fs.order[1:len[i]]
    error[i] <- classifier(dat.tr[, sel, drop = F], cl.tr,
      dat.te[, sel, drop = F], cl.te,
      method = cl.method, ...
    )$err
  }

  ret <- list(
    cl.method = cl.method, fs.len = len, error = error,
    fs.method = fs.method, fs.order = fs.order, fs.rank = fs.rank
  )

  return(ret)
}

#' =========================================================================
#' wll-29-04-2008: Wrapper function for validation of feature selection by
#'                 classification.
#' wll-29-10-2008: Give a logical value to validate all fs or not.
#' lwc-07-10-2011: use get.fs.len again but remove the last one.
#' Note: 1. Similar but more complete function is frankvali
#'       2. should change cl.method as method
fs.cl <- function(dat, cl, fs.order = colnames(dat), fs.len = 1:ncol(dat),
                  cl.method = "svm", pars = valipars(), all.fs = FALSE, ...) {
  len <- get.fs.len(ncol(dat), fs.len = fs.len)
  if (!all.fs) {
    len <- len[1:(length(len) - 1)] #' remove the last one
  }
  nlen <- length(len)

  res <- sapply(1:nlen, function(i) {
    id <- fs.order[1:len[i]] #' extract index of selected features
    #' cat("\n--Feature Length = :",i,"\n"); flush.console()
    res <- aam.cl(dat[, id, drop = F], cl, method = cl.method, pars = pars, ...)
  })
  res <- t(res)
  rownames(res) <- len
  res
}

#' =========================================================================
#' lwc-27-06-2011: Wrapper function for validation of feature selection by
#'                classification.
#' Note: This function evaluate all features given by user either in
#' individual feature or aggregated features.
fs.cl.1 <- function(dat, cl, fs.order = colnames(dat), cl.method = "svm",
                    pars = valipars(), agg_f = FALSE, ...) {
  len <- length(fs.order)
  if (agg_f) {
    res <- sapply(1:len, function(i) {
      id <- fs.order[1:i] #' aggregation of features
      res <- aam.cl(dat[, id, drop = F], cl, method = cl.method, pars = pars, ...)
    })
  } else {
    res <- sapply(1:len, function(i) {
      id <- fs.order[i] #' individual feature
      res <- aam.cl(dat[, id, drop = F], cl, method = cl.method, pars = pars, ...)
    })
  }
  res <- t(res)
  rownames(res) <- 1:len
  res
}

#' =========================================================================
#' lwc-27-06-2011: Wrapper function for validation of feature selection by
#'                 classification.
#' lwc-19-05-2012: replace aam.cl with accest in order to get more results.
#' lwc-22-10-2012: To get aam, call perf.aam.
#' lwc-21-01-2014: To get other outcome such as SE and CI, need to provide
#' extra code scripts. Refer to frankvali. Usages:
#' usages
#'  data(abr1)
#'  dat <- abr1$pos
#'  x <- preproc(dat[, 110:500], method = "log10")
#'  y <- factor(abr1$fact$class)
#'  dat <- dat.sel(x, y, choices = c("1", "2"))
#'  x.1 <- dat[[1]]$dat
#'  y.1 <- dat[[1]]$cls
#'  pars <- valipars(sampling = "cv", niter = 4, nreps = 4)
#'
#'  #' multi-classes
#'  fs <- fs.rf(x, y)
#'  ord <- fs$fs.order[1:50]
#'  res <- fs.cl.2(x, y,
#'    fs.order = ord, cl.method = "svm", pars = pars,
#'    agg_f = TRUE
#'  )
#'  perf.aam(res)
#'
#'  #' two-classes
#'  fs <- fs.rf(x.1, y.1)
#'  ord <- fs$fs.order[1:50]
#'  res.1 <- fs.cl.2(x.1, y.1,
#'    fs.order = ord, cl.method = "svm", pars = pars,
#'    agg_f = TRUE
#'  )
#'  perf.aam(res.1)
#'
fs.cl.2 <- function(dat, cl, fs.order = colnames(dat), cl.method = "svm",
                    pars = valipars(), agg_f = FALSE, ...) {
  len <- length(fs.order)
  if (agg_f) {
    res <- lapply(1:len, function(i) {
      id <- fs.order[1:i] #' aggregation of features
      res <- accest(dat[, id, drop = F], cl,
        method = cl.method,
        pars = pars, ...
      )
      #' res <- aam.cl(dat[,id, drop=F],cl, method=cl.method, pars=pars,...)
    })
  } else {
    res <- lapply(1:len, function(i) {
      id <- fs.order[i] #' individual feature
      res <- accest(dat[, id, drop = F], cl,
        method = cl.method,
        pars = pars, ...
      )
      #' res <- aam.cl(dat[,id, drop=F],cl, method=cl.method, pars=pars,...)
    })
  }
  #' res <- t(res)
  #' rownames(res) <- 1:len
  names(res) <- 1:len
  res
}

#' =======================================================================
#' lwc-21-01-2014: get average of acc, auc and mar from outcome of fs.cl.2
perf.aam <- function(res) {
  tmp <- sapply(res, function(x) { #'    x = res[[1]]
    acc <- x$acc
    auc <- ifelse(!is.null(x$auc), x$auc, NA)
    mar <- ifelse(!is.null(x$mar), x$mar, NA)
    res <- c(acc = acc, auc = auc, mar = mar)
  })
  return(t(tmp))
}

#' =======================================================================
#' lwc-24-05-2012: Wrapper function for perf.aam.
perf <- function(res) {
  perf <- lapply(res, function(x) perf.aam(x))
}

#' ======================================================================
#' Generate a sequence of feature number
#' History:
#'  25-10-2006: commence
#'  31-10-2006: add user defined sequence
#'  15-11-2006: fix a bug in returning a decreasing vector
#' Usages:
#' get.fs.len(10,fs.len=c(1,5,3,11.2,7.8,23,1,0))
#' get.fs.len(200,fs.len="half")
get.fs.len <- function(p, fs.len = c("power2")) {
  if (!is.character(fs.len)) {
    fs.len <- as.vector(fs.len)
    fs.len <- as.integer(fs.len)
    fs.len <- fs.len[fs.len <= p & fs.len > 0]
    fs.len <- c(p, fs.len)
    x <- unique(fs.len)
  } else {
    fs.len <- match.arg(fs.len, c("full", "half", "power2"))

    if (fs.len == "full") {
      x <- seq(p, 1)
    } else if (fs.len == "half") {
      x <- tmp <- p
      while (tmp > 1) {
        tmp <- trunc(tmp / 2)
        x <- c(x, tmp)
      }
    } else {
      n <- ceiling(log2(p))
      x <- 2^(n:0)
      x[1] <- p
    }
  }

  #' x <- sort(x,decreasing = T) #' must be decreasing order for SVM-RFE
  x <- sort(x, decreasing = F)
  return(x)
}

#' ======================================================================
#' lwc-02-02-2007: convert feature rank to feature order
#' NOTE: The vector of fs.rank should have variable names.
#' Usages:
#'  load("D:\R_lwc\data\R-W2-GC\31_01_2007\class_rfrankvali_auto.RData")
#'  fs.rank.list <- do.call("cbind",fs$rfe$fs.list)
#'  tmp     <- mt:::rank2order(fs.rank.list[,1])
#'  fs.order.list <- sapply(1:ncol(fs.rank.list),
#'                          function(x) mt:::rank2order(fs.rank.list[,x]))
#' Internal function.
rank2order <- function(fs.rank) {
  fs.order <- order(fs.rank)
  tmp <- names(fs.rank[fs.order])
  if (!is.null(tmp)) fs.order <- tmp

  return(fs.order)
}

#' =======================================================================
#' wll-12-03-2007: convert feature order to feature rank
#' Internal function.
order2rank <- function(fs.order) {
  fs.rank <- order(fs.order)
  names(fs.rank) <- fs.order[fs.rank]
  return(fs.rank)
}

#'  1) feat.freq
#'  2) feat.cons
#'  3) feat.mfs
#'  4) feat.mfs.stab
#'  5) feat.mfs.stats
#'  6) feat.agg
#'  7) feat.rank.re
#'  8) fs.plsvip
#'  9) fs.plsvip.1
#' 10) fs.plsvip.2
#' 11) fs.pls
#' 12) fs.pca
#' 13) fs.bw
#' 14) fs.relief
#' 15) fs.rfe
#' 16) fs.snr
#' 17) fs.auc
#' 18) fs.rf
#' 19) fs.rf.1
#' 20) fs.welch
#' 21) fs.welch.1
#' 22) fs.wilcox
#' 23) fs.anova
#' 24) fs.kruskal
#' 25) frankvali.default
#' 26) frankvali
#' 27) frankvali.formula
#' 28) boxplot.frankvali
#' 29) print.frankvali
#' 30) summary.frankvali
#' 31) print.summary.frankvali
#' 32) frank.err
#' 33) fs.cl
#' 34) fs.cl.1
#' 35) fs.cl.2
#' 36) perf.aam
#' 37) perf
#' 38) get.fs.len
#' 39) rank2order
#' 40) order2rank
