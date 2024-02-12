#' ========================================================================
#' Calculate classification accuracy.
#' History:
#'   31-07-06: Commence
#'   15-09-06: Not restrict types of classification. Just check the validity
#'             of models and the output of predict method.
#'   10-01-07: some minor changes.
#'   30-03-07: revise
#'   30-06-07: Add posterior, margin and AUC
#'   03-07-07: For loocv, no AUC are available. (Should check in code)
#'   14-12-07: method and pred.dunc can be a function or a character string
#'             naming function.
#'   11-05-12: Re-write calculation of err, auc and mar. Keep res.all.
accest.default <- function(dat, cl, method,
                           pred.func = predict, pars = valipars(),
                           tr.idx = NULL, ...) {
  #' validity checking
  if (missing(dat) || missing(cl)) {
    stop("data set or class are missing")
  }
  if (missing(method)) { #' 01-06-2008:  you can use: if(is.na(method))
    stop("'method' is missing")
  }

  #' if(length(dim(dat)) != 2)
  #'   stop("'mat' must be a matrix or data frame")

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
  n <- nrow(dat)
  rownames(dat) <- NULL # strip off the row names
  #' lwc-15-09-2006: avoid the warnings messages when calling SVM with
  #' bootstrap methods. The warnings: 'row.names duplicated' come from:
  #'        df <- na.action(data.frame(y, x))
  #' in svm.default's '#' subsetting and na-handling for matrices' segment.

  #' wll-14-12-2007: either a function or a character string naming the
  #' function to be pred.func.
  if (is.function(method)) {
    method <- deparse(substitute(method))
  }
  pred.func <-
    if (is.function(pred.func)) {
      pred.func
    } else if (is.character(pred.func)) {
      get(pred.func)
    } else {
      eval(pred.func)
    }

  #' construct index of train data
  if (is.null(tr.idx)) {
    if (pars$sampling == "cv" && pars$nreps > n) {
      pars$sampling <- "loocv"
    }
    tr.idx <- trainind(cl, pars = pars)
  }
  pars$niter <- length(tr.idx)
  pars$nreps <- length(tr.idx[[1]])

  #' Main computation
  res.all <- list()
  cat("\nACCEST Iteration (", method, ", ", pars$niter, "):", sep = "")
  for (i in 1:pars$niter) {
    cat(" ", i, sep = "")
    flush.console() #' for Windows
    train.ind <- tr.idx[[i]]
    res <- list()
    for (j in 1:length(train.ind)) {
      dat.tr <- dat[train.ind[[j]], , drop = F]
      cl.tr <- cl[train.ind[[j]]]
      dat.te <- dat[-train.ind[[j]], , drop = F]
      cl.te <- cl[-train.ind[[j]]]
      res[[j]] <- classifier(dat.tr, cl.tr, dat.te, cl.te,
        method = method, pred.func = pred.func
      )
    }
    names(res) <- paste("reps_", seq(1, pars$nreps), sep = "")
    res.all[[i]] <- res
  }
  cat("\n")
  names(res.all) <- paste("Iter_", seq(1, pars$niter), sep = "")

  #' calculate error rate
  err.all <- lapply(res.all, function(x) {
    func <- function(x) sapply(x, function(y) y$err)
    func(x)
  })
  err.all <- do.call(rbind, err.all)
  err.iter <- rowMeans(err.all) #' apply(err.all,1,mean)
  err <- mean(err.iter)

  #' calculate mean of margin
  if (!is.null(res.all[[1]][[1]]$margin)) {
    mar.all <- lapply(res.all, function(x) {
      func <- function(x) sapply(x, function(y) mean(y$margin))
      func(x)
    })
    mar.all <- do.call(rbind, mar.all)
    mar.iter <- rowMeans(mar.all)
    mar <- mean(mar.iter)
  } else {
    mar.all <- NULL
    mar.iter <- NULL
    mar <- NULL
  }

  #' calculate AUC
  if (!is.null(res.all[[1]][[1]]$auc) && pars$sampling != "loocv") {
    auc.all <- lapply(res.all, function(x) {
      func <- function(x) sapply(x, function(y) y$auc)
      func(x)
    })
    auc.all <- do.call(rbind, auc.all)
    auc.iter <- rowMeans(auc.all)
    auc <- mean(auc.iter)
  } else {
    auc.all <- NULL
    auc.iter <- NULL
    auc <- NULL
  }

  #' calculate accuracy rate
  acc.all <- 1 - err.all
  acc.iter <- 1 - err.iter
  acc <- 1 - err
  #' if (pars$niter > 1) acc.std  <- sd(acc.iter) else  acc.std <- NULL

  #' process bootstrap acc
  if (pars$sampling == "boot") {
    resub <- classifier(dat, cl,
      dat.te = NULL, cl.te = NULL,
      method = method, ...
    )
    err.boot <- lapply(err.iter, function(x) boot.err(x, resub))
    #' reshape
    err.boot <- t(sapply(err.boot, function(x) do.call("c", x)))
    acc.boot <- 1 - err.boot
  }

  #' overall confusion matrix
  #' lwc-01-06-2007: Do not use sapply here because sometimes get non-equal
  #' fold.
  all.cl <- lapply(res.all, function(x) {
    foo <- function(x) lapply(x, function(y) y$cl)
    foo(x)
  })
  all.cl <- unlist(unlist(all.cl, recursive = F, use.names = F))

  all.pred <- lapply(res.all, function(x) {
    foo <- function(x) lapply(x, function(y) y$pred)
    foo(x)
  })
  all.pred <- unlist(unlist(all.pred, recursive = F, use.names = F))

  #' overall confusion matrix
  conf <- table(all.cl, all.pred)

  #' construct output
  ret <- list(
    method = method,
    acc = acc,
    acc.iter = acc.iter,
    acc.all = acc.all,
    auc = auc,
    auc.iter = auc.iter,
    auc.all = auc.all,
    mar = mar,
    mar.iter = mar.iter,
    mar.all = mar.all,
    err = err,
    err.iter = err.iter,
    err.all = err.all,
    sampling = switch(pars$sampling,
      "loocv"  = "leave-one-out cross-validation",
      "cv"     = "cross validation",
      "boot"   = "bootstrap",
      "rand"   = "randomised validation (holdout)"
    ),
    niter = pars$niter,
    nreps = pars$nreps,
    conf = conf,
    res.all = res.all #' lwc-10-05-2012: keep all results.
  )
  if (pars$sampling == "boot") {
    ret$acc.boot <- acc.boot
  }

  class(ret) <- "accest"
  return(ret)
}

#' ========================================================================
#' wll-01-07-2007: add AUC and margin
print.accest <- function(x, digits = 3, ...) {
  #'  cat("\nCall:\n", deparse(x$call), "\n")
  cat("\nMethod:\t\t\t", x$method)
  cat("\nAccuracy:\t\t", round(x$acc, digits))
  if (!is.null(x$auc)) {
    cat("\nAUC:\t\t\t", round(x$auc, digits))
  }
  if (!is.null(x$mar)) {
    cat("\nMargin:\t\t\t", round(x$mar, digits))
  }

  cat("\n\nNo. of iteration:\t", x$niter)
  cat("\nSampling:\t\t", x$sampling)
  cat("\nNo. of replications:\t", x$nreps)

  cat("\n\nOverall confusion matrix of training data:\n")
  print(x$conf)

  invisible(x)
}

#' ========================================================================
summary.accest <- function(object, ...) {
  structure(object, class = "summary.accest")
}

#' ========================================================================
print.summary.accest <- function(x, digits = 3, ...) {
  print.accest(x)
  cat("\nAccuracy on each iteration:\n")
  print(round(x$acc.iter, digits))
  invisible(x)
}

#' =======================================================================
plot.accest <- function(x, main = NULL, xlab = NULL, ylab = NULL, ...) {
  if (x$niter == 1) {
    stop("Number of iteration (niter) must be greater than 1")
  }

  if (is.null(main)) {
    main <- paste("Performance of `", x$method, "'",
      sep = " ", "(",
      x$sampling, ")"
    )
  }

  if (is.null(xlab)) xlab <- "Index of niter"
  if (is.null(ylab)) ylab <- "Accuracies"

  plot(x$acc.iter,
    type = "b", main = main, xlab = xlab, ylab = ylab,
    col = "blue", ...
  )
}

#' ========================================================================
accest <- function(dat, ...) UseMethod("accest")

#' ========================================================================
accest.formula <- function(formula, data = NULL, ..., subset,
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

  ret <- accest.default(x, y, ..., na.action = na.action)

  ret$call <- call
  ret$call[[1]] <- as.name("accest")
  ret$terms <- Terms
  if (!is.null(attr(m, "na.action"))) {
    ret$na.action <- attr(m, "na.action")
  }
  class(ret) <- c("accest.formula", class(ret))
  return(ret)
}

#' =======================================================================
#' Estimates the accuracy of pairwise comparison.
#' History:
#'   10-10-2006: Create
#'   03-12-2006: process fold of cv or scv is greater than number of data
#'   row
#'   10-01-2007: vectorised
binest <- function(dat, cl, choices = NULL, method, pars = valipars(), ...) {

  #' construct the pairwise comparison data set
  dat.sub <- .dat.sel(dat, cl, choices = choices)

  len <- sapply(dat.sub$cl, length)

  if (pars$sampling == "cv" && pars$nreps > min(len)) {
    pars$sampling <- "loocv"
  }
  if (pars$sampling == "loocv") pars$niter <- 1

  n <- nrow(dat.sub$com)
  acc <- sapply(c(1:n), function(x) {
    accest(dat.sub$dat[[x]], dat.sub$cl[[x]],
      method = method,
      pars = pars, ...
    )$acc
  })

  com <- apply(dat.sub$com, 1, paste, collapse = " ~ ")
  names(acc) <- com

  ret <- list(
    com = dat.sub$com, acc = acc, method = method,
    niter = pars$niter, sampling = pars$sampling
  )
  if (pars$sampling != "loocv") ret$nreps <- pars$nreps
  return(ret)
}

#' =======================================================================
#' wll-29-04-2008: Wrapper function for re-sampling based classification
#' with multiple classifiers
aam.mcl <- function(x, y, method, pars = valipars(), ...) {
  res <- lapply(method, function(m) {
    cat("\nClassifier = :", m, "\n")
    flush.console()
    aam.cl(x, y, method = m, pars, ...)
  })
  names(res) <- method
  res <- do.call(rbind, res)
}

#' ========================================================================
#' wll-29-04-2008: Calculate acc, auc and mar by calling accest
aam.cl <- function(x, y, method, pars = valipars(), ...) {
  val <- accest(x, y, method = method, pars = pars, ...)
  acc <- val$acc
  auc <- ifelse(!is.null(val$auc), val$auc, NA)
  mar <- ifelse(!is.null(val$mar), val$mar, NA)
  res <- c(acc = acc, auc = auc, mar = mar)
  return(round(res, digits = 3))
}

#' ========================================================================
#' lwc-05-07-2006: Calculate classification rate.
#' lwc-16-09-2006: Arguments checking.
#' lwc-30-10-2006: Convert arguments as vectors
#' lwc-05-05-2010: Change function name from class.rate
#' lwc-19-05-2-12: Add kappa and change class.error to class.acc.
#' Usage:
#'  a = c(1,2,1,2,1,2,1,2,1,2)
#'  b = c(2,2,2,1,1,2,1,2,1,2)
#'  c = rbind(a,b)
#'  d = c(1,1,1,1,1,1,1,1,1,1)
#'  cl.rate(a,b)
#'  cl.rate(a,d)
#'  cl.rate(a,c)
cl.rate <- function(obs, pre) {
  observed <- as.vector(obs)
  predicted <- as.vector(pre)

  if (length(observed) != length(predicted)) {
    stop(" 'observed' and 'predicted' must have the same length.")
  }

  acc <- sum(observed == predicted) / length(observed)

  #' NOTE: If arguments are factor, it will trigger an error when they
  #'       have different levels.
  err <- (1 - acc)
  con <- table(observed = observed, predicted = predicted)
  kappa <- classAgreement(con)$kappa #' from e1071
  con <- cbind(con, class.acc = diag(con) / rowSums(con))
  res <- list(acc = acc, err = err, con.mat = con, kappa = kappa)
  return(res)
}

#' =======================================================================
#' lwc-28-04-2010: Classification performance.
#' lwc-19-05-2012: add kappa, con.mat and positive in the returned list.
#' Note: Also see cl.auc and cl.roc for ROC assessment.
#'   True positive rate:  tpr = TP/P = TP/TP+FN
#'   False positive rate: fpr = FP/N = FP/FP+TN
#'   Accuracy:            acc = P*tpr + N*(1-fpr)
#' References: Tom Fawcett, ROC Graphs: Notes and Practical Considerations
#'             for Data Mining Researchers, January 7, 2003.
#' Arguments:
#'   obs 	Factor or vector of observed class.
#'   pre 	Factor or vector of predicted class.
cl.perf <- function(obs, pre, pos = levels(as.factor(obs))[2]) {
  obs <- factor(obs)
  pre <- factor(pre)
  pos <- match.arg(pos, levels(obs))

  obs <- relevel(obs, ref = pos) #' put pos in the 1st position
  pre <- relevel(pre, ref = pos) #' put pos in the 1st position

  #' confusion matrix
  con <- table(obs, pre)
  kappa <- classAgreement(con)$kappa
  con <- cbind(con, class.acc = diag(con) / rowSums(con))

  #' change levels as pos(Positive) and neg (Negative)
  levels(obs) <- list("neg" = levels(obs)[2], "pos" = levels(obs)[1])
  levels(pre) <- list("neg" = levels(pre)[2], "pos" = levels(pre)[1])
  #' levels(obs) <- levels(pre) <- c("pos", "neg")

  #' Temp function for number of True or False (positive or negative)
  #' (Not used here)
  tf <- function(x) {
    x.obs <- grepl(x, obs)
    x.pre <- grepl(x, pre)
    Tr <- sum(x.pre & x.obs)
    Fa <- sum(x.pre) - Tr
    list(Tr = Tr, Fa = Fa)
  }

  #' pos.tf <- tf("pos")
  pos.obs <- grepl("pos", obs)
  pos.pre <- grepl("pos", pre)
  TP <- sum(pos.pre & pos.obs)
  Pos <- sum(pos.pre)
  FP <- Pos - TP

  #' neg.tf <- tf("neg")
  neg.obs <- grepl("neg", obs)
  neg.pre <- grepl("neg", pre)
  TN <- sum(neg.pre & neg.obs)
  Neg <- sum(neg.pre)
  FN <- Neg - TN

  P <- TP + FN
  N <- FP + TN

  tpr <- TP / (TP + FN)
  fpr <- FP / (FP + TN)
  acc <- (TP + TN) / (P + N)

  spec <- 1 - fpr #' Specificity
  sens <- tpr #' Sensitivity or recall

  perf <- list(
    acc = acc, tpr = tpr, fpr = fpr, sens = sens, spec = spec,
    con.mat = con, kappa = kappa, positive = pos
  )
  return(perf)
}

#' =====================================================================
#' lwc-11-03-2007: Area under Receiver Operating Curve (ROC).
#' lwc-05-05-2010: Major modification.
#' lwc-01-05-2012: adjust auc
#' Note: 1.) AUC varies between 0.5 and 1.0 for sensible models; the higher
#'           the better.
#'       2.) AUC should be adjusted as:
#'           auc[auc < 0.5] <- 1 - auc[auc < 0.5]
#'       3.) This is the implementation of cutoff-independent AUC.
#'       4.) The AUC has an important statistical property: the AUC of a
#'           classifier is equivalent to the probability that the classifier
#'           will rank a randomly chosen positive instance higher than a
#'           randomly chosen negative instance. This is equivalent to the
#'           Wilcoxon test of ranks.
#'              -- An introduction to ROC analysis by Tom Fawcett
#'       5.) Codes come from auRROC function in package limma.
cl.auc <- function(stat, label, pos = levels(as.factor(label))[2]) {
  if (missing(label) || missing(stat)) stop("arguments miss")
  if (length(label) != length(stat)) stop("lengths differ")

  #' convert label as factor
  label <- factor(label)
  if (nlevels(label) != 2) stop("'label' must be two categorical data")

  #' sort out which level is pos and convert it to "1".
  pos <- match.arg(pos, levels(label))
  label <- relevel(label, ref = pos) #' put pos in the 1st position
  levels(label) <- list("0" = levels(label)[2], "1" = levels(label)[1])

  label <- label[order(stat, decreasing = T)]
  label <- as.numeric(levels(label))[as.integer(label)]
  tmp <- cumsum(label) / sum(label)
  auc <- mean(tmp[label == 0])

  #' if (auc < 0.5) auc <- 1 - auc            #' lwc-08-06-2010: fix it
  auc[auc < 0.5] <- 1 - auc[auc < 0.5] #' lwc-01-05-2012: uncommented.

  return(auc)
}

#' =====================================================================
#' lwc-11-03-2007:Area under Receiver Operating Curve (ROC)
#' Note: 1. Internal function with no documents.
#'       2. It is old version of cl.auc.
#'       3. It is modified from auROC function in package limma.
auc <- function(stat, label) {
  if (missing(label) || missing(stat)) stop("arguments miss")
  if (length(label) != length(stat)) stop("lengths differ")

  if (is.factor(label)) {
    label <- as.numeric(levels(label))[as.integer(label)]
  }
  if (!all(sort(unique(label)) == c(0, 1))) {
    stop("'label' must take values 0 or 1")
  }

  label <- label[order(stat, decreasing = T)]
  tmp <- cumsum(label) / sum(label)
  auc <- mean(tmp[label == 0])
  return(auc)
}

#' =======================================================================
#' lwc-12-03-2007: ROC
#' lwc-05-05-2010: Major modification.
#' lwc-15-05-2012: Major changes:
#'                 1.) The cutoff points are sorted stat
#'                 2.) The rule is >=, not >.
#' Note: 1. This function is partly taken from ROCfuns.R in package ROC.
#' Note: sensitivity and specificity are estimated as:
#'   sensitivity (true positive rate) =
#'               Positives correctly classified / Total positives
#'   specificity (true negative rate) =
#'               Negatives correctly classified / Total negatives
cl.roc <- function(stat, label, pos = levels(as.factor(label))[2],
                   plot = TRUE, ...) {
  roc_rule <- function(thres, x) ifelse(x >= thres, 1, 0)

  if (missing(label) || missing(stat)) stop("arguments miss")
  if (length(label) != length(stat)) stop("lengths differ")

  #' convert label as factor
  label <- factor(label)
  if (nlevels(label) != 2) stop("'label' must be two categorical data")

  #' sort out which level is Positive
  pos <- match.arg(pos, levels(label))
  label <- relevel(label, ref = pos) #' Reorder Levels of Factor
  levels(label) <- list("0" = levels(label)[2], "1" = levels(label)[1])

  #' thresholds
  thres <- unique(sort(stat))

  #' prediction based on thresholds
  pred <- lapply(thres, roc_rule, x = stat)

  perf <- lapply(pred, function(x) {
    sens <- mean(x[label == 1])
    spec <- mean(1 - x[label == 0])
    fpr <- 1 - spec #' false positive rate: 1 - spec
    tpr <- sens #' true positive rate:  sens
    acc <- sum(label == x) / length(label)
    c(acc = acc, tpr = tpr, fpr = fpr, sens = sens, spec = spec)
  })
  perf <- do.call(rbind, perf)
  perf <- cbind(cutoff = thres, perf)
  perf <- data.frame(perf)

  auc <- cl.auc(stat, label)

  if (plot) { #' plot ROC curve
    main <- "ROC Curve"
    xlab <- "False positive Rate"
    ylab <- "True positive Rate"
    #' xlab <- "1 - Specificity"
    #' ylab <- "Sensitivity"

    plot(perf$fpr, perf$tpr,
      type = "n", xlim = c(0, 1), ylim = c(0, 1),
      main = main, xlab = xlab, ylab = ylab, ...
    )
    points(perf$fpr, perf$tpr, col = "red", type = "b", pch = 19)
    abline(
      h = seq(0, 1, by = .1), v = seq(0, 1, by = .1), lty = 3,
      lwd = 0.5, col = "grey"
    )
    abline(0, 1)
  }

  return(list(perf = perf, auc = auc, positive = pos))
}

#' =======================================================================
#' lwc-15-05-2012: Another version of ROC which call cl.perf.
#' lwc-19-05-2012: minor changes.
#' Note: 1.) The cutoff points are sorted stat
#'       2.) The rule is >=, not >.
#'       3.) call cl.perf to get acc, tpr, fpr, sens and spec.
cl.roc.1 <- function(stat, label, pos = levels(as.factor(label))[2],
                     plot = TRUE, ...) {
  if (missing(label) || missing(stat)) stop("arguments miss")
  if (length(label) != length(stat)) stop("lengths differ")

  label <- factor(label)
  if (nlevels(label) != 2) stop("'label' must be two categorical data")

  #' sort out which level is Positive
  pos <- match.arg(pos, levels(label))
  label <- relevel(label, ref = pos) #' The 1st level is positive

  #' thresholds
  thres <- unique(sort(stat))

  #' prediction based on thresholds
  pred <- lapply(thres, function(x) {
    tmp <- ifelse(stat >= x, levels(label)[1], levels(label)[2])
    tmp <- factor(tmp, levels = levels(label))
  })
  #' tidy up perf
  perf <- lapply(pred, function(x) {
    cl.perf(label, x, pos = pos)[c("acc", "tpr", "fpr", "sens", "spec")]
  })
  perf <- do.call(rbind, perf)
  perf <- cbind(cutoff = thres, perf)
  perf <- data.frame(perf)

  #' cutoff - independent AUC
  auc <- cl.auc(stat, label)

  if (plot) { #' plot ROC curve
    main <- "ROC Curve"
    xlab <- "False positive Rate"
    ylab <- "True positive Rate"
    #' xlab <- "1 - Specificity"
    #' ylab <- "Sensitivity"

    plot(perf$fpr, perf$tpr,
      type = "n", xlim = c(0, 1), ylim = c(0, 1),
      main = main, xlab = xlab, ylab = ylab, ...
    )
    points(perf$fpr, perf$tpr, col = "red", type = "b", pch = 19)
    abline(
      h = seq(0, 1, by = .1), v = seq(0, 1, by = .1), lty = 3,
      lwd = 0.5, col = "grey"
    )
    abline(0, 1)
  }

  return(list(perf = perf, auc = auc, positive = pos))
}

#' =======================================================================
#' lwc-01-12-06: Wrapper function for single classifier
#' History:
#'  15-09-06: check predict's output
#'  24-03-07: add dat.te and cl.te as NULLs.
#'  30-06-07: Add the posterior, margin and AUC. I doubt using average of
#'            margin as an assessment factor like AUC for classification.
#'  02-07-07: Deal with constant variables which have zero SD. Any
#'            possible bugs or mistakes in both programming and algorithm?
#'  04-07-07: check validity of auc
#'  06-07-07: Fix a bug
#'  03-10-07: Restore user defined predict function.
#'  14-12-07: method and pred.dunc can be a function or a character string
#'            naming function.
#'  09-10-08: Spot a potential bug in dealing with collinearity
#'  19-05-12: remove names of margin.
#'  22-05-12: minor changes for dots processing for svm
classifier <- function(dat.tr, cl.tr, dat.te = NULL, cl.te = NULL, method,
                       pred.func = predict, ...) {

  #' lwc-14-12-2007: either a function or a character string naming the
  #' function to be pred.func.
  if (is.function(method)) {
    method <- deparse(substitute(method))
  }
  pred.func <-
    if (is.function(pred.func)) {
      pred.func
    } else if (is.character(pred.func)) {
      get(pred.func)
    } else {
      eval(pred.func)
    }

  #' 05-10-2007: re-write this wrapper to avoid multiple arguments
  if (method == "knn") {
    method <- c("knn_wrap")
    knn_wrap <- function(train, cl, ...) list(train = train, cl = cl, ...)
    pred.func <- function(object, newdata, k = 1, ...) {
      #' knn(train=object$train, test=newdata,cl=object$cl, k=k,...)
      knn(train = object$train, test = newdata, cl = object$cl, k = k)
    }
  }

  if (is.null(dat.te) || is.null(cl.te)) {
    dat.te <- dat.tr
    cl.te <- cl.tr
  }

  #' Especially for SVM
  idx <- which(apply(dat.tr, 2, sd) > .Machine$double.eps)
  dat.tr <- dat.tr[, idx, drop = F]
  dat.te <- dat.te[, idx, drop = F]
  #' 09-10-08: Potential bug here. If idx is invalid, dat.tr will be NA. See
  #'           preproc.sd for fixing it. It happens in the situation where
  #'           the data set is a single column matrix and all values are
  #'           same.

  nte <- length(cl.te) #' number of rows in test data
  dat.all <- rbind(dat.te, dat.tr)
  #' cl.all  <- factor(c(cl.te,cl.tr))
  cl.all <- factor(c(as.character(cl.te), as.character(cl.tr)))
  #' Note-04-07-2007: This is proper way to merge factors.

  #' lwc-22-05-2012: minor changes
  if (method != "svm") {
    model <- do.call(method, c(list(dat.tr, cl.tr), list(...)))
  } else {
    #' pre-process only for SVM
    dots <- list(...)
    #' wl-08-11-2021, Mon: Use this one 
    if (hasArg("probability")) dots$probability <- NULL
    # if (hasArg(probability)) dots$probability <- NULL
    model <- do.call(method, c(list(dat.tr, cl.tr), probability = T, dots))
    #' Note-30-06-2007: Using probability = T only for SVM.
  }

  #' pred  <- pred.func(model, dat.all,...)     # predict the entire data set
  pred <- pred.func(model, data.frame(dat.all), ...) #' lwc-22-05-2012:

  if (!is.list(pred)) {
    pred.te <- pred[1:nte]
    if (method == "svm") {
      prob <- attr(predict(model, dat.all, probability = TRUE), "probabilities")
      prob.te <- prob[1:nte, , drop = F] #' test prob
    } else if (method == "randomForest") {
      prob <- predict(model, dat.all, type = "prob")
      prob.te <- prob[1:nte, , drop = F] #' test prob
    } else {
      prob.te <- NULL
    }
  } else {
    if (!is.null(pred$class)) { # for 'lda' and 'qda'
      pred.te <- pred$class[1:nte]
      prob.te <- pred$posterior[1:nte, , drop = F]
    } else {
      stop("predict does not return a list with component 'class'.")
    }
  }

  #' calculate error rate
  err <- sum(cl.te != pred.te) / length(cl.te)
  #' err <- cl.rate(cl.te, pred.te)$err.rate

  #' Sort the prob by the column names. (for AUC calculation purpose)
  #' if (!is.null(prob.te)) {
  #'   dfnames <- colnames(prob.te)
  #'   dfnames <- sort(dfnames)
  #'   prob.te <- prob.te[,dfnames,drop=F]
  #' }
  #' Note-06-07-07: Above code lines have bugs, e.g. if colnames are like
  #' 3,2,10. After sorting, colnames are 10,2,3 not 2,3,10.
  prob.te <- prob.te[, levels(cl.te), drop = F] #' lwc-06-07-07: for AUC purpose

  #' calculate margin
  if (!is.null(prob.te)) {
    margin <- .marg(prob.te, cl.te)
    names(margin) <- NULL #' lwc-19-05-2012:
    if (length(levels(cl.te)) == 2 && length(cl.te) > 1) {
      #' calculate AUC if two-class problem
      cl.tmp <- cl.te
      levels(cl.tmp) <- c(0, 1)
      stat <- prob.te[, 2]
      auc <- cl.auc(stat, cl.tmp)
      #' auc <- auc(stat,cl.tmp)
    } else {
      auc <- NULL
    }
  } else {
    margin <- NULL
    auc <- NULL
  }

  ret <- list(
    err = err, cl = cl.te, pred = pred.te, posterior = prob.te,
    acc = 1 - err, margin = margin, auc = auc
  )
  return(ret)
}

#' ========================================================================
#' lwc-30-06-2007: calculate the margin of a classifier based on the
#' posterior
#' NOTE: 1. This function hacked from package 'randomForest'. For more
#'       description, see package 'randomForest'.
#'       2. Internal function
.marg <- function(prob, observed) {
  if (missing(observed) || missing(prob)) stop("arguments miss")
  if (length(observed) != nrow(prob)) stop("lengths differ")

  if (any(prob > 1)) {
    prob <- sweep(prob, 1, rowSums(prob), "/")
  }
  observed <- as.factor(observed)

  mat <- data.frame(prob, observed)
  names(mat) <- c(dimnames(prob)[[2]], "observed")
  #' NOTE-lwc: Ater data.frame() operation, the colnames may be changed (if
  #' the colnames are numbers). The above line is to restore the original
  #' colnames.

  nlev <- length(levels(observed))

  ans <- apply(mat, 1, function(x) {
    pos <- match(x[nlev + 1], names(x))
    t1 <- as.numeric(x[pos])
    t2 <- max(as.numeric(x[-c(pos, nlev + 1)]))
    t1 - t2
  })

  names(ans) <- observed
  ans
}

#' =========================================================================
#' lwc-16-08-2006: Calculate bootstrap, .632 bootstrap and .632 plus
#' bootstrap error rate
#' History:
#'   16-08-2006: truncate r
#'   23-03-2007: major changes in principles
boot.err <- function(err, resub) {

  #' apparent error rate/re-substitution rate ( not training error rate)
  err.ae <- resub$err
  cl <- resub$cl
  pred <- resub$pred

  #' .632 bootstrap error
  err.b632 <- 0.368 * err.ae + 0.632 * err

  gamma <-
    sum(outer(cl, pred, function(x, y) ifelse(x == y, 0, 1))) / (length(cl)^2)
  r <- (err - err.ae) / (gamma - err.ae)
  r <- ifelse(err > err.ae & gamma > err.ae, r, 0)
  #' lwc-16-08-2006: if r still falls outside of [0,1], truncate it to 0.
  #' if (r > 1 || r < 0) r <- 0

  errprime <- min(err, gamma)
  #' weight <- .632/(1-.368*r)
  #' err.b632p <- (1-weight)*err.ae + weight*err

  err.b632p <-
    err.b632 + (errprime - err.ae) * (0.368 * 0.632 * r) / (1 - 0.368 * r)

  ret <- list(ae = err.ae, boot = err, b632 = err.b632, b632p = err.b632p)

  return(ret)
}

#' =========================================================================
#' Control parameters for validation using in estimation of classification
#' and feature selection.
#' lwc-03-08-2006: commence
#' lwc-21-03-2007: revise
valipars <- function(sampling = "cv", niter = 10, nreps = 10,
                     strat = FALSE, div = 2 / 3) {
  sampling <- match.arg(sampling, c("loocv", "cv", "boot", "rand"))

  if (sampling == "cv" && nreps < 2) {
    stop("Number of fold (nreps) for cv must be greater than 1.")
  }

  if (sampling == "loocv") {
    res <- list(sampling = sampling, niter = 1)
  } else if (sampling == "rand") {
    res <- list(
      sampling = sampling, niter = niter, nreps = nreps,
      strat = strat, div = div
    )
  } else {
    res <- list(
      sampling = sampling, niter = niter, nreps = nreps,
      strat = strat
    )
  }

  class(res) <- "valipars"
  return(res)
}

#' =======================================================================
#' Generate index for training or feature ranking
#'   lwc-27-10-2006: commence
#'   lwc-29-10-2006: empty class checking
#'   lwc-04-12-2006: support iteration information
#'   lwc-13-02-2007: Fix a bug (set.seed(71))
#'   lwc-22-03-2007: revise
trainind <- function(cl, pars = valipars()) {
  if (!inherits(pars, "valipars")) {
    stop("pars not of class valipars")
  }

  cl <- factor(cl) #' lwc-17-09-2007: drop the factor levels

  idx.func <- function(cl, pars = valipars()) {
    n <- length(cl)

    if (pars$sampling == "loocv") {
      train.ind <- lapply(1:n, function(i) seq(1, n)[-i])
    } else {
      #' 1) each class must have at least 2 samples in training index
      #' 2) each class must have at least 1 sample in test index
      emp_f <- T
      while (emp_f) {
        emp_f <- !emp_f
        switch(pars$sampling,
          "cv" = {
            train.ind <- cv.idx(cl, pars$nreps, strat = pars$strat)
          },
          "rand" = {
            train.ind <- rand.idx(cl, pars$nreps,
              strat = pars$strat,
              div = pars$div
            )
          },
          "boot" = {
            train.ind <- boot.idx(cl, pars$nreps, strat = pars$strat)
          }
        )
        for (i in 1:length(train.ind)) {
          if (any(table(cl[train.ind[[i]]]) < 2) ||
            any(table(cl[-train.ind[[i]]]) < 1)) {
            emp_f <- !emp_f
            break
          }
        } #' end of for
      } #' end of while
    } #' end of else

    return(train.ind)
  }

  tr.idx <- list()
  for (i in 1:pars$niter) {
    tr.idx[[i]] <- idx.func(cl, pars = pars)
  }
  names(tr.idx) <- paste("Iter_", 1:pars$niter, sep = "")

  return(tr.idx)
}

#' ========================================================================
boot.idx <- function(x, nreps, strat = FALSE) {
  n <- length(x)

  if (strat) {
    x <- factor(x) #' drops the levels that do not occur
    idx <- sample(1:n, n, replace = F)
    #' shuffle the original x, #' idx  <- c(1:n)
    x <- x[idx]

    v <- length(levels(x))
    #' index of each factor
    #' s.idx <- lapply(1:v, function(i) idx[which(x == levels(x)[i])])
    #' lwc-05-10-2010: bug fixed
    s.idx <- lapply(1:v, function(i) which(x == levels(x)[i]))

    train.ind <- lapply(1:nreps, function(y) { #' y is not used.
      tmp <- lapply(s.idx, function(x) sample(x, length(x), replace = T))
      do.call("c", tmp)
    })
    #' shuffle the results
    train.ind <- lapply(train.ind, function(x) sample(x, length(x), replace = F))
  } else {
    train.ind <- lapply(1:nreps, function(x) sample(1:n, n, replace = T))
  }
  return(train.ind)
}

#' ========================================================================
rand.idx <- function(x, nreps, strat = FALSE, div = 2 / 3) {
  n <- length(x)

  if (strat) {
    x <- factor(x) #' drops the levels that do not occur
    idx <- sample(1:n, n, replace = F)
    #' shuffl the original x, #' idx  <- c(1:n)
    x <- x[idx]

    v <- length(levels(x))
    #' index of each factor
    #' s.idx <- lapply(1:v, function(i) idx[which(x == levels(x)[i])])
    #' lwc-05-10-2010: bug fixed
    s.idx <- lapply(1:v, function(i) which(x == levels(x)[i]))

    train.ind <-
      lapply(1:nreps, function(y) { #' y is not used.
        tmp <- lapply(
          s.idx,
          function(x) sample(x, trunc(length(x) * div), replace = F)
        )
        do.call("c", tmp)
      })
    #' shuffl the results
    train.ind <-
      lapply(train.ind, function(x) sample(x, length(x), replace = F))
  } else {
    train.ind <-
      lapply(1:nreps, function(x) sample(1:n, trunc(n * div), replace = F))
  }

  return(train.ind)
}

#' ========================================================================
cv.idx <- function(x, nreps, strat = FALSE) {

  #' One change has been made to get different results for each calling.
  #' from package ipred
  ssubset <- function(y, k, strat = TRUE) {
    #' if (!is.factor(y)) stop("y is not of class factor") #' lwc-14-04-2007
    if (!is.factor(y)) y <- as.factor(y)

    N <- length(y)
    nlevel <- table(y)
    nindx <- list()

    #' lwc-29-10-2006: changed
    indx <- sample(1:N, N, replace = F)
    y <- y[indx]
    #'  indx <- 1:N
    outindx <- list()
    if (strat) {
      for (j in 1:length(nlevel)) {
        nindx <- c(nindx, list(indx[which(y == levels(y)[j])]))
      }
      kmat <- kfoldcv(k, N, nlevel)
      for (i in 1:k) {
        sset <- kmat[, i]
        kindx <- c()
        for (j in 1:length(nlevel)) {
          if (i > 1) {
            kindx <- c(kindx, nindx[[j]][(sum(kmat[
              j,
              1:(i - 1)
            ]) + 1):sum(kmat[j, 1:i])])
          } else {
            kindx <- c(kindx, nindx[[j]][1:kmat[j, 1]])
          }
        }
        kindx <- kindx[!is.na(kindx)]
        outindx <- c(outindx, list(kindx))
      }
      return(outindx)
    } else {
      kmat <- kfoldcv(k, N)
      nindx <- indx
      for (i in 1:k) {
        if (i > 1) {
          outindx <- c(
            outindx,
            list(nindx[(sum(kmat[1:(i - 1)]) + 1):sum(kmat[1:i])])
          )
        } else {
          outindx <- c(outindx, list(nindx[1:kmat[1]]))
        }
      }
    }
    return(outindx)
  }

  #' from package ipred
  kfoldcv <- function(k, N, nlevel = NULL) {
    if (is.null(nlevel)) { # no stratification
      if (k > N) {
        return(c(rep(1, N), rep(0, k - N)))
      }
      fl <- floor(N / k)
      ce <- ceiling(N / k)
      if (fl == ce) {
        return(rep(fl, k))
      } else {
        return(c(rep(ce, round((N / k - fl) * k)), rep(fl, round((1 - (N / k -
          fl)) * k))))
      }
    } else { # stratification
      # if (!is.integer(nlevel)) stop("nlevel is not a vector if integers")
      kmat <- matrix(0, ncol = k, nrow = length(nlevel))
      for (i in 1:length(nlevel)) {
        kmat[i, ] <- kfoldcv(k, nlevel[i])
      }
      return(kmat)
    }
  }

  n <- length(x)
  #' get index of test
  test.ind <- ssubset(x, nreps, strat = strat)
  #' get index of training
  train.ind <- lapply(1:nreps, function(i) seq(1, n)[-test.ind[[i]]])
  #' shuffle the results
  train.ind <- lapply(train.ind, function(x) sample(x, length(x), replace = F))
  return(train.ind)
}

#'  1) accest.default
#'  2) print.accest
#'  3) summary.accest
#'  4) print.summary.accest
#'  5) plot.accest
#'  6) accest
#'  7) accest.formula
#'  8) binest
#'  9) aam.mcl
#' 10) aam.cl
#' 11) cl.rate
#' 12) cl.perf
#' 13) cl.auc
#' 14) auc
#' 15) cl.roc
#' 16) cl.roc.1
#' 17) classifier
#' 18) .marg
#' 19) boot.err
#' 20) valipars
#' 21) trainind
#' 22) boot.idx
#' 23) rand.idx
#' 24) cv.idx
