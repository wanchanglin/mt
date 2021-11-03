#' lwc-07-09-2010, Tue: Functions not documented

#' ========================================================================
#' Generates Class Indicator Matrix from a Factor. 
#' A matrix which is zero except for the column corresponding to the class.
#' Internal function.  From package NNET
class.ind <- function(cl) {
  n <- length(cl)
  cl <- as.factor(cl)
  x <- matrix(0, n, length(levels(cl)))
  x[(1:n) + n * (unclass(cl) - 1)] <- 1
  dimnames(x) <- list(names(cl), levels(cl))
  x
}

#' ========================================================================
#' lwc-07-07-2011: batch shifting: remove mean withing each batch/block
#' Internal function.
#' Arguments:
#'  x - data matrix
#'  y - categorical data for batch/block information
#' References:
#'   Silvia Wagner, et.al, Tools in Metabonomics: An Integrated Validation
#'   Approach for LC-MS Metabolic Profiling of Mercapturic Acids in Human
#'   Urine Anal. Chem., 2007, 79 (7), pp 2918-2926, DOI: 10.1021/ac062153w
batch.shift <- function(x, y, type = "mean") {
  x <- as.data.frame(x)

  g.mean <- sapply(x, function(x) tapply(x, y, type, na.rm = T))
  g.mean <- sapply(1:ncol(x), function(i) g.mean[, i][y])
  x <- x - g.mean

  return(x)
}

#' ========================================================================
#' lwc-02-06-2011: Relative standard deviation of data in column
#' Internal function.
rsd <- function(x) {
  mn <- colMeans(x, na.rm = TRUE)
  std <- apply(x, 2, sd, na.rm = TRUE)
  #' std <- sd(x,na.rm=TRUE)           #' sd(<data.frame>) is deprecated.
  res <- 100 * std / mn
  return(res)
}

#' =========================================================================
#' tic() and toc() functions.
#' Internal function.
#' Modifed from package MATLAB by David Enot
tic <- function(gcFirst = FALSE) {
  if (gcFirst == TRUE) {
    gc(verbose = FALSE)
  }
  assign("savedTime", proc.time()[3], envir = .GlobalEnv)
  invisible()
}

toc <- function(echo = TRUE) {
  prevTime <- get("savedTime", envir = .GlobalEnv)
  diffTimeSecs <- proc.time()[3] - prevTime
  if (echo) {
    cat(sprintf("elapsed time is %f seconds", diffTimeSecs), "\n")
    return(invisible())
  }
  else {
    return(diffTimeSecs)
  }
}

#' ========================================================================
#' lwc-17-05-2011: Write list of data matrix to an Excel's XLS file.
#' Internal function.
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
#' lwc-06-04-2011: write data frames to an Excel file using Perl
#' Internal function. 
#' Slight modification of WriteXLS.R (2.1.0, 2010-09-18) 
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

#' =======================================================================
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

#'  1) class.ind
#'  2) batch.shift
#'  3) rsd
#'  4) tic
#'  5) toc
#'  6) list2xls
#'  7) WriteXLS
#'  8) sub.samp
#'  9) panel.qqconf
#' 10) LogReg
#' 11) predict.LogReg
#' 12) LogRegAnn
#' 13) predict.LogRegAnn
#' 14) pred.pls
