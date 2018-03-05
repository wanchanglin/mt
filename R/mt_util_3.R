## wll-25-11-2015: all functions are not documented.
## wll-25-11-2015: What is overlap's purpose?

## =========================================================================
## lwc-17-05-2011: slight modification of colwise in package plyr.
##   Prevent data frame with string as factors.
## Usages:
##   a <- data.frame(str=letters[1:5],num=c(1:5))
##   a[,1]        ## factors
##   a[,2]        ## numeric
##   b <- colwise.1(as.character)(a)
##   b[,1]        ## chracters
## --------------------------------------------------------------------------
colwise.1 <- function(.fun, .cols = true) {
  if (!is.function(.cols)) {
    .cols <- as.quoted(.cols)
    filter <- function(df) as.data.frame(eval.quoted(.cols, df))
  } else {
    filter <- function(df) Filter(.cols, df)
  }

  function(df, ...) {
    stopifnot(is.data.frame(df))
    filtered <- filter(df)
    if (ncol(filtered) == 0) return(data.frame())

    df <- as.data.frame(lapply(filtered, .fun, ...), stringsAsFactors=F)
    ## df <- as.data.frame(lapply(filtered, .fun, ...))
    names(df) <- names(filtered)
    df
  }
}


## ======================================================================
## wll-18-02-2008: calculate the overlap.
## Usage:
##   x <- c(2.225350, 2.743615) ## lower and upper
##   y <- c(3.157918, 3.264999) ## lower and upper
##   overlap(x,y)
## Return:
##   A scalar values. Possitive indicates no overlap and width of boundary
##   and negative for the depth of overlap.
overlap <- function(x,y){
  tmp   <- c(x,y)
  tmp.s <- sort(tmp)

  ## sum of two differences
  dif <- sum(diff(tmp)[c(1,3)])

  ## difference between max and min of all 4 values
  dif1 <- diff(tmp.s[c(1,4)])

  overlap <- dif1 - dif
  names(overlap) <- NULL
  overlap
}

## ======================================================================
## wll-18-02-2008: Wrapper function for pairwise.lap
compare.lap <- function(i, j, x) {
  ## i = "Ailsa"
  ## j = "Brodick"
  xi <- x[i,]
  xj <- x[j,]
  overlap(xi,xj)
}

## ===========================================================================
## wll-18-02-2008:Creates table of lap values for pairwise comparisons
## Note: This function is modified from pairwise.table in stats package
## Arguments:
##   compare.lap - Function to compute (raw) p value given indices i and j
##   level.names - Names of the group levels
##
pairwise.lap <- function(compare.lap, level.names, x) {
  ix <- seq(along=level.names)
  names(ix) <- level.names
  pp <- outer(ix[-1], ix[-length(ix)],function(ivec, jvec)
    sapply(seq(along=ivec), function(k) {
             i<-ivec[k]
             j<-jvec[k]
             if (i > j) compare.lap(names(i), names(j), x) else NA ## wll-18-02-2008: add names
           }))

  lap   <- pp[lower.tri(pp,T)]
  ## sort out the pairwise names
  dname <- dimnames(pp)
  tmp   <- outer(dname[[1]], dname[[2]],paste, sep="~")
  names(lap) <- tmp[lower.tri(tmp,T)]

  lap
  ##  pp
}

## =======================================================================
## wll-02-06-2008: pairwise differences
## Usages:
##   x <- rnorm(6)
##   names(x) = paste("V", 1:length(x), sep="")
##   p.diff(x)
p.diff <- function(x){
  n <- length(x)
  if (is.null(names(x))) names(x) <- 1:n

  prs <- cbind(rep(1:n, each = n), 1:n)
  com <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
  res <- x[com[, 2]] - x[com[, 1]]
  names(res) <- paste(names(x)[com[, 2]], names(x)[com[, 1]], sep = "~")
  res
}

## ========================================================================
## wll-02-06-2008: Marc Schwartz's pairwise differences
## Usage:
##   mat <- matrix(rnorm(15), nrow=5)
##   colnames(mat) <- paste("C",1:3, sep="")
##   rownames(mat) <- paste("R",1:5, sep="")
##   pairwise.diffs(mat)
pairwise.diffs <- function(x) {
  stopifnot(is.matrix(x))

                                        # create column combination pairs
  prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
  col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]

                                        # do pairwise differences
  result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]

                                        # set colnames
  if(is.null(colnames(x)))
    colnames(x) <- 1:ncol(x)

  colnames(result) <- paste(colnames(x)[col.diffs[, 1]], ".vs.",
                            colnames(x)[col.diffs[, 2]], sep = "")
  result
}

## ==========================================================================
## wll-18-06-2008: Nice way to print pairwise comparison matrix such as of
##   the coefficient of correlation analysis and pairwise test
##   NOTE: x must be a pairwise matrix of correlation analysis or p-values.
pairwise.print <- function(x){
  x[upper.tri(x)] <- NA
  diag(x) <- NA
  x <- x[-1,-ncol(x)]
  print(x, na.print="")      ## debug
  ## return(x)
}

## ========================================================================
## lwc-28-11-2013: Geometric mean.
## Note 1: The geometric mean is defined as the nth root (where n is the
##   count of numbers) of the product of the numbers.
## Note 2: Geometric mean is usually restricted to positive inputs, because
##   otherwise the answer can have an imaginary component. If you really want
##   the geometric mean of negative inputs, use the second method but convert
##   the input to be a complex number first.
##
## comp.x <- as.complex(c(-5,-4,4,5))
## geom.mean2(comp.x)
## [1] 0+4.472136i
## -------------------------------------------------------------------------
## geo.mean1 <- function(x) prod(x)^(1/length(x))
## geo.mean2 <- function(x) exp(mean(log(x)))
geo.mean <- function(x, na.rm = FALSE){
	if(na.rm) x <- x[!is.na(x)]
	n <- length(x)
	prod(x)^(1/n)
}

## ======================================================================
## lwc-05-08-2010: Round 0.5 upwards by  Duncan Murdoch.
## Note: Excel always round up numbers ending in 5 but R's round dosen't.
roundup <- function(x, digits=0) {
  factor <- 10^digits
  trunc(x*factor + 0.5)/factor
}

## =========================================================================
## lwc-27-07-2010: memory management
## R is an in-memory application, so every new object you create takes up
## RAM. (Yes, there are ways around that, but that's a topic for another
## article.) If you're working on a small machine (say, a 32-bit Windows
## system with 1Gb of RAM or less) you might need to be careful with the
## object you create. This StackOverflow question offers some useful tips
## for managing objects in RAM, including code for the function lsos to list
## objects sorted by size: improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

## shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

## =========================================================================
## lwc-08-07-2014: As there are zero counts as well, we use a convenience
##  function clog() providing a continuity-corrected logarithm. Taken from
##  "Regression Models for Count Data in R".
clog <- function(x) log(x + 0.5)

## =========================================================================
## lwc-05-09-2013: Get all objects size
my.ls <- function (pos = 1, sorted = FALSE, envir = as.environment(pos)) {
  .result <- sapply(ls(envir = envir, all.names = TRUE), function(..x)
    object.size(eval(as.symbol(..x), envir = envir))
                    )
  if (sorted) {
    .result <- rev(sort(.result))
  }
  .ls <- as.data.frame(rbind(as.matrix(.result), `**Total` = sum(.result)))
  names(.ls) <- "Size"
  .ls$Size <- formatC(.ls$Size, big.mark = ",", digits = 0, format = "f")
  .ls$Mode <- c(unlist(lapply(rownames(.ls)[-nrow(.ls)], function(x)
    mode(eval(as.symbol(x), envir = envir)))), "-------")
  .ls
}


## ==================================================
## TOC on 25-11-2015
## ==================================================
## (1). colwise.1
## (4). overlap
## (5). compare.lap
## (6). pairwise.lap
## (7). p.diff
## (8). pairwise.diffs
## (9). pairwise.print
## (12). geo.mean
## (13). roundup
## (14). .ls.objects
## (16). lsos
## (17). clog
## (18). my.ls
