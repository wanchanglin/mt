%% lwc-21-03-2007: commence
%% lwc-16-02-2010: add fs.rf.1
%% lwc-18-10-2011: re-write examples. Will use these exaples for all
%% feature selection.
%% wll-08-12-2015: add codes for comparison of 'fs.rf' and 'fs.rf.1'.

\name{fs.rf}
\alias{fs.rf}
\alias{fs.rf.1}

\title{ Feature Selection Using Random Forests (RF) }

\description{ Feature selection using Random Forests (RF). }

\usage{
  fs.rf(x,y,\dots)
  fs.rf.1(x,y,fs.len="power2",\dots)
}

% --------------------------------------------------------------------
\arguments{

  \item{x}{ A data frame or matrix of data set. }
  \item{y}{ A factor or vector of class. }
  \item{fs.len}{ Method or numeric sequence for feature lengths. For
    details, see \code{\link{get.fs.len}} }
  \item{\dots}{ Arguments to pass to \code{randomForests}. }
}

\value{
  A list with components:
  \item{fs.rank}{A vector of feature ranking scores.}
  \item{fs.order}{A vector of feature order from best to worst.}
  \item{stats}{A vector of measurements. For \code{fs.rf}, it is Random
    Forest important score. For \code{fs.rf.1}, it is a dummy variable
    (current ignored). }
}

\details{
  \code{fs.rf.1} select features based on successively eliminating the least 
  important variables.
}

\author{ Wanchang Lin \email{wll@aber.ac.uk} }

% ----------------------------------------------------------------------
\examples{
## prepare data set
data(abr1)
cls <- factor(abr1$fact$class)
dat <- abr1$pos
## dat <- abr1$pos[,110:1930]

## fill zeros with NAs
dat <- mv.zene(dat)

## missing values summary
mv <- mv.stats(dat, grp=cls) 
mv    ## View the missing value pattern

## filter missing value variables
## dim(dat)
dat <- dat[,mv$mv.var < 0.15]
## dim(dat)

## fill NAs with mean
dat <- mv.fill(dat,method="mean")

## log transformation
dat <- preproc(dat, method="log10")

## select class "1" and "2" for feature ranking
ind <- grepl("1|2", cls)
mat <- dat[ind,,drop=FALSE] 
mat <- as.matrix(mat)
grp <- cls[ind, drop=TRUE]   

## apply random forests for feature selection/ranking
res   <- fs.rf(mat,grp)
res.1 <- fs.rf.1(mat,grp)

## compare the results
fs <- cbind(fs.rf=res$fs.order, fs.rf.1=res.1$fs.order)

## plot the important score of 'fs.rf' (not 'fs.rf.1')
score <- res$stats
score <- sort(score, decreasing=T)
plot(score)

## --------------------------------------------------------------------
## Since Random Forest is an random method, the results are different in
## different runs.
set.seed(1234)

## apply random forests for feature selection/ranking
tmp <- lapply(1:10, function(x){
  res   <- fs.rf(mat,grp)
  res.1 <- fs.rf.1(mat,grp)
  
  list(order = res$fs.order, order.1 = res.1$fs.order,
       rank  = res$fs.rank,  rank.1  = res.1$fs.rank)
})

## compare the results:
## 1.) feature orders
order <- lapply(tmp, function(x) x$order)
order <- do.call(cbind, order)

order.1 <- lapply(tmp, function(x) x$order.1)
order.1 <- do.call(cbind, order.1)

## 2.) feature ranks
rank <- lapply(tmp, function(x) x$rank)
rank <- do.call(cbind, rank)

rank.1 <- lapply(tmp, function(x) x$rank.1)
rank.1 <- do.call(cbind, rank.1)

}

\keyword{classif}
