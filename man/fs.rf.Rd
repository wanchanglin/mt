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

\author{ Wanchang Lin  }

% ----------------------------------------------------------------------
\examples{
data(abr1)
cls <- factor(abr1$fact$class)
dat <- abr1$pos

## fill zeros with NAs
dat <- mv.zene(dat)

## missing values summary
mv <- mv.stats(dat, grp=cls) 
mv    ## View the missing value pattern

## filter missing value variables
dat <- dat[,mv$mv.var < 0.15]

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
score <- sort(score, decreasing = TRUE)
plot(score)

}
\keyword{classif}
