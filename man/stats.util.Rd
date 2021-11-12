%% lwc-09-01-2014: commence
%% wll-01-12-2015: add an argument for 'fc'
\name{stats.util}
\alias{stats.mat}
\alias{stats.vec}

\title{
  Statistical Summary Utilities for Two-Classes Data
}
\description{
  Functions to summarise two-group data.
}
\usage{
stats.vec(x,y, method="mean",test.method = "wilcox.test",fc=TRUE,\dots)
stats.mat(x,y, method="mean",test.method = "wilcox.test",
          padj.method= "fdr",fc=TRUE,\dots)
}

% --------------------------------------------------------------------
\arguments{

  \item{x}{ A data frame or matrix of data set for \code{stats.mat} or
    data vector for \code{stats.vec}. }
  
  \item{y}{ A factor or vector of class. Two classes only. }

  \item{method}{method for group center such as \code{\link{mean}} or
    \code{\link{median}}. }

  \item{test.method}{method for p-values from parametric test such as
    \code{\link{t.test}} or non-parametric test such as
    \code{\link{wilcox.test}}. }

  \item{padj.method}{method for p-values correction. Can be one in
    \code{\link{p.adjust.methods}}: "holm", "hochberg", "hommel",
    "bonferroni", "BH", "BY", "fdr" and "none". }

  \item{fc}{a flag for fold-change.} 

  \item{\dots}{ Additional parameters. }
  
}

% ----------------------------------------------------------------------
\value{
  
  \code{stats.vec} returns an vector of center, group center, fold change, 
  log2 fold change, AUC and  p-value. 

  \code{stats.mat}, which is an wrapper function of \code{stats.vec},
  returns an data frame of center, group center, fold change, log2 
  fold-change, AUC, p-value and adjusted p-values.
  
}

\note{
  If \code{x} has negative values, the results of fold-change and
  log2 fold-change are not reliable.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------
\examples{
data(iris)
sel <- c("setosa", "versicolor")
grp <- iris[,5]
idx <- match(iris[,5],sel,nomatch = 0) > 0

cls <- factor(iris[idx,5])
dat <- iris[idx,1:4]

## stats for the single vector
stats.vec(dat[[1]],cls, method = "median",test.method = "wilcox.test")

## use matrix format
tab <- stats.mat(dat,cls, method = "median",test.method = "wilcox.test",
                 padj.method = "BH")
sapply(tab, class)

}

\keyword{util}
