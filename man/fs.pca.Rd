% wll-04-02-2008

\name{fs.pca}
\alias{fs.pca}

\title{
  Feature Selection by PCA  
}
\description{
  Feature selection using PCA loadings.
}
\usage{
  fs.pca(x,thres=0.8, \dots)
}

% --------------------------------------------------------------------
\arguments{
  \item{x}{
  A data frame or matrix of data set. 
  }
  \item{thres}{
    The threshold of the cumulative percentage of PC's explained variances.
  }
  \item{\dots}{
    Additional arguments to \code{\link{prcomp}}. 
  }
}

% ----------------------------------------------------------------------------
\value{
  A list with components:
  \item{fs.rank}{A vector of feature ranking scores.}
  \item{fs.order}{A vector of feature order from best to worst.}
  \item{stats}{A vector of measurements.}
}

\details{
  Since PCA loadings is a matrix with respect to PCs, the Mahalanobis distance of 
  loadings is applied to select the features. (Other ways, for example, the sum 
  of absolute values of loadings, or squared root of loadings, can be used.)
  
  It should be noticed that this feature selection method is unsupervised.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

\seealso{
  \code{\link{feat.rank.re}}
}

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

## feature selection by PCA
res <- fs.pca(dat)
names(res)

}

\keyword{classif}
