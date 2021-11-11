% lwc-12-04-2007

\name{fs.bw}
\alias{fs.bw}
\title{
  Feature Selection Using Between-Group to Within-Group (BW) Ratio
}

\description{
  Feature selection using ratio of between-group to within-group sums of squares 
  (BW).
}
\usage{
  fs.bw(x,y,\dots)
}

% --------------------------------------------------------------------
\arguments{
  \item{x}{
  A data frame or matrix of data set. 
  }
  \item{y}{
  A factor or vector of class.
  }
  \item{\dots}{
  Arguments to pass(current ignored).
  }
}


% ----------------------------------------------------------------------------
\value{
  A list with components:
  \item{fs.rank}{A vector of feature ranking scores.}
  \item{fs.order}{A vector of feature order from best to worst.}
  \item{stats}{A vector of measurements.}
}

\references{
  Dudoit, S., Fridlyand, J. and Speed, T.P. Comparison of discrimination methods 
  for classification of tumours using gene expression data. 
  \emph{Journal of the American Statistical Association}. Vol.97, No.457, 77-87.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
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

## apply BW ratio method for feature selection/ranking
res <- fs.bw(mat,grp)
names(res)

}

\keyword{classif}
