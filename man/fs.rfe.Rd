% lwc-21-03-2007

\name{fs.rfe}
\alias{fs.rfe}
\title{
  Feature Selection Using SVM-RFE 
}
\description{
  Feature selection using Support Vector Machine based on Recursive Feature Elimination (SVM-RFE) 
}
\usage{
  fs.rfe(x,y,fs.len="power2",\dots)
}

% --------------------------------------------------------------------
\arguments{
  \item{x}{
  A data frame or matrix of data set. 
  }
  \item{y}{
  A factor or vector of class.
  }
  \item{fs.len}{
  Method for feature lengths used in SVM-RFE computation.  
  For details, see \code{\link{get.fs.len}}. 

  }
  \item{\dots}{
  Arguments to pass to \code{svm}.
  }
}

% ----------------------------------------------------------------------------
\value{
  A list with components:
  \item{fs.rank}{A vector of feature ranking scroes.}
  \item{fs.order}{A vector of feature order from best to worst.}
}


% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

\seealso{
  \code{\link{feat.rank.re}}, \code{\link{get.fs.len}}
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

## apply RFE method for feature selection/ranking
res <- fs.rfe(mat,grp)
names(res)

}

\keyword{classif}
