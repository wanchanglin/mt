% wll-22-11-2007

\name{fs.pls}
\alias{fs.pls}
\alias{fs.plsvip}
\alias{fs.plsvip.1}
\alias{fs.plsvip.2}

\title{
  Feature Selection Using PLS  
}
\description{
  Feature selection using coefficient of regression and VIP values of PLS.
}
\usage{
  fs.pls(x,y, pls="simpls",ncomp=10,\dots)
  fs.plsvip(x,y, ncomp=10,\dots)
  fs.plsvip.1(x,y, ncomp=10,\dots)
  fs.plsvip.2(x,y, ncomp=10,\dots)
}

% --------------------------------------------------------------------
\arguments{
  \item{x}{
  A data frame or matrix of data set. 
  }
  \item{y}{
  A factor or vector of class.
  }
  \item{pls}{
  A method for calculating PLS scores and loadings. The following methods are supported:
    \itemize{
      \item \code{simpls:} SIMPLS algorithm. 
      \item \code{kernelpls:} kernel algorithm.
      \item \code{oscorespls:} orthogonal scores algorithm. 
    }
  For details, see \code{\link[pls]{simpls.fit}}, \code{\link[pls]{kernelpls.fit}} and
  \code{\link[pls]{oscorespls.fit}} in package \pkg{pls}.
  }

  \item{ncomp}{
    The number of components to be used.  
  }
  \item{\dots}{
    Arguments passed to or from other methods.
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
  \code{fs.pls} ranks the features by regression coefficient of PLS. Since the 
  coefficient is a matrix due to the dummy multiple response variables designed 
  for the classification (category) problem, the Mahalanobis distance of 
  coefficient is applied to select the features. (Other ways, for example, the sum 
  of absolute values of coefficient, or squared root of coefficient, can be used.)
  
  \code{fs.plsvip} and \code{fs.plsvip.1} carry out feature selection based on the 
  the Mahalanobis distance and absolute values of PLS's VIP, respectively. 
  
  \code{fs.plsvip.2} is similar to \code{fs.plsvip} and \code{fs.plsvip.1}, but 
  the category response is not treated as dummy multiple response matrix. 
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

## apply PLS methods for feature selection
res.pls      <- fs.pls(mat,grp, ncomp=4)
res.plsvip   <- fs.plsvip(mat,grp, ncomp=4)
res.plsvip.1 <- fs.plsvip.1(mat,grp, ncomp=4)
res.plsvip.2 <- fs.plsvip.2(mat,grp, ncomp=4)

## check differences among these methods
fs.order <- data.frame(pls      = res.pls$fs.order,
                       plsvip   = res.plsvip$fs.order,
                       plsvip.1 = res.plsvip.1$fs.order,
                       plsvip.2 = res.plsvip.2$fs.order)
head(fs.order, 20)

}

\keyword{classif}
