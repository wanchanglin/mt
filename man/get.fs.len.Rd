% lwc-01-11-2006

\name{get.fs.len}
\alias{get.fs.len}
\title{
  Get Length of Feature Subset for Validation
}
\description{
  Get feature lengths for feature selection validation by classification.
}
\usage{
  get.fs.len(p,fs.len=c("power2"))
}

% --------------------------------------------------------------------
\arguments{
  \item{p}{
  Number of features in the data set. 
  }
  \item{fs.len}{
    Method or numeric sequence for feature lengths. It can be a numeric vector as
    user-defined feature lengths, or methods:
    \itemize{
      \item \code{full}. The feature lengths are \code{p, \dots, 2, 1}. This is an exhaustive 
      method. If \code{p} is too large, it will consume a lot of time and hence it is not practical.
      \item \code{half}. The feature lengths are the sequence of 
      \code{p, p/2, p/2/2, \dots, 1}.
      \item \code{power2}. The feature lengths are the sequence of 
      \code{p, 2^(log2(p)-1), \dots, 2^1, 2^0}. 
    }
  }
}

% ----------------------------------------------------------------------------
\value{
  An descending order numeric vector of feature lengths. 
}

\details{
	The generation of feature length is used in the validation of feature subsets 
	by classification. The feature length decide the lengths of feature subset starting 
	from top of the full feature order list.
}

\note{
  The last length of feature returned is always \code{p}.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

\seealso{
  \code{\link{fs.rfe}}, \code{\link{frank.err}},  \code{\link{frankvali}}
}

% ----------------------------------------------------------------------
\examples{

data(abr1)
dat <- abr1$pos

## number of featres
p <- ncol(dat)

## predefined feature lengths. The returned will be descending order 
## vector with the first one is p.
(vec <- get.fs.len(p, fs.len=c(1,2,3,4,5,6,7,8,9,10)))

## use all features as feature lengths
(vec.full <- get.fs.len(p, fs.len="full"))

## use "half"
(vec.half <- get.fs.len(p, fs.len="half"))

## use "power2"
(vec.power2 <- get.fs.len(p, fs.len="power2"))

}

\keyword{manip}
