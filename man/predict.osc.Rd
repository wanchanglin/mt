% wll-05-06-2007: 
%
\name{predict.osc}
\alias{predict.osc}
\title{
  Predict Method for Class 'osc'
}
\description{
  Pre-processing of new data by \code{osc}.
}
\usage{
\method{predict}{osc}(object, newdata,\dots)
}
% ----------------------------------------------------------------------------
\arguments{
  \item{object}{
  Object of class \code{osc}.
  }
  \item{newdata}{
  A matrix or data frame of cases to be corrected by OSC.  
  }
  \item{\dots}{
  Arguments based from or to other methods.
  }
}
% ----------------------------------------------------------------------------
\details{
This function is a method for the generic function \code{predict()} for
class \code{osc}. If \code{newdata} is omitted, the corrected data set used in model of 
\code{osc} will be returned.
}
% ----------------------------------------------------------------------------
\value{
  A list containing the following components:
  \item{x}{
    A matrix of OSC corrected data set.
  }
  \item{Q2}{
  The fraction of variation in X after OSC correction for the new data.
  }
}
% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}
% ----------------------------------------------------------------------------
\seealso{
  \code{\link{osc}}, \code{\link{osc_wold}}, \code{\link{osc_sjoblom}},
  \code{\link{osc_wise}}
}
% ----------------------------------------------------------------------------
\examples{
data(abr1)
cl   <- factor(abr1$fact$class)
dat  <- abr1$pos

## divide data as training and test data
idx <- sample(1:nrow(dat), round((2/3)*nrow(dat)), replace=FALSE) 

## construct train and test data 
train.dat  <- dat[idx,]
train.t    <- cl[idx]
test.dat   <- dat[-idx,]        
test.t     <- cl[-idx] 

## build OSC model based on the training data
res <- osc(train.dat, train.t, method="wold",osc.ncomp=2, pls.ncomp=4)
names(res)
res
summary(res)

## pre-process test data by OSC
test <- predict(res,test.dat)
test.dat.1 <- test$x
}

\keyword{manip}
