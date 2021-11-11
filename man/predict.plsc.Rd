% wll-23-05-2007: 
% wll-03-10-2007: add plsc stuff 
%
\name{predict.plsc}
\alias{predict.plsc}
\alias{predict.plslda}

\title{
  Predict Method for Class 'plsc' or 'plslda'
}
\description{
  Prediction of test data using \code{plsc} or \code{plslda}.
}
\usage{
\method{predict}{plsc}(object, newdata,\dots)
\method{predict}{plslda}(object, newdata,\dots)
}
% ----------------------------------------------------------------------------
\arguments{
  \item{object}{
  Object of class \code{plsc} or \code{plslda}.
  }
  \item{newdata}{
  A matrix or data frame of cases to be classified.  
  }
  \item{\dots}{
  Arguments based from or to other methods.
  }
}
% ----------------------------------------------------------------------------
\details{
Two functions are methods for the generic function \code{predict()} for
class \code{plsc} or \code{plslda}. If \code{newdata} is omitted, the results of 
training data in \code{plsc} or \code{plslda} object will be returned.
}
% ----------------------------------------------------------------------------
\value{
  A list with components:
  \item{class}{
  The predicted class (a factor).
  }
  \item{posterior}{
  The posterior probabilities for the predicted classes.
  }
  \item{x}{
  The rotated test data by the projection matrix of PLS.
  }
}
% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}
% ----------------------------------------------------------------------------
\seealso{
\code{\link{plsc}}, \code{\link{plot.plsc}},\code{\link{plslda}}, \code{\link{plot.plslda}}
}
% ----------------------------------------------------------------------------
\examples{
data(iris3)

tr    <- sample(1:50, 25)
train <- rbind(iris3[tr,,1], iris3[tr,,2], iris3[tr,,3])
test  <- rbind(iris3[-tr,,1], iris3[-tr,,2], iris3[-tr,,3])
cl    <- factor(c(rep("s",25), rep("c",25), rep("v",25)))

## model fit using plsc and plslda without tuning of ncomp
(z.plsc       <- plsc(train, cl))
(z.plslda     <- plslda(train, cl))
## predict for test data
pred.plsc    <- predict(z.plsc, test)
pred.plslda  <- predict(z.plslda, test)

## plot the projected test data.
grpplot(pred.plsc$x, pred.plsc$class, main="PLSC: Iris") 
grpplot(pred.plslda$x, pred.plslda$class, main="PLSLDA: Iris") 
}

\keyword{classif}
