% wll-02-07-2007: 
%
\name{predict.pcalda}
\alias{predict.pcalda}
\title{
  Predict Method for Class 'pcalda'
}
\description{
  Prediction of test data using \code{pcalda}.
}
\usage{
\method{predict}{pcalda}(object, newdata,\dots)
}
% ----------------------------------------------------------------------------
\arguments{
  \item{object}{
  Object of class \code{pcalda}.
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
This function is a method for the generic function \code{predict()} for
class \code{pcalda}. If \code{newdata} is omitted, the results of training data 
in \code{pcalda} object will be returned.
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
  The rotated test data by the projection matrix of LDA.
  }
}
% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}
% ----------------------------------------------------------------------------
\seealso{
\code{\link{pcalda}}, \code{\link{plot.pcalda}}
}
% ----------------------------------------------------------------------------
\examples{
data(iris3)

tr    <- sample(1:50, 25)
train <- rbind(iris3[tr,,1], iris3[tr,,2], iris3[tr,,3])
test  <- rbind(iris3[-tr,,1], iris3[-tr,,2], iris3[-tr,,3])
cl    <- factor(c(rep("s",25), rep("c",25), rep("v",25)))

z     <- pcalda(train, cl)
pred  <- predict(z, test)

## plot the projected data.
grpplot(pred$x, pred$class, main="PCALDA: Iris") 
}

\keyword{classif}
