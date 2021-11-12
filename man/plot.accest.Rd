% lwc-10-10-2006: First draft
%
\name{plot.accest}
\alias{plot.accest}
\title{
  Plot Method for Class 'accest'
}
\description{
  Plot accuracy rate of each iteration.
}

\usage{
\method{plot}{accest}(x, main = NULL, xlab = NULL, ylab = NULL, \dots)
}

% ----------------------------------------------------------------------------
\arguments{
  \item{x}{
  An object of class \code{accest}.
  }
  \item{main}{
  An overall title for the plot.  
  }
  \item{xlab}{
   A title for the x axis.
  }
  \item{ylab}{
  A title for the y axis.
  }
  \item{\dots}{
  Additional arguments to the plot. 
  }
}

\value{
  Returns plot of class \code{accest}.
}

% ----------------------------------------------------------------------------
\details{
This function is a method for the generic function \code{plot()} for class 
\code{accest}. It plots the accuracy rate against the index of iterations. 
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}


\seealso{
  \code{\link{accest}}
}

\examples{
# Iris data
data(iris)
# Stratified cross-validation of PCALDA for Iris data
pars <- valipars(sampling="cv", niter=10, nreps=10, strat=TRUE)
acc  <- accest(Species~., data = iris, method = "pcalda", pars = pars)
              
acc
summary(acc)
plot(acc)
}

\keyword{plot}
