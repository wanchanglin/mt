% lwc-30-01-2007: First draft
%
\name{plot.maccest}
\alias{plot.maccest}
\title{
  Plot Method for Class 'maccest'
}
\description{
  Plot accuracy rate with standard derivation of each classifier.
}

\usage{
\method{plot}{maccest}(x, main = NULL, xlab = NULL, ylab = NULL, \dots)
}

% ----------------------------------------------------------------------------
\arguments{
  \item{x}{
  An object of class \code{maccest}.
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

% ----------------------------------------------------------------------------
\details{
This function is a method for the generic function \code{plot()} for class 
\code{maccest}. It plots the accuracy rate with standard derivation against the 
classifiers. 
}

\value{
  Returns plot of class \code{maccest}.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}


\seealso{
  \code{\link{maccest}}, \code{\link{boxplot.maccest}}
}

\examples{
# Iris data
data(iris)
x      <- subset(iris, select = -Species)
y      <- iris$Species

method <- c("randomForest","svm","pcalda","knn")
pars   <- valipars(sampling="boot", niter = 10, nreps=4)
res    <- maccest(x, y, method=method, pars=pars,
                  comp="anova",kernel="linear") 

res
plot(res)
}

\keyword{plot}
