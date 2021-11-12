% lwc-30-01-2007: First draft
%
\name{boxplot.maccest}
\alias{boxplot.maccest}
\title{
  Boxplot Method for Class 'maccest'
}
\description{
  Boxplot method for the accuracy rate of each classifier.
}

\usage{
\method{boxplot}{maccest}(x,  \dots)
}

% ----------------------------------------------------------------------------
\arguments{
  \item{x}{
  An object of class \code{maccest}.
  }
  \item{\dots}{
  Additional arguments to the plot, such as \code{main}, \code{xlab} and 
  \code{ylab}. 
  }
}

% ----------------------------------------------------------------------------
\details{
This function is a method for the generic function \code{boxplot()} for class 
\code{maccest}. It plots the accurary rate for each classifier.
}

\value{
  Returns boxplot of class \code{maccest}.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

\seealso{
  \code{\link{maccest}}, \code{\link{plot.maccest}}
}

\examples{
# Iris data
data(iris)
x      <- subset(iris, select = -Species)
y      <- iris$Species

method <- c("randomForest","svm","knn")
pars   <- valipars(sampling="cv", niter = 2, nreps=5)
tr.idx <- trainind(y, pars=pars)
res    <- maccest(x, y, method=method, pars=pars,
                  comp="anova",kernel="linear") 

res
boxplot(res)
}

\keyword{plot}
