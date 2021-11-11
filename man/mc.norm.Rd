% lwc-30-01-2007:

\name{mc.norm}
\alias{mc.norm}
\title{
  Normality Test by Shapiro-Wilk Test
}
\description{
  Perform Shapiro-Wilk normality test by \code{shapiro.test} and plot the 
  density function and boxplot.
}

\usage{
mc.norm(x, \dots)
}

% --------------------------------------------------------------------
\arguments{
  \item{x}{
    A matrix or data frame to be tested. 
  }
  \item{\dots}{Additional arguments pass to \code{shapiro.test}.}
}

% ----------------------------------------------------------------------
\value{
  Object of \code{shapiro.test}.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

\seealso{
  \code{\link{maccest}}, \code{\link{mc.anova}}
}

% ----------------------------------------------------------------------
\examples{

# Iris data
data(iris)
x      <- subset(iris, select = -Species)
y      <- iris$Species

method <- c("randomForest","svm","pcalda","knn")
pars   <- valipars(sampling="boot", niter = 10, nreps=10)
res    <- maccest(x, y, method=method, pars=pars,
                  comp="anova",kernel="linear") 
res
norm.t <- mc.norm(res$acc.iter)
}

\keyword{classif}
