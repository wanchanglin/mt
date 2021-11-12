% lwc-30-01-2007:

\name{mc.anova}
\alias{mc.anova}
\title{
  Multiple Comparison by 'ANOVA' and Pairwise Comparison by 'HSDTukey Test'
}
\description{
  Performs multiple comparison by \code{ANOVA} and pairwise comparison by 
  \code{HSDTukey Test}.
}

\usage{
mc.anova(x, \dots)
}

% --------------------------------------------------------------------
\arguments{
  \item{x}{
    A matrix or data frame to be tested. 
  }
  \item{\dots}{Additional arguments pass to \code{anova} or \code{HSDTukey test}.}
}

% ----------------------------------------------------------------------
\value{
  A list with components:
  \item{anova}{Hypothesis test results of \code{anova}.}
  \item{tukey}{Hypothesis test results of \code{HSDTukey.test}.}
  \item{gl.pval}{Global or overall p value returned by \code{anova}.}
  \item{mc.pval}{Pairwise p value returned by \code{HSDTukey.test}.}
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

\seealso{
  \code{\link{maccest}}, \code{\link{mc.fried}}
}

% ----------------------------------------------------------------------
\examples{

# Iris data
data(iris)
x      <- subset(iris, select = -Species)
y      <- iris$Species

method <- c("randomForest","svm","pcalda","knn")
pars   <- valipars(sampling="boot", niter = 10, nreps=4)
res    <- maccest(x, y, method=method, pars=pars, comp="anova") 

res
htest <- mc.anova(res$acc.iter)

oldpar <- par(mar = c(5,10,4,2) + 0.1)
plot(htest$tukey,las=1)   ## plot the tukey results
par(oldpar)
}

\keyword{classif}
