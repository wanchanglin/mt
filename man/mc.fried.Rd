% lwc-30-01-2007:

\name{mc.fried}
\alias{mc.fried}
\title{
  Multiple Comparison by 'Friedman Test' and Pairwise Comparison by 'Wilcoxon Test'
}
\description{
  Performs multiple comparison by \code{Friedman test} and pairwise comparison by 
  \code{Wilcoxon Test}.
}

\usage{
mc.fried(x, p.adjust.method = p.adjust.methods,\dots)
}

% --------------------------------------------------------------------
\arguments{
  \item{x}{
    A matrix or data frame to be tested. 
  }
  \item{p.adjust.method}{Method for adjusting p values (see \code{\link{p.adjust}}).}  
  \item{\dots}{Additional arguments pass to \code{friedman.test} or \code{pairwise.wilcox.test}.}
}

% ----------------------------------------------------------------------
\value{
  A list with components:
  \item{fried}{Hypothesis test results of \code{friedman.test}.}
  \item{wilcox}{Hypothesis test results of \code{pairwise.wilcox.test}.}
  \item{gl.pval}{Global or overall p value returned by \code{friedman.test}.}
  \item{mc.pval}{Pairwise p value returned by \code{pairwise.wilcox.test}.}
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
pars   <- valipars(sampling="cv", niter = 10, nreps=4)
res    <- maccest(x, y, method=method, pars=pars,
                  comp="fried",kernel="linear") 

res

htest <- mc.fried(res$acc.iter)

}

\keyword{classif}
