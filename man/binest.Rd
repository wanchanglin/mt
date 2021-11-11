% lwc-11-10-2006
% lwc-31-10-2006: minor changes

\name{binest}
\alias{binest}
\title{
  Binary Classification
}
\description{
  Binary classification.
}

\usage{
  binest(dat, cl, choices = NULL, method, pars=valipars(),\dots) 
}

% --------------------------------------------------------------------
\arguments{
  \item{dat}{
    A matrix or data frame containing the explanatory variables.
  }
  \item{cl}{
    A factor specifying the class for each observation.
  }
  \item{choices}{
    The vector or list of class labels to be chosen for binary
    classification. For details, see \code{\link{dat.sel}}.
  }
  \item{method}{
    Classification method to be used. For details, see
    \code{\link{accest}}.
  }             
  \item{pars}{
    A list of parameters of the resampling method. For details, see
    \code{\link{valipars}}.
  }
  \item{\dots}{Additional parameters to \code{method}.}
}

% ----------------------------------------------------------------------
\value{
  A list with components:
  \item{com}{A matrix of combination of the binary class labels.}
  \item{acc}{A table of classification accuracy for the binary
    combination in each iteration.
  }
  \item{method}{Classification method used.}
  \item{sampling}{Sampling scheme used.} 
  \item{niter}{Number of iterations.}
  \item{nreps}{Number of replications in each iteration if resampling is
     not \code{loocv}. }
}

% -----------------------------------------------------------------------
\author{
  Wanchang Lin
}

\seealso{
  \code{\link{accest}}, \code{\link{valipars}}, \code{\link{dat.sel}}
}

% ----------------------------------------------------------------------
\examples{
# iris data set
data(iris)
dat <- subset(iris, select = -Species)
cl  <- iris$Species

## PCALDA with cross-validation
pars    <- valipars(sampling="cv",niter = 6, nreps = 5)
binpcalda <- binest(dat,cl,choices=c("setosa"), method="pcalda", pars = pars)

## SVM with leave-one-out cross-validation. SVM kernel is 'linear'.
pars   <- valipars(sampling="loocv")
binsvm <- binest(dat,cl,choices=c("setosa","virginica"), method="svm",
                 pars = pars, kernel="linear")

## randomForest with bootstrap
pars  <- valipars(sampling="boot",niter = 5, nreps = 5)
binrf <- binest(dat,cl,choices=c("setosa","virginica"), 
                method="randomForest", pars = pars)

## KNN with randomised validation. The number of neighbours is 3.
pars   <- valipars(sampling="rand",niter = 5, nreps = 5)
binknn <- binest(dat,cl,choices = list(c("setosa","virginica"),
                                       c("virginica","versicolor")), 
                 method="knn",pars = pars, k = 3)

}

\keyword{classif}
