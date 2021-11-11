% wll-30-01-2007:commence
% wll-03-07-2007:add auc and margin

\name{mbinest}
\alias{mbinest}
\title{
  Binary Classification by Multiple Classifier
}
\description{
  Binary classification by multiple classifier.
}

\usage{
mbinest(dat, cl, choices = NULL, method, pars=valipars(),\dots) 
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
  The vector or list of class labels to be chosen for binary classification.
  For details, see \code{\link{dat.sel}}.
  }
  \item{method}{
    Multiple classification methods to be used. For details, see \code{\link{maccest}}.
  }             
  \item{pars}{
    A list of parameters of the resampling method. See \code{\link{valipars}} for details.
  }
  \item{\dots}{Additional parameters to \code{method}.}
}

% ----------------------------------------------------------------------
\value{
  A list with components:
  \item{all}{All results of classification.}
  \item{com}{A matrix of the combinations of the binary class labels.}
  \item{acc}{A table of classification accuracy for the binary combination.}
  \item{mar}{Prediction margin.}
  \item{auc}{The area under receiver operating curve (AUC).}
  \item{method}{Classification methods used.}
  \item{niter}{Number of iterations.}
  \item{sampling}{Sampling scheme used.} 
  \item{nreps}{Number of replications in each iteration if sampling is not \code{loocv}. }
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

\seealso{
  \code{\link{maccest}}, \code{\link{maccest}},\code{\link{valipars}}, \code{\link{dat.sel}}
}

% ----------------------------------------------------------------------
\examples{
## iris data set
data(iris)
dat    <- subset(iris, select = -Species)
cl     <- iris$Species
method <- c("svm","pcalda")

pars  <- valipars(sampling="cv",niter = 10, nreps = 5)
res   <- mbinest(dat,cl,choices=c("setosa"), method=method,
                  pars = pars, kernel="linear")

## combine prediction accuracy, AUC and margin 
z      <- round(cbind(res$acc,res$auc,res$mar),digits=3)
colnames(z) <- c(paste(method,".acc", sep=""),paste(method,".auc", sep=""),
                 paste(method,".mar", sep=""))

}

\keyword{classif}
