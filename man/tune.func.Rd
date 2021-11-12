% wll-03-10-2007: add plsc stuff 
%

\name{tune.func}
\alias{tune.func}
\alias{tune.plsc}
\alias{tune.plslda}
\alias{tune.pcalda}

\title{Functions for Tuning Appropriate Number of Components}
\description{
  Tune appropriate number of components (\code{ncomp}) for \code{plsc}, 
  \code{plslda} or \code{pcalda}.
}

\usage{
tune.plsc(x,y, pls="simpls",ncomp=10, tune.pars,\dots)

tune.plslda(x,y, pls="simpls",ncomp=10, tune.pars,\dots)

tune.pcalda(x,y, ncomp=NULL, tune.pars,\dots)
}

\arguments{
  \item{x}{
    A matrix or data frame containing the explanatory variables if no formula is
    given as the principal argument.
  }
  \item{y}{
    A factor specifying the class for each observation if no formula principal 
    argument is given.
  }
  \item{pls}{
  A method for calculating PLS scores and loadings. The following methods are supported:
    \itemize{
      \item \code{simpls:} SIMPLS algorithm. 
      \item \code{kernelpls:} kernel algorithm.
      \item \code{oscorespls:} orthogonal scores algorithm. 
    }
  For details, see \code{\link[pls]{simpls.fit}}, \code{\link[pls]{kernelpls.fit}} and
  \code{\link[pls]{oscorespls.fit}} in package \pkg{pls}.
  }
  \item{ncomp}{
    The number of components to be used in the classification.  
  }
  \item{tune.pars}{
    A list of parameters using by the resampling method. 
  	See \code{\link{valipars}} for details.
  }
  \item{\dots}{Further parameters passed to \code{tune}.}
}

\value{
  A list including: 
  \item{ncomp}{The best number of components.}
  \item{acc.tune}{Accuracy rate of components.}
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

\seealso{
  \code{\link{plsc}}, \code{\link{plslda}}, \code{\link{pcalda}},\code{\link{valipars}}
}

% ----------------------------------------------------------------------------
\examples{
\dontrun{
data(abr1)
cl   <- factor(abr1$fact$class)
dat  <- preproc(abr1$pos , y=cl, method=c("log10"),add=1)[,110:500]

## divide data as training and test data
idx <- sample(1:nrow(dat), round((2/3)*nrow(dat)), replace=FALSE) 

## construct train and test data 
train.dat  <- dat[idx,]
train.t    <- cl[idx]
test.dat   <- dat[-idx,]        
test.t     <- cl[-idx] 

## tune the best number of components
ncomp.plsc   <- tune.plsc(dat,cl, pls="simpls",ncomp=20)
ncomp.plslda <- tune.plslda(dat,cl, pls="simpls",ncomp=20)
ncomp.pcalda <- tune.pcalda(dat,cl, ncomp=60)

## model fit
(z.plsc   <- plsc(train.dat,train.t, ncomp=ncomp.plsc$ncomp))
(z.plslda <- plslda(train.dat,train.t, ncomp=ncomp.plslda$ncomp))
(z.pcalda <- pcalda(train.dat,train.t, ncomp=ncomp.pcalda$ncomp))

## or indirect use tune function in model fit
z.plsc   <- plsc(train.dat,train.t, ncomp=20, tune=TRUE)
z.plslda <- plslda(train.dat,train.t, ncomp=20, tune=TRUE)
z.pcalda <- pcalda(train.dat,train.t, ncomp=60, tune=TRUE)

## predict test data
pred.plsc   <- predict(z.plsc, test.dat)$class
pred.plslda <- predict(z.plslda, test.dat)$class
pred.pcalda <- predict(z.pcalda, test.dat)$class

## classification rate and confusion matrix
cl.rate(test.t, pred.plsc)
cl.rate(test.t, pred.plslda)
cl.rate(test.t, pred.pcalda)

}
}

\keyword{models}
