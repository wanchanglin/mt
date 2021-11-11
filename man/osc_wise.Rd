% wll-05-06-2007: 
%
\name{osc_wise}
\alias{osc_wise}

\title{
  Orthogonal Signal Correction (OSC) Approach by Wise and Gallagher.
}
\description{
  Orthogonal signal correction (OSC) approach by Wise and Gallagher.
}

\usage{
  osc_wise(x, y, center=TRUE,osc.ncomp=4,pls.ncomp=10,
           tol=1e-3,iter=20,\dots)
}

% ----------------------------------------------------------------------------
\arguments{
  \item{x}{
    A numeric data frame or matrix to be pre-processed. 
  }
  \item{y}{
    A vector or factor specifying the class for each observation.
  }
  \item{center}{
  A logical value indicating whether the data set should be centred by column-wise.
  }
  \item{osc.ncomp}{
    The number of components to be used in the OSC calculation.  
  }
  \item{pls.ncomp}{
    The number of components to be used in the PLS calculation.  
  }
  \item{tol}{
  A scalar value of tolerance for OSC computation.
  }
  \item{iter}{
  The number of iteration used in OSC calculation.
  }
  \item{\dots}{
    Arguments passed to or from other methods.
  }
}

% ----------------------------------------------------------------------------
\value{
  A list containing the following components:
  \item{x}{
    A matrix of OSC corrected data set.
  }
  \item{R2}{
  R2 statistics. It is calculated as the fraction of variation in X after OSC correction.
  }
  \item{angle}{
  An angle used for checking if scores \code{t} is orthogonal to \code{y}. An 
  angle close to 90 degree means that orthogonality is achieved in the correction process.
  }
  \item{w}{
  A matrix of OSC weights.
  }
  \item{p}{
   A matrix of OSC loadings.
  }
  \item{t}{
  A matrix of OSC scores.
  }
  \item{center}{
  A logical value indicating whether the data set has been centred by column-wise.
  }
}

% ----------------------------------------------------------------------------
\references{
Westerhuis, J. A., de Jong, S., Smilde, A, K. (2001). Direct orthogonal signal correction. 
\emph{Chemometrics Intell. Lab. Syst.}, 56: 13-25.

Wise, B. M. and Gallagher, N.B. \emph{http://www.eigenvector.com/MATLAB/OSC.html}.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------------
\seealso{
  \code{\link{osc}}, \code{\link{predict.osc}}, \code{\link{osc_sjoblom}},
  \code{\link{osc_wold}}
}

% ----------------------------------------------------------------------------
\examples{
data(abr1)
cl   <- factor(abr1$fact$class)
dat  <- abr1$pos

## divide data as training and test data
idx <- sample(1:nrow(dat), round((2/3)*nrow(dat)), replace=FALSE) 

## construct train and test data 
train.dat  <- dat[idx,]
train.t    <- cl[idx]
test.dat   <- dat[-idx,]        
test.t     <- cl[-idx] 

## build OSC model based on the training data
res <- osc_wise(train.dat, train.t)
names(res)

## pre-process test data by OSC
test.dat.1 <- predict.osc(res,test.dat)$x

}

\keyword{manip}
