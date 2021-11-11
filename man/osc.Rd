% wll-05-06-2007: 
%
\name{osc}
\alias{osc}
\alias{osc.default}
\alias{osc.formula}
\alias{print.osc}
\alias{summary.osc}
\alias{print.summary.osc}

\title{
  Orthogonal Signal Correction (OSC)
}
\description{
  Data pre-processing by orthogonal signal correction (OSC).
}

\usage{
osc(x, \dots)

\method{osc}{default}(x, y, method="wold",center=TRUE,osc.ncomp=4,pls.ncomp=10,
   tol=1e-3, iter=20,\dots)

\method{osc}{formula}(formula, data = NULL, \dots, subset, na.action = na.omit)
}

% ----------------------------------------------------------------------------
\arguments{
  \item{formula}{
    A formula of the form \code{groups ~ x1 + x2 + \dots}  That is, the
    response is the grouping factor and the right hand side specifies
    the (non-factor) discriminators.
  }
  \item{data}{
    Data frame from which variables specified in \code{formula} are
    preferentially to be taken.
  }
  \item{x}{
    A matrix or data frame containing the explanatory variables if no formula is
    given as the principal argument.
  }
  \item{y}{
    A factor specifying the class for each observation if no formula principal 
    argument is given.
  }
  \item{method}{
  A method for calculating OSC weights, loadings and scores. The following methods are supported:
    \itemize{
      \item \code{wold:} Original Wold et al approach. 
      \item \code{sjoblom:} Sjoblom et al approach.
      \item \code{wise:} Wise and Gallagher approach. 
    }
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
  \item{subset}{
    An index vector specifying the cases to be used in the training
    sample.  
  }
  \item{na.action}{
    A function to specify the action to be taken if \code{NA}s are found. The 
    default action is \code{na.omit}, which leads to rejection of cases with 
    missing values on any required variable. An alternative is \code{na.fail}, 
    which causes an error if \code{NA} cases are found. 
  }	
}

% ----------------------------------------------------------------------------
\value{
  An object of class \code{osc} containing the following components:
  \item{x}{
    A matrix of OSC corrected data set.
  }
  \item{R2}{
  R2 statistics. It is calculated as the fraction of variation in X after OSC 
  correction for the calibration (training) data.
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
  \item{call}{
    The (matched) function call.
  }
  \item{center}{
  A logical value indicating whether the data set has been centred by column-wise.
  }

  \item{osc.ncomp}{
  The number of component used in OSC computation.
  }
  \item{pls.ncomp}{
  The number of component used in PLS computation.
  }
  \item{method}{
  The OSC algorithm used.
  }
}

% ----------------------------------------------------------------------------
\note{
  This function may be called giving either a formula and
  optional data frame, or a matrix and grouping factor as the first
  two arguments. 
}

% ----------------------------------------------------------------------------
\references{
Wold, S., Antti, H., Lindgren, F., Ohman, J.(1998). Orthogonal signal correction of near infrared
spectra. \emph{Chemometrics Intell. Lab. Syst.}, 44: 175-185.

Westerhuis, J. A., de Jong, S., Smilde, A, K. (2001). Direct orthogonal signal correction. 
\emph{Chemometrics Intell. Lab. Syst.}, 56: 13-25.

Sjoblom. J., Svensson, O., Josefson, M., Kullberg, H., Wold, S. (1998). An evaluation of
orthogonal signal correction applied to calibration transfer of near infrared
spectra. \emph{Chemometrics Intell. Lab. Syst.},44: 229-244.

Svensson, O., Kourti, T. and MacGregor, J.F. (2002). An investigation of orthogonal 
correction algorithms and their characteristics. \emph{Journal of Chemometrics}, 16:176-188.

Wise, B. M. and Gallagher, N.B. \emph{http://www.eigenvector.com/MATLAB/OSC.html}.
}


% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------------
\seealso{
  \code{\link{predict.osc}}, \code{\link{osc_wold}}, \code{\link{osc_sjoblom}},
  \code{\link{osc_wise}}
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
res <- osc(train.dat, train.t, method="wise", osc.ncomp=2, pls.ncomp=4)
names(res)
res
summary(res)

## pre-process test data by OSC
test.dat.1 <- predict(res,test.dat)$x

}

\keyword{manip}
