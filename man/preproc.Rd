% lwc-01-12-2006:
% lwc-22-02-2010: Add other functions

\name{preproc}
\alias{preproc}
\alias{preproc.sd}
\alias{preproc.const}

\title{
  Pre-process Data Set
}
\description{
  Pre-process a data frame or matrix by different methods.
}
\usage{
preproc (x, y=NULL,method="log",add=1)

preproc.sd(x, y=NULL, na.rm = FALSE)

preproc.const(x, y, tol = 1.0e-4)
}
\arguments{
  \item{x}{
    A numeric data frame or matrix to be pre-processed. 
  }
  \item{y}{
    A factor specifying the group. It is only used by the 
    method \code{TICnorm} in \code{preproc}.
  }
  \item{method}{
    A method used to pre-process the data set. The following methods are supported:
    \itemize{
      \item \code{center:} Centering
      \item \code{auto:} Auto scaling
      \item \code{range:} Range scaling
      \item \code{pareto:} Pareto scaling
      \item \code{vast:} Vast scaling
      \item \code{level:} Level scaling
      \item \code{log:} Log transformation (default)
      \item \code{log10:} Log 10 transformation
      \item \code{sqrt:} Square root transformation
      \item \code{asinh:} Inverse hyperbolic sine transformation
      \item \code{TICnorm:} TIC normalisation
    }
  }
  \item{na.rm}{A logical value indicating whether NA values should be stripped 
       before the computation proceeds.
  }
  \item{add}{
  A numeric value for addition used in the logarmath transformation \code{log}
  and \code{log10}.
  }
  \item{tol}{	A tolerance to decide if a matrix is singular; it will reject 
   variables and linear combinations of unit-variance variables whose variance 
   is less than tol^2.
   } 
}

\value{
  A pre-processed data set.
}

\details{
  \code{preproc} transforms data set by provided \code{method}.  

  \code{preproc.sd} removes variables which have (near) zero S.D with or without 
   respect to class/grouped information.

  \code{preproc.const} removes variables appears to be constant within groups / classes.
}

\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------
\references{
Berg, R., Hoefsloot, H., Westerhuis, J., Smilde, A. and Werf, M. (2006),
Centering, scaling, and transformations: improving the biological
information content of metabolomics data, \emph{BMC Genomics}, 7:142
}

\examples{
data(abr1)
cl   <- factor(abr1$fact$class)
dat  <- abr1$pos

## normalise data set using "TICnorm"
z.1  <- preproc(dat, y=cl, method="TICnorm")

## scale data set using "log10"
z.2 <- preproc(dat,method="log10", add=1)

## or scale data set using "log" and "TICnorm" sequentially
z.3 <- preproc(dat,method=c("log","TICnorm"), add=0.1)

}

\keyword{manip}
