% wll-23-05-2007: commence
% wll-03-10-2007: add plsc stuff
%
\name{plsc}
\alias{plsc}
\alias{plsc.default}
\alias{plsc.formula}
\alias{print.plsc}
\alias{summary.plsc}
\alias{print.summary.plsc}
\alias{plslda}
\alias{plslda.default}
\alias{plslda.formula}
\alias{print.plslda}
\alias{summary.plslda}
\alias{print.summary.plslda}

\title{
  Classification with PLSDA
}
\description{
  Classification with partial least squares (PLS) or PLS plus linear discriminant 
  analysis (LDA).
}

\usage{
plsc(x, \dots)

plslda(x, \dots)

\method{plsc}{default}(x, y, pls="simpls",ncomp=10, tune=FALSE,\dots)

\method{plsc}{formula}(formula, data = NULL, \dots, subset, na.action = na.omit)

\method{plslda}{default}(x, y, pls="simpls",ncomp=10, tune=FALSE,\dots)

\method{plslda}{formula}(formula, data = NULL, \dots, subset, na.action = na.omit)
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
  \item{tune}{
    A logical value indicating whether the best number of components should be tuned.
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
\details{
\code{plcs} implements PLS for classification. In details, the categorical response 
vector \code{y} is converted into a numeric matrix for regression by PLS and the 
output of PLS is convert to posteriors by \code{softmax} method. 
The classification results are obtained based on the posteriors. \code{plslda} 
combines PLS and LDA for classification, in which, PLS is for dimension reduction
and LDA is for classification based on the data transformed by PLS. 

Three PLS functions,\code{\link[pls]{simpls.fit}}, 
\code{\link[pls]{kernelpls.fit}} and  \code{\link[pls]{oscorespls.fit}}, are  
implemented in package \pkg{pls}.
}

% ----------------------------------------------------------------------------
\value{
  An object of class \code{plsc} or \code{plslda} containing the following components:
  \item{x}{
    A matrix of the latent components or scores from PLS.
  }
  \item{cl}{
  The observed class labels of training data. 
  }
  \item{pred}{
   The predicted class labels of training data. 
  }
  \item{conf}{
  The confusion matrix based on training data. 
  }
  \item{acc}{
  The accuracy rate of training data. 
  }
  \item{posterior}{
  The posterior probabilities for the predicted classes. 
  }
  \item{ncomp}{
  The number of latent component used for classification. 
  }
  \item{pls.method}{
  The PLS algorithm used.
  }
  \item{pls.out}{
  The output of PLS.
  }
  \item{lda.out}{
  The output of LDA used only by \code{plslda}.
  }
  \item{call}{
    The (matched) function call.
  }
}

% ----------------------------------------------------------------------------
\note{
  Two functions may be called giving either a formula and
  optional data frame, or a matrix and grouping factor as the first
  two arguments. 
}

% ----------------------------------------------------------------------------
\references{
  Martens, H. and Nas, T. (1989) \emph{Multivariate calibration.}
  John Wiley & Sons.
}


% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------------
\seealso{
  \code{\link[pls]{kernelpls.fit}}, \code{\link[pls]{simpls.fit}}, 
  \code{\link[pls]{oscorespls.fit}}, \code{\link{predict.plsc}},
  \code{\link{plot.plsc}}, \code{\link{tune.func}}
}

% ----------------------------------------------------------------------------
\examples{
library(pls)  
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

## apply plsc and plslda
(res   <- plsc(train.dat,train.t, ncomp = 20, tune = FALSE))
## Estimate the mean squared error of prediction (MSEP), root mean squared error
## of prediction (RMSEP) and R^2 (coefficient of multiple determination) for 
## fitted PLSR model 
MSEP(res$pls.out)
RMSEP(res$pls.out)
R2(res$pls.out)

(res.1  <- plslda(train.dat,train.t, ncomp = 20, tune = FALSE))
## Estimate the mean squared error of prediction (MSEP), root mean squared error
## of prediction (RMSEP) and R^2 (coefficient of multiple determination) for 
## fitted PLSR model 
MSEP(res.1$pls.out)
RMSEP(res.1$pls.out)
R2(res.1$pls.out)

\dontrun{
## with function of tuning component numbers
(z.plsc   <- plsc(train.dat,train.t, ncomp = 20, tune = TRUE))
(z.plslda <- plslda(train.dat,train.t, ncomp = 20, tune = TRUE))

## check nomp tuning results
z.plsc$ncomp
plot(z.plsc$acc.tune)
z.plslda$ncomp
plot(z.plslda$acc.tune)

## plot
plot(z.plsc,dimen=c(1,2,3),main = "Training data",abbrev = TRUE)
plot(z.plslda,dimen=c(1,2,3),main = "Training data",abbrev = TRUE)

## predict test data
pred.plsc   <- predict(z.plsc, test.dat)$class
pred.plslda <- predict(z.plslda, test.dat)$class

## classification rate and confusion matrix
cl.rate(test.t, pred.plsc)
cl.rate(test.t, pred.plslda)

}
}

\keyword{classif}
