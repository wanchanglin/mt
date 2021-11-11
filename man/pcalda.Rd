% wll-02-06-2007
%
\name{pcalda}
\alias{pcalda}
\alias{pcalda.default}
\alias{pcalda.formula}
\alias{print.pcalda}
\alias{summary.pcalda}
\alias{print.summary.pcalda}

\title{
  Classification with PCADA
}
\description{
  Classification with combination of principal component analysis (PCA) and linear discriminant 
  analysis (LDA).
}

\usage{
pcalda(x, \dots)

\method{pcalda}{default}(x, y, center = TRUE, scale. = FALSE, ncomp = NULL,
       tune=FALSE,\dots)

\method{pcalda}{formula}(formula, data = NULL, \dots, subset, na.action = na.omit)
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
  \item{center}{
    A logical value indicating whether \code{x} should be shifted to zero 
    centred by column-wise.
  }
  \item{scale.}{
    A logical value indicating whether \code{x} should be scaled to have unit 
    variance by column-wise before the analysis takes place. 
  }
  \item{ncomp}{
    The number of principal components to be used in the classification. If
    \code{NULL} and \code{tune=TRUE}, it is the row number of \code{x} minus the 
    number of class indicating in \code{y}. If \code{NULL} and \code{tune=FALSE}, 
    it is the half of row number of \code{x}.
    
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
A critical issue of applying linear discriminant analysis (LDA) is both the
singularity and instability of the within-class scatter matrix. In practice, 
there are often a large number of features available, but the total number of 
training patterns is limited and commonly less than the dimension of the feature 
space. To tackle this issue, \code{pcalda} combines PCA and LDA for 
classification. It uses PCA for dimension reduction. The rotated data resulted 
from PCA will be the input variable to LDA for classification. 
}

% ----------------------------------------------------------------------------
\value{
  An object of class \code{pcalda} containing the following components:
  \item{x}{
    The rotated data on discriminant variables.
  }
  \item{cl}{
  The observed class labels of training data. 
  }
  \item{pred}{
   The predicted class labels of training data. 
  }
  \item{posterior}{
  The posterior probabilities for the predicted classes. 
  }
  \item{conf}{
  The confusion matrix based on training data. 
  }
  \item{acc}{
  The accuracy rate of training data. 
  }
  \item{ncomp}{
  The number of principal components used for classification. 
  }
  \item{pca.out}{
  The output of PCA.
  }
  \item{lda.out}{
  The output of LDA.
  }
  \item{call}{
    The (matched) function call.
  }
}

% ----------------------------------------------------------------------------
\note{
  This function may be called giving either a formula and
  optional data frame, or a matrix and grouping factor as the first
  two arguments. 
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------------
\seealso{
  \code{\link{predict.pcalda}}, \code{\link{plot.pcalda}}, \code{\link{tune.func}}
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

## apply pcalda
model    <- pcalda(train.dat,train.t)
model
summary(model)

## plot
plot(model,dimen=c(1,2),main = "Training data",abbrev = TRUE)
plot(model,main = "Training data",abbrev = TRUE)

## confusion matrix
pred.te  <- predict(model, test.dat)$class
table(test.t,pred.te)

}

\keyword{classif}
