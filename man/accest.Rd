% lwc-14-09-2006: First draft
% lwc-09-10-2006: Modify the method description and note.
%
\name{accest}
\alias{accest}
\alias{aam.cl}
\alias{aam.mcl}
\alias{accest.formula}
\alias{accest.default}
\alias{print.accest}
\alias{summary.accest}
\alias{print.summary.accest}

\title{Estimate Classification Accuracy By Resampling Method}

% ------------------------------------------------------------------------
\description{
  Estimate classification accuracy rate by resampling method. 
}

% ------------------------------------------------------------------------
\usage{
accest(dat, \dots)

\method{accest}{default}(dat, cl, method, pred.func=predict,pars = valipars(), 
       tr.idx = NULL, \dots) 

\method{accest}{formula}(formula, data = NULL, \dots, subset, na.action = na.omit)

aam.cl(x,y,method, pars = valipars(),\dots)

aam.mcl(x,y,method, pars = valipars(),\dots)

}


% ------------------------------------------------------------------------
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
  \item{dat,x}{
    A matrix or data frame containing the explanatory variables if no
    formula is given as the principal argument.
  }
  \item{cl,y}{
    A factor specifying the class for each observation if no formula
    principal argument is given.
  }
  \item{method}{
    Classification method whose accuracy rate is to be estimated, such as 
    \code{randomForest}, \code{svm}, \code{knn} and \code{lda}. 
    For details, see \code{note} below. Either a function or a character 
    string naming the function to be called.
 
  }   
  \item{pred.func}{
    Predict method (default is \code{predict}). Either a function or a
    character string naming the function to be called.
  }
  \item{pars}{
    A list of parameters using by the resampling method such as 
    \emph{Leave-one-out cross-validation}, \emph{Cross-validation}, 
    \emph{Bootstrap} and \emph{Randomised validation (holdout)}.
  	See \code{\link{valipars}} for details.
  }
  \item{tr.idx}{
    User defined index of training samples. Can be generated by \code{trainind}.
  }
  \item{\dots}{Additional parameters to \code{method}.}
  \item{subset}{
    Optional vector, specifying a subset of observations to be used.
  }
  \item{na.action}{
    Function which indicates what should happen when the data
    contains \code{NA}'s, defaults to \code{\link{na.omit}}.
  } 
}

% ----------------------------------------------------------------------------
\details{
  The accuracy rates of classification are estimated by techniques such
  as \emph{Random Forest}, \emph{Support Vector Machine},
  \emph{k-Nearest Neighbour Classification} and \emph{Linear
  Discriminant Analysis} based on resampling methods, including
  \emph{Leave-one-out cross-validation}, \emph{Cross-validation},
  \emph{Bootstrap} and \emph{Randomised validation (holdout)}.
}
% ----------------------------------------------------------------------------
\note{
  The \code{accest} can take any classification models if their argument
  format is \code{model(formula, data, subset, na.action, ...)} and
  their corresponding method \code{predict.model(object, newdata, ...)}
  can either return the only predicted class label or a list with a
  component called \code{class}, such as \code{lda} and \code{pcalda}.
  
  If classifier \code{method} provides posterior probabilities, the
  prediction margin \code{mar} will be generated, otherwise \code{NULL}.
  
  If classifier \code{method} provides posterior probabilities and the
  classification is for two-class problem, \code{auc} will be generated,
  otherwise \code{NULL}.
  
  \code{aam.cl} is a wrapper function of \code{accest}, returning
  accuracy rate, AUC and classification margin. \code{aam.mcl} accepts
  multiple classifiers in one running.
}

% ----------------------------------------------------------------------------
\value{
  \code{accest} returns an object including the components:
  \item{method}{Classification method used.}
  \item{acc}{Overall accuracy rate.}
  \item{acc.iter}{Average accuracy rate for each iteration.}
  \item{acc.all}{Accuracy rate for each iteration and replication.}

  \item{auc}{Overall area under receiver operating curve (AUC).}
  \item{auc.iter}{Average AUC for each iteration.}
  \item{auc.all}{AUC for each iteration and replication.}
  
  \item{mar}{Overall prediction margin.}
  \item{mar.iter}{Average prediction margin for each iteration.}
  \item{mar.all}{Prediction margin for each iteration and replication.}

  \item{err}{Overall error rate.}
  \item{err.iter}{Average error rate for each iteration.}
  \item{err.all}{Error rate for each iteration and replication.}

  \item{sampling}{Sampling scheme used.} 
  \item{niter}{Number of iteration.}
  \item{nreps}{Number of replications in each iteration if resampling is 
     not \code{loocv}. }
  \item{conf}{Overall confusion matrix.}
  \item{res.all}{All results which can be further processed.}      
  \item{acc.boot}{ A list of bootstrap accuracy such as \code{.632} and
    \code{.632+} if the resampling method is bootstrap.
  } 
  

	\code{aam.cl} returns a vector with \code{acc} (accuracy), 
	\code{auc}(area under ROC curve) and \code{mar}(class margin).


	\code{aam.mcl} returns a matrix with columns of \code{acc} (accuracy), 
	\code{auc}(area under ROC curve) and \code{mar}(class margin).

}

% ------------------------------------------------------------------------
\author{
  Wanchang Lin
}

\seealso{
  \code{\link{binest}}, \code{\link{maccest}}, \code{\link{valipars}}, 
  \code{\link{trainind}}, \code{\link{classifier}}
}

% ------------------------------------------------------------------------
\examples{
# Iris data
data(iris)
# Use KNN classifier and bootstrap for resampling
acc <- accest(Species~., data = iris, method = "knn",
              pars = valipars(sampling = "boot",niter = 2, nreps=5))
acc
summary(acc)
acc$acc.boot

# alternatively the traditional interface:
x <- subset(iris, select = -Species)
y <- iris$Species

## -----------------------------------------------------------------------
# Random Forest with 5-fold stratified cv 
pars   <- valipars(sampling = "cv",niter = 4, nreps=5, strat=TRUE)
tr.idx <- trainind(y,pars=pars)
acc1   <- accest(x, y, method = "randomForest", pars = pars, tr.idx=tr.idx)
acc1
summary(acc1)
# plot the accuracy in each iteration
plot(acc1)

## -----------------------------------------------------------------------
# Forensic Glass data in chap.12 of MASS
data(fgl, package = "MASS")    # in MASS package
# Randomised validation (holdout) of SVM for fgl data
acc2 <- accest(type~., data = fgl, method = "svm", cost = 100, gamma = 1, 
              pars = valipars(sampling = "rand",niter = 10, nreps=4,div = 2/3) )
              
acc2
summary(acc2)
# plot the accuracy in each iteration
plot(acc2)

## -----------------------------------------------------------------------
## Examples of amm.cl and aam.mcl
aam.1 <- aam.cl(x,y,method="svm",pars=pars)
aam.2 <- aam.mcl(x,y,method=c("svm","randomForest"),pars=pars)

## If use two classes, AUC will be calculated
idx <- (y == "setosa")
aam.3 <- aam.cl(x[!idx,],factor(y[!idx]),method="svm",pars=pars)
aam.4 <- aam.mcl(x[!idx,],factor(y[!idx]),method=c("svm","randomForest"),pars=pars)

}

\keyword{classif}
