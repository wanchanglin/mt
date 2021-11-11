% lwc-16-09-2006
% lwc-05-05-2010:
\name{cl.rate}
\alias{cl.rate}
\alias{cl.perf}
\alias{cl.roc}
\alias{cl.auc}

\title{
  Assess Classification Performances
}
\description{
  Assess classification performances.
}
\usage{
  cl.rate(obs, pre)
  cl.perf(obs, pre, pos=levels(as.factor(obs))[2])
  cl.roc(stat, label, pos=levels(as.factor(label))[2], plot=TRUE, \dots)
  cl.auc(stat, label, pos=levels(as.factor(label))[2])
}
\arguments{
  \item{obs}{
  Factor or vector of observed class.
  }
  \item{pre}{
  Factor or vector of predicted class.
  }
  \item{stat}{
  Factor or vector of statistics for positives/cases. 
  }
  \item{label}{
  Factor or vector of label for categorical data.
  }
  \item{pos}{
  Characteristic string for positive.
  }
  \item{plot}{
   Logical flag indicating whether ROC should be plotted.  
  }
  \item{\dots}{ Further arguments for plotting.  }
}

% ----------------------------------------------------------------------------
\details{
  \code{cl.perf} gets the classification performances such as accuracy rate and 
  false positive rate. \code{cl.roc} computes receiver operating characteristics 
  (ROC).  \code{cl.auc} calculates area under ROC curve. Three functions are only 
  for binary class problems.  
}

\note{
   AUC varies between 0.5 and 1.0 for sensible models; the higher the better. If 
   it is less than 0.5, it should be corrected by \code{1 - AUC}. Or re-run it by 
   using \code{1 - stat}. 
}
% ----------------------------------------------------------------------
\value{
  \code{cl.rate} returns a list with components:
  \item{acc}{ Accuracy rate of classification.}
  \item{err}{ Error rate of classification.}
  \item{con.mat}{ Confusion matrix. }
  \item{kappa}{ Kappa Statistics.}
  
  \code{cl.perf} returns a list with components:
  \item{acc}{ Accuracy rate}
  \item{tpr}{ True positive rate}
  \item{fpr}{ False positive rate}
  \item{sens}{ Sensitivity}
  \item{spec}{ Specificity}
  \item{con.mat}{ Confusion matrix. }
  \item{kappa}{ Kappa Statistics.}
  \item{positive}{ Positive level.}

  \code{cl.roc} returns a list with components:
  \item{perf}{A data frame of \code{acc}, \code{tpr},\code{fpr},\code{sens},
              \code{spec} and \code{cutoff} (thresholds).}
  \item{auc}{ Area under ROC curve}
  \item{positive}{ Positive level.}


  \code{cl.auc} returns a scalar value of AUC.
  
}

\references{
  Fawcett, F. (2006) \emph{An introduction to ROC analysis}. 
  \emph{Pattern Recognition Letters}.  vol. 27, 861-874.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------
\examples{

## Measurements of Forensic Glass Fragments
library(MASS)
data(fgl, package = "MASS")    # in MASS package
dat <- subset(fgl, grepl("WinF|WinNF",type))
## dat <- subset(fgl, type %in% c("WinF", "WinNF"))
x   <- subset(dat, select = -type)
y   <- factor(dat$type)

## construct train and test data 
idx   <- sample(1:nrow(x), round((2/3)*nrow(x)), replace = FALSE) 
tr.x  <- x[idx,]
tr.y  <- y[idx]
te.x  <- x[-idx,]        
te.y  <- y[-idx] 

model <- lda(tr.x, tr.y)

## predict the test data results
pred  <- predict(model, te.x)

## classification performances
obs <- te.y
pre <- pred$class   
cl.rate(obs, pre)
cl.perf(obs, pre, pos="WinNF")
## change positive as "WinF"
cl.perf(obs, pre, pos="WinF")

## ROC and AUC
pos  <- "WinNF"            ## or "WinF"
stat <- pred$posterior[,pos]
## levels(obs) <- c(0,1)

cl.auc (stat,obs, pos=pos)
cl.roc (stat,obs, pos=pos)

## test examples for ROC and AUC
label <- rbinom(30,size=1,prob=0.2)
stat  <- rnorm(30)
cl.roc(stat,label, pos=levels(factor(label))[2],plot = TRUE)
cl.auc(stat,label,pos=levels(factor(label))[2])

## if auc is less than 0.5, it should be adjusted by 1 - auc. 
## Or re-run them:
cl.roc(1 - stat,label, pos=levels(factor(label))[2],plot = TRUE)
cl.auc(1 - stat,label,pos=levels(factor(label))[2])

}

\keyword{classif}
