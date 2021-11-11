% wll-11-01-2007
% wll-02-07-2007: add some stuff of posteriot, marin and auc

\name{classifier}
\alias{classifier}
\title{
  Wrapper Function for Classifiers
}
\description{
  Wrapper function for classifiers. The classification model is built up on the 
  training data and error estimation is performed on the test data.
}
\usage{
classifier(dat.tr, cl.tr, dat.te=NULL, cl.te=NULL, method,
           pred.func=predict,\dots)
}

% --------------------------------------------------------------------
\arguments{
  \item{dat.tr}{
  A data frame or matrix of training data. The classification model are built 
  on it.
  }
  \item{cl.tr}{
  A factor or vector of training class.
  }
  \item{dat.te}{
  A data frame or matrix of test data. Error rates are calculated on this data set.
  }
  \item{cl.te}{
  A factor or vector of test class.
  }
  \item{method}{
    Classification method to be used. Any classification methods can be employed 
    if they have method \code{predict} (except \code{knn}) with output of predicted class
    label or one component with name of \code{class} in the returned list, such as 
    \code{randomForest}, \code{svm}, \code{knn} and \code{lda}. 
    Either a function or a character string naming the function to be called
  }
  \item{pred.func}{
  Predict method (default is \code{predict}). Either a function or a character 
  string naming the function to be called. 
  }
  \item{\dots}{Additional parameters to \code{method}.}
}

% ----------------------------------------------------------------------------
\value{
  A list including components:
  \item{err}{Error rate of test data.}
  \item{cl}{The original class of test data.}
  \item{pred}{The predicted class of test data.}
  \item{posterior}{
    Posterior probabilities for the classes if \code{method} provides posterior 
    output.    
  }
  \item{acc}{ Accuracy rate of classification.}
  \item{margin}{
   The margin of predictions from classifier \code{method} if it provides posterior 
   output. 
    
    The margin of a data point is defined as the proportion of probability for the 
    correct class minus maximum proportion of probabilities for the other classes. 
    Positive margin means correct classification, and vice versa. This idea come
    from package \pkg{randomForest}. For more details, see 
    \code{\link[randomForest]{margin}}.
  }
  \item{auc}{
   The area under receiver operating curve (AUC) if classifier \code{method} 
   produces posterior probabilities and the classification is for two-class 
   problem. 
  }
}

% ----------------------------------------------------------------------------
\note{
  The definition of margin is based on the posterior probabilities. Classifiers,
  such as \code{\link[randomForest]{randomForest}}, \code{\link[e1071]{svm}},
  \code{\link[MASS]{lda}}, \code{\link[MASS]{qda}}, \code{\link{pcalda}} and
  \code{\link{plslda}}, do output posterior probabilities. But 
  \code{\link[class]{knn}} does not.   
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

\seealso{
  \code{\link{accest}}, \code{\link{maccest}}
}

% ----------------------------------------------------------------------
\examples{
data(abr1)
dat <- preproc(abr1$pos[,110:500], method="log10")  
cls <- factor(abr1$fact$class)        

## tmp <- dat.sel(dat, cls, choices=c("1","2"))
## dat <- tmp[[1]]$dat
## cls <- tmp[[1]]$cls

idx <- sample(1:nrow(dat), round((2/3)*nrow(dat)), replace = FALSE) 
## constrcuct train and test data 
train.dat  <- dat[idx,]
train.cl   <- cls[idx]
test.dat   <- dat[-idx,]       
test.cl    <- cls[-idx] 

## estimates accuracy
res <- classifier(train.dat, train.cl, test.dat, test.cl, 
                  method="randomForest")
res
## get confusion matrix
cl.rate(obs=res$cl, res$pred)   ## same as: cl.rate(obs=test.cl, res$pred)

## Measurements of Forensic Glass Fragments
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

res.1 <- classifier(tr.x, tr.y, te.x, te.y, method="svm")
res.1
cl.rate(obs=res.1$cl, res.1$pred) 

## classification performance for the two-class case.
pos <- "WinF"                              ## select positive level
cl.perf(obs=res.1$cl, pre=res.1$pred, pos=pos)
## ROC and AUC
cl.roc(stat=res.1$posterior[,pos],label=res.1$cl, pos=pos)

}

\keyword{classif}
