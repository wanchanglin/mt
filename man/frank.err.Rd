% lwc-13-11-2006

\name{frank.err}
\alias{frank.err}
\title{
  Feature Ranking and Validation on Feature Subset
}
\description{
  Get feature ranking on the training data and validate selected feature subsets by
  estimating their classification error rate.
}
\usage{
frank.err(dat.tr, cl.tr, dat.te, cl.te, cl.method="svm",
          fs.method="fs.auc", fs.order=NULL, fs.len="power2", \dots)
}

% --------------------------------------------------------------------
\arguments{
  \item{dat.tr}{
  A data frame or matrix of training data. Feature ranking and classification model 
  are carried on this data set.
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
  \item{cl.method}{
    Classification method to be used. Any classification methods can be employed 
    if they have method \code{predict} (except \code{knn}) with output of predicted class
    label or one component with name of \code{class} in the returned list, such as 
    \code{randomForest}, \code{svm}, \code{knn} and \code{lda}. 
  }
  \item{fs.method}{
    Feature ranking method. If \code{fs.order} is not \code{NULL}, it is ignored.
  }             
  \item{fs.order}{
   A vector of feature order. Default is \code{NULL} and then the feature selection will be 
   performed on the training data. 
  }
  \item{fs.len}{
    The lengths of feature subsets used for validation. For details, see \code{\link{get.fs.len}}. 
  }
  \item{\dots}{Additional parameters to \code{fs.method} or \code{cl.method}.}
}

% ----------------------------------------------------------------------------
\value{
  A list with components:
  \item{cl.method}{Classification method used.}
  \item{fs.len}{The lengths of feature subsets used for validation.}
  \item{error}{Error rate for each feature length.}
  \item{fs.method}{Feature ranking method used.}
  \item{fs.order}{Feature order vector.}
  \item{fs.rank}{Feature ranking score vector.}
}


% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

\seealso{
  \code{\link{frankvali}}, \code{\link{get.fs.len}}
}

% ----------------------------------------------------------------------
\examples{
data(abr1)
dat <- abr1$pos
x   <- preproc(dat[,110:500], method="log10")  
y   <- factor(abr1$fact$class)        

dat <- dat.sel(x, y, choices=c("1","6"))
x.1 <- dat[[1]]$dat
y.1 <- dat[[1]]$cls

idx <- sample(1:nrow(x.1), round((2/3)*nrow(x.1)), replace=FALSE) 
## construct train and test data 
train.dat  <- x.1[idx,]
train.cl   <- y.1[idx]
test.dat   <- x.1[-idx,]   
test.cl    <- y.1[-idx] 

## validate feature selection on some feature subsets
res <- frank.err(train.dat, train.cl, test.dat, test.cl, 
                 cl.method="knn", fs.method="fs.auc",  
                 fs.len="power2")
names(res)
## full feature order list
res$fs.order

## validation on subsets of feature order 
res$error

## or first apply feature selection
fs <- fs.auc(train.dat,train.cl)
## then apply error estimation for each selected feature subset
res.1 <- frank.err(train.dat, train.cl, test.dat, test.cl, 
                   cl.method="knn", fs.order=fs$fs.order,  
                   fs.len="power2")

res.1$error

}

\keyword{classif}
