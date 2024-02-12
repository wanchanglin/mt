% lwc-16-09-2006
% lwc-15-02-2010
% lwc-15-04-2010: add combn.pw
\name{dat.sel}
\alias{dat.sel}
\alias{combn.pw}

\title{
  Generate Pairwise Data Set
}
\description{
  Generate index or data set of pairwise combination based on class labels.
}
\usage{
  combn.pw(cls, choices = NULL) 
  
  dat.sel(dat, cls, choices = NULL)
}

% --------------------------------------------------------------------
\arguments{
  \item{dat}{
  A data frame or matrix of data set. 
  }
  \item{cls}{
  A factor or vector of class labels or categorical data.
  }
  \item{choices}{
  The vector or list of class labels to be chosen for binary combination.
  }
}

\details{
  If \code{choices} is \code{NULL}, all binary combinations will be computed.
  If \code{choices} has one class label, the comparisons between this one and 
  any other classes will be calculated. If \code{choices} has more than two 
  classes, all binary combinations in \code{choices} will be generated. 
  For details, see \code{examples} below.
}

% ----------------------------------------------------------------------
\value{
  \code{combn.pw} returns a data frame of index (logical values).
  
  \code{dat.set} returns a list of list with components:
  \item{dat}{
  Pairwise data set.
  }
  \item{cls}{
  Pairwise class label.
  }
}

\seealso{
  Applications of \code{dat.sel} in \code{\link{pca_plot_wrap}}, 
  \code{\link{lda_plot_wrap}} and \code{\link{pls_plot_wrap}}.

}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------
\examples{
data(iris)
x <- subset(iris, select = -Species)
y <- iris$Species
 
## generate data set with class "setosa" and "virginica"
binmat.1 <- dat.sel(x,y,choices=c("setosa","virginica"))
names(binmat.1)

## generate data sets for "setosa" vs other classes. These are: 
## "setosa" and "versicolor", "setosa" and "virginica".
binmat.2 <- dat.sel(x,y,choices=c("setosa"))
names(binmat.2)

## generate data set with combination of each class. These are:  
## "setosa" and "versicolor", "setosa" and "virginica",  
## "versicolor" and "virginica" 
binmat.3 <- dat.sel(x,y,choices= NULL)
names(binmat.3)

data(abr1)
cls <- factor(abr1$fact$class)
dat <- preproc(abr1$pos, method="log")

## There are some examples of 'choices'
choices <- c("2")
choices <- c("2","3","4")
choices <- list(c("2","3"),c("4","5"))
choices <- NULL
idx <- combn.pw(cls,choices=choices)  

dat.pw <- dat.sel(dat, cls,choices=choices)

}

\keyword{manip}
