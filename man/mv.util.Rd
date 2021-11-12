% lwc-19-10-2011: commence
\name{mv.util}
\alias{mv.stats}
\alias{mv.fill}
\alias{mv.zene}

\title{
  Missing Value Utilities
}
\description{
  Functions to handle missing values of data set.
}
\usage{
mv.stats(dat,grp=NULL,\dots) 
  
mv.fill(dat,method="mean",ze_ne = FALSE)

mv.zene(dat)
  
}

% --------------------------------------------------------------------
\arguments{
  \item{dat}{
  A data frame or matrix of data set. 
   }
  \item{grp}{
  A factor or vector of class.
  }
  \item{method}{Univariate imputation method for missing value. For details, see 
  examples below.
  }
  \item{ze_ne}{ 
     A logical value indicating whether the zeros or negatives should be 
     treated as missing values.
  }
  \item{\dots}{
    Additional parameters to \code{mv.stats} for plotting using \pkg{lattice}.
  }
}

% ----------------------------------------------------------------------
\value{
   \code{mv.fill} returns an imputed data frame. 

   \code{mv.zene} returns an NA-filled data frame.
   
   
   \code{mv.stats} returns a list including the components:
    \itemize{ 
     \item \code{mv.overall}: Overall missng value rate.
     \item \code{mv.var}: Missing value rate per variable (column).
     \item \code{mv.grp}: A matrix of missing value rate for different groups
     if argument \code{grp} is given.
     \item \code{mv.grp.plot}:  An object of class \code{trellis} for plotting  
     of \code{mv.grp} if argument \code{grp} is given.
    }
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------
\examples{
data(abr1)
dat <- abr1$pos[,1970:1980]
cls <- factor(abr1$fact$class)

## fill zeros with NAs
dat <- mv.zene(dat)

## missing values summary
mv <- mv.stats(dat, grp=cls) 
plot(mv$mv.grp.plot)

## fill NAs with mean
dat.mean <- mv.fill(dat,method="mean")

## fill NAs with median
dat.median <- mv.fill(dat,method="median")

## -----------------------------------------------------------------------
## fill NAs with user-defined methods: two examples given here.
## a.) Random imputation function:
rand <- function(x,...) sample(x[!is.na(x)], sum(is.na(x)), replace=TRUE)

## test this function:
(tmp <- dat[,1])        ## an vector with NAs
## get the randomised values for NAs
rand(tmp)

## fill NAs with method "rand"
dat.rand <- mv.fill(dat,method="rand")

## b.) "Low" imputation function:
"low" <- function(x, ...) {
  max(mean(x,...) - 3 * sd(x,...), min(x, ...)/2)
}
## fill NAs with method "low"
dat.low <- mv.fill(dat, method="low") 

## summary of imputed data set
df.summ(dat.mean)

}

\keyword{util}
