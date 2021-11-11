% lwc-11-10-2011: commence
\name{df.util}
\alias{df.summ}
\alias{vec.summ}
\alias{vec.summ.1}

\title{
  Summary Utilities
}
\description{
  Functions to summarise data.
}
\usage{
df.summ(dat, method=vec.summ,\dots)
  
vec.summ(x)

vec.summ.1(x)
}

% --------------------------------------------------------------------
\arguments{
  \item{dat}{
  A data frame or matrix of data set. 
  }
  \item{x}{
  A vector value.
  } 
  \item{method}{Summary method such as \code{vec.summ} and \code{vec.summ.1}.
  For user-defined methods, see examples below.
  }
  \item{\dots}{
    Additional parameters to \code{method} function. 
  }
}

% ----------------------------------------------------------------------
\value{
   \code{df.summ} returns a summarised data frame. 

   \code{vec.summ} returns an vector of number of variables (exclusing NAs),
    minimum, mean, median, maximum and standard derivation.
    
   \code{vec.summ.1} returns an vector of number of variables (exclusing NAs),
   mean, median, 95\% confidence interval of median, IQR and standard derivation.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------
\examples{
data(abr1)
dat <- (abr1$pos)[,110:150]
cls <- factor(abr1$fact$class)

## sort out missing value
dat <- mv.zene(dat)

## summary of an individual column
vec.summ(dat[,2])
vec.summ.1(dat[,2])

## summary of data frame
summ   <- df.summ(dat)                       ## default: vec.summ
summ.1 <- df.summ(dat, method=vec.summ.1)

## summary by groups
by(dat, list(cls=cls), df.summ)

## User-defined summary function: 
vec.segment <- function(x, bar=c("SD", "SE", "CI"))
{  
  bar <- match.arg(bar)

  centre <- mean(x, na.rm = TRUE)

  if (bar == "SD") {
    stderr <- sd(x, na.rm = TRUE)        ## Standard derivation (SD)
    lower  <- centre - stderr
    upper  <- centre + stderr
  } else if (bar == "SE") {      ## Standard error(SE) of mean
    stderr <- sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))
    ## stderr <- sqrt(var(x, na.rm = TRUE)/length(x[complete.cases(x)]))
    lower  <- centre - stderr
    upper  <- centre + stderr
  } else if (bar == "CI") {      ## Confidence interval (CI), here 95%.
    conf   <- t.test(x)$conf.int
    lower  <- conf[1]
    upper  <- conf[2]
  } else {
    stop("'method' invalid")
  }

  res <- c(lower=lower, centre=centre, upper=upper)
  return(res)
}

## test it
vec.segment(dat[,2])
summ.2 <- df.summ(dat, method=vec.segment, bar="SE")

## ----------------------------------------------------------
#' iris data
df.summ(iris)

#' Group summary
## library(plyr)
## ddply(iris, .(Species), df.summ)
## (tmp <- dlply(iris, .(Species), df.summ, method=vec.segment))
##do.call("rbind", tmp)

#' or you can use summarise to get the group summary for single variable:
## ddply(iris, .(Species), summarise, 
##      mean=mean(Sepal.Length), std=sd(Sepal.Length))

}

\keyword{util}
