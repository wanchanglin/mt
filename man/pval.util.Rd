% lwc-21-03-2007

\name{pval.util}
\alias{pval.test}
\alias{pval.reject}

\title{
  P-values Utilities
}
\description{
  Functions to handle p-values of data set.
}
\usage{
  pval.test(x,y, method="oneway.test",\dots)

  pval.reject(adjp,alpha)

}

% --------------------------------------------------------------------
\arguments{
  \item{x}{
  A data frame or matrix of data set. 
  }
  \item{y}{
  A factor or vector of class.
  }
  \item{method}{Hypothesis test such as \code{t.test} and \code{wilcox.test}.}
  \item{adjp}{A matrix-like p-values of simultaneously testing.}
  \item{alpha}{A vector of cutoff of p-values or Type I error rate.}
  \item{\dots}{
  Arguments to pass.
  }
}

% ----------------------------------------------------------------------------
\value{
  \code{pval.test} returns a vector of p-values.

  \code{pval.reject} returns a matrix indicating rejected number according to 
  cutoff.
}


% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------
\examples{
library(lattice)

## Example for pval.test and pval.reject
## prepare data set
data(abr1)
cls <- factor(abr1$fact$class)
dat <- abr1$pos[,200:500]
dat <- preproc(dat, method="log")

## select class "1" and "2" for feature ranking
ind <- grepl("1|2", cls)
dat <- dat[ind,,drop=FALSE] 
cls <- cls[ind, drop=TRUE]   

## univariate p-values and its adjusted p-values
pval <- sort(pval.test(dat, cls, method="t.test"))

## adjust p-values
pval.ad <- sapply(c("fdr","bonferroni","BY"), function(y){
  p.adjust(pval, method=y)
})
pval.ad <- cbind(raw=pval, pval.ad)
pval.reject(pval.ad,c(0.005, 0.01, 0.05))

## plot the all p-values
tmp <- cbind(pval.ad, idx=1:nrow(pval.ad))
tmp <- data.frame(tmp)

# pval_long <- melt(tmp, id="idx")
pval_long <- data.frame(idx = tmp$idx, stack(tmp, select = -idx))
pval_long <- pval_long[c("idx", "ind", "values")]
names(pval_long) <- c("idx", "variable", "value")

pval.p <- xyplot(value~idx, data=pval_long, groups=variable,
                par.settings = list(superpose.line = list(lty=c(1:7))),
                as.table = TRUE, type="l", 
                par.strip.text = list(cex=0.65), ylim=c(-0.005, 1.0),
                ylab="P-values", xlab="Index of variables",
                main="p-values",
                auto.key = list(lines=TRUE, points = FALSE,space="right"),
                panel = function(x, y,...) {
                  panel.xyplot(x, y, ...)
                  panel.abline(h = 0.05, col = "red",lty =2)
                })
pval.p

}

\keyword{util}
