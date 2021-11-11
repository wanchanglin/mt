% lwc-15-02-2010
% lwc-25-02-2010: give example

\name{feat.agg}
\alias{feat.agg}

\title{ Rank aggregation by Borda count algorithm }

\description{ Use Borda count to get the final feature order. }

\usage{ feat.agg(fs.rank.list) }

\arguments{
  \item{fs.rank.list}{ A data frame of feature orders by different
    feature selectors. }
}

\value{
  A list with components:
  \item{fs.order}{Final feature order. }
  \item{fs.rank}{Aggregated rank list by Borda count. }
}

\author{ Wanchang Lin  }

\seealso{ \code{\link{feat.rank.re}}, \code{\link{feat.mfs}} }

% ----------------------------------------------------------------------
\examples{
data(abr1)
dat <- preproc(abr1$pos[,200:400], method="log10")  
cls <- factor(abr1$fact$class)

## feature selection without resampling
fs <- feat.mfs(dat, cls, method=c("fs.anova","fs.rf","fs.rfe"), 
               is.resam=FALSE)
## rank aggregation 
fs.1 <- feat.agg(fs$fs.rank)
names(fs.1)
}

\keyword{classif}
