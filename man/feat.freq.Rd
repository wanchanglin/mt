%% lwc-15-03-2007:
%% wll-04-12-2015:tidy up 
\name{feat.freq}
\alias{feat.freq}

\title{Frequency and Stability of Feature Selection }

\description{Frequency and stability of feature selection. }

\usage{
  feat.freq(x,rank.cutoff=50,freq.cutoff=0.5)
}

% --------------------------------------------------------------------
\arguments{
  \item{x}{ A matrix or data frame of feature orders. }
  \item{rank.cutoff}{A numeric value for cutoff of top features.}
  \item{freq.cutoff}{A numeric value for cutoff of feature frequency.}
}

% -------------------------------------------------------------------------
\value{
  A list with components:
  \item{freq.all}{Feature frequencies. }
  \item{freq}{Feature frequencies larger than \code{freq.cutoff}. }
  \item{stability}{Stability rate of feature ranking.}
  \item{rank.cutoff}{ Top feature order cut-off used. }
  \item{freq.cutoff}{ Feature frequency cut-off used. }
}

\references{
  Davis, C. A., et al., (2006) Reliable gene signatures for microarray
  classification: assessment of stability and performance.
  \emph{Bioinformatics}, vol.22, no.19, 2356 - 2363.

  Michiels, S., et al., (2005) Prediction of cancer outcome with
  microarrays: a multiple random validation strategy. \emph{Lancet},
  vol.365, 488 - 492.
}

\author{ Wanchang Lin  }

\seealso{ \code{\link{feat.rank.re}} }

% ----------------------------------------------------------------------
\examples{
## prepare data set
data(abr1)
cls <- factor(abr1$fact$class)
dat <- abr1$pos
## dat <- abr1$pos[,110:1930]

## fill zeros with NAs
dat <- mv.zene(dat)

## missing values summary
mv <- mv.stats(dat, grp=cls) 
mv    ## View the missing value pattern

## filter missing value variables
## dim(dat)
dat <- dat[,mv$mv.var < 0.15]
## dim(dat)

## fill NAs with mean
dat <- mv.fill(dat,method="mean")

## log transformation
dat <- preproc(dat, method="log10")

## select class "1" and "2" for feature ranking
ind <- grepl("1|2", cls)
mat <- dat[ind,,drop=FALSE] 
mat <- as.matrix(mat)
grp <- cls[ind, drop=TRUE]   

## use resampling method of bootstrap 
pars   <- valipars(sampling="boot",niter=10, nreps=5)
z      <- feat.rank.re(mat,grp,method="fs.plsvip",pars = pars)

## compute the frequency and stability of feature selection 
freq <- feat.freq(z$order.list,rank.cutoff=50,freq.cutoff=0.5)

}

\keyword{classif}
