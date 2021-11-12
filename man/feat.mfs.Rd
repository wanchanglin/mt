% wll-12-12-2007:
% lwc-16-02-2010:
%
\name{feat.mfs}
\alias{feat.mfs}
\alias{feat.mfs.stab}
\alias{feat.mfs.stats}

\title{
  Multiple Feature Selection
}

\description{ Multiple feature selection with or without resampling
procedures. }

\usage{
feat.mfs(x,y,method,pars = valipars(),is.resam = TRUE, \dots)
         
feat.mfs.stab(fs.res,rank.cutoff = 20,freq.cutoff = 0.5)

feat.mfs.stats(fs.stats,cumu.plot=FALSE, main="Stats Plot", 
               ylab="Values", xlab="Index of variable", \dots)

}

\arguments{

  \item{x}{ A matrix or data frame containing the explanatory variables.}

  \item{y}{ A factor specifying the class for each observation. }

  \item{method}{ Multiple feature selection/ranking method to be used. }

  \item{pars}{ A list of resampling scheme. See \code{\link{valipars}}
    for details. }
  
  \item{is.resam}{ A logical value indicating whether the resampling
    should be applied. }
  
  \item{fs.res}{ A list obtained by running \code{feat.mfs} .} 

  \item{rank.cutoff}{Cutoff of top features for frequency calculating.}

  \item{freq.cutoff}{ Cutoff of feature frequency. }

  \item{fs.stats}{A matrix of feature statistics or values outputted by
    \code{feat.mfs} }
  
  \item{cumu.plot}{ A logical value indicating the cumulative scores
  should be plotted. }

  \item{main,xlab,ylab}{Plot parameters}

  \item{\dots}{Additional parameters.}
}

% -------------------------------------------------------------------------
\value{

  \code{feat.mfs} returns a list with components:

  \item{fs.order}{A data frame of feature order from best to worst.}
  \item{fs.rank}{A matrix of feature ranking scores.}
  \item{fs.stats}{A matrix of feature statistics or values.}
  \item{all}{A list of output of \code{\link{feat.rank.re}} for each feature 
    selection method.}
  
  \code{feat.mfs.stab} returns a list with components:

  \item{fs.freq}{Feature frequencies larger than \code{freq.cutoff}. }
  \item{fs.subs}{Feature with frequencies larger than \code{freq.cutoff}. }
  \item{fs.stab}{Stability rate of feature ranking.}
  \item{fs.cons}{A matrix of feature consensus table based on feature
    frequency.}
  
  \code{feat.mfs.stats} returns a list with components:

  \item{stats.tab}{A statistical values with their corresponding names.}
  \item{stats.long}{Long-format of statistical values for plotting.} 
  \item{stats.p}{ An object of class "trellis".}

}

\details{

  \code{feat.mfs.stab} summarises multiple feature selection only when
  resampling strategy is employed (i.e. \code{is.resam} is \code{TRUE}
  when calling \code{feat.mfs}). It obtains these results based on
  \code{feat.mfs}'s returned value called \code{all}.
  
  \code{feat.mfs.stats} handles the statistical values or scores. Its
  purpose is to provide a guidance in selecting the best number of
  features by spotting the elbow point. This method should work in
  conjunction with plotting of p-values and their corresponding adjusted
  values such as FDR and Bonferroni in the multiple hypothesis test.
}

\note{

  The feature order can be computed directly from the overall statistics
  \code{fs.stats}. It is, however, slightly different from
  \code{fs.order} obtained by rank aggregation when resampling is
  employed.
  
  The \code{fs.cons} and \code{fs.freq} are computed based on
  \code{fs.order}.
  
}

\author{ Wanchang Lin  }

\seealso{ \code{\link{feat.rank.re}}, \code{\link{feat.freq}} }

% -------------------------------------------------------------------------
\examples{
\dontrun{
library(lattice)	
data(abr1)
dat <- preproc(abr1$pos[,200:400], method="log10")  
cls <- factor(abr1$fact$class)

tmp <- dat.sel(dat, cls, choices=c("1","2"))
x   <- tmp[[1]]$dat
y   <- tmp[[1]]$cls

fs.method <- c("fs.anova","fs.rf","fs.rfe")
fs.pars   <- valipars(sampling="cv",niter=10,nreps=5)
fs <- feat.mfs(x, y, fs.method, fs.pars)   ## with resampling
names(fs)

## frequency, consensus and stabilities of feature selection 
fs.stab <- feat.mfs.stab(fs)
print(fs.stab$fs.cons,digits=2,na.print="")

## plot feature selection frequency
freq <- fs.stab$fs.freq
dotplot(freq$fs.anova, type="o", main="Feature Selection Frequencies")
barchart(freq$fs.anova)

## rank aggregation 
fs.agg <- feat.agg(fs$fs.rank)

## stats table and plotting
fs.stats <- fs$fs.stats
tmp <- feat.mfs.stats(fs.stats, cumu.plot = TRUE)
tmp$stats.p
fs.tab <- tmp$stats.tab
## convert to matrix
fs.tab <- list2df(un.list(fs.tab))

## without resampling
fs.1 <- feat.mfs(x, y, method=fs.method, is.resam = FALSE)
}
}
\keyword{classif}
