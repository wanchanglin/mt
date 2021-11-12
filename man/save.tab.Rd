% lwc-15-02-2010: Previous name is my.save.tab
%
\name{save.tab}
\alias{save.tab}
\title{
  Save List of Data Frame or Matrix into CSV File
}
\description{
  Save a list of data frame or matrix into a CSV file
}
\usage{
  save.tab(x, filename="temp.csv", firstline="\n")
}

% --------------------------------------------------------------------
\arguments{
  \item{x}{
  A list of data frame or matrix. 
  }
  \item{filename}{
  A character string for saved file name.
  }
  \item{firstline}{
  A string giving some description of the saved file.
  }
}

\details{
  This function gives a quick option to save a set of data frame or matrix into a
  single table file.
}

\value{
  No return value, called for side effects.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

\seealso{
  \code{\link{write.table}}
}

% ----------------------------------------------------------------------
\examples{
\dontrun{
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
fs <- fs[1:3]

## save consistency of feature selection
filename  <- "fs.csv"
firstline <- paste('\nResults of feature selection', sep='')

save.tab(fs, filename, firstline)
}
}

\keyword{manip}
