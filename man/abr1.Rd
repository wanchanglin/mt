% lwc-12-10-2006: First draft
% lwc-20-02-2007: minor changes

\name{abr1}
\alias{abr1}

\title{abr1 Data}
\usage{
data(abr1)
}

\description{
  An FIE-MS data.
}

\details{
  \code{abr1} is an FIE-MS data matrices developed from analysis of
  samples representing a time course of pathogen attack in a model plant
  species (Brachypodium distachyon). The data was developed in a single
  batch with all samples randomised using a Thermo LTQ linear ion trap.
  Both positive and negative ion mode are given (\code{abr1$pos} and
  \code{abr1$neg}).
}

\value{
  A list with the following elements:
  \item{fact}{A data frame containing experimental meta-data.}
  \item{pos}{A data frame for positive data with 120 observations and
    2000 variables.}
  \item{neg}{A data frame for negative data with 120 observations and
    2000 variables.}
}

\examples{
# Load data set
data(abr1)

# Select data set
dat <- abr1$neg

# number of observations and variables
dim(dat)

# Transform data
dat.log   <- preproc(dat, method = "log")
dat.sqrt  <- preproc(dat, method = "sqrt")
dat.asinh <- preproc(dat, method = "asinh")

op <- par(mfrow=c(2,2), pch=16)
matplot(t(dat),main="Original",type="l",col="blue",
     ylab="Intensity")
matplot(t(dat.log),main="Log",type="l",col="green",
     ylab="Intensity")
matplot(t(dat.sqrt),main="Sqrt",type="l",col="red",
     ylab="Intensity")
matplot(t(dat.asinh),main="ArcSinh)",type="l",col="black",
     ylab="Intensity")
par(op)
mtext("Data set", line=2.5, font=3, cex=1.5)

}    
\keyword{datasets}
