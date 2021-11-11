% wll-23-05-2007: commence
% wll-03-10-2007: add plsc stuff 
% wll-13-12-2007: major changes

%
\name{plot.plsc}
\alias{plot.plsc}
\alias{plot.plslda}

\title{
  Plot Method for Class 'plsc' or 'plslda'
}
\description{
  Plot latent components of \code{plsc} or \code{plslda}.
}

\usage{
\method{plot}{plsc}(x, dimen, \dots)

\method{plot}{plslda}(x, dimen, \dots)
}

% ----------------------------------------------------------------------------
\arguments{
  \item{x}{ An object of class \code{plsc} or \code{plslda}. }
  \item{dimen}{ The index of latent components to be used for the plot. }
  \item{\dots}{ Further arguments.  See corresponding entry in
    \code{\link{xyplot}} for non-trivial details. One argument is
    \code{ep}: an integer for plotting ellipse. \code{1} and \code{2} for 
    plotting overall and group ellipse, respectively. Otherwise, none.
    For details, see \code{\link{panel.elli.1}}. }
}

% ----------------------------------------------------------------------------
\details{
  Two functions are methods for the generic function \code{plot()} of
  class \code{plsc} and \code{plslda}. 

  If the length of \code{dimen} is greater than 2, a pairs plot is used.
  If the length of \code{dimen} is equal to 2, a scatter plot is drawn.
  Otherwise, the dot plot is drawn for the single component.
}

\value{An object of class \code{"trellis"}.}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}


\seealso{
  \code{\link{plsc}}, \code{\link{predict.plsc}},\code{\link{plslda}}, 
  \code{\link{predict.plslda}}, \code{\link{pls.plot.wrap}},
  \code{\link{panel.elli.1}}.
}

\examples{
data(abr1)
cl   <- factor(abr1$fact$class)
dat  <- abr1$pos

mod.plsc    <- plsc(dat,cl,ncomp=4)
mod.plslda  <- plslda(dat,cl,ncomp=4)

## Second component versus first
plot(mod.plsc,dimen=c(1,2),main = "Training data", ep = 2)
plot(mod.plslda,dimen=c(1,2),main = "Training data", ep = 2)

## Pairwise scatterplots of several components 
plot(mod.plsc, main = "Training data", ep = 1)
plot(mod.plslda, main = "Training data", ep = 1)

## single component
plot(mod.plsc,dimen=c(1),main = "Training data")
plot(mod.plslda,dimen=c(1),main = "Training data")

}

\keyword{plot}
