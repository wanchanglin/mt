% wll-02-07-2007: 
%
\name{plot.pcalda}
\alias{plot.pcalda}
\title{
  Plot Method for Class 'pcalda'
}
\description{
  Plot linear discriminants of \code{pcalda}.
}

\usage{
\method{plot}{pcalda}(x, dimen, \dots)
}

% ----------------------------------------------------------------------------
\arguments{
  \item{x}{An object of class \code{pcalda}. }
  \item{dimen}{ The index of linear discriminants to be used for the plot. }
  \item{\dots}{ Further arguments.  See corresponding entry in
    \code{\link{xyplot}} for non-trivial details. One argument is
    \code{ep}: an integer for plotting ellipse. \code{1} and \code{2} for 
    plotting overall and group ellipse, respectively. Otherwise, none.
    For details, see \code{\link{panel.elli.1}}. }
}

% ----------------------------------------------------------------------------
\details{
  This function is a method for the generic function \code{plot()} for
  class \code{pcalda}. If the length of \code{dimen} is greater
  than 2, a pairs plot is used. If the length of \code{dimen} is equal
  to 2, a scatter plot is drawn. Otherwise, the dot plot is drawn for
  the single component.
}

\value{An object of class \code{"trellis"}.}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}


\seealso{
  \code{\link{pcalda}}, \code{\link{predict.pcalda}},
  \code{\link{lda_plot_wrap}},\code{\link{panel.elli.1}}.
}

\examples{
data(abr1)
cl   <- factor(abr1$fact$class)
dat  <- abr1$pos

model <- pcalda(dat,cl)

## Second component versus first
plot(model,dimen=c(1,2),main = "Training data",ep=2)
## Pairwise scatterplots of several components 
plot(model,main = "Training data",ep=1)

## The first component
plot(model,dimen=c(1),main = "Training data")
}

\keyword{plot}
