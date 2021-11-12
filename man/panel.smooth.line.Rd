% lwc-14-02-2010:
%
\name{panel.smooth.line}
\alias{panel.smooth.line}
\title{ Panel Function for Plotting Regression Line}
\description{
  \pkg{lattice} panel function for plotting regression line with red colour. 
}
\usage{
panel.smooth.line(x, y,\dots) 
}
\arguments{
  \item{x, y}{ Variables to be plotted.}
  \item{\dots}{Further arguments. See corresponding entry in
    \code{\link{xyplot}} for non-trivial details.  }
}

\value{
	An object of class \code{"trellis"}.
}

\author{Wanchang Lin }

\examples{
library(lattice) 
data(iris)
splom(~iris[,1:4], varname.cex = 1.0, pscales = 0, panel = panel.smooth.line)
}


\keyword{plot}
