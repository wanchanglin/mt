% wll-10-06-2015: copy from pcaplot.Rd and change it in Vim.

\name{mdsplot}
\alias{mdsplot}

\title{Plot Classical Multidimensional Scaling}
\description{
  Plot metric MDS with categorical information.
}

\usage{
mdsplot(x, y, method = "euclidean", dimen = c(1,2), \dots)
}

\arguments{
  \item{x}{A matrix or data frame to be plotted. }
  \item{y}{A factor or vector giving group information of columns of \code{x}.}

  \item{method}{ 
  The distance measure to be used. This must be one of "euclidean", "maximum", 
  "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring
   can be given. It is only for \code{mds_plot_wrap}.
  }

  \item{dimen}{
  A vector of index of dimentonal to be plotted. Only two dimentions are 
  are allowed.}
  \item{\dots}{
  Further arguments to \code{\link{prcomp}} or \code{lattice}. See 
  corresponding entry in \code{\link{xyplot}} for non-trivial details of 
  \code{lattice}. For \code{pcaplot}, one argument is
    \code{ep}: an integer for plotting 95\% ellipse. \code{1} and \code{2} for 
    plotting overall and group ellipse, respectively. Otherwise, none.
  For details, see \code{\link{panel.elli.1}}. 
  }
}

\value{
 \code{mdsplot} returns an object of class \code{"trellis"}.
}


\author{
  Wanchang Lin 
}

\seealso{
  \code{\link{grpplot}}, \code{\link{panel.elli}}, \code{\link{mds_plot_wrap}},
  \code{\link{pcaplot}}
}


% ----------------------------------------------------------------------------
\examples{
## examples of 'mdsplot'

data(iris)
x <- iris[,1:4]
y <- iris[,5]
mdsplot(x,y, dimen=c(1,2),ep=2)
mdsplot(x,y, dimen=c(2,1),ep=1)

tmp <- mdsplot(x,y, ep=2, conf.level = 0.9)
tmp

## change symbol's color, type and size
mdsplot(x, y, main="IRIS DATA", cex=1.2,
  auto.key=list(space="right", col=c("black","blue","red"), cex=1.2),
  par.settings = list(superpose.symbol = list(col=c("black","blue","red"),
                                              pch=c(1:3))))
}

\keyword{plot}

