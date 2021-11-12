% wll-13-12-2007:
% wll-15-07-2015: move 'ep' for 'pcaplot' into \dots.
\name{pcaplot}
\alias{pcaplot}
\alias{pca.plot}
\alias{pca.comp}

\title{Plot Function for PCA with Grouped Values}
\description{
  Plot function for PCA with grouped values.
}

\usage{
pcaplot(x, y, scale = TRUE, pcs = 1:2, \dots)

pca.plot(x, y, scale=TRUE, abbrev = FALSE, ep.plot=FALSE,\dots)

pca.comp(x, scale=FALSE, pcs=1:2,\dots)
}

\arguments{

  \item{x}{A matrix or data frame to be plotted. }

  \item{y}{A factor or vector giving group information of columns of
    \code{x}.}

  \item{scale}{ A logical value indicating whether the data set \code{x}
    should be scaled. }

  \item{pcs}{A vector of index of PCs to be plotted.}

  \item{ep.plot}{ A logical value indicating whether the ellipse should
    be plotted. }

  \item{abbrev}{ Whether the group labels are abbreviated on the plots.
    If \code{abbrev > 0} this gives \code{minlength} in the call to
    \code{abbreviate}. }

  \item{\dots}{ Further arguments to \code{\link{prcomp}} or
    \code{lattice}. See corresponding entry in \code{\link{xyplot}} for
    non-trivial details of \code{lattice}. For \code{pcaplot}, one
    argument is \code{ep}: an integer for plotting ellipse. \code{1} and
    \code{2} for plotting overall and group ellipse, respectively.
    Otherwise, none. For details, see \code{\link{panel.elli.1}}. }

}

\value{
 \code{pcaplot} returns an object of class \code{"trellis"}.

 \code{pca.comp} returns a list with components:
 \item{scores}{ PCA scores}
 \item{vars}{Proportion of variance}
 \item{varsn}{A vector of string indicating the percentage of variance.}

}

\note{
  Number of columns of \code{x} must be larger than 1. \code{pcaplot} uses
  \code{lattice} to plot PCA while \code{pca.plot} uses the basic graphics 
  to do so. \code{pca.plot} plots PC1 and PC2 only.
}  

\author{
  Wanchang Lin 
}

\seealso{ \code{\link{grpplot}}, \code{\link{panel.elli.1}},
  \code{\link{pca.plot.wrap}} }


% ----------------------------------------------------------------------------
\examples{
## examples of 'pcaplot'
data(iris)
pcaplot(iris[,1:4], iris[,5],pcs=c(2,1),ep=2)
## change confidence interval (see 'panel.elli.1')
pcaplot(iris[,1:4], iris[,5],pcs=c(1,2),ep=2, conf.level = 0.9)

pcaplot(iris[,1:4], iris[,5],pcs=c(2,1),ep=1,
        auto.key=list(space="top", columns=3))
pcaplot(iris[,1:4], iris[,5],pcs=c(1,3,4))
tmp <- pcaplot(iris[,1:4], iris[,5],pcs=1:3,ep=2)
tmp

## change symbol's color, type and size
pcaplot(iris[,1:4], iris[,5],pcs=c(2,1),main="IRIS DATA", cex=1.2,
  auto.key=list(space="right", col=c("black","blue","red"), cex=1.2),
  par.settings = list(superpose.symbol = list(col=c("black","blue","red"),
                                              pch=c(1:3))))

## compare pcaplot and pca.plot. 
pcaplot(iris[,1:4], iris[,5],pcs=c(1,2),ep=2)
pca.plot(iris[,1:4], iris[,5], ep.plot = TRUE)

## an example of 'pca.comp'
pca.comp(iris[,1:4], scale = TRUE, pcs=1:3)

}

\keyword{plot}

