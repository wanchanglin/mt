% lwc-19-12-2007:
% lwc-21-02-2010: change name as grpplot
% lwc-15-07-2015: remove 'ep' from function argument and add the link of
%                 panel.elli.1.

\name{grpplot}
\alias{grpplot}
\title{Plot Matrix-Like Object by Group}
\description{Plot matrix-like object by group}

\usage{
  grpplot(x, y, plot = "pairs", \dots)
}

\arguments{
  \item{x}{A matrix or data frame to be plotted. }
  \item{y}{A factor or vector giving group information of columns of \code{x}.}
  \item{plot}{
    One of plot types: \code{strip}, \code{box}, \code{density} and
    \code{pairs}. }
  \item{\dots}{Further arguments. See corresponding entry in
    \code{\link{xyplot}} for non-trivial details. One argument is
    \code{ep}: an integer for plotting ellipse. \code{1} and
    \code{2} for plotting overall and group ellipse, respectively.
    Otherwise, none. For details, see \code{\link{panel.elli.1}}. }
}

\value{An object of class \code{"trellis"}.}

% ----------------------------------------------------------------------------
\seealso{
  \code{\link{panel.elli.1}}, \code{\link{pcaplot}},
  \code{\link{pca_plot_wrap}},  \code{\link{lda_plot_wrap}}, 
  \code{\link{pls_plot_wrap}}.
}

\author{Wanchang Lin }

\examples{
  data(iris)
  grpplot(iris[,1:4], iris[,5],plot="strip", main="IRIS DATA")
  grpplot(iris[,1:4], iris[,5],plot="box", main="IRIS DATA")
  grpplot(iris[,1:4], iris[,5],plot="density", main="IRIS DATA")
  grpplot(iris[,1:4], iris[,5],plot="pairs",main="IRIS DATA",ep=2)

  ## returns an object of  class "trellis".
  tmp <- grpplot(iris[,c(2,1)], iris[,5],main="IRIS DATA",ep=2)
  tmp
  
  ## change symbol's color, type and size
  grpplot(iris[,c(2,1)], iris[,5],main="IRIS DATA", cex=1.5,
         auto.key=list(space="right", col=c("black","blue","red")),
         par.settings = list(superpose.symbol = list(col=c("black","blue","red"),
                                                     pch=c(1:3))))

}

\keyword{plot}
