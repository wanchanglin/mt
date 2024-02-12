% lwc-14-02-2010:
% lwc-15-07-2015: major changes

\name{panel.elli}
\alias{panel.elli}
\alias{panel.outl}
\alias{panel.elli.1}
\title{ Panel Function for Plotting Ellipse and outlier }

\description{
  \pkg{lattice} panel functions for plotting grouped ellipse and outlier
}

\usage{
panel.elli(x, y, groups = NULL,conf.level = 0.975, \dots)
panel.elli.1(x, y, subscripts, groups=NULL, conf.level = 0.975,
             ep=0, com.grp=NULL, no.grp=NULL, ell.grp=NULL, \dots)
panel.outl(x, y, subscripts, groups=NULL, conf.level = 0.975, labs, \dots)
}

\arguments{
  \item{x, y}{Variables to be plotted.}
  \item{conf.level}{Confident level for ellipse}
  \item{groups, subscripts}{ Internal parameters for Lattice.}
  \item{labs}{Labels for potential outliers.}
  \item{ep}{ An integer for plotting ellipse. \code{1} and \code{2} for 
    plotting overall and group ellipse, respectively. Otherwise, none. }
  \item{com.grp}{A list of characters to select which combination of
    groups to be plotted.}
  \item{no.grp}{A list of characters to select which individual group
    not to be plotted. Note it will be overridden by \code{com.grp}.
    If no \code{com.grp} and \code{no.grp}, ellipses of each individual
    group  will be plotted.}
  \item{ell.grp}{Another categorical vector used for plotting ellipse.
    If provided, \code{ep},\code{com.grp} and \code{no.grp} will be
    ignored. It should be different from default \code{groups}, but has
    the same length of \code{groups}. For details, see example below.}
  \item{\dots}{Further arguments. See corresponding entry in
    \code{\link{xyplot}} for non-trivial details.  }
}

\value{
  Retuns objects of class \code{"trellis"}.
}

\details{
  \code{panel.elli} is modified from function
  \code{\link{panel.ellipse}} in package \pkg{latticeExtra}.

  \code{panel.elli.1} gives more control on how to plot ellipse for the
  current group. It also provides an option to plot ellipse based on
  another user-defined groups.

  \code{panel.outl} plots the labels of data points outside the ellipse.
  These data points can be treated as potential outliers.
}

\note{ \code{panel.elli.1} can be called by functions \code{grpplot},
  \code{pcaplot}, \code{mdsplot}, \code{pca_plot_wrap},
  \code{mds_plot_wrap}, \code{pls_plot_wrap} and \code{lda_plot_wrap} by
  passing argument of \code{ep}. See examples of these function for
  details. }

\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------------
\seealso{
  \code{\link{grpplot}}, \code{\link{pcaplot}}, \code{\link{mdsplot}}.
}

\examples{
library(lattice) 
data(iris)

## =====================================================================
## Examples of calling 'panel.elli' and 'panel.outl'
xyplot(Sepal.Length ~ Petal.Length, data = iris, groups=Species,
        par.settings = list(superpose.symbol = list(pch=c(15:17)),
                            superpose.line = list(lwd=2, lty=1:3)),
        panel = function(x, y, ...) {
            panel.xyplot(x, y, ...)
            panel.elli(x, y, ..., type="l",lwd=2)
            panel.outl(x,y, ...)
        },
        auto.key = list(x = .1, y = .8, corner = c(0, 0)),
        labs=rownames(iris), conf.level=0.9,adj = -0.5)

## Without groups
xyplot(Sepal.Length ~ Petal.Length, data = iris,
        par.settings = list(plot.symbol = list(cex = 1.1, pch=16)),
        panel = function(x, y, ...) {
            panel.xyplot(x, y, ...)
            panel.elli(x, y, ..., type="l", lwd = 2)
            panel.outl(x,y, ...)
        },
        auto.key = list(x = .1, y = .8, corner = c(0, 0)),
        labs=rownames(iris), conf.level=0.9,adj = -0.5)


## With conditioning
xyplot(Sepal.Length ~ Petal.Length|Species, data = iris, 
        par.settings = list(plot.symbol = list(cex = 1.1, pch=16)),
        layout=c(2,2),
        panel = function(x, y, ...) {
            panel.xyplot(x, y, ...)
            panel.elli(x, y, ..., type="l", lwd = 2)
            panel.outl(x,y, ...)
        },
        auto.key = list(x = .6, y = .8, corner = c(0, 0)),
        adj = 0,labs=rownames(iris), conf.level=0.95)

## =====================================================================
## Examples of 'panel.elli.1'
xyplot(Sepal.Length ~ Petal.Length, data = iris, groups=Species,
        ## ---------------------------------------------------
        ## Select what to be plotted.
        ep=2,
        ## com.grp = list(a="setosa",b=c("versicolor", "virginica")),
        ## no.grp = "setosa", ## Not draw ellipse for "setosa"
        ## ---------------------------------------------------
        par.settings = list(superpose.symbol = list(pch=c(15:17)),
                            superpose.line = list(lwd=2, lty=1:3)),
        panel = function(x, y, ...) {
            panel.xyplot(x, y, ...)
            panel.elli.1(x, y, ...)
            panel.outl(x,y, ...)
        },
        auto.key = list(points = TRUE, rectangles = FALSE, space = "right"),
        adj = 0,labs=rownames(iris), conf.level=0.95)

xyplot(Sepal.Length ~ Petal.Length, data = iris, groups=Species,
        ## ---------------------------------------------------
        ## Select what to be plotted.
        ep=2,
        ## com.grp = list(a="setosa",b=c("versicolor", "virginica")),
        no.grp = c("setosa","versicolor"),## Only draw "virginica"
        ## ---------------------------------------------------
        par.settings = list(superpose.symbol = list(pch=c(15:17)),
                            superpose.line = list(lwd=2, lty=1:3)),
        panel = function(x, y, ...) {
            panel.xyplot(x, y, ...)
            panel.elli.1(x, y, ...)
        },
        auto.key = list(x = .1, y = .8, corner = c(0, 0)))

xyplot(Sepal.Length ~ Petal.Length, data = iris, groups=Species,
        ## ---------------------------------------------------
        ## Select what to be plotted.
        ep=2,
        com.grp = list(a="setosa",b=c("versicolor", "virginica")),
        ## no.grp = "setosa", ## Not draw ellipse for "setosa"
        ## ---------------------------------------------------
        par.settings = list(superpose.symbol = list(pch=c(15:17)),
                            superpose.line = list(lwd=2, lty=1:3)),
        panel = function(x, y, ...) {
            panel.xyplot(x, y, ...)
            panel.elli.1(x, y, ...)
        },
        auto.key = list(x = .1, y = .8, corner = c(0, 0)))

  xyplot(Sepal.Length ~ Petal.Length, data = iris, groups=Species, ep=1,
        par.settings = list(superpose.symbol = list(pch=c(15:17)),
                            superpose.line = list(lwd=2, lty=1:3)),
        panel = function(x, y, ...) {
            panel.xyplot(x, y, ...)
            panel.elli.1(x, y, ...)
        },
        auto.key = list(points = TRUE, rectangles = FALSE, space = "right"))

## =====================================================================
## Another data set from package MASS
require(MASS)
data(Cars93)

## Plot ellipse based on original groups: DriveTrain
xyplot(Price~EngineSize, data=Cars93, groups=DriveTrain, ep=2,
        par.settings = list(superpose.symbol = list(pch=c(15:17)),
                            superpose.line = list(lwd=2, lty=1:3)),
        panel = function(x, y, ...) {
            panel.xyplot(x, y, ...)
            panel.elli.1(x, y, ...)
        },
        auto.key = list(points = TRUE, rectangles = FALSE, space = "right"))

## But we want to plot ellipse using AirBags
xyplot(Price~EngineSize, data=Cars93, groups=DriveTrain, 
        ell.grp=Cars93$AirBags,
        par.settings = list(superpose.symbol = list(pch=c(15:17)),
                            superpose.line = list(lwd=2, lty=1:3)),
        panel = function(x, y, ...) {
            panel.xyplot(x, y, ...)
            panel.elli.1(x, y, ...)
        },
        auto.key = list(points = TRUE, rectangles = FALSE, space = "right"))

}

\keyword{plot}
