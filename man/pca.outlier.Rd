%% lwc-01-09-2010
%% wll-29-11-2015: give more information about ellipse and outlier

\name{pca.outlier}
\alias{pca.outlier}
\alias{pca.outlier.1}

\title{ Outlier detection by PCA }

\description{
  Outlier detection by the Mahalanobis distances of PC1 and PC2. Also
  plot PC1 and PC2 with its confidence ellipse.
}

\usage{
  pca.outlier(x, center = TRUE, scale=TRUE,conf.level = 0.975,...) 

  pca.outlier.1(x, center = TRUE, scale=TRUE, conf.level = 0.975, 
              group=NULL, main = "PCA", cex=0.7,...)
}

% --------------------------------------------------------------------
\arguments{

  \item{x}{ A data frame or matrix. }

  \item{center}{ A logical value indicating whether the variables should
  	be shifted to be zero centred before PCA analysis takes place. }

  \item{scale}{ A logical value indicating whether the variables should
  	be scaled to have unit variance before PCA analysis takes place. }
   
  \item{conf.level}{ The confidence level for controlling the cutoff of
    the Mahalanobis distances. }

  \item{group}{ A string character or factor indicating group
  information of row of \code{x}. It is used only for plotting. }

  \item{main}{An overall title for PCA plot.}

  \item{cex}{ A numerical value giving the amount by which plotting text
    and symbols should be magnified relative to the default. }

  \item{\dots}{Further arguments for plotting}
  
}

% ----------------------------------------------------------------------------
\value{
  A list with components:

  \item{plot}{plot object of class \code{"trellis"} by
    \code{pca.outlier} only. }
  \item{outlier}{Outliers detected.}
  \item{conf.level}{Confidence level used.}
  \item{mah.dist}{Mahalanobis distances of each data sample.}
  \item{cutoff}{ Cutoff of Mahalanobis distances used for outlier detection.}
}

\note{ Examples of \code{\link{panel.elli}} and \code{\link{panel.outl}}
  give more general information about ellipses and outliers. If you
  ONLY want to plot outliers based on PCA in a general way, for
  example, outliers in different groups or in conditional panel, you can
  write an wrapper function combining with \code{\link{pca.comp}},
  \code{\link{panel.elli}} and \code{\link{panel.outl}}. It is quite
  similiar to the implementation of \code{\link{pca_plot_wrap}}.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------------
\seealso{\code{\link{pcaplot}}, \code{\link{grpplot}},
  \code{\link{panel.outl}},\code{\link{panel.elli}},
  \code{\link{pca_plot_wrap}} }


% ----------------------------------------------------------------------
\examples{
  data(iris)

  ## call lattice version
  pca.outlier(iris[,1:4], adj=-0.5)
  ## plot group
  pca.outlier(iris[,1:4], adj=-0.5,groups=iris[,5])
  ## more information about groups
  pca.outlier(iris[,1:4],groups=iris[,5],adj = -0.5, xlim=c(-5, 5),
                auto.key = list(x = .05, y = .9, corner = c(0, 0)),
                par.settings = list(superpose.symbol=list(pch=rep(1:25))))

  ## call basic graphic version
  pca.outlier.1(iris[,1:4])
  ## plot group
  pca.outlier.1(iris[,1:4], group=iris[,5])

}

\keyword{plot}
