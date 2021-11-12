% lwc-16-02-2010: commence
% lwc-02-09-2010: add corrgram
% lwc-14-09-2010: add cor.heat.gram.
\name{cor.util}
\alias{cor.cut}
\alias{cor.hcl}
\alias{cor.heat}
\alias{corrgram.circle}
\alias{corrgram.ellipse}
\alias{cor.heat.gram}
\alias{hm.cols}

\title{
  Correlation Analysis Utilities
}
\description{
  Functions to handle correlation analysis on data set.
}
\usage{
cor.cut(mat,cutoff=0.75,abs.f = FALSE,
        use = "pairwise.complete.obs", method = "pearson",\dots)
  
cor.hcl(mat, cutoff=0.75, use = "pairwise.complete.obs", 
        method = "pearson",fig.f=TRUE, hang=-1,
        horiz = FALSE, main = "Cluster Dendrogram", 
        ylab = ifelse(!horiz, "1 - correlation",""), 
        xlab = ifelse(horiz, "1 - correlation",""),\dots)
  
cor.heat(mat, use = "pairwise.complete.obs", method = "pearson",
         dend = c("right", "top", "none"),\dots)
  
corrgram.circle(co, 
                col.regions = colorRampPalette(c("red", "white", "blue")),
                scales = list(x = list(rot = 90)), \dots) 

corrgram.ellipse(co,label=FALSE,
                 col.regions = colorRampPalette(c("red", "white", "blue")),
                 scales = list(x = list(rot = 90)), \dots) 

cor.heat.gram(mat.1, mat.2, use = "pairwise.complete.obs", 
              method = "pearson", main="Heatmap of correlation", 
              cex=0.75, \dots)
              
hm.cols(low = "green", high = "red", n = 123)

}

% --------------------------------------------------------------------
\arguments{
  \item{mat, mat.1, mat.2}{
    A data frame or matrix. It should be noticed that \code{mat.1} and \code{mat.2}
    must have the same number of row.
   }
  \item{cutoff}{A scalar value of threshold.}
  \item{abs.f}{
   Logical flag indicating whether the absolute values should be used.
  }
  \item{fig.f}{
   Logical flag indicating whether the dendrogram of correlation matrix should be plotted.
  }
  \item{hang}{The fraction of the plot height by which labels should hang below 
    the rest of the plot. A negative value will cause the labels to hang down 
    from 0. See \code{\link{plot.hclust}}.}

  \item{horiz}{Logical indicating if the dendrogram should be drawn
    \emph{horizontally} or not.}

  \item{main, xlab, ylab}{Graphical parameters, see \code{\link{plot.default}}.}
  
	\item{dend}{Character string indicating whether to draw 'right', 'top' or 
	'none' dendrograms}.
  
	\item{use}{Argument for \code{\link{cor}}. An optional character string giving a method for computing 
    covariances in the presence of missing values.  This must be (an 
    abbreviation of) one of the strings \code{"everything"}, \code{"all.obs"}, 
    \code{"complete.obs"},  \code{"na.or.complete"}, or \code{"pairwise.complete.obs"}.}
  
  \item{method}{Argument for \code{\link{cor}}. A character string indicating which correlation coefficient (or 
    covariance) is to be computed.  One of  \code{"pearson"}, \code{"kendall"}, or \code{"spearman"},
    can be abbreviated.}

  \item{co}{Correlation matrix}
  \item{label}{ 
     A logical value indicating whether the correlation coefficient should be 
     plotted.
  }
  \item{\dots}{
    Additional parameters to \pkg{lattice}.
  }
  \item{col.regions}{Color vector to be used}
  \item{scales}{A list determining how the x- and y-axes (tick marks and labels) 
     are drawn. More details, see \code{\link{xyplot}}.
  }
  \item{cex}{A numeric multiplier to control character sizes for axis labels}.
  \item{low}{ Colour for low value}
  \item{high}{ Colour for high value}
  \item{n}{The number of colors (>= 1) to be in the palette} 
}

% ----------------------------------------------------------------------
\value{
   \code{cor.cut} returns a data frame with three columns, in which the first and second columns
   are variable names and their correlation (lager than cutoff) are
   given in the third column.
   
   \code{cor.hcl} returns a list with components of each cluster group and all correlation
   coefficients. 
   
   \code{cor.heat} returns an object of class "trellis".
   
   \code{corrgram.circle} returns an object of class "trellis".
   
   \code{corrgram.ellipse} returns an object of class "trellis".
   
   \code{cor.heat.gram} returns a list including the components:
    \itemize{ 
     \item \code{cor.heat}: An object of class  "trellis" for correlation heatmap 
       ordered by the hierarchical clustering.
     \item \code{cor.gram}: An object of class  "trellis" for corrgrams with 
      circle ordered by the hierarchical clustering.
     \item \code{cor.short}: A matrix of correlation coefficient in short format.
     \item \code{cor.long}: A matrix of correlation coefficient in long format.
    }
}

\details{
  \code{cor.cut} returns the pairs with correlation coefficient larger than \code{cutoff}.
  
  \code{cor.hcl} computes hierarchical cluster analysis based on correlation 
    coefficient. For other graphical parameters, see \code{\link{plot.dendrogram}}. 
  
  \code{cor.heat} display correlation heatmap using \pkg{lattice}.

  \code{corrgram.circle} and \code{corrgram.ellipse} display corrgrams with
  circle and ellipse. The functions are modified from codes given in 
  Deepayan Sarkar's \code{Lattice: Multivariate Data Visualization with R, 
  13.3.3 Corrgrams as customized level plots, pp:238-241}.
  
  \code{cor.heat.gram} handles the correlation of two data sets which have the 
  same row number. The best application is correlation between MS data 
  (metabolites) and meta/clinical data.
  
  \code{hm.cols} creates a vector of n contiguous colors for heat map.  
}

% ----------------------------------------------------------------------------
\references{
  Michael Friendly (2002). 
  \emph{Corrgrams: Exploratory displays for correlation matrices}.
  The American Statistician, 56, 316--324.

  D.J. Murdoch, E.D. Chow (1996).
  \emph{A graphical display of large correlation matrices}.
  The American Statistician, 50, 178--180.

  Deepayan Sarkar (2008).  
  \emph{Lattice: Multivariate Data Visualization with R}. Springer.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------
\examples{
data(iris)
cor.cut(iris[,1:4],cutoff=0.8, use="pairwise.complete.obs")
cor.hcl(iris[,1:4],cutoff=0.75,fig.f = TRUE)
ph <- cor.heat(iris[,1:4], dend="top")
ph
update(ph, scales = list(x = list(rot = 45)))

## change heatmap color scheme
cor.heat(iris[,1:4], dend="right", xlab="", ylab="",
  col.regions = colorRampPalette(c("green", "black", "red")))

## or use hm.cols
cor.heat(iris[,1:4], dend="right", xlab="", ylab="", col.regions = hm.cols())

## prepare data set
data(abr1)
cls <- factor(abr1$fact$class)
dat <- preproc(abr1$pos[,110:1930], method="log10")  

## feature selection
res <- fs.rf(dat,cls)
## take top 20 features
fs  <- res$fs.order[1:20]
## construct the data set for correlation analysis
mat <- dat[,fs]

cor.cut(mat,cutoff=0.9)
ch <- cor.hcl(mat,cutoff=0.75,fig.f = TRUE, xlab="Peaks")
## plot dendrogram horizontally with coloured labels.
ch <- cor.hcl(mat,cutoff=0.75,fig.f = TRUE, horiz=TRUE,center=TRUE, 
              nodePar = list(lab.cex = 0.6, lab.col = "forest green", pch = NA),
              xlim=c(2,0))

names(ch)
cor.heat(mat,dend="right")
cor.heat(mat,dend="right",col.regions = colorRampPalette(c("yellow", "red")))

## use corrgram with order by the hierarchical clustering
co  <- cor(mat, use="pairwise.complete.obs")
ord <- order.dendrogram(as.dendrogram(hclust(as.dist(1-co)))) 
corrgram.circle(co[ord,ord], main="Corrgrams with circle")
corrgram.ellipse(co[ord,ord], label = TRUE, main = "Corrgrams with circle",
                 col.regions = hm.cols())

## if without ordering
corrgram.circle(co, main="Corrgrams with circle")

## example of cor.heat.gram 
fs.1  <- res$fs.order[21:50]
mat.1 <- dat[,fs.1]

res.cor <- 
  cor.heat.gram(mat, mat.1, main="Heatmap of correlation between mat.1 and mat.2")
names(res.cor)
res.cor$cor.heat
res.cor$cor.gram

}

\keyword{util}
