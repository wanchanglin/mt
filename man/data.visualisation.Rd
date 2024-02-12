% lwc-20-10-2010:
% lwc-10-06-2015: add mds_plot_wrap
% lwc-15-07-2015: more explaination for 'ep'.
% lwc-12-02-2024: change dot as underscore to avoid:
%  S3 methods shown with full name in Rd file 'data.visualisation.Rd':
%    ‘lda.plot.wrap.1’
%  The \usage entries for S3 methods should use the \method markup and not
%  their full name.
\name{data.visualisation}
\alias{pca_plot_wrap}
\alias{mds_plot_wrap}
\alias{lda_plot_wrap}
\alias{lda_plot_wrap.1}
\alias{pls_plot_wrap}
\title{
  Grouped Data Visualisation by PCA, MDS, PCADA and PLSDA
}
\description{
  Grouped data visualisation by PCA, MDS, PCADA and PLSDA.
}
\usage{
  pca_plot_wrap(data.list,title="plotting",\dots) 

  mds_plot_wrap(data.list,method="euclidean",title="plotting",\dots) 

  lda_plot_wrap(data.list,title="plotting",\dots) 

  lda_plot_wrap.1(data.list,title="plotting",\dots) 

  pls_plot_wrap(data.list,title="plotting",\dots) 
}

\arguments{

  \item{data.list}{ A two-layer list structure, in which the second
    layer include a data frame called \code{dat} and a factor of class
    label called \code{cls}. Noticed that names of the first layer of
    \code{data.list} should be given. \code{data.list} can be produced
    by \code{\link{dat.sel}}. }

  \item{method}{ The distance measure to be used. This must be one of
    "euclidean", "maximum", "manhattan", "canberra", "binary" or
    "minkowski". Any unambiguous substring can be given. It is only for
    \code{mds_plot_wrap}. }

  \item{title}{A part of title string for plotting.}
  
  \item{\dots}{Further arguments to \code{lattice}. See corresponding
    entry in \code{\link{xyplot}} for non-trivial details of
    \code{lattice}. One argument is \code{ep}: an integer flag for
    ellipse. \code{1} and \code{2} for plotting overall and group
    ellipse, respectively. Otherwise, none. For details, see
    \code{\link{panel.elli.1}}. } }

\value{
  \code{mds_plot_wrap} returns a handle for MDS plot.

  All other four functions return a list with components: the first one
  is an object of class \code{"trellis"} for data visualisation; the
  second one is also an object of class \code{"trellis"} but plotting
  the corresponding variables, PCs (principal components), LDs (linear
  discrimniants) and LCs (latent components); and the third one is a
  matrix of these variables.

}

\note{

  There is a slight differences between \code{lda_plot_wrap.1} and
  \code{lda_plot_wrap}. The former plots the two-class grouped data,
  which has one linear discriminant (LD1), with strip plot. The later
  plots the two-class data by LD1 vs LD2 which is identical to LD1.
  Hence \code{lda_plot_wrap} is more general and can be applied to
  fusion of two and more class data sets.

}

%% -----------------------------------------------------------------------
\author{
  Wanchang Lin 
}


%% -----------------------------------------------------------------------
\seealso{
  \code{\link{pcaplot}}, \code{\link{mdsplot}}, \code{\link{plot.pcalda}}, 
  \code{\link{plot.plsc}}, \code{\link{dat.sel}}, \code{\link{grpplot}},
  \code{\link{panel.elli.1}}. 
}


%% ------------------------------------------------------------------------
\examples{
data(iris)
x <- subset(iris, select = -Species)
y <- iris$Species
## generate data list by dat.sel
iris.pw <- dat.sel(x,y,choices=NULL)
names(iris.pw)

pca.p <- pca_plot_wrap(iris.pw, ep=2)
pca.p[[1]]     ## visualised by PCA
pca.p[[2]]     ## plot PCA variables
pca.p[[3]]     ## matrix of PCA variables

mds.p <- mds_plot_wrap(iris.pw)
mds.p

pls.p  <- pls_plot_wrap(iris.pw)
pls.p[[1]]
pls.p[[2]]
pls.p[[3]]

lda.p <- lda_plot_wrap.1(iris.pw)
lda.p[[1]]
lda.p[[2]]
lda.p[[3]]
lda_plot_wrap(iris.pw)$lda.p

## only plot iris data
ph <- pca_plot_wrap(list(list(dat=x, cls=y)))$pca.p  
## Not given data names
ph
update(ph, strip=FALSE)       ## strip is an argument of lattice

tmp <- list(iris.dat=list(dat=x, cls=y))
pca_plot_wrap(tmp)$pca.p
pca_plot_wrap(tmp,strip=FALSE)$pca.p
pls_plot_wrap(tmp,strip=FALSE)$pls.p
lda_plot_wrap(tmp,strip=FALSE)$lda.p

data(abr1)
cls <- factor(abr1$fact$class)
dat <- preproc(abr1$pos, method="log")
## pair-wise data set
dat.pw <- dat.sel(dat, cls,choices=c("2","3","4"))

## add mult-class
idx <- grep("2|3|4",cls)
cls.234 <- factor(cls[idx])
dat.234 <- dat[idx,,drop = FALSE]

## combine all
dat.tmp <- c(dat.pw, 
             "2~3~4"=list(list(dat=dat.234,cls=cls.234)),
             all=list(list(dat=dat, cls=cls)))

## PCA
ph <- pca_plot_wrap(dat.tmp, title="abr1", par.strip.text = list(cex=0.75), 
                     scales=list(cex =.75,relation="free"), ep=2) 
## See function grpplot for usage of ep.
ph[[1]]
ph[[2]]                     

##PLSDA
ph <- pls_plot_wrap(dat.tmp, title="abr1", par.strip.text = list(cex=0.75), 
                     scales=list(cex =.75,relation="free"), ep=2) 
ph[[1]]
ph[[2]]                     

## PCADA
ph <- lda_plot_wrap(dat.tmp, title="abr1", par.strip.text = list(cex=0.75), 
                     scales=list(cex =.75,relation="free")) 
ph[[1]]
ph[[2]]    
}

\keyword{plot}

