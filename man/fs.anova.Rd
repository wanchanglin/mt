% lwc-21-03-2007

\name{fs.anova}
\alias{fs.anova}
\title{
  Feature Selection Using ANOVA
}
\description{
  Feature selection using ANOVA.
}
\usage{
  fs.anova(x,y,\dots)
}

% --------------------------------------------------------------------
\arguments{
  \item{x}{
  A data frame or matrix of data set. 
  }
  \item{y}{
  A factor or vector of class.
  }
  \item{\dots}{
  Arguments to pass to method.
  }
}


% ----------------------------------------------------------------------------
\value{
  A list with components:
  \item{fs.rank}{A vector of feature ranking scores.}
  \item{fs.order}{A vector of feature order from best to worst.}
  \item{stats}{A vector of statistics.}
  \item{pval}{A vector of p values.}
}


% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------
\examples{
## prepare data set
data(abr1)
cls <- factor(abr1$fact$class)
dat <- abr1$pos
## dat <- abr1$pos[,110:1930]

## fill zeros with NAs
dat <- mv.zene(dat)

## missing values summary
mv <- mv.stats(dat, grp=cls) 
mv    ## View the missing value pattern

## filter missing value variables
## dim(dat)
dat <- dat[,mv$mv.var < 0.15]
## dim(dat)

## fill NAs with mean
dat <- mv.fill(dat,method="mean")

## log transformation
dat <- preproc(dat, method="log10")

## select class "1" and "2" for feature ranking
ind <- grepl("1|2", cls)
mat <- dat[ind,,drop=FALSE] 
mat <- as.matrix(mat)
grp <- cls[ind, drop=TRUE]   

## apply ANOVA method for feature selection/ranking
res <- fs.anova(mat,grp)
names(res)

}

\keyword{classif}
