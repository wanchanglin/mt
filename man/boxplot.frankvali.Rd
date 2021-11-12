% lwc-13-11-2006: First draft
%
\name{boxplot.frankvali}
\alias{boxplot.frankvali}
\title{
  Boxplot Method for Class 'frankvali'
}
\description{
  Boxplot method for error rate of each feature subset.
}

\usage{
\method{boxplot}{frankvali}(x,  \dots)
}

% ----------------------------------------------------------------------------
\arguments{
  \item{x}{
  An object of class \code{frankvali}.
  }
  \item{\dots}{
  Additional arguments to the plot, such as \code{main}, \code{xlab} and 
  \code{ylab}. 
  }
}

% ----------------------------------------------------------------------------
\details{
This function is a method for the generic function \code{boxplot()} for class 
\code{frankvali}. It plots the error rate of each feature subset.
}

\value{
Returns boxplot of class \code{frankvali}.
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

\seealso{
  \code{\link{frankvali}}
}

\examples{
data(abr1)
dat <- abr1$pos[,110:500]

x   <- preproc(dat, method="log10")  
y   <- factor(abr1$fact$class)        

dat <- dat.sel(x, y, choices=c("1","2"))
x.1 <- dat[[1]]$dat
y.1 <- dat[[1]]$cls

pars <- valipars(sampling="cv",niter=2,nreps=4)
res  <- frankvali(x.1,y.1,fs.method = "fs.rfe",fs.len = "power2", 
                  cl.method = "knn",pars = pars)
res
summary(res)
boxplot(res)                  
}

\keyword{plot}
