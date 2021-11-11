% lwc-12-04-2007

\name{fs.relief}
\alias{fs.relief}
\title{
  Feature Selection Using RELIEF Method
}
\description{
  Feature selection using RELIEF method.
}
\usage{
  fs.relief(x,y, m=NULL, k=10, \dots)
}

% --------------------------------------------------------------------
\arguments{
  \item{x}{
  A data frame or matrix of data set. 
  }
  \item{y}{
  A factor or vector of class.
  }
  \item{m}{
  Number of instances to sample without replacement. Default is \code{NULL}
  which takes all instances for computation.
  }
  \item{k}{
  Number of nearest neighbours used to estimate feature relevance. 
  }
  \item{\dots}{
  Arguments to pass to method (current ignore).    
  }
}

\details{
This function implements the \bold{Relief} algorithm's extension called 
\bold{ReliefF}, which applies to multi-class problem and searches for \code{k} of its 
nearest neighbours from the same class, called \emph{hits}, and also \code{k} 
nearest neighbours from each of the different classes, called \emph{misses}. 
}

% ----------------------------------------------------------------------------
\value{
  A list with components:
  \item{fs.rank}{A vector of feature ranking scores.}
  \item{fs.order}{A vector of feature order from best to worst.}
  \item{stats}{A vector of measurements.}
}

% ----------------------------------------------------------------------------
\references{
  Kira, K. and Rendel, L. (1992). 
  The Feature Selection Problem: Traditional Methods and a new algorithm. 
  \emph{Proc. Tenth National Conference on Artificial Intelligence}, MIT Press, 
  129 - 134. 
  
  Kononenko, I., Simes, E., and Robnik-Sikonja, M. (1997). 
  Overcoming the Myopia of Induction Learning Algorithms with RELIEFF. 
  \emph{Applied Intelligence}, Vol.7, 1, 39-55. 

  Kononenko, I. (1994) Estimating Attributes: Analysis and Extensions of RELIEF,
  \emph{European Conference on Machine Learning}, Ed. Francesco Bergadano and
   Luc De Raedt, 171-182, Springer

  Robnik-Sikonja, M. and Kononenko, I. (2003) Theoretical and Empirical Analysis
    of ReliefF and RReliefF, \emph{Machine Learning}, 53, 23 - 69.
}


% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------
\examples{
data(iris)
x <- subset(iris, select = -Species)
y <- iris$Species

fs <- fs.relief(x, y, m=20,k=10)

}

\keyword{classif}
