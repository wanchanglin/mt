\name{trainind}
\alias{trainind}
\title{ Generate Index of Training Samples }

\description{
  Generate index of training samples. The sampling scheme includes leave-one-out 
  cross-validation (\code{loocv}), cross-validation (\code{cv}), randomised 
  validation (\code{random}) and bootstrap (\code{boot}).
}

\usage{
trainind(cl, pars = valipars())
}
\arguments{
  \item{cl}{
  A factor or vector of class.
  }
  \item{pars}{A list of sampling parameters for generating training index. It has 
  the same structure as the output of \code{valipars}. See \code{valipars} for details.}
}

\value{
  Returns a list of training index. 
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

\seealso{\code{\link{valipars}}}

\examples{
## A trivia example
x <- as.factor(sample(c("a","b"), 20, replace=TRUE))
table(x)
pars <- valipars(sampling="rand", niter=2, nreps=4, strat=TRUE,div=2/3)
(temp <- trainind(x,pars=pars))
(tmp  <- temp[[1]])
x[tmp[[1]]];table(x[tmp[[1]]])     ## train idx
x[tmp[[2]]];table(x[tmp[[2]]])
x[tmp[[3]]];table(x[tmp[[3]]])
x[tmp[[4]]];table(x[tmp[[4]]])

x[-tmp[[1]]];table(x[-tmp[[1]]])   ## test idx
x[-tmp[[2]]];table(x[-tmp[[2]]])
x[-tmp[[3]]];table(x[-tmp[[3]]])
x[-tmp[[4]]];table(x[-tmp[[4]]])

# iris data set
data(iris)
dat <- subset(iris, select = -Species)
cl  <- iris$Species

## generate 5-fold cross-validation samples
cv.idx <- trainind(cl, pars = valipars(sampling="cv", niter=2, nreps=5))

## generate leave-one-out cross-validation samples
loocv.idx <- trainind(cl, pars = valipars(sampling = "loocv"))

## generate bootstrap samples with 25 replications
boot.idx <- trainind(cl, pars = valipars(sampling = "boot", niter=2,
                                           nreps=25))

## generate randomised samples with 1/4 division and 10 replications. 
rand.idx <- trainind(cl, pars = valipars(sampling = "rand", niter=2, 
                                           nreps=10, div = 1/4))


}

\keyword{manip}
