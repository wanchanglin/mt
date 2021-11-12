% lwc-15-02-2010: 
%
\name{list.util}
\alias{list2df}
\alias{un.list}
\alias{shrink.list}

\title{
  List Manipulation Utilities 
}
\description{
  Functions to handle manipulation of list. 
}
\usage{
list2df(x)

un.list(x, y="")

shrink.list(x)
}

% --------------------------------------------------------------------
\arguments{
  \item{x}{A list to be manipulated. }
  \item{y}{A character or string of separator.}
}

\value{
  \code{list2df} returns a data frame.
  \code{un.list} returns a list.
	\code{shrink.list} retuns a list.
}

\details{
  \code{list2df} converts a list with components of vector to a data 
  frame. Shorter vectors will be filled with 
  NA. It is useful to convert rugged vectors into a data frame which can 
  be written to an Excel file.
  
  \code{un.list} collapses higher-depths list to 1-depth list. 
       This function uses recursive programming skill to tackle any depths 
       of list.

  \code{shrink.list} removes all NULL or NA entries from a list.
}

\seealso{
  \code{\link{feat.mfs}}
}

% ----------------------------------------------------------------------------
\author{
  Wanchang Lin 
}

% ----------------------------------------------------------------------
\examples{
 ## See examples of function feat.mfs for the usages of list2df and un.list.
a <- list(x=1, y=NA, z=NULL)
b <- list(x=1, y=NA)
c <- list(x=1, z=NULL)

shrink.list(a)
shrink.list(b)
shrink.list(c)

}

\keyword{util}
