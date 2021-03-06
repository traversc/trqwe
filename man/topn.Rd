% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/my_functions.r
\name{topn}
\alias{topn}
\title{Highest elements in a vector.}
\usage{
topn(x, n = 100, value = F, lowest = F)
}
\arguments{
\item{x}{A numeric vector.}

\item{n}{The number of top elements to return.}

\item{value}{If TRUE, returns the values of the top elements.  If FALSE, returns the indices.}

\item{lowest}{If TRUE, returns the lowest elements instead of the highest.}
}
\value{
A vector containing the indices or the values of the top elements.
}
\description{
Finds the top elements in a vector very quickly.  Equivalent ot \code{-sort(-x, partial=1:n)}
}
\examples{
naive_top <- function(x, n) {
   -sort(-x, partial=1:n)
}
x <- runif(1e7)
microbenchmark(naive_top(x,100), topn(x,100,value=T), times=10)

Unit: milliseconds
                   expr       min        lq     mean    median        uq
      naive_top(x, 100) 1070.0180 1071.5951 1075.964 1072.3520 1073.9989
topn(x, 100, value = T)  433.6682  433.8882  434.771  434.4986  435.6029
}
\seealso{
\url{http://stackoverflow.com/questions/18450778/}
}
