\name{eq8.wire}
\alias{eq8.wire}
\title{Calculate the lamda values}
\description{
This function calculates the lamda values in equation 2.8.  
}
\usage{
eq8.wire(m,g,nb,X,W,U,V,sigma.e,sigma.b,sigma.c,nh,contrast)
}
\arguments{
\item{m}{The number of columns in the data}
\item{g}{The number of components in the mixture model}
\item{nb}{The number of columns in the design matrix X}
\item{X}{The design matrix X}
\item{W}{The design matrix W}
\item{U}{The design matrix U}
\item{V}{The design matrix V}
\item{sigma.e}{The covariance parameters of error item e}
\item{sigma.b}{The covariance parameters of random item b}
\item{sigma.c}{The covariance parameters of random item c}
\item{nh}{A vector with each element as the number of genes in each cluster}
\item{contrast}{The contrast vector, for example, c(1/2,1/2,-1) or c(1,0,-1), etc. }
}
\details{
The default contrast for two classes is c(1,-1), and three classes c(1/2,1/2,-1).
And it is called by function \code{scores.wire} and \code{wj2.permuted}.
}

\value{
A vector of lamda.
}
\seealso{
\code{\link{wj2.permuted}}}


\keyword{cluster}
