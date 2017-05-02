\name{wire.init.fit}
\alias{wire.init.fit}
\alias{wire.init.reg}
\alias{tau.estep.wire}
\title{Get the initial values}
\description{
The routinnes to  fit mixture models to the data and to obtain
initial partition or initial values. 
}
\usage{
wire.init.fit(dat,X,qe,n,m,g,nkmeans,nrandom=0)
wire.init.reg(dat,X,qe,n,m,g,cluster)
tau.estep.wire(dat,pro,mu,sigma,n,m,g)
}
\arguments{
\item{dat}{The dataset, an n by m numeric matrix, where n is number of observations and m the dimension of data.}
\item{X}{The design matrix}
\item{n}{The number of observations}
\item{m}{The number of variables}
\item{g}{The number of components in the mixture model}
\item{qe}{The number of columns of design matrix W}
\item{cluster}{A vector of integers specifying the initial partitions of the data}	
\item{nkmeans}{An integer to specify the number of KMEANS partitions to be used to find the best initial values.}
\item{nrandom}{An integer to specify the number of random partitions to be used to find the best initial values; the default value is 0.}
\item{pro}{A vector of mixing proportions pi}
\item{mu}{A numeric matrix with each column corresponding to the mean.}
\item{sigma}{The covaraince m by m by g array.}
}
\details{
These functions are called internally.
}
\value{
  \item{pro}{A vector of mixing proportions pi}
  \item{beta}{A numeric matrix with each column corresponding to the mean.}
  \item{sigma.e}{The covaraince of error}
  \item{cluster}{A vector of final partition}
  \item{loglik}{The loglikelihood at convergence}
  \item{lk}{A vector of loglikelihood at each EM iteration}
}
\seealso{\code{\link{emmixwire}}}

\keyword{cluster}
\keyword{datasets}
