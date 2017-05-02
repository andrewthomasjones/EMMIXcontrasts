\name{wj2.permuted}
\alias{wj2.permuted}
\alias{pvalue.wire}
\title{The null distribution and p-value}
\description{
This function caculates the null distribution and p-values of the weighted contrast W_j.  
}
\usage{
wj2.permuted(data,ret,nB=99,contrast=NULL) 
pvalue.wire(wj,wj0)	
}
\arguments{
\item{data}{The dataset, an n by m numeric matrix, where n is number of observations and m the dimension of data}
\item{ret}{The return list of function emmixwire}
\item{nB}{The number of permutations}
\item{contrast}{A two- or three- dimensional vector the contrast(s) for the class differences}
\item{wj}{An n-dimensional vector containing the value of the statistic W_j for each gene}
\item{wj0}{An n by nB matrix with its columns as the statistic W_j for each permutation}
%%%\item{mu}{An adjustment number, default value is zero.}
}
\details{
The number of classes of samples is either two or three, and
the default contrast for two classes is c(1,-1), and three classes c(1,0,-1).
}

\value{
An n by nB matrix with its columns as the statistic Wj for each permutation.
}
\seealso{
\code{\link{emmixwire}},\code{\link{scores.wire}}.
}
\examples{

\dontrun{

dat <- read.table("GSE36703_37628_col.txt",header=FALSE,sep='\t')
rownames(dat) <- 1:nrow(dat)
set.seed(12345)
ret <-emmixwire(dat,g=3,ncov=3,nvcov=1,n1=5,n2=6,n3=3,
		debug=1,itmax=1000,epsilon=1e-5)

###calculate the W_j
wj <- scores.wire(ret,contrast=c(0.5,0.5,-1))

### the null distribution of W_j
wj0 <- wj2.permuted(dat,ret,nB=19)
### the p-values of W_j
pv  <- pvalue.wire(wj,wj0)

}
}
\keyword{cluster}
\keyword{datasets}

