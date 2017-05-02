\name{scores.wire}
\alias{scores.wire}
\title{The weighted contrast}
\description{
This function caculates the weighted contrast W_j in order to find out the most significant dfferentially expressed genes.  
}
\usage{
scores.wire(obj,contrast=NULL)
}
\arguments{
\item{obj}{The return list of function emmixwire}
\item{contrast}{The vector of the specified contrast}
}
\details{
The number of classes of samples is either two or three.
}

\value{
The vector of the statistic Wj
}
\seealso{
\code{\link{wire.init.fit}},\code{\link{scores.wire}}.
}
\examples{

\dontrun{

dat <- read.table("GSE36703_37628_col.txt",header=FALSE,sep='\t')

rownames(dat) <- 1:nrow(dat)

###normalize the rows
x <- DoRows(dat)

set.seed(12345)

ret <-emmixwire(x,g=3,ncov=3,nvcov=1,n1=5,n2=6,n3=3,
		debug=1,itmax=1000,epsilon=1e-5)

###calculate the W_j
wj <- scores.wire(ret,contrast=c(0.5,0.5,-1))

}
}
\keyword{cluster}
\keyword{datasets}

