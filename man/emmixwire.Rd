\name{emmixwire}
\alias{emmixwire}
\title{The EMMIX model with random effects}
\description{
The function emmixwire fits the data with the specified EMMIX-WIRE model.  
}
\usage{
emmixwire(dat,g=1,ncov=3,nvcov=1,n1=0,n2=0,n3=0,X=NULL,W=NULL,U=NULL,V=NULL,
cluster=NULL,init=NULL,debug=0,itmax=1000,epsilon=1e-5,nkmeans=5,nrandom=0)}
\arguments{
\item{dat}{The dataset, an n by m numeric matrix, where n is number of observations and m the dimension of data}
\item{g}{The number of components in the mixture model}
\item{ncov}{A small integer indicating the type of covariance structure of random item b }
\item{nvcov}{0 or 1, indicating whether or not to include random item c in the model}
\item{n1,n2,n3}{The number of samples in each class}
\item{X}{The design matrix X}
\item{W}{The design matrix W}
\item{U}{The design matrix U}
\item{V}{The design matrix V}
\item{cluster}{ A vector of integers specifying the initial partitions of the data; the default is NULL}	
\item{init}{A list containing the initial parameters for the mixture model. See details. The default value is NULL}	
\item{itmax}{ A big integer specifying the maximum number of iterations to apply; the default value is 1000}
\item{epsilon}{A small number used to stop the EM algorithm loop when the relative difference between log-likelihood at each iteration become sufficient small; the default value is 1e-5}
\item{nkmeans}{An integer to specify the number of KMEANS partitions to be used to find the best initial values}
\item{nrandom}{An integer to specify the number of random partitions to be used to find the best initial values}
\item{debug}{A logical value, if it is TRUE, the output will be printed out; FALSE silent; the default value is TRUE}
}
\details{

The combination of ncov and nvcov defines the covariance structure 
of the random effect item b and c respectively. For example,when both ncov and nvcov are zeros,
it is a mixtue of normal regression model. Specifically, when ncov=0,
random effect b is ignored; when ncov=1 or 2, all components share a common
covariance of random item b; when ncov=2, the common covariance matrix is diagonal;
when ncov=3 or 4, each component has their own covariance matrix of item b; 
when ncov=4, the covariance matrices of item b are diagonal;
when ncov=5, the covariance matrices of item b are sigma(h)*I(m)(
diagonal scale parameter matrix with same identical diagonal element values)

}

\value{
  \item{error}{Error code, 0 = normal exit;  1 = did not converge within \code{itmax} iterations}
  \item{loglik}{The log likelihood at convergence}
  \item{np}{The total number of parameters}
  \item{AIC}{Akaike Information Criterion (AIC) }
  \item{BIC}{Bayes Information Criterion (BIC)}
  \item{pro}{A vector of mixing proportions}
  \item{beta}{A numeric matrix with each column corresponding to the location parameter.}
  \item{sigma.e}{The covariance parameters of error item e}
  \item{sigma.b}{The covariance parameters of random item b}
  \item{sigma.c}{The covariance parameters of random item c}
  \item{cluster}{A vector of final partition}
  \item{eb}{The conditional expectation of random item b}
  \item{ec}{The conditional expectation of random item c}
  \item{lk}{A vector of log likelihood at each EM iteration}
  \item{tau}{An n by g matrix of posterior probability for each data point}
\item{X}{The design matrix X}
\item{W}{The design matrix W}
\item{U}{The design matrix U}
\item{V}{The design matrix V}
\item{g}{The number of components in the mixture model}

}
\seealso{
\code{\link{scores.wire}}}

\examples{

\dontrun{
data(expr.norm)
data(mapping.unique)

dat = expr.norm
map = mapping.unique
#map= 1: S/C=1, 2535;
#map=-1: "empties", 10131;
#map=-2: "mixed", ambiguous, 13, excluded;
#map> 1: 1331 nonnulls,  S/C > 1.

# in summary, 
#1331 (9.5 percent) nonnulls;
#and 
#12,666 (90.5 percent) true nulls.
#---------------------
dat = as.matrix(dat[map>=-1,])
map = map[map>=-1]
#---------------------

y  <- log(dat)
set.seed(123)
ret <- emmixwire(y,g=3,ncov=3,nvcov=1,n1=3,n2=3,n3=0,
                 debug=0,itmax=1000,epsilon=1e-5,nkmeans=5)

### alternatively, 
#X <- U <- cbind(c(1,1,1,0,0,0),c(0,0,0,1,1,1))
#m<-6   # m is number of columns
#V<-diag(m)
#W <-rep(1,m)
#ret <- emmixwire(y,g=3,ncov=3,nvcov=1,X=X,W=W,U=U,V=V,
#       debug=0,itmax=1000,epsilon=1e-5,nkmeans=5)

###calculate the weighted contrast W_j
wj <- scores.wire(ret)
names(wj) <- names(map)
###top 1000 genes
wire.1000 <- names(map)[order(abs(wj),decreasing=TRUE)][1:1000]
###the number of false non-nulls in the top 1000 genes 
sum(map[wire.1000]==1) + sum( map[wire.1000]==-1)
#119

##alternatively
### the null distribution of W_j
wj0 <- wj2.permuted(y,ret,nB=19)
pv  <- pvalue.wire(wj,wj0)
wire.1000 <- names(map)[order(pv,decreasing=0)][1:1000]
###the number of false non-nulls in the top 1000 genes 
sum(map[wire.1000]==1) + sum( map[wire.1000]==-1)
#119
hist(pv,50)

}


}

\keyword{cluster}
