\name{hedenlc}
\alias{hedenlc}
\docType{data}
\title{
The Hendenfalk data
}
\description{
These data are column normalized and used as an example.
}
\usage{data(hedenlc)}
\format{
This is a 3226 by 15 matrix. The first seven columns are from class I and last eight columns are from class II.  
}

\examples{
\dontrun{
data(hedenlc)
###
set.seed(123456)
obj<-emmixwire(hedenlc,g=5,ncov=3,nvcov=1,n1=7,n2=8,
      debug=1,itmax=1000,epsilon=1e-4) 
### to save the estimation results for later use
#save(obj,file="ret3.rbin")
###reload the results
#load("ret3.rbin")
### calculate the weighted contrasts W_j
Wj  <- scores.wire(obj)
}
}
\keyword{datasets}
