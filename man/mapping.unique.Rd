\name{mapping.unique}
\alias{mapping.unique}
\docType{data}
\title{
The map of goldenspike data
}
\description{
The map of genes to their S/C values
}
\usage{data(mapping.unique)}
\details{
It is a vector with names = probe set names and values = fold changes, 
with "-1" depicting probe sets whose target RNAs were not spiked in ("empty"), 
 and "-2" marking "mixed" probe sets. In summary,  
 1:  S/C=1, 2535;
-1:  "empties", 10131;
-2:  "mixed", ambiguous, 13, excluded;
>1: 1331 nonnulls.
}
\source{
http://www2.ccr.buffalo.edu/halfon/spike/spikedownloads.html
}
\examples{
data(mapping.unique)
}
\keyword{datasets}
