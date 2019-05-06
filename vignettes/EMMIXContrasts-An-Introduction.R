## ---- echo=FALSE, cache=TRUE---------------------------------------------
library(EMMIXcontrasts2)
library(Biobase)
knitr::opts_chunk$set(fig.width=8, fig.height=6) 
set.seed(123)

## ---- echo=FALSE, cache=TRUE---------------------------------------------
library(repmis)
repmis::source_data("https://github.com/andrewthomasjones/goldenspikeData/blob/master/goldenspike.Rdata?raw=true")

## ---- echo=FALSE, cache=TRUE---------------------------------------------
#library(erccdashboard)
#data(SEQC.Example)

## ---- echo=FALSE,include=FALSE, cache=TRUE-------------------------------
#add extra column to goldenspike make clearer which are null genes
goldenspike$isNull<-abs(goldenspike$fold)==1
gold2<-as.matrix(goldenspike[,seq_len(6)])

Emmix<-emmixwire(gold2,g='BIC', maxg=3, ncov=3,nvcov=3,n1=3,n2=3, debug=0,itmax=100,epsilon=1e-4)
M<-nrow(gold2)
scores<-abs(scores.wire(Emmix))
pTrue<-cumsum(goldenspike$isNull[order(scores, decreasing = TRUE)])/(seq_len(M))

## ---- echo=FALSE, fig.show='hold', cache=TRUE----------------------------
library(ggplot2)
library(reshape2)
combined <- data.frame(n = seq_len(M), pTrue = pTrue)
combined2<-melt(combined, id='n')
combined3<-subset(combined2, n<=1000)
p<-ggplot(data=combined3, aes(x=n, y=value, colour=variable)) + geom_line()
p<-p+scale_colour_discrete(guide=FALSE)+scale_x_continuous(name = "Gene Rank")+ scale_y_continuous(name = "FPR")+theme_classic() 
p<-p+ggtitle(paste("Golden Spike - false positive rate among top ranked genes"))
print(p)

## ---- echo=FALSE, cache=TRUE---------------------------------------------
library(vsn)
library(colonCA)
data(colonCA, package = 'colonCA')
alon_data<-(exprs(vsn2(colonCA)))

## ---- echo=FALSE, cache=TRUE---------------------------------------------
library(fission)
data("fission")

