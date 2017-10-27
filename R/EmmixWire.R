#imports
#'@importFrom stats coef kmeans lm resid
#'@importFrom MASS ginv
NULL

#E-step
.tau.estep.wire<-function(dat, pro, mu, sigma, n, m, g)
{
    #calls C++ code
    ret<-estep(dat, n, m, g, pro, mu, sigma)
    return(ret)        
}

#Sets up matrix M
.matM.wire<-function(B, C, tau, g, m)
{
    M<-array(0, dim=c(m, m, g))
    for(h in seq_len(g)){
        M[, , h]<-solve(B[, , h]+C[, , h]*sum(tau[, h]))
    }
    return(M)
}
#--------------------------------------------

#E-Step for B
.eb.estep<-function(tau, y, mu, DU, U, B, C, M, g, n, m, qb)
{
    eb1<-array(0, dim=c(n, qb, g))
    eb2<-array(0, c(g, qb, qb))
    invB<-array(0, c(m, m, g))
        
        for(h in seq_len(g))
        {
        invB[, , h]<-solve(B[, , h])
        # E(bi|y)
        eb1[, , h]<-t( DU[, , h]%*%t(U)%*%invB[, , h]%*%((t(y)-mu[, h])-
            c(C[, , h]%*%M[, , h]%*%colSums(t(t(y)-mu[, h])*tau[, h])) ))
        #E(bi%*%bi^t|y)
        eb2[h, , ]<-(sum(tau[, h])*DU[, , h]-(sum(tau[, h])-1)*
            DU[, , h]%*%t(U)%*%invB[, , h]%*%U%*%DU[, , h]-    
            DU[, , h]%*%t(U)%*%M[, , h]%*%U%*%DU[, , h])
        DU[, , h]<-(eb2[h, , ]+ t(eb1[, , h]*tau[, h])%*%eb1[, , h])/
            sum(tau[, h])
        }
    list(DU=DU, eb1=eb1, invB=invB)
}

#E-Step for C
.ec.estep<-function(tau, y, mu, sigma.c, V, M, g, n, qc)
{
    ec1<-array(0, c(qc, g))
    ec2<-kc<-rep(0, g)
    #
    for(h in seq_len(g))
    {
        ec1[, h]<-sigma.c[h]*t(V)%*%M[, , h]%*%colSums(t(t(y)-mu[, h])*tau[, h])
        kc[h] <-sigma.c[h]*qc-sum(diag(t(V)%*%M[, , h]%*%V))*
            (sigma.c[h]^2*sum(tau[, h]))
        ec2[h]<-c(t(ec1[, h])%*%ec1[, h]+kc[h])/qc
    }    
    #
    list(ec2=ec2, ec1=ec1)
}

#E-Step for Errors
.ee.estep<-function(y, mu, tau, U, V, W, A, invB, M, g, n, m, qe, dw, eb1, ec1)
{
    ae<-ee<-array(0, dim=c(n, m, g))
    ke<-rep(0, g)
    thet<-matrix(0, ncol=g, nrow=qe)
    mi<-diag(t(W)%*%W)
    for(h in seq_len(g))
    {
        ae[, , h]<-t(mu[, h]+U%*%t(eb1[, , h])+c(V%*%ec1[, h]))            
        ee[, , h]<- (y-ae[, , h])
        for(id in seq_len(ncol(W)))
        {
            AL<-A[, , h]*W[, id]
            ke[h]<-(sum(tau[, h])*sum(diag(AL))-sum(diag(t(AL)%*%M[, , h]%*%AL))
                -(sum(tau[, h])-1)*sum(diag( t(AL)%*%invB[, , h]%*%AL)))
            
            thet[id, h]<-( sum(c(t((t(ee[, , h])*W[, id])^2)*tau[, h]))+ke[h])/
                (mi[id]*sum(tau[, h]))
        }# end of loop
    } #end of loop
    list(sigma.e=thet, ae=ae)
}

#returns which cluster allocation (MAP) from tau matrix
.tau2cluster<-function(tau)
{
    apply(tau, FUN=which.max, MARGIN=1)
}

#returns covarainces for each group h (in g)
.getcov <-function(msigma, sumtau, n, m, g, ncov)
{
    sigma<-array(0, c(m, m))
    
    if( (ncov == 1)|(ncov == 2))
    {
        for(h in seq_len(g)){
            sigma<-sigma+sumtau[h]*msigma[, , h]
        }
        sigma<-as.matrix(sigma/n)
        
        if(ncov == 2){
            sigma<-diag(c(diag(sigma)), m)
            for(h in seq_len(g)){
                msigma[, , h]=sigma
            }
        }
    }
    
    if(m > 1)
    {
        if(ncov == 4){
            for(h in seq_len(g)){
                msigma[, , h]<-diag(c(diag(msigma[, , h])), m)
            }
        }
        
        if(ncov == 5){
            for(h in seq_len(g)){
                msigma[, , h]<-diag(sum(diag(msigma[, , h]))/m, m)
            }
        }
    }
    
    return(msigma)
}

# do EMMIX-WIRE analysis from initial values
.fit.emmix.wire<-function(dat, X, W, U, V, pro, beta, 
                    sigma.e, sigma.b, sigma.c, 
                    n, m, g, nb, qb, qc, qe, 
                    debug, ncov, nvcov, itmax, epsilon, log=TRUE)
{
    # controls
    flag<-error<-0
    lk<-rep(0, itmax)
    
    oldpro<-pro
    nbeta<-beta
    
    
    VV<-V%*%t(V)
    dw<-diag(t(W)%*%W)
    xxx<-ginv(t(X)%*%X)%*%t(X)
    mu<-matrix(0, ncol=g, nrow=m)
    
    #main EM loop
    
    for(i in seq_len(itmax))
    {
        mobj<-.mat.ABC.wire(U, VV, W, sigma.e, sigma.b, sigma.c, g, m)
        A<-mobj$A
        B<-mobj$B
        C<-mobj$C
        BC<-mobj$BC
        
        for(h in seq_len(g)) {
            mu[, h] <- as.vector(X%*%beta[, h])
        }
        
        # E-step
        eobj<-.tau.estep.wire(dat, oldpro, mu, BC, n, m, g)
        pro<-eobj$pro
        tau<-eobj$tau
        lk[i] <-eobj$loglik
        sumtau<- colSums(tau)
        M<-.matM.wire(B, C, tau, g, m)
        obj1 <- .eb.estep(tau, dat, mu, sigma.b, U, B, C, M, g, n, m, qb)
        obj2 <- .ec.estep(tau, dat, mu, sigma.c, V, M, g, n, qc)
        obj3 <- .ee.estep(dat, mu, tau, U, V, W, A, 
                obj1$invB, M, g, n, m, qe, dw, obj1$eb1, obj2$ec1)
        # M-step
        
        if(ncov > 0){
            sigma.b<-obj1$DU
        }else{
            for(h in seq_len(g)){
                sigma.b[, , h]<-diag(0, qb)
            }
        }
        
        if( (ncov > 0) & (ncov!=3) & (ncov!="AR") ){
            sigma.b <- .getcov(sigma.b, sumtau, n, qb, g, ncov)
        }
        
        if(nvcov > 0){
            sigma.c<-obj2$ec2
        }else{
            sigma.c<-rep(0, g)
        }
        
        sigma.e<-obj3$sigma.e
        
        #--------------------------------
        
        for(h in seq_len(g)){
            nbeta[, h]<-(beta[, h]+xxx%*%M[, , h]%*%A[, , h]%*%
                    colSums(t(t(dat)-mu[, h])*tau[, h])/sumtau[h])
        }
        #--------------------------------
        #
        loglik <- lk[i]
        
        if(debug){
            message('\n', i, 'th, loglik=', loglik)
        }
        
        beta   <- nbeta;
        oldpro <- pro
        
        if(i <= 10) next            
        if(log){
            if(abs(lk[i]-lk[i-10])<epsilon*abs(lk[i-10])) {flag<-1;break}
        } else {
            if(max(abs(c(nbeta-beta, pro-oldpro))) < epsilon) {flag<-1;break}
        }
        
    } # end of loop 
    
    if(flag == 0){
        error <- 1
    }
    
    #get the the final partition
    
    for(h in seq_len(g)) {
        mu[, h]<-as.vector(X%*%beta[, h]+V%*%obj2$ec1[, h])
    }
    
    eobj2<-.tau.estep.wire(dat, pro, mu, B, n, m, g)
    cluster<-.tau2cluster(eobj2$tau)
    
    # 
    # BIC & AIC
    nu<-switch(paste(ncov), 
            '0'=  0,            # without u
            '1'=  qb*(1+qb)/2,     #common covariance
            '2'=  qb,              #common diagonal covariance
            '3'=  g*(qb*(1+qb)/2), #general covariance
            '4'=  g*qb,            #general diagonal covariance
            '5'=  g  )            #sigma(h)*I_m
    
    nv<-ifelse(nvcov > 0, g, 0)
    np<-(g-1)+nb*g + nu + nv+ g*qe
    #
    BIC<--2*loglik+np*log(n)
    AIC<--2*loglik+np*2
    
    if(debug){
        message('\n', g, "BIC=", BIC, "AIC=", AIC, 
        "\nloglik=", loglik, "np=", np)
    }
    
    #return values
    
    ret <- list(error=error, loglik=loglik, np=np, 
                BIC=BIC, AIC=AIC, cluster=cluster, 
                pro=pro, beta=beta, sigma.e=  sigma.e)
    ret$lk=lk[lk!=0]
    ret$tau=eobj2$tau
    
    if(ncov == 1||ncov == 2||ncov == 3||ncov == 4||ncov == 5)
    {
        ret$sigma.b <- sigma.b
        ret$eb      <- obj1$eb1
    }
    
    if(nvcov > 0)
    {
        ret$sigma.c <- sigma.c
        ret$ec      <- obj2$ec1
    }
    #######################
    return(ret)
}

#'@name wire.init.fit
#'@title Get the initial values
#'@description The routinnes to  fit mixture models to the data and
#'to obtain  initial partition or initial values. 
#'@param dat The dataset,  an n by m numeric matrix,  
#'where n is number of observations and m the dimension of data.
#'@param X The design matrix.
#'@param n The number of observations.
#'@param m The number of variables.
#'@param g The number of components in the mixture model.
#'@param qe The number of columns of design matrix W.
#'@param nkmeans An integer to specify the number of KMEANS partitions
#'to be used to find the best initial values.
#'@param nrandom An integer to specify the number of random partitions
#'to be used to find the best initial values; the default value is 0.
#'@details These functions are called internally.
#'@return A list containing
#'\item{pro}{A vector of mixing proportions pi}
#'\item{beta}{A numeric matrix with each column corresponding to the mean.}
#'\item{sigma.e}{The covaraince of error}
#'\item{cluster}{A vector of final partition}
#'\item{loglik}{The loglikelihood at convergence}
#'\item{lk}{A vector of loglikelihood at each EM iteration}
#'@seealso \code{\link{emmixwire}}
#'@keywords cluster datasets
#'@export
wire.init.fit<-function(dat, X, qe, n, m, g, nkmeans, nrandom=0)
{
    wire.init.km<-function(dat, X, qe, n, m, g)
    {
        cluster<-rep(1, n)        
        if(g > 1){
            cluster<- kmeans(dat, g, nstart=5)$cluster
        }
        wire.init.reg(dat, X, qe, n, m, g, cluster)
    }
    
    wire.init.rd<-function(dat, X, qe, n, m, g)
    {
        cluster<-rep(1, n)        
        if(g > 1){
            cluster<- sample(seq_len(g), n, replace=TRUE)
        }
        wire.init.reg(dat, X, qe, n, m, g, cluster)
    }
    
    found<-NULL
    found$loglik<- -Inf
    
    if(nkmeans > 0) {
        for(j in seq_len(nkmeans))
        {    
            initobj<-try(wire.init.km(dat, X, qe, n, m, g))    
            if(class(initobj)!="try-error"){
                if(initobj$loglik > found$loglik){
                    found<-initobj
                }
            }
        }
        
    }
    if(nrandom > 0) {
        for(j in seq_len(nrandom))
        {
            initobj<-try(wire.init.rd(dat, X, qe, n, m, g))
            if(class(initobj)!="try-error"){
                if(initobj$loglik > found$loglik){
                    found<-initobj
                }
            }
        
        }
    }
    return(found)
}

#'@name wire.init.reg
#'@title Get the initial values
#'@description The routinnes to  fit mixture models to the data and
#'to obtain  initial partition or initial values. 
#'@param dat The dataset,  an n by m numeric matrix,
#'where n is number of observations and m the dimension of data.
#'@param X The design matrix.
#'@param n The number of observations.
#'@param m The number of variables.
#'@param g The number of components in the mixture model.
#'@param qe The number of columns of design matrix W.
#'@param cluster A vector of integers specifying the initial partitions 
#'of the data.
#'@details These functions are called internally.
#'@return A list containing
#'\item{pro}{A vector of mixing proportions pi}
#'\item{beta}{A numeric matrix with each column corresponding to the mean.}
#'\item{sigma.e}{The covaraince of error}
#'\item{cluster}{A vector of final partition}
#'\item{loglik}{The loglikelihood at convergence}
#'\item{lk}{A vector of loglikelihood at each EM iteration}
#'@seealso \code{\link{emmixwire}}
#'@keywords cluster datasets
#'@export
wire.init.reg<-function(dat, X, qe, n, m, g, cluster)
{
    # set x for regression
    xx<-as.matrix(X)
    beta<-matrix(0, nrow=ncol(xx), ncol=g)
    sigma<-pro<-rep(0, g)
    mu     <- array(0, c(m, g))
    msigma <- array(0, c(m, m, g))
    lk <- rep(0, g)
    for( ij in seq_len(g)){
        ni<-sum(cluster == ij)
        if(ni == 0){
            warning("empty cluster found!")
            next
        }
        nn     <- ni
        #pile up y
        y      <- c(t(dat[cluster == ij, ]))    
        if(nn > 1)
        {
            for(i in 2:nn){
                xx <- rbind(xx, X)
            }
        }

        obj<-lm(y~0+xx)
        beta[, ij]<-coef(obj)
        sigma[ij]<-sum((resid(obj))^2)/(nn*m)    
        pro[ij]<-ni/n
        xx<-as.matrix(X) # reset x for next iteration

    }
    
    # loglikelihood
    
    for(h in seq_len(g)) {
        mu[, h] <-c(X%*%beta[, h])
        msigma[, , h] <- diag(sigma[h], m)
        ni <- sum(cluster == h)
        if(ni == 0){
            warning("empty cluster found!")
            next
        }
    }     
    
    ooo <- .tau.estep.wire(dat, pro, mu, msigma, n, m, g)
    loglik <- ooo$loglik
    sigma.e<-matrix(0, ncol=g, nrow=qe)
    for(i in seq_len(qe)){
        sigma.e[i, ]<-sigma
    }
    return(list(beta=beta, sigma.e=sigma.e, pro=pro, loglik=loglik, lk=lk))
}


#'@title The EMMIX model with random effects
#'@description The function emmixwire fits the data
#'with the specified EMMIX-WIRE model as in [1].
#'@param dat The dataset,  an n by m numeric matrix, 
#'where n is number of observations and m the dimension of data.
#'@param g The number of components in the mixture model.
#'@param ncov A small integer indicating the type of covariance structure of 
#'item b .
#'@param nvcov 0 or 1,  indicating whether or not to include
#'random item c in the model.
#'@param n1 The number of samples in class 1.
#'@param n2 The number of samples in class 2.
#'@param n3 The number of samples in class 3.
#'@param X The design matrix X
#'@param W The design matrix W
#'@param U The design matrix U
#'@param V The design matrix V
#'@param cluster  A vector of integers specifying the initial 
#'partitions of the data;
#'the default is NULL.
#'@param init A list containing the initial parameters for the mixture model. 
#'See details. The default value is NULL.
#'@param itmax  A big integer specifying the maximum number of iterations to 
#'apply;
#'the default value is 1000.
#'@param epsilon A small number used to stop the EM algorithm loop
#'when the relative difference
#'between log-likelihood at each iteration become sufficient small; 
#'the default value is 1e-5.
#'@param nkmeans An integer to specify the number of KMEANS partitions
#'to be used to find the best initial values.
#'@param nrandom An integer to specify the number of random partitions
#'to be used to find the best initial values.
#'@param debug A logical value,  if it is TRUE,  the output will be printed out;
#'FALSE silent; the default value is FALSE.

#'@details  The combination of ncov and nvcov defines the covariance structure 
#'of 
#'the random effect item b and c respectively.#'For example, when both ncov
#'and nvcov are zeros, #'it is a mixtue of normal regression model. 
#'Specifically,  when ncov=0, random effect b is ignored; when ncov=1 or 2,
#'all components share a common covariance of random item b; when ncov=2,  
#'the common covariance matrix is diagonal;
#'when ncov=3 or 4,  each component has their own covariance matrix of item b; 
#'when ncov=4,  the covariance matrices of item b are diagonal;  when ncov=5,
#'the covariance matrices of item b are 
#'sigma(h)*I(m)(diagonal scale parameter matrix with same
#'identical diagonal element values).

#'@return A list containing the following:
#'\item{error}{Error code,  0=  normal exit;  1=  did not converge 
#'within \code{itmax} iterations}
#'\item{loglik}{The log likelihood at convergence}
#'\item{np}{The total number of parameters}
#'\item{AIC}{Akaike Information Criterion (AIC) }
#'\item{BIC}{Bayes Information Criterion (BIC)}
#'\item{pro}{A vector of mixing proportions}
#'\item{beta}{A numeric matrix with each column corresponding to 
#'the location parameter.}
#'\item{sigma.e}{The covariance parameters of error item e}
#'\item{sigma.b}{The covariance parameters of random item b}
#'\item{sigma.c}{The covariance parameters of random item c}
#'\item{cluster}{A vector of final partition}
#'\item{eb}{The conditional expectation of random item b}
#'\item{ec}{The conditional expectation of random item c}
#'\item{lk}{A vector of log likelihood at each EM iteration}
#'\item{tau}{An n by g matrix of posterior probability for each data point}
#'\item{X}{The design matrix X}
#'\item{W}{The design matrix W}
#'\item{U}{The design matrix U}
#'\item{V}{The design matrix V}
#'\item{g}{The number of components in the mixture model}

#'@seealso \code{\link{scores.wire}}
#'@references [1] Ng, S. K., McLachlan, G. J., Wang, K., Nagymanyoki, Z., Liu,
#'S., & Ng, S. W. (2014). Inference on differences between classes using
#'cluster-specific contrasts of mixed effects. Biostatistics, 16(1), 98-112.
#'@examples 
#'\dontrun{
#'data(expr.norm)
#'data(mapping.unique)
#'
#'dat=  expr.norm
#'map=  mapping.unique
#'#map= 1: S/C=1,  2535;
#'#map=-1: "empties",  10131;
#'#map=-2: "mixed",  ambiguous,  13,  excluded;
#'#map> 1: 1331 nonnulls,   S/C > 1.
#'
#'# in summary,  
#'#1331 (9.5 percent) nonnulls;
#'#and 
#'#12, 666 (90.5 percent) true nulls.
#'#---------------------
#'dat=  as.matrix(dat[map>=-1, ])
#'map=  map[map>=-1]
#'#---------------------
#'
#'y  <- log(dat)
#'set.seed(123)
#'ret <- emmixwire(y, g=3, ncov=3, nvcov=1, n1=3, n2=3, n3=0, 
#'         debug=0, itmax=1000, epsilon=1e-5, nkmeans=5)
#'
#'### alternatively,  
#'#X <- U <- cbind(c(1, 1, 1, 0, 0, 0), c(0, 0, 0, 1, 1, 1))
#'#m<-6   # m is number of columns
#'#V<-diag(m)
#'#W <-rep(1, m)
#'#ret <- emmixwire(y, g=3, ncov=3, nvcov=1, X=X, W=W, U=U, V=V, 
#'#    debug=0, itmax=1000, epsilon=1e-5, nkmeans=5)
#'
#'###calculate the weighted contrast W_j
#'wj <- scores.wire(ret)
#'names(wj) <- names(map)
#'###top 1000 genes
#'wire.1000 <- names(map)[order(abs(wj), decreasing=TRUE)][1:1000]
#'###the number of false non-nulls in the top 1000 genes 
#'sum(map[wire.1000] == 1) + sum( map[wire.1000] == -1)
#'#119
#'
#'##alternatively
#'### the null distribution of W_j
#'wj0 <- wj2.permuted(y, ret, nB=19)
#'pv  <- pvalue.wire(wj, wj0)
#'wire.1000 <- names(map)[order(pv, decreasing=0)][1:1000]
#'###the number of false non-nulls in the top 1000 genes 
#'sum(map[wire.1000] == 1) + sum( map[wire.1000] == -1)
#'#119
#'hist(pv, 50)
#'}


#'@export
emmixwire<-function(dat, g = 1, ncov = 3, nvcov = 0, n1 = 0, n2 = 0, n3 = 0, 
                    X = NULL, W = NULL, U = NULL, V = NULL, 
                    cluster = NULL, init = NULL, debug = FALSE, 
                    itmax = 1000, epsilon = 1e-5, nkmeans = 5, nrandom = 0)
{
    
    
    ###############    
    
    dat<-as.matrix(dat)
    n<-nrow(dat)
    m<-ncol(dat)
    
    #
    if(n1 > 0 && n2 > 0 ){
        if(n3 == 0){
            X<-U<-cbind(rep(c(1, 0), c(n1, n2)), rep(c(0, 1), c(n1, n2)))
    }else{
        X<-U<-cbind(rep(c(1, 0, 0), c(n1, n2, n3)), 
                    rep(c(0, 1, 0), c(n1, n2, n3)),
                    rep(c(0, 0, 1), c(n1, n2, n3)))
    }
        W <-rep(1, m)
        V<-diag(m)
    }else{
    # check the matrix U
        if(ncov == 0){
            U<- diag(m)
        }else{
            if(is.null(U)){
                U <- cbind(rep(1, m))
            }
        }

        # check the matrix V
        if(is.null(V) || nvcov == 0){
            V  <- diag(m)
        }
        # check the matrix W
        if(is.null(W)){
            W <- cbind(rep(1, m))
        }

    }
    
    #
    # check the matrix X
    if(is.null(X)){ stop("X must be specified")}
    
    X<-as.matrix(X);U<-as.matrix(U)
    V<-as.matrix(V);W<-as.matrix(W)
    qb<-ncol(U);qc<-ncol(V);qe<-ncol(W);nb<-ncol(X)
    ##################################
    #  some variables
    #
    tuv <- 0.2
    # initialize the sigma_b and sigma_c
    sigma.b<-array(0, c(qb, qb, g))
    for(h in seq_len(g)){
        if(qb > 1){
            diag(sigma.b[, , h])<-rep(tuv, qb)
        }else{
            sigma.b[1, 1, h] <- tuv
        }
    }
    sigma.c<-rep(tuv, g)
    
    message("initializing ...")
    #part 1: initial values
    if(!is.null(init)){
        found<-init
    } else {

        if(is.null(cluster))
        {
            found <- wire.init.fit(dat, X, qe, n, m, g, nkmeans, nrandom)
        }else{
            found <- wire.init.reg(dat, X, qe, n, m, g, cluster)
        }
    }
    
    if(length(found)<4){
        stop("not found inital values")
    }
    
    #part 2: call the main estimate procedure
    message("Fitting mixture model with EM algorithm ...")
    ret<-.fit.emmix.wire(dat, X, W, U, V, 
                found$pro, found$beta, found$sigma.e,
                sigma.b, sigma.c, 
                n, m, g, nb, qb, qc, qe, 
                debug, ncov, nvcov, itmax, epsilon)
    
    message(" done.")
    
    if(qb == m && (ncov == 4 || ncov == 2)){
        tmp <- array(0, c(m, g))
        for(h in seq_len(g)){
            tmp[, h] <- diag(ret$sigma.b[, , h])
        }
        ret$sigma.b <- tmp
    }
    
    if(ncov == 5){
        tmp <- rep(0, g)
        for(h in seq_len(g)){
            tmp[h] <- diag(ret$sigma.b[, , h])[1]
        }
        ret$sigma.b <- tmp
    }
    
    ret$g<-g
    ret$m<-m
    ret$nb<-nb
    ret$X<-X
    ret$W<-W
    ret$U<-U
    ret$V<-V
    
    return(ret)
}

#'@title Calculate the lambda values
#'@description This function calculates the lamda values
#'in equation 2.8.
#'@param m The number of columns in the data
#'@param g The number of components in the mixture model
#'@param nb The number of columns in the design matrix X
#'@param X The design matrix X
#'@param W The design matrix W
#'@param U The design matrix U
#'@param V The design matrix V
#'@param sigma.e The covariance parameters of error item e
#'@param sigma.b The covariance parameters of random item b
#'@param sigma.c The covariance parameters of random item c
#'@param nh A vector with each element as the number of genes in each cluster
#'@param contrast The contrast vector,  for example, 
#'c(1/2, 1/2, -1) or c(1, 0, -1),  etc. 
#'@details The default contrast for two classes is c(1, -1), 
#'and three classes c(1/2, 1/2, -1). 
#'And it is called by function \code{scores.wire} and \code{wj2.permuted}.
#'@return A vector of lamda values.
#'@seealso \code{\link{wj2.permuted}} \code{\link{scores.wire}}
#'@export
eq8.wire <-function(m, g, nb, X, W, U, V, 
                    sigma.e, sigma.b, sigma.c, nh, contrast)
{
    
    omega <- rep(0, g)
    if(is.null(contrast)){
        if(nb == 2){
            K1<-c(1, -1, rep(0, m))
            K2<-c(1, -1)
        }else{
            if(nb == 3){ 
                K1<-c(1, 0, -1, rep(0, m))
                K2<-c(1, 0, -1)
            }
        }
    }else{
        if(nb == 2) {
            if(length(contrast)!=2){
                stop("contrast should be a vector with length of 2")
            }
        K1<-c(contrast, rep(0, m))
        K2<-c(contrast)
        } else {
            if(nb == 3)
            {
                if(length(contrast)!=3){
                    stop("contrast should be a vector with length of 3")
                }
            K1<-c(contrast, rep(0, m))
            K2<-c(contrast)
            }
        }
    }
    XV <- cbind(X, V)

    for(h in seq_len(g))
    {
        # A
        if(ncol(W) > 1)
            A <- diag(1/c(W%*%c(sigma.e[, h])))
        else
            A <- diag(1/c(W[, 1]*sigma.e[1, h]))
        
        # B    
        B <- solve(sigma.b[, , h])
        
        # C is nvcov=  0 no C
        if(!is.null(sigma.c)){
            C <- (1/sigma.c[h]) * diag(ncol(V))
        }else{
            C <- 0 * diag(ncol(V))
        }
        # E
        E <- solve(t(U) %*% A %*% U + B)
        #XVAU
        XVAU <- t(XV)%*%(A%*%U)
        
        
        # P
        P <- t(XV) %*% A %*% XV
        
        P[2+seq_len(m), 2+seq_len(m)] <- P[2+seq_len(m), 2+seq_len(m)] + C
        
        P <- P - XVAU %*% E %*% t(XVAU)
        #P sometimes singular
        P <- ginv(P)/nh[h]
        
        #KOK
        XVAUEK <- XVAU %*% E %*% K2
        KOK <- (t(K1)-t(XVAUEK)) %*% P %*% K1
        KOK <- KOK - t(K1)%*% P %*%  XVAUEK + K2 %*% E %*% K2
        KOK <- KOK + t(XVAUEK)  %*% P %*%  XVAUEK
        omega[h] <- c(KOK)
    }
    sqrt(omega)    
}

#'@title The weighted contrast
#'@description This function caculates the weighted contrast W_j 
#'in order to find out the most significant dfferentially expressed genes.  
#'
#'@param obj The return list of function emmixwire.
#'@param contrast The vector of the specified contrast.
#'@param useZ use the latent variable allocation Z_i (default,  TRUE)
#'or the posterior probability (FALSE).
#'@details The number of classes of samples is either two or three.
#'@return The vector of the statistic Wj
#'@seealso \code{\link{wire.init.fit}} \code{\link{scores.wire}}
#'@examples
#'
#'\dontrun{
#'
#'dat <- read.table("GSE36703_37628_col.txt", header=FALSE, sep='\t')
#'
#'rownames(dat) <- seq_len(nrow(dat))
#'
#'###normalize the rows
#'x <- DoRows(dat)
#'
#'set.seed(12345)
#'
#'ret <-emmixwire(x, g = 3, ncov = 3, nvcov = 1, n1 = 5, n2 = 6, n3 = 3, 
#'             debug = 1, itmax = 1000, epsilon = 1e-5)
#'
#'###calculate the W_j
#'wj <- scores.wire(ret, contrast = c(0.5, 0.5, -1))
#'
#'}
#'
#'@keywords cluster datasets
#'@export
scores.wire <-function(obj,  contrast = NULL,  useZ = TRUE) 
{
    if(obj$nb !=2 && obj$nb !=3) {
        stop("only two or three classes can be compared.")
    }
    
    # get the sqrt KOK 
    ooo <- eq8.wire(obj$m, obj$g, obj$nb,  obj$X,  obj$W,  obj$U,  obj$V,  
            obj$sigma.e,  obj$sigma.b, 
            obj$sigma.c, colSums( obj$tau), contrast)
    # contrasts
    if(obj$nb == 2){
        if(is.null(contrast)){
            d1d2<-(t(( obj$eb[, 1, ]- obj$eb[, 2, ]))+ 
                    ( obj$beta[1, ]- obj$beta[2, ]))/ooo
        }else{ 
            d1d2<-(t(( obj$eb[,1,]*contrast[1]+ obj$eb[,2,] *contrast[2]))
                +(obj$beta[1,]*contrast[1] + obj$beta[2,]*contrast[2]))/ooo  
        }
    }else{ 

        if(is.null(contrast)){
            d1d2<-(t(( obj$eb[, 1, ]- obj$eb[, 3, ]))+
                    (obj$beta[1, ]- obj$beta[3, ]))/ooo
        }else{ 
            d1d2<-(t(( obj$eb[, 1, ] * contrast[1] + obj$eb[, 2, ]  *
                    contrast[2]  + obj$eb[, 3, ]  * contrast[3] ))
                    + ( obj$beta[1, ] * contrast[1] + obj$beta[2, ] *
                    contrast[2]  + obj$beta[3, ] * contrast[3] ))/ooo
        }
    }
    
    
    # equation 5
    #
    
    n<-length(obj$cluster)
    ret<-array(0, n)
    
    if(useZ){
        for(i in seq_len(n)){
            ret[i]<-t(d1d2)[i, obj$cluster[i]]
        }
    }else{
        ret<-c(rowSums(  obj$tau * t(d1d2)))
    }
    
    return(ret)
    
}

# permutation and null distribution
# B=99 permutations for class labels
# when calculate the W_j,  only re-do the numerator, 
# but keep the denominator same!
NULL
#'@title The null distribution
#'@description This function caculates the null distribution 
#'of the weighted contrast W_j. 
#'@param data The dataset,  an n by m numeric matrix,  
#'where n is number of observations and m the dimension of data
#'@param ret The return list of function emmixwire
#'@param nB The number of permutations
#'@param contrast A two- or three- dimensional vector the contrast(s) 
#'for the class differences
#'@param seed random seed for the permutations. 
#'@details The number of classes of samples is either two or three,  and 
#'the default contrast for two classes is c(1, -1),  
#'and three classes c(1, 0, -1).
#'@return An n by nB matrix with its columns as the statistic Wj
#'for each permutation.
#'@seealso \code{\link{emmixwire}} \code{\link{scores.wire}}.
#'@examples
#'\dontrun{
#'
#'dat <- read.table("GSE36703_37628_col.txt", header=FALSE, sep='\t')
#'rownames(dat) <- seq_len(nrow(dat))
#'set.seed(12345)
#'ret <-emmixwire(dat, g=3, ncov=3, nvcov=1, n1=5, n2=6, n3=3, 
#'         debug=1, itmax=1000, epsilon=1e-5)
#'
#'###calculate the W_j
#'wj <- scores.wire(ret, contrast=c(0.5, 0.5, -1))
#'
#'### the null distribution of W_j
#'wj0 <- wj2.permuted(dat, ret, nB=19)
#'### the p-values of W_j
#'pv  <- pvalue.wire(wj, wj0)
#'
#'}
#'@keywords cluster datasets
#'@export
wj2.permuted <- function(data, ret, nB=99, contrast=NULL,  seed=1234) {
    g  <- ret$g
    X<-as.matrix(ret$X)    
    U<-as.matrix(ret$U)
    V<-as.matrix(ret$V)
    W<-as.matrix(ret$W)
    
    qb<-ncol(U)
    qc<-ncol(V)
    qe<-ncol(W)
    
    VV<-V%*%t(V)
    
    data<-as.matrix(data)
    
    n<-nrow(data)
    m<-ncol(data)
    mu<-matrix(0, ncol=g, nrow=m)
    
    for(h in seq_len(g)){ 
        mu[, h] <- as.vector(X %*% ret$beta[, h])
    }

    mobj <- .mat.ABC.wire(U, VV, W, ret$sigma.e, ret$sigma.b, ret$sigma.c, g, m)
    A <- mobj$A
    B <- mobj$B
    C <- mobj$C
    BC<- mobj$BC    
    #-------------------------------------
    # get the sqrt KOK 
    
    ooo <- eq8.wire(m, g, ret$nb, X, W, U, V, ret$sigma.e, ret$sigma.b,  
        ret$sigma.c, colSums(ret$tau), contrast)
    
    ooo[is.na(ooo) ] <- 1e+10
    
    #-------------------------------------    
    set.seed(seed)
    wj0 <- array(0, c(n, nB))
    for(b in seq_len(nB)) { #do B permutation
        da <- data[, sample(seq_len(m), m, replace=FALSE)]    
        #   get tau for da        
        tau  <- .tau.estep.wire(da, ret$pro, mu, BC, n, m, g)$tau
        M <- .matM.wire(B, C, tau, g, m)
        ###########################################
        # update eb,             
        eb <- .eb.estep(tau, da, mu, ret$sigma.b, U, B, C, M, g, n, m, qb)$eb1
        # ec is fixed
        ec <- ret$ec
        
        ###########################################
        
        # get new tau        
        mu2    <- mu    
        for(h in seq_len(g)){
            if(!is.null(ec[, h])){
                mu2[, h]<- c(X%*%ret$beta[, h]+V%*%ec[, h])
            }else{
                mu2[, h]<- c(X%*%ret$beta[, h])
            }
        }
        tau    <- .tau.estep.wire(da, ret$pro, mu2, B, n, m, g)$tau        
        
        ###########################################
        
        # contrasts
        if(ret$nb == 2){
            if(is.null(contrast)){
                d1d2<-(t(( eb[,1,]- eb[,2,]))+ 
                    ( ret$beta[1,]- ret$beta[2,]))/ooo
            }else{ 
                d1d2<-(t(( eb[,1,]*contrast[1]+eb[,2,]*contrast[2]))+
                    (ret$beta[1,]*contrast[1]+ret$beta[2,]*contrast[2]))/ooo  
            }
        }else{ 
            if(is.null(contrast)){
                d1d2<-(t((eb[,1,]-eb[,3,]))+(ret$beta[1,]-ret$beta[3,]))/ooo
            }else{
                d1d2<-(t(( eb[,1,]*contrast[1]+eb[,2,]*contrast[2]+
                    eb[,3,]*contrast[3]))+(ret$beta[1,]*contrast[1]+
                    ret$beta[2,]*contrast[2]+ret$beta[3,]*contrast[3]))/ooo
            }
        }
        
        # equation 5
        wj0[, b] <- (rowSums(  tau * t(d1d2)))
        
    } #end of B permutations    
    return(wj0)
}
# permutation and null distribution
# B=99 permutations for class labels
# when calculate the W_j,  only re-do the numerator, 
# but keep the denominator same!
NULL
#'@title The p-values for the weighted contrast W_j. 
#'@description This function caculates the p-value of the weighted contrast,
#'W_j. 
#'@param wj The weighted contrast W_j,  output of \code{\link{scores.wire}}.
#'@param wj0 The null distribution for W_j,  
#'output of \code{\link{wj2.permuted}}.
#'@details The p-values relate to the tests set by the
#'contrast in \code{\link{scores.wire}}.
#'@return A vector of p-values.
#'@seealso \code{\link{emmixwire}} \code{\link{scores.wire}}
#'\code{\link{wj2.permuted}}.
#'@examples
#'\dontrun{
#'
#'dat <- read.table("GSE36703_37628_col.txt", header=FALSE, sep='\t')
#'rownames(dat) <- seq_len(nrow(dat))
#'set.seed(12345)
#'ret <-emmixwire(dat, g=3, ncov=3, nvcov=1, n1=5, n2=6, n3=3, 
#'         debug=1, itmax=1000, epsilon=1e-5)
#'
#'###calculate the W_j
#'wj <- scores.wire(ret, contrast=c(0.5, 0.5, -1))
#'
#'### the null distribution of W_j
#'wj0 <- wj2.permuted(dat, ret, nB=19)
#'### the p-values of W_j
#'pv  <- pvalue.wire(wj, wj0)
#'
#'}
#'@keywords cluster datasets
#'@export
pvalue.wire <- function(wj, wj0){    
    n0 <- length(wj)
    nn <- length(wj0)
    pv <- rep(0, n0)
    for(j in seq_len(n0)) {
        pv[j] <- sum( abs(c(wj, wj0)) >= abs(wj[j])) /(nn+n0)
    }
    return(pv)    
}


.mat.ABC.wire<-function(U, VV, W, sigma.e, DU, sigma.c, g, m)
{
    A<-B<-C<-BC<-array(0, dim=c(m, m, g))
    for(h in seq_len(g)){
        if(ncol(W) > 1){
            A[, , h]<-diag(as.vector(W%*%sigma.e[, h]))
        }else{
            A[, , h]<-diag(c(W[, 1]*sigma.e[1, h]))
        }
        B[, , h]<-A[, , h]+U%*%DU[, , h]%*%t(U)
        if(!is.null(sigma.c[h])){
            C[, , h]<-sigma.c[h]*VV
            BC[, , h]<-B[, , h]+C[, , h]
        }else{
            BC[, , h]<-B[, , h]
        }
    }
    return(list(A=A, B=B, C=C, BC=BC))
}
NULL
#'@title The goldenspike gene expression dataset
#'
#'@description A dataset containing the goldenspike data "10a.dat".
#'
#'@docType data
#'@keywords datasets
#'@name expr.norm
#'@usage data(expr.norm)
#'@source http://www2.ccr.buffalo.edu/halfon/spike/spikedownloads.html
#'@format A data frame with 14010 rows (genes) and 6 variables (samples).
NULL

#'@title The Hendenfalk breast cancer dataset.
#'@description These data are column normalized 
#'@docType data
#'@keywords datasets
#'@name hedenlc
#'@usage data(hedenlc)
#'@format This is a 3226 by 15 matrix. The first seven columns are from class I
#'and last eight columns are from class II.  
#'@examples
#'\dontrun{
#'data(hedenlc)
#'###
#'set.seed(123456)
#'obj<-emmixwire(hedenlc, g=5, ncov=3, nvcov=1, n1=7, n2=8, 
#'             debug=1, itmax=1000, epsilon=1e-4) 
#'### to save the estimation results for later use
#'#save(obj, file="ret3.rbin")
#'###reload the results
#'#load("ret3.rbin")
#'### calculate the weighted contrasts W_j
#'Wj  <- scores.wire(obj)
#'}
NULL
#'@name mapping.unique
#'@docType data
#'@title The map of goldenspike data
#'@description The map of genes to their S/C values
#'@usage data(mapping.unique)
#'@format It is a vector with names=  probe set names and values=
#'fold changes,  with "-1" depicting probe sets whose target RNAs were not 
#'spiked in ("empty"),  and "-2" marking "mixed" probe sets. In summary,   
#'1:  S/C=1,  2535;
#'-1:  "empties",  10131;
#'-2:  "mixed",  ambiguous,  13,  excluded;
#'>1: 1331 nonnulls.
#'@source  http://www2.ccr.buffalo.edu/halfon/spike/spikedownloads.html
#'@examples 
#'data(mapping.unique)
#'@keywords datasets
NULL
