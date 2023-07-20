### this file contains all util functions 
### @ Christopher Rieser
library(igraph)
library(Rcpp)
library(RcppArmadillo)
library(abind)
library(robustbase)
library(ICSNP)
library(robustbase)


Sys.setenv(PATH="%PATH%;C:/Rtools/gcc-4.6.3/bin;c:/Rtools/bin")
sourceCpp(paste0(here::here(),"/src/mcdcov.cpp"))


# edgewise mcd function
# this is the main function
# X is the response observations as a data matrix
# Z is the covariate data matrix
# Lp is the Laplacian matrix
# edgesmat is a |E|x2 matrix with the edge indices in each row (no edge twice)
# type = const -> no covariates; type = "reg" regression case where Z is not only vector of ones
# typeinit  are initial estimates for covariances 
edgewiseDetMcd <- function(X, Z = NULL, Lp, edgesmat,
                           type = "const", typeinit = c("tan","spm","msc","sps"), alpha = 0.5, 
                           maxit = 20, maxerr = 1e-5){
  
  
  # function to turn corr matrix back into cov matrix 
  covrev <- function(A,ll,rr){ sweep(sweep(A,2,rr,"*"),1,ll,"*") }
  
  # init
  X  = as.matrix(X)
  Z  = as.matrix(Z)
  p = dim(X)[2]
  n = dim(X)[1]
  Li = gpower(Lp,-1)
  
  # replace edge matrix with one directonal one
  W = matrix(0,n,n); W[edgesmat] = 1; diag(W) = 0; W = (W+t(W)); W[W!=0] = 1
  graphObj = graph_from_adjacency_matrix(W, mode="undirected") 
  edgesmat = as_edgelist(graphObj, names = FALSE)
  
  # check 
  if(is.null(Z) & type=="const"){ message("Center is moot") }
  if(is.null(Z) & type=="reg"){  stop("Z must be given if type equal to reg") }
  
  # initial estimators 
  init  = mcdinit(X=X,Z=Z,Lp=Lp,edgesmat=edgesmat,typeinit=typeinit,type=type)
  ninit = dim(init$Scorrest)[3]
  
  # solve with edgewise mcd
  sol = detmcdBest(X,
                   Z,
                   init$Sest,
                   init$thest,
                   init$xscales,
                   dim(init$Sest)[3],
                   Lp,
                   edgesmat,
                   Li,
                   alpha,
                   maxit,
                   maxerr)

  
  # get outlying edges 
  inds = which(sol$outlIndicator == 1)
  if(length(inds)!=0){
    sol$outledges = edgesmat[inds,]
  }else{
    sol$outledges = NA
  }
  
  # compute center
  if(type=="const"){
    sol$mu = matrix(1,nrow=n) %*% apply(X,2,median)#ICSNP::spatial.median(X)
  }else if(type == "reg"){
    sol$mu = Z %*% sol$theta
  }

  
  return(sol)
}


# function to compute power of symmetric matrix 
# generalized power of A
gpower <- function(A, power = -1){
  
  # init tol
  tol   = sqrt(.Machine$double.eps)

  # eigen decomposition
  eig.decompt = eigen(A, symmetric = TRUE)
  eig.values  = eig.decompt$values
  thresh      =  max(tol * eig.values[1L], 0)
  positives   = eig.values > thresh

  ddiag = eig.values[positives]^power
  val   = eig.decompt$vectors[, positives] %*% (t(eig.decompt$vectors[, positives]) * ddiag)

  return(val)  
}



# function for the initial covariance estimator
covinit <- function(X, Y, type = "spm"){
  

  #init
  px = dim(X)[2]
  py = dim(Y)[2]
  n  = dim(X)[1]
  corOut = matrix(0,px,py)
  covOut = matrix(0,px,py)

  # set scaling function
  sc <- function(u){Qn(u)}

  # function to turn corr matrix back into cov matrix 
  covrev <- function(A,ll,rr){ 
    sweep(sweep(A,2,rr,"*"),1,ll,"*") 
    }
  
  # cov2cor
  cov2cor <- function(A){
    dd = diag(A)
    dai = dd; dai[dai!=0] = 1/sqrt(dai[dai!=0])
    return(sweep(sweep(A,2,dai,"*"),1,dai,"*") )
  }
  
  # function to turn calculate correlation matrix 
  corcust <- function(x,y){ 
      # compute correlation matrix 
      n = dim(x)[1]
      x = as.matrix(x); y = as.matrix(y)
      xc = colMeans(x)
      yc = colMeans(y)
      xs = apply(x,2,sd); xsi = xs; xsi[xsi!=0] = 1/xsi[xsi!=0]
      ys = apply(y,2,sd); ysi = ys; ysi[ysi!=0] = 1/ysi[ysi!=0]
      xcs = sweep(x,2,xc,"-")
      ycs = sweep(y,2,yc,"-")
      xcs = sweep(xcs,2,xsi,"*")
      ycs = sweep(ycs,2,ysi,"*")
      cm = 1/(n-1)*t(xcs)%*%ycs

      return(cm)
  }
  
  # columnwise scales
  Xscs  = apply(X, 2, function(u){ sc(u) })
  Yscs  = apply(Y, 2, function(u){ sc(u) })
  Xscsi = Xscs; Xscsi[Xscsi!=0] = 1/Xscsi[Xscsi!=0]
  Yscsi = Yscs; Yscsi[Yscsi!=0] = 1/Yscsi[Yscsi!=0]

  # scale
  Xsc    = sweep(X,2,Xscsi,"*")#scale(X, center = F, scale = Xscs)
  Ysc    = sweep(Y,2,Yscsi,"*")#scale(Y, center = F, scale = Yscs)
  
  # compute scores
  if(type == "tan"){
    
    Xz     = apply(Xsc, 2, function(u){ tanh(u) })
    Yz     = apply(Ysc, 2, function(u){ tanh(u) })
    corOut = corcust(Xz,Yz)

  }else if(type == "spm"){

    Xz     = apply(Xsc, 2, function(u){ uu = rank(u)})
    Yz     = apply(Ysc, 2, function(u){ uu = rank(u)})
    corOut = corcust(Xz,Yz)

  }else if(type == "msc"){
    
    Xz     = apply(Xsc, 2, function(u){ uu = (rank(u)-1/3)/(n+1/3); qnorm(uu,0,1)})
    Yz     = apply(Ysc, 2, function(u){ uu = (rank(u)-1/3)/(n+1/3); qnorm(uu,0,1) })
    corOut = corcust(Xz,Yz)

  }else if(type == "sps"){
    
    Xz     = t(apply(Xsc, 1, function(u){ if(all(u==0)){return(u)}else{ return(u / sqrt(sum(u^2)))} }))
    Yz     = t(apply(Ysc, 1, function(u){ if(all(u==0)){return(u)}else{ return(u / sqrt(sum(u^2)))}  }))
    corOut = corcust(Xz,Yz)

  }
      
  # adjust eigenvalues
  svddec  = svd(corOut)
  adjscx  = apply(t(t(svddec$u) %*% t(Xsc)), 2, function(u){ sc(u) })
  adjscy  = apply(t(t(svddec$v) %*% t(Ysc)), 2, function(u){ sc(u) })
  corOut  = sweep(svddec$u,2,adjscx,"*") %*% (sweep(t(svddec$v),1,adjscy,"*")) 
  if(all(dim(X)==dim(Y))){ corOut = cov2cor(corOut)}
  
  # turn into covariance matrix
  covOut = covrev(corOut, Xscs, Yscs)

  # return correlation matrix and scales 
  return(list(corOut = corOut,covOut = covOut,Xscs = Xscs, Yscs = Yscs))
}


# function to compute initial estimates 
mcdinit <- function(X,Z,Lp,edgesmat,typeinit = c("tan","spm","msc","sps"),type = "const"){
  
  # init 
  n = dim(X)[1]
  p = dim(X)[2]
  q = dim(Z)[2]
  nedges = dim(edgesmat)[1]
  
  
  # edge differences
  Xedges  = edgediff(X,edgesmat,Lp)
  Zedges  = edgediff(Z,edgesmat,Lp)
  
  # initial estimators 
  if(type == "const"){
     thest = array(0,dim = c(1,p,length(typeinit)))
     inits    = lapply(typeinit,function(tp){covinit(Xedges, Xedges, type = tp)})
     Sest     = abind( lapply(inits,function(init){init$covOut}), along = 3)
     xscales  = abind( lapply(inits,function(init){init$Xscs}), along = 2)
  }else if(type == "reg"){
    # compute initial theta
    initsth = lapply(typeinit,function(tp){
          zzi = gpower(covinit(Zedges, Zedges, type = tp)$covOut,-1) 
          zx  = covinit(Zedges, Xedges, type = tp)$covOut
          return(zzi %*% zx)
      })
    thest = abind( initsth , along=3 )

    # initial correlation matrices 
    inits = lapply(seq_along(typeinit),function(idx){
          thest   = thest[,,idx]; tp = typeinit[idx]
          Xcedges = edgediff(X -  Z %*% thest,edgesmat,Lp)
          covinit(Xcedges, Xcedges, type = tp)
    })
    Sest     = abind(lapply(inits,function(init){init$covOut}), along=3 )
    xscales  = abind(lapply(inits,function(init){init$Xscs})  , along=2)
  }

  return(list(Xedges = Xedges, Zedges = Zedges, Sest = Sest, thest = thest, xscales = xscales))
}


# function to calculate mcdL
mcdL <- function(X,Z,Lp){ 
  
    # compute robust fit 
    Lp2 = gpower(Lp,1/2)
    coeffs  = lapply(1:dim(X)[2],function(pidx){ robustbase::ltsReg(x = Lp2 %*% Z,y = Lp2 %*% X[,pidx], intercept = FALSE)$coefficients })
    coeffs  = do.call("cbind",coeffs)
    muest   = Z %*% coeffs
    mod     = robustbase::covMcd(giso(X-muest,Lp)) 
    res = list(Sigma = mod$cov, theta = coeffs)
  
  return(res)
}

mcdE <- function(X,Z,Lp,edgesmat){
      # init
      p = dim(X)[1]; N = dim(Lp)[1]

      # estimate mu
      Lp2 = gpower(Lp,1/2)
      coeffs  = lapply(1:dim(X)[2],function(pidx){ robustbase::ltsReg(x = Lp2 %*% Z,y =  Lp2 %*% X[,pidx], intercept = FALSE)$coefficients })
      coeffs  = do.call("cbind",coeffs)
      muest   = Z %*% coeffs
  
      # solve edgewise with MCD implementation
      Xedges    = edgediff(X-muest,edgesmat,Lp)
      soll      = DetMCD::DetMCD(Xedges,alpha = 0.5)
      res = list(Sigma = soll$cov, theta = coeffs)
      res
      
}


