### this file contains all the functions for simulations 
### @ Christopher Rieser

# load libraries
library(igraph)
library(cccd)  # for graph weights
source(paste0(here::here(),"/R/utils.R"))


# error measurements
errs <- function(cov0,theta0,cov1,theta1,Lp,X,Z,Zcorr,edgesmat){
  
  # init
  p     = dim(cov0)[1]
  N     = dim(Lp)[1]
  q     = dim(theta0)[2]
  covih = gpower(cov0,-1/2)
  ci    = covih%*%covih
  mu0   = Zcorr %*% theta0   # Zcorr as used for the true edge outliers 
  mu1   = Zcorr %*% theta1
  Lpi   = gpower(Lp,-1)
  quant = qchisq(0.975,p)


  # err for theta 
  errtheta = p * norm(theta0-theta1,"F") /   norm(theta0,"F")

  # error of determinants
  d0      = det(cov0)

  # KL
  mur0 = Z %*% theta0 
  mur1  = Z %*% theta1
  errKL = sum(diag((mur1-mur0) %*% ci %*% t(mur1-mur0))) 
  errKL = errKL + q*(-p + log(d0 / det(cov1)) + sum(diag(ci %*% cov1)))
  errKL = 0.5 / q * errKL
  
  # prf scores
  thresh   = edgethresh(edgesmat,Lp, Lpi,quant)
  mhds     = edgemd2(X, mu0, cov0, edgesmat, Lp)
  mhdsEst  = edgemd2(X, mu1, cov1, edgesmat, Lp)
  outedges    =  edgesmat[which(mhds > thresh),]
  outedgesEst = edgesmat[which(mhdsEst > thresh),]
  prfs = prf(outedges, outedgesEst, N)
  
  
  return(c(errKL = errKL,errtheta = errtheta,prfs = prfs))
}

# function to calculate f scores
prf <- function(out.edges, out.edges.est, N){
  
  
  #adjacency matrices
  Wot = matrix(0,N,N); Woest = matrix(0,N,N)
  Wot[out.edges] = 1; Wot = Wot + t(Wot); diag(Wot) = 0; Wot[Wot!=0] =1
  Woest[out.edges.est] = 1; Woest = Woest + t(Woest); diag(Woest) = 0; Woest[Woest!=0] =1


  # calculate precision recall
  prec = if(sum(Wot)!= 0){ 
    sum(Woest * Wot) / sum(Wot) 
  }else if(sum(Wot)== 0 & sum(Woest)==0){ 1 }else{ 0 }
  rec  = if(sum(Woest)!=0){ 
    sum(Woest * Wot) / sum(Woest) 
  }else if(sum(Wot)== 0 & sum(Woest)==0){ 1 }else{ 0 }
  fs   = if(prec!=0 | rec!=0){ 
    2 * prec * rec / (prec + rec) 
  }else{ 0 }
  
  return(list(precision = prec, recall = rec, fscore = fs))
}



# function to generate random graph 
randgraph <- function(N, type = "knn"){
  
  # set some fixed parameters for random graphs
  rconn = sqrt(log(N) / (pi*N)) # fro geometric rg
  prob  = 0.05 # for erdos 
  sc1  = 1  # for scale free graph 
  sc2  = 2  # for scale free graph 
  nnbr = 5  # for knn graph 
  
  # knn graph 
  if(type == "knn"){
    x = runif(N,0,1)
    y = runif(N,0,1)
    locations = cbind(x,y)
    graphObj = nng(locations, k = nnbr,mutual = TRUE)
  }
  
  # geometric graph 
  if(type == "geom"){
    x = runif(N,0,1)
    y = runif(N,0,1)
    locations = cbind(x,y)
    W = (as.matrix(dist(locations)) <=  rconn)  
    W[W]     = 1
    diag(W)  = 0
    graphObj = graph_from_adjacency_matrix(W, mode="undirected")  
  }
  
  # erdos graph  
  if(type == "erdos"){ graphObj = erdos.renyi.game(N,prob,type ="gnp") }
  
  # scale free graph 
  if(type == "scalefree"){ graphObj = sample_pa(n = N,power = sc1, m = sc2) }
   
  
  # symmetrize
  edgesmat = as_edgelist(graphObj, names = FALSE)
  W = matrix(0,N,N); W[edgesmat] = 1; diag(W) = 0; W = W + t(W); W[W!=0] = 1
  
  # get L, edgelist and graph object
  graphObj = graph_from_adjacency_matrix(W, mode="undirected") 
  W        = as.matrix(as_adjacency_matrix(graphObj))
  
  # # generate random weights 
  weights  = rep(0,N*(N-1)/2)
  inds     = (W[lower.tri(W)]!=0)
  weights[inds]   = runif(sum(inds),0,1)
  W               = matrix(0,N,N)
  W[lower.tri(W)] = weights
  W = t(W)+W
  
  # laplacian
  L        = diag(rowSums(W)) - W
  edgesmat = as_edgelist(graphObj, names = FALSE)
  
  return(list(edgesmat = edgesmat, L = L,W = W, graphobj = graphObj))
}



# function to simulate outliers 
outlsim <- function(N, p, q , Lp, edgesmat, 
                    corrper, corrtype = "multiply", gamma = 1, zeta ){
  
  
  # function to get incident edges
  indcid.edges <- function(n.idx){
    aa = edgesmat[edgesmat[,1] == n.idx,]
    bb = edgesmat[edgesmat[,2] == n.idx,]
    return((rbind(aa,bb)))
  }
  # function to get incident nodes
  indcid.nodes <- function(n.idx){
    aa = edgesmat[edgesmat[,1] == n.idx,2]
    bb = edgesmat[edgesmat[,2] == n.idx,1]
    return(unique(c(aa,bb)))
  }
  # function to check if in edge set
  inedgeset <- function(edge,edges){
    edge.bool = (edge[1] == edges[,1] & edge[2] == edges[,2])
    edge.bool = edge.bool | (edge[2] == edges[,1] & edge[1] == edges[,2])
    return(any(edge.bool))
  }
  
  
  # init 
  Lpi2 = gpower(Lp,-1/2)
  Lpi  = Lpi2 %*% Lpi2
  mu   = matrix(0,N,p) 
  X    = matrix(0,N,p)
  nedges = dim(edgesmat)[1]
  cutoff = 0.975 # cut off for chi square quantile 

  # random center 
  Z        = matrix(runif(N*q,-1,1),nrow = N,ncol = q)
  theta    = matrix(rnorm(q*p,0,1),nrow = q, ncol = p)
  mu       = Z %*% theta
  
  # simulate random correlation matrix 
  S  = matrix(rnorm(p^2),p,p)   # covariance matrix
  US = eigen(S %*% t(S))$vectors
  S  = US %*% diag(runif(p,1,50)) %*% t(US)
  S  = cov2cor(S)
  
  # simulate data
  Si   = solve(S)
  Xcen = matrix(rnorm(N*p,0,1),N,p)  #data matrix
  Xcen = Xcen %*% gpower(S,1/2)
  Xcen = Lpi2 %*% Xcen
  X    = Xcen + mu
  
  # corrupt data 
  if(corrper == 0 ){
    
      # nothing to do here
      Xcorr    = X
      Zcorr    = Z 
      ind.corr = NA
      corr.edges = NA
      

  }else if(corrper > 0 & corrper <= 0.5){
    
     # init
     Xcorr  = X
     Zcorr  = Z
     thresh  = edgethresh(edgesmat, Lp, Lpi, qchisq(cutoff,p))


     # probability of choosing a node
     di             = diag(Lp)
     node.pr        = rep(1,N)
     #node.pr[di!=0] = 1/di[di!=0]
     node.pr        = node.pr / sum(node.pr)

     # get nodes to corrupt 
     maxloop = nedges #ceiling(nedges*0.2)
     edgecorrnbr= 0  # nbr of corrupted edges
     nodes.cnt  = 0  # nbr of corrupted nodes
     nodes      = 1:N # non-corrupted nodes
     nodes.corr = numeric(0)  
     counter    = 1 

     # generate outliers
     # spatial permutation outliers
     pcscores = prcomp(Xcen, center = FALSE)$x[,1] # (X-mu) %*% u1        # first pc scores
     pcorddec = order(pcscores, decreasing = TRUE)  # decreasing order of Pcs
     pcordinc = order(pcscores, decreasing = FALSE) # increasing order of Pcs

     # nbr of incident nodes per nodes changed
     nbrincedges = sapply(1:floor(N/2), function(idx){
                           incident.edges = lapply(c(pcorddec[1:idx],pcorddec[(N-idx+1):(N)]),indcid.edges)
                           if(length(incident.edges)==0){return(0)}
                           incident.edges = unique(do.call("rbind",incident.edges))
                           dim(incident.edges)[1]
                          })   
     smw  = (nbrincedges / nedges) > corrper
     if(all(smw)){ ncor = 1 }else{
        ncor = min(min(which((nbrincedges / nedges) > corrper)),floor(N/2))  
     }
     
     # corrupt data by exchanging quantiles
     upperidx = pcorddec[1:ncor]  # sample(1:N,ncor)#
     loweridx = pcordinc[1:ncor] # sample(setdiff(1:N,upperidx),ncor)#
     
     # corrupt center if regression case 
     nodes.corr          = unique(c(loweridx,upperidx))
     Zcorr[nodes.corr,]  = matrix(runif(q*length(nodes.corr),-zeta,zeta),length(nodes.corr),q)
     #mucorr[nodes.corr,] = Zcorr[nodes.corr,] %*% theta
     #Xcorr[nodes.corr,]  = Xcorr[nodes.corr,] + mucorr[nodes.corr,]
 
     
     ll = if(length(loweridx)==1){matrix(Xcorr[loweridx,],nrow = 1)}else{  Xcorr[loweridx,] }
     uu = if(length(upperidx)==1){matrix(Xcorr[upperidx,],nrow = 1)}else{  Xcorr[upperidx,] }
     Xcorr[loweridx,] =  uu[sample(1:length(loweridx),length(loweridx)),]
     Xcorr[upperidx,] =  ll[sample(1:length(upperidx),length(upperidx)),]
     
  }
  
  # outlying edges
  mhds     = edgemd2(Xcorr, mu, S, edgesmat, Lp)
  thresh   = edgethresh(edgesmat,Lp, Lpi, qchisq(cutoff,p))
  outedges = edgesmat[which(mhds > thresh),]

  return(list(Xcorr = Xcorr, X = X, S = S, mu = mu, outedges = outedges, Z = Z, Zcorr = Zcorr, theta = theta))
}



