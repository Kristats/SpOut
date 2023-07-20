library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
source(paste0(here::here(),"/R/utils.R"))
source(paste0(here::here(),"/R/simfuns.R"))


# run simulation for different graph size and growing N
nreps      = 250
graphtype  = as.data.frame(c("knn","scalefree","erdos","geom"))  # "erdos", "geom"
perouts    = as.data.frame(c(0,0.05,0.1,0.2,0.3)) # 0.2
ps         = c(3,10)
ns         = c(50,100,200,300)#seq(50,300,50)
nps        = expand.grid(ps,ns)
simparams  = merge(merge(perouts,graphtype),nps)
colnames(simparams) = c("perout","graphtype","p","N")


# set dimension for regression and outlyingness
q = 7
gamma = -1
zeta = 10

# start simulation loop
for(sidx in 1:dim(simparams)[1]){ 
  
  #init
  errors.big = array(0,dim = c(5,3,nreps))

  # set current paramters 
  p         = simparams[sidx,]$p
  N         = simparams[sidx,]$N
  corrtype  = simparams[sidx,]$corrtype
  perout    = simparams[sidx,]$perout
  graphtype = simparams[sidx,]$graphtype
  
  # print current setup
  print(paste0("loopidx:",sidx," p:",p," N:",N," perout: ",
               perout," grtype:",graphtype))
 
  # print index
  print(sidx)                                          
                                  
  # run simumlation
  scores = lapply(1:nreps, function(iii){
      # generate graph                                 
      rg = randgraph(N, type = graphtype)
      Lp = rg$L
      edgesmat = rg$edgesmat

      # simulate data
      outl   = outlsim(N = N, p = p,q = q, Lp = Lp, edgesmat = edgesmat, corrper = perout, corrtype = corrtype, gamma = gamma, zeta = zeta)
      Xcorr  = outl$Xcorr
      cov0   = outl$S  
      mu0    = outl$mu
      Z      = outl$Z
      Zcorr  = outl$Zcorr
      covi0  = solve(cov0) 
      Lgest = gpower(Lp,-1)
      theta0 = outl$theta
      
      # solve with different methods
      sol =  tryCatch(edgewiseDetMcd(Xcorr, Zcorr, Lp, edgesmat,type = "reg",typeinit = c("tan","spm","msc","sps"), 
                           alpha = 0.5, maxit = 40, maxerr = 1e-5), 
                           error =  function(e) return(NA),
                           warning = function(w) return(NA))
      mcdEd=  tryCatch(mcdE(Xcorr,Zcorr,Lp,edgesmat), 
                  error =  function(e) return(NA),
                  warning = function(w) return(NA))
      Sest = tryCatch({
                         thest = gpower(t(Zcorr) %*% Lp %*%  Zcorr,-1) %*% t(Zcorr) %*% Lp  %*% Xcorr
                         Sest  = 1/dim(Xcorr)[1]*t(Xcorr-Zcorr%*%thest) %*% Lp %*% (Xcorr-Zcorr%*%thest) 
                         list(Sigma = Sest, theta = thest)
                       },
                  error =  function(e) return(NA),
                  warning = function(w) return(NA))

      # scores
      errs1 = if(!any(is.na(sol))){errs(cov0,theta0, sol$Sigma,sol$theta,Lp,Xcorr,Z,Zcorr,edgesmat)  }else{ rep(NA,5) }
      errs2 = if(!any(is.na(mcdEd))){ errs(cov0,theta0, mcdEd$Sigma,mcdEd$theta,Lp,Xcorr,Z,Zcorr,edgesmat) }else{  rep(NA,5) }
      errs3  = if(!any(is.na(Sest))){ errs(cov0,theta0, Sest$Sigma,Sest$theta,Lp,Xcorr,Z,Zcorr,edgesmat) }else{  rep(NA,5) }
      err.df = cbind(unlist(errs1),unlist(errs2),unlist(errs3))
      colnames(err.df) = c("edgemcd","mcdEd","std")      

      # set precision/recall/Fscore to NA as no meaning in this setting
      if(perout == 0){ err.df[3:5,] = NA }
      
      errors.big[,,iii] = err.df                
      return(err.df)
  })

  # reformate scores
  errors.big = abind(scores, along = 3)
  colnames(errors.big) = c("edgemcd","mcdEd","std")
  rownames(errors.big) = c("kl","theta","prec","rec","fs")
  print(apply(errors.big,1:2,function(A){mean(A,na.rm = TRUE)}))

  # safe simulation result
  path = paste0(here::here(),"/simres/",sidx,"_",p,"_",N,"_",
               perout,"_",graphtype,".rds")
  saveRDS(errors.big,file = path)
  #readRDS(path)
}  


 

#### make plots
# read all files 
simPath  = paste0(getwd(),"/simres/")
allfiles = list.files(simPath)


# make data frame 
plotdf = lapply(c("kl","theta","prec","rec","fs"), function(errortype){
  
          plotdf.loc = lapply(allfiles, function(path){
          
          # gets the parameters of file 
          params = strsplit(path,"_")  
          p = params[[1]][2] ; N = params[[1]][3]; perout = params[[1]][4]; 
          graphtype = strsplit(as.character(params[[1]][5]),split = ".rds")[[1]][1]
          params = c(p, N, perout, graphtype)
        
          #get data
          df           = readRDS(paste0(simPath,path))
          df           = t(df[errortype,,]) #apply(df,1:2,function(A){mean(A)})
          colnames(df) = c("edgemcd","mcd","std")
          dfnames      = colnames(df)
        
          # return
          df = as.data.frame(df)
          df = cbind(df,do.call("rbind",lapply(1:nrow(df),function(.){params})), rep(errortype,nrow (df)))
          colnames(df) = c(dfnames,"p","N", "perout", "graphtype","errortype")
          return(df)
          
        })
        plotdf.loc = do.call("rbind",plotdf.loc)
        return(plotdf.loc)
})
plotdf = do.call("rbind",plotdf)
plotdf = as.data.frame(plotdf)

# name
plotdf        = tibble(plotdf)
plotdf$p      = as.numeric(plotdf$p)
plotdf$N      = as.numeric(plotdf$N)
plotdf$perout = as.numeric(plotdf$perout)


# delete for plot
plotdf = plotdf %>% filter(graphtype != "geom")  # omit permute outliers

## make all plot in pdf

#pdf("allplots.pdf")
iters = unique(plotdf$errortype)

for(iter.idx in 1:length(iters)){
  
  print(iter.idx)
  errortype = iters[iter.idx]

  if(errortype == "rec" | errortype == "prec"){ next }
  
  # filter data
  inds = (which(plotdf$errortype == errortype))
  plotdf.filt = plotdf[inds,] %>%  # outltype == "row", filter(graphtype != "erdos")
    group_by(p,perout,graphtype,N,errortype) %>%
    summarize(edgemcd = mean(edgemcd,na.rm=TRUE), 
              #mcdL = mean(mcdL,na.rm=TRUE),
              mcd  = mean(mcd,na.rm=TRUE),
              std = mean(std,na.rm=TRUE)) 
  plotdf.filt = plotdf.filt %>% gather(edgemcd,mcd,std,key="method",value = "mean.error")%>% 
    arrange((perout)) 
  plotdf.filt$perout = 100 * as.numeric(plotdf.filt$perout)
  
  if(errortype == "kl" | errortype == "fro" | errortype == "theta"){
    plotdf.filt$mean.error = log(plotdf.filt$mean.error)
  }
  
  
  #safe plot 
  png(filename = paste0(here::here(),"/plots/ ",errortype,corrtype,".png"),width = 6000, height = 6000,res = 500)
  
  
   if(errortype == "kl"){ gtt = ggtitle(paste("Log(KL)")) }
   if(errortype == "theta"){gtt = ggtitle(paste("Log(RD)")) }
   if(errortype == "fs"){ gtt = ggtitle(paste("F-score"))  }
              
  if(errortype == "fs"){ 
    pl =  ggplot(plotdf.filt %>% filter(perout != 0), aes(x = N, y = mean.error)) + 
    geom_line(aes(group = interaction(method,p),linetype = as.factor(p),colour = method),size = 1) + 
    ylab("Performance measure")+
    xlab("N") +
    facet_grid(graphtype ~ perout,scales = "free",
               labeller = labeller(perout = c("5" = "5%", "10" = "10%","20" = "20%","30" = "30%"))) +
    gtt + ylim(c(0,1))
  }else{
    pl =  ggplot(plotdf.filt, aes(x = N, y = mean.error)) + 
    geom_line(aes(group = interaction(method,p),linetype = as.factor(p),colour = method),size = 1) + 
    ylab("Performance measure")+
    xlab("N") +
    facet_grid(graphtype ~ perout,scales = "free",
               labeller = labeller(perout = c("0" = "0%","5" = "5%", "10" = "10%","20" = "20%","30" = "30%"))) +
    gtt 
    
  }
  
  pl = pl + theme(strip.text = element_text(size = 13, colour = "black")) + 
                scale_linetype_manual("Dimension",values = c("solid", "dashed")) +
                scale_colour_manual("Method", values = c("black", "green", "red")) +
                theme_linedraw(base_size = 25) + 
    theme(strip.background = element_rect(fill="grey"),
          strip.text = element_text(colour = "black"),
          legend.position="bottom", legend.key.size = unit(1.2, "cm")) 
              

  
  print(pl)
  dev.off()
  
}
#dev.off()


