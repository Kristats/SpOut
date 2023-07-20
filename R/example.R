# This file contains all the computations and plots for the example of French departements
# for the 2015 elections 
### @ Christopher Rieser

library(xtable)
library(patchwork)
library(ggplot2)
library(sf)
library(tidyr)
library(cowplot)
library(robCompositions)
library(dplyr)
library(MASS)
library(gridExtra)
source(paste0(here::here(),"/R/utils.R"))
source(paste0(here::here(),"/R/simfuns.R"))



# read data 
df = read_sf(paste0(here::here(),"/data/fulldf.shp"))
df = df[-c(46,78),]  # remove corsica
df = dplyr::select(df,-"age_min",-"ag_2540",-"ag_5564")  # remove certain ages 
df = df %>% mutate(centroid = st_centroid(geometry))

# log transform data - later on turned into ilr (!)
df = df  %>% mutate(ag_1839 = log(ag_1839), age_65 = log(age_65),ag_4064 = log(ag_4064))
df = df  %>% mutate(N_CAPBE = log(N_CAPBE), bac = log(bac), dplm_sp = log(dplm_sp))
df = df  %>% mutate(  AZ = log(AZ),BE = log(BE), FZ = log(FZ), GU = log(GU), OQ = log(OQ))
df = df  %>% mutate(left = log(left), right = log(right), others = log(others))

# use department data
centers = df$centroid
centers = data.frame(st_coordinates(centers))
N  = dim(df)[1]

# get edge set
edges = st_intersects(df, df,remove_self = T)
edgesmat = lapply(seq_len(length(edges)), function(i){ 
  cbind(rep(i,length(edges[[i]])),edges[[i]]) 
})
edgesmat = do.call("rbind",edgesmat)
edgesmat = edgesmat[order(edgesmat[,1]),]


# make weight matrix
W = matrix(0,N,N); W[edgesmat] = 1; diag(W) = 0; W = 0.5*(t(W) + W)
graphObj = graph_from_adjacency_matrix(W, mode="undirected") 
edgesmat = as_edgelist(graphObj, names = FALSE)

# Laplacian
Lp  = diag(rowSums(W))-W

### make basic plot
pl = ggplot() +
      geom_sf(data = df) +
      geom_point(data = centers,aes(x = X,y = Y),col = "red")

for(eidx in 1:dim(edgesmat)[1]){
  
  edge = edgesmat[eidx,]
  pl = pl +
    geom_line(data = centers[edge,],aes(x = X,y = Y))+
    xlab("")+
    ylab("")
  
}
pl = pl + theme(#panel.background = element_rect(fill = "white"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_line(color = "white", size = 0.8))

ggsave(paste0(here::here(),"/plots/","network.png"),pl)


# get data
#names
xnames = colnames(df)[5:7]
znames = colnames(df)[8:23]

# get data
X = as.matrix(df[,xnames])[,1:length(xnames)]
Z = as.matrix(df[,znames])[,1:length(znames)]
X = matrix(unlist(X),ncol = dim(X)[2])
Z = matrix(unlist(Z),ncol = dim(Z)[2])
colnames(X) = xnames
colnames(Z) = znames
p = dim(X)[2]
q = dim(Z)[2]

# scale Z
#Z = scale(Z)
Zilr = matrix(0,nrow = N,ncol = q-3)
ages = znames[1:3]
educs = znames[4:6]
secs  = znames[7:11]
Vzages  = orthbasis(length(ages))$V
Vzeduc = orthbasis(length(educs))$V
Vzsecs = orthbasis(length(secs))$V

# turn into clr
X = t(apply(X,1,function(u){u-mean(u)}))
Z[,ages] = t(apply(Z[,ages],1,function(u){u-mean(u)}))
Z[,educs] = t(apply(Z[,educs],1,function(u){u-mean(u)}))
Z[,secs] = t(apply(Z[,secs],1,function(u){u-mean(u)}))

# turn into ilr
Zilr = cbind(Z[,ages] %*% Vzages,
             Z[,educs] %*% Vzeduc,
             Z[,secs] %*% Vzsecs,
             Z[,12:q])

Vx   = orthbasis(p)$V
Xilr = as.matrix(X %*% Vx)

# plot 
#combn(1:dim(X)[2],2)
vars = rbind(1:3,1:3)
apply(vars,2,function(xidx){
        # each variable with another one at neighbor position   
        plotdf = t(apply(rbind(edgesmat,cbind(edgesmat[,2],edgesmat[,1])),1,function(edge){
                                       c(X[edge[1],xidx[1]],X[edge[2],xidx[2]])}))
        plotdf = as.data.frame(plotdf)
        colnames(plotdf) = c("position","neighbor")
        ggplot(plotdf,aes(x = position, y = neighbor))+ geom_point() + geom_smooth(method='lm', formula= y~x)
})


# solve compositional model
model = edgewiseDetMcd(X = Xilr, Z = Zilr, Lp, edgesmat,
                       type = "reg", typeinit =  c("tan","spm","msc","sps"), alpha = 0.5,
                       maxit = 50, maxerr = 1e-8)

# transform back 
model$Sigma = Vx %*% model$Sigma %*% t(Vx)
model$mu    = model$mu%*%t(Vx)
model$theta = rbind(Vzages %*% model$theta[1:2,] ,
                     Vzeduc %*%model$theta[3:4,],
                     Vzsecs %*% model$theta[5:8,],
                     model$theta[9:13,]) 
model$theta = model$theta %*% t(Vx)

# draw solution
# get outlying edges for given quant
Lgest        = gpower(Lp,-1)
prec         = gpower(model$Sigma,-1)
mhds.sol     = edgemd2(X, model$mu, model$Sigma, edgesmat, Lp)
thresh.sol   = edgethresh(edgesmat,Lp, Lgest, qchisq(0.975,p-1))
outedges.sol = edgesmat[which(mhds.sol > thresh.sol),]
mhds.sol.sc = mhds.sol / thresh.sol

# outlying rate 
dim(outedges.sol)[1] / dim(edgesmat)[1]
plot(sqrt(mhds.sol / thresh.sol), xlab = "Edge Index", ylab = "Edgewise Mahalanobis distance")
abline(h = sqrt(1))

# plot edges that are outlying 
pp = ggplot(data.frame(x = seq.int(dim(mhds.sol)[1]),y = sqrt(mhds.sol / edgethresh(edgesmat,Lp, Lgest, 1))),aes(x,y))+ 
              geom_point()+
              geom_hline(yintercept=sqrt(qchisq(0.975,p-1)))+
              ylab(expression(sqrt(hat(Delta)[ij]^clr~"/"~(w[ij]~(l[ii]+l[jj]-2~l[ij]))))) + xlab("Edge Index") +
              theme_gray(base_size = 16)

ggsave(paste0(here::here(),"/plots/","edgeoutls.png"),pp)

#plot (mhds.sol / thresh.sol) on edges


# align for node outliers 
nmhd2  = diag((X - model$mu) %*% prec %*% t(X - model$mu))
nodesc = median(nmhd2 / diag(Lgest) ) / qchisq(0.5,p-1)
nodemhd2 = sqrt(nmhd2 /  (nodesc * qchisq(0.975,p-1) * diag(Lgest)))

# make plot of France 
astgmap = theme(
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_line(color = "white", size = 0.6))
plout = ggplot()+ geom_sf(data = cbind(df,val = as.factor(nodemhd2 > 1)), 
                          colour = "black",
                          aes(fill = val))+
                  scale_colour_gradient(low = "white", high = "black",
                          guide = guide_colorbar(label = TRUE,
                                          draw.ulim = TRUE, 
                                          draw.llim = TRUE,
                                          # here comes the code change:
                                          frame.colour = "black",
                                          title.vjust = 1,
                                          title.hjust = 1.3))+
        scale_fill_manual(values = c("limegreen", "red"))+ 
        astgmap + ylab("")+xlab("")


# add edge outliers
for(eidx in 1:dim(outedges.sol)[1]){
  edge  = outedges.sol[eidx,]
  idx   = which((edge[1] == edgesmat[,1] & edge[2] == edgesmat[,2]) | (edge[1] == edgesmat[,2] & edge[2] == edgesmat[,1]))
  plout = plout + geom_line(data = data.frame(centers[edge,],delta =  sqrt(mhds.sol.sc[idx]*qchisq(0.975,p-1)) ),
                            aes(x = X,y = Y,color = delta), size = 0.6)
}
plout = plout + guides(fill = "none") + 
                labs(title = "Edgewise outliers\nDepartments of France", 
                     colour = expression(sqrt(hat(Delta)[ij]^clr~"/"~(w[ij]~(l[ii]+l[jj]-2~l[ij]))))) +
                theme(legend.position="bottom") 
# save
ggsave(paste0(here::here(),"/plots/","FranceEdgeOutl.png"),plout)

# make plot that zooms into Paris
ileidx = c(86,30,14,50,47,87,23)  # indices for ile de france
df$NOM_DEP[ileidx]  # ile de france 
ilenames = df$NOM_DEP[ileidx]
dd  = X[ileidx,]
rownames(dd) = df$NOM_DEP[ileidx]
dd
pplout = ggplot()+ geom_sf(data = cbind(df[ileidx,],val = as.factor(nodemhd2[ileidx] > 1)), 
                          colour = "black",
                          aes(fill = val))+
                  scale_colour_gradient(low = "white", high = "black",
                          guide = guide_colorbar(label = TRUE,
                                          draw.ulim = TRUE, 
                                          draw.llim = TRUE,
                                          # here comes the code change:
                                          frame.colour = "black",
                                          title.vjust = 1))+
                  scale_fill_manual(values = c("limegreen", "red"))+ 
                  astgmap + ylab("")+xlab("")

# outlying edges for paris 
poutedges.sol = outedges.sol[outedges.sol[,1] %in% ileidx & outedges.sol[,2] %in% ileidx,]

for(eidx in 1:dim(poutedges.sol)[1]){
  edge  = poutedges.sol[eidx,]
  idx   = which((edge[1] == edgesmat[,1] & edge[2] == edgesmat[,2]) | (edge[1] == edgesmat[,2] & edge[2] == edgesmat[,1]))
  pplout = pplout + geom_line(data = data.frame(centers[edge,],delta =  sqrt(mhds.sol.sc[idx]*qchisq(0.975,p-1)) ),
                            aes(x = X,y = Y,color = delta), size = 1)
}
pplout = pplout + guides(fill = "none",color = "none") + 
                  labs(title = "Edgewise outliers\nDepartments of Ile de France", colour = expression(tilde(Delta)[ij])) +
                  theme(legend.position="bottom")

# safe
ggsave(paste0(here::here(),"/plots/","ParisEdgeOutl.png"),pplout)

# safe both 
com = plout + pplout & theme(legend.position = "bottom")
com = com + plot_layout(guides = "collect")
ggsave(paste0(here::here(),"/plots/","FPEdgeOutl.png"),com)

# back to percentages
Xp = t(apply(exp(X),1,function(u){u / sum(u)}))


# function to plot departements log ratios
lrmaps <- function(df, Xp, idx, title = TRUE,...){
        xlr = log(Xp[,"left"] /Xp[,"right"]); xlo = log(Xp[,"left"]/Xp[,"others"]); xro = log(Xp[,"right"]/Xp[,"others"])
        psc = max(abs(xlr[idx]),abs(xlo[idx]),abs(xro[idx]))
        plr = ggplot()+ geom_sf(data = cbind(df[idx,], val = xlr[idx]), colour = "black",aes(fill = val)) +
                          scale_fill_gradient2(limits = c(-psc,psc))+
                          geom_sf_label(data = df[idx,], aes(label = NOM_DEP), color = "black", ...)+ 
                          astgmap + ylab("")+xlab("")+ guides(fill="none")
        plo = ggplot()+ geom_sf(data = cbind(df[idx,], val = xlo[idx]), colour = "black",aes(fill = val))+
                          scale_fill_gradient2(limits = c(-psc,psc))+
                          geom_sf_label(data = df[idx,], aes(label = NOM_DEP), color = "black", ...)+ 
                          astgmap + ylab("")+xlab("")+ guides(fill="none")
        pro = ggplot()+ geom_sf(data = cbind(df[idx,], val = xro[idx]), colour = "black",aes(fill = val))+
                          scale_fill_gradient2(limits = c(-psc,psc))+
                          geom_sf_label(data = df[idx,], aes(label = NOM_DEP), color = "black", ...)+ 
                          astgmap + ylab("")+xlab("")+ guides(fill="none")
        if(title == TRUE){
            plr = plr + labs(title = "log(left/right)")
            plo = plo + labs(title = "log(left/others)")
            pro = pro + labs(title = "log(right/others)") 
        }
        
        com = plr + plo + pro & theme(legend.position = "bottom")
        com = com + plot_layout(guides = "collect")
        return(com)
}


# plot for log ratios of voter data
com = lrmaps(df, Xp, seq.int(dim(Xp)[1]), title = TRUE, alpha = 0, size = 0)

# for Lot, Correze, Cantal show data
lccidx = c(71,38,16)
lccidx = c(lccidx,11, 27,39,59,18,77,32,74,36) # adding also surrinding departements
df$NOM_DEP[lccidx]
xsub = Xp[lccidx,]
rownames(xsub) = df$NOM_DEP[lccidx]
zsub = as.matrix(Z[lccidx,znames])[,seq.int(znames)]
rownames(zsub) = df$NOM_DEP[lccidx]
print(xsub)
print(zsub)

# plot
comlcc = lrmaps(df, Xp, lccidx,  title = TRUE, alpha = 0.3, size = 1.5,label.size =NA)
comlcc = comlcc + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ggsave(paste0(here::here(),"/plots/","lrlcc.png"),comlcc)


# for Arieges and Pyrenees-Orientales
apidx = c(22,46)
apidx = c(apidx,49,72) # ass also surrounding departements
df$NOM_DEP[apidx]
xsub = Xp[apidx,]
rownames(xsub) = df$NOM_DEP[apidx]
zsub = as.matrix(Z[apidx,znames])[,seq.int(znames)]
rownames(zsub) = df$NOM_DEP[apidx]
print(xsub)
print(zsub)


comap = lrmaps(df, Xp, apidx, title = FALSE, alpha = 0.3, size = 2,label.size =NA)
comap = comap + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ggsave(paste0(here::here(),"/plots/","lrap.png"),comap)


# for Seine-Saint-Denis, Hauts-de-Seine, 
# Seine-et-Marne, Val-de-Marne and Val d’Oise
shsvvidx = c(86,30,14,50,47,87,23)
df$NOM_DEP[shsvvidx]
xsub = Xp[shsvvidx,]
rownames(xsub) = df$NOM_DEP[shsvvidx]
zsub = as.matrix(Z[shsvvidx,znames])[,seq.int(znames)]
rownames(zsub) = df$NOM_DEP[shsvvidx]
print(xsub)
print(zsub)

comshsvv = lrmaps(df, Xp, shsvvidx, title = FALSE, alpha = 0.3, size = 1.5,label.size =NA) 
comshsvv = comshsvv + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ggsave(paste0(here::here(),"/plots/","lrparis.png"),comshsvv)


comall = plot_grid(comlcc, comap, comshsvv, nrow = 3, align = "hv")
ggsave(paste0(here::here(),"/plots/","lrall.png"),comall)


# table of coefficients 
sds    = apply(Z,2,function(u){sd(u)}) #  scaling 
coeffs = sweep(model$theta,1,sds,"*")
coeffs = sapply(1:p, function(pidx){ coeffs[,pidx]})
coeffs = round(coeffs,4)
rownames(coeffs) =  colnames(Z)
colnames(coeffs) = c("left", "right", "others")
xtable(t(coeffs), type = "latex")

# or
xtable(t(coeffs)[,1:10], type = "latex")
xtable(t(coeffs)[,11:16], type = "latex")


# plot diagnostics 
resids = gpower(Lp,1/2) %*% (X-model$mu) %*% gpower(model$Sigma,-1/2)
rc1 = max(abs(resids[,1]))
rc2 = max(abs(resids[,2]))
rc3 = max(abs(resids[,3]))
rc = max(rc1,rc2,rc3); rc1 = rc; rc2 = rc; rc3 = rc  # if same scale for all 
pr1 = ggplot() + geom_sf(data = cbind(df,resids = resids[,1]), aes(fill = resids)) +
                 scale_fill_gradient2(midpoint = 0,  limits = c(-rc1, rc1)) +
                 astgmap + ylab("")+xlab("") + labs(title = "Residual Left") + 
                 theme(plot.title = element_text(size=17))+
                 guides(fill = "none")
pr2 = ggplot() + geom_sf(data = cbind(df,resids = resids[,2]), aes(fill = resids)) +
                 scale_fill_gradient2(midpoint = 0,  limits = c(-rc2, rc2))+
                 astgmap + ylab("")+xlab("")+ labs(title = "Residual Right")+ 
                 theme(plot.title = element_text(size=17))+
                 guides(fill = "none")
pr3 = ggplot() + geom_sf(data = cbind(df,resids = resids[,3]), aes(fill = resids)) +
                 scale_fill_gradient2(midpoint = 0,  limits = c(-rc3, rc3)) +
                 astgmap + ylab("")+xlab("")+ labs(title = "Residual Others")+ 
                 theme(plot.title = element_text(size=17))+
                 guides(fill = "none")

prr1 = ggplot(data.frame(x = 1:N,y = resids[,1]),aes(x=x,y=y)) + 
        geom_point() + coord_equal(ratio = 6)+theme_gray(base_size = 16) +
        ylab("") + xlab("Node index") + ylim(c(-rc1, rc1))
prr2 = ggplot(data.frame(x = 1:N,y = resids[,2]),aes(x=x,y=y)) + 
        geom_point() +  coord_equal(ratio = 6)+theme_gray(base_size = 16) +
        ylab("") + xlab("Node index") + ylim(c(-rc2, rc2))
prr3 = ggplot(data.frame(x = 1:N,y = resids[,3]),aes(x=x,y=y)) + 
        geom_point() +  coord_equal(ratio = 6)+theme_gray(base_size = 16) +
        ylab("") + xlab("Node index") + ylim(c(-rc3, rc3))


pp = grid.arrange(arrangeGrob(pr1,prr1,ncol = 2,widths=c(2,0.5)),
                  arrangeGrob(pr2,prr2,ncol = 2,widths=c(2,0.5)),
                  arrangeGrob(pr3,prr3,ncol = 2,widths=c(2,0.5)),nrow = 3)
grid.arrange(pr1,prr1,pr2,prr2,pr3,prr3)

pp1 = grid.arrange(pr1,prr1,nrow = 1)
pp2 = grid.arrange(pr2,prr2,nrow = 1)
pp3 = grid.arrange(pr3,prr3,nrow = 1)

#pp = grid.arrange(p1,p2,p3,nrow = 3)

ggsave(paste0(here::here(),"/plots/","residualplots1.png"),pp1)
ggsave(paste0(here::here(),"/plots/","residualplots2.png"),pp2)
ggsave(paste0(here::here(),"/plots/","residualplots3.png"),pp3)


# residual plots
plot(model$mu[,1],X[,1])
plot(model$mu[,2],X[,2])
plot(model$mu[,3],X[,3])


# ternary plot 
expResids = exp(resids)
expResids = t(apply(expResids,1,function(u){u/sum(u)}))
colnames(expResids) = c("left","right","others")
pter = robCompositions::ternaryDiag(expResids)
pter = compositions::plot.acomp(expResids)
ggsave(paste0(here::here(),"/plots/","ternary.png"),pter)




# maybe also plot positions vs neighbour for the resiudals
vars = rbind(1:3,1:3)
apply(vars,2,function(xidx){
        plotdf = t(apply(rbind(edgesmat,cbind(edgesmat[,2],edgesmat[,1])),1,function(edge){
                                       c(resids[edge[1],xidx[1]],resids[edge[2],xidx[2]])}))
        plotdf = as.data.frame(plotdf)
        colnames(plotdf) = c("position","neighbor")
        ggplot(plotdf,aes(x = position, y = neighbor))+ geom_point() #+ geom_smooth(method='lm', formula= y~x)
})



# residual maps for Paris 
rc1 = max(abs(resids[ileidx,1]))
rc2 = max(abs(resids[ileidx,2]))
rc3 = max(abs(resids[ileidx,3]))
#rc = max(rc1,rc2,rc3); rc1 = rc; rc2 = rc; rc3 = rc  # if same scale for all 
pr1 = ggplot() + geom_sf(data = cbind(df[ileidx,],resids = resids[ileidx,1]), aes(fill = resids)) +
                 scale_fill_gradient2(midpoint = 0,  limits = c(-rc1, rc1)) +
                 astgmap + ylab("")+xlab("")
pr2 = ggplot() + geom_sf(data = cbind(df[ileidx,],resids = resids[ileidx,2]), aes(fill = resids)) +
                 scale_fill_gradient2(midpoint = 0,  limits = c(-rc2, rc2))+
                 astgmap + ylab("")+xlab("")
pr3 = ggplot() + geom_sf(data = cbind(df[ileidx,],resids = resids[ileidx,3]), aes(fill = resids)) +
                 scale_fill_gradient2(midpoint = 0,  limits = c(-rc3, rc3)) +
                 astgmap + ylab("")+xlab("")

grid.arrange(pr1,pr2,pr3, nrow = 1)


# also calculate model without any robustness
thetnr = ginv(t(Z)%*% Lp %*% Z) %*% t(Z) %*% Lp %*% X
munr   = Z %*% thetnr
Signr = t(X - munr) %*% Lp %*% (X - munr)





