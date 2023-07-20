// 
// Code by Christopher Rieser 
// 


#define ARMA_WARN_LEVEL 0
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;



// function to compute edgewise difference data frame 
//[[Rcpp::export]]
mat edgediff(const mat &x,const mat &edgesmat, const mat &Lp){
  
  // init 
  const int p=x.n_cols,ned = edgesmat.n_rows;
  mat xedges(ned,p);
  int e1, e2;
  
  for(int eidx=0;eidx<ned;eidx++){
    e1 = edgesmat(eidx,0)-1;
    e2 = edgesmat(eidx,1)-1;
    xedges.row(eidx) = (x.row(e1)-x.row(e2))*sqrt(abs(Lp(e1,e2)));
  }
  return(xedges);
}

// function to compute reduced set of 
// scaled eigenvectors of L 
//[[Rcpp::export]]
mat eigL(const mat &L){

  // init
  vec eigval;
  mat eigvec;
  double thresh;
  const int n=L.n_cols;

  // eigen decomp
  eig_sym(eigval, eigvec, L);
  thresh = std::max(eigval(n-1)*(1e-8),double(0));

    
  // get indices
  mat eigvecred = eigvec.cols(find(eigval >= thresh));
  mat eigvalred = eigval(find(eigval >= thresh));
  const int nr = eigvecred.n_cols;
  
  // multiply with eigenvalues
  for(int i=0;i<nr;i++){ eigvecred.col(i) = eigvecred.col(i)*sqrt(abs(eigvalred(i))); }
  
  return(eigvecred);
}

// function to compute the isometric map
//[[Rcpp::export]]
mat giso(const mat &x,const mat &L){
  
  // get scaled eigenvecs and multiply
  return(eigL(L).t()*x);
}



// function to compute one c step 
// returns the new estimate of Sigma and theta
//[[Rcpp::export]]
void cstep(const uvec &mdord,
          const mat &x,
          const mat &z,
          const int &h,
          mat &Sig,
          mat &theta){

  
  // initialize 
  const int n=x.n_rows,p=x.n_cols,q=z.n_cols;
  mat xSub(h,p);
  mat zSub(h,q);
	mat Data2(n,p);
	vec md2(n);
	
  // copy data
	for(int i=0;i<h;i++){
	  xSub.row(i)=x.row(mdord(i));
	  zSub.row(i)=z.row(mdord(i));
	} 
  
  // compute center of current subset
  mat qm,rm;
  qr_econ(qm,rm,zSub);                  
  theta = solve(rm,(trans(qm)*xSub));
  
  // new sigma 
  xSub -= zSub*theta;
  Sig  = (1/(double)(h-1))*xSub.t()*xSub; // new covariance 

}

// function to compute mahalanobis distances - without wij 
//[[Rcpp::export]]
vec mds(const mat &x, 
        const mat &z,
        const mat &Sig,
        const mat &theta){
  
  // init
  const int n=x.n_rows,p=x.n_cols;
  mat Sigi(p,p);
  mat u(1,p);
  vec mdss(n);
  
  // inverse of covariance
  Sigi = pinv(Sig);

  // mahalanobis distances
  for(int i=0;i<n;i++){
    u       = x.row(i)-z.row(i)*theta;
    mdss(i) = sqrt(abs(as_scalar(u*Sigi*u.t())));  // distances 
  }
  return(mdss);
}


// function computing edgewise deltas_ij
//[[Rcpp::export]]
vec edgemd2(const mat &xo, 
            const mat &mu,
            const mat &Sig,
            const mat &edgesmat,
            const mat &L){
  
    // init
  const int p=xo.n_cols,nedges = edgesmat.n_rows; 
  mat Sigi(p,p);
  mat u1(1,p);
  mat u2(1,p);
  mat u(1,p);
  vec mdsw(nedges);
  int iidx,jidx;
  double wij;
  
  // inverse of covariance
  Sigi = pinv(Sig);

   for(int i=0;i<nedges;i++){
     // get edge indices
     iidx = edgesmat(i,0)-1;
     jidx = edgesmat(i,1)-1;
     
     // calculate edge difference 
     u1 = xo.row(iidx)-mu.row(iidx);
     u2 = xo.row(jidx)-mu.row(jidx);
     u  = u1 - u2;
     
     // compute 
     wij = -L(iidx,jidx);
     mdsw(i) = (abs(as_scalar(u*Sigi*u.t())))*wij;
   }
  
  return(mdsw);
}


// function to compute edgewise thresholds
//[[Rcpp::export]]
vec edgethresh(const mat &edgesmat,
               const mat &L,
               const mat &Li,
               const double &quant){
   // init 
   const int nedges=edgesmat.n_rows;
   vec trs(nedges);
   int iidx,jidx;
   double lii,ljj,lij, wij;
    
   for(int i=0;i<nedges;i++){
     // get edge indices
     iidx = edgesmat(i,0)-1;
     jidx = edgesmat(i,1)-1;
     
     // compute 
     lii = Li(iidx,iidx);
     ljj = Li(jidx,jidx);
     lij = Li(iidx,jidx);
     wij = -L(iidx,jidx);
     trs(i) = wij*(lii+ljj-2*lij)*quant;
   }
  
  return(trs);
}


// function to do reweighting step
void reweightStep(const mat &xedges,
                  const mat &zedges,
                  mat &Sigma, 
                  mat &theta, 
                  const vec &thress, 
                  const vec &mdss,
                  uvec &outl){
  
  // init 
  const int nedges=thress.n_elem;
  //uvec outl(nedges);
  uvec mdord(nedges);
  mdord.zeros();
  int h; 
  int counter=0;
  
  // reweighting step
  for(int i=0;i<nedges;i++){ 
    if(pow(mdss(i),2)<(thress(i))){
         outl(i) = 0;
         mdord(counter) = i;
         counter++;
    }else{ outl(i) = 1;}
  }
  // recompute solution 
  mdord.resize(counter);
  h  = counter;//mdord.n_elem;
  //Rcpp::Rcout << h; // print 
  //Rcpp::Rcout << mdord.n_elem; // print 

   // reweighting step
  cstep(mdord,xedges,zedges,h,Sigma,theta);
  
}


// function to compute solution given initial values 
void detmcd(const mat &xedges,
          const mat &zedges,
          mat &Sig0,
          mat &theta0,
          const mat &L,
          const mat &edgesmat,
          const double &alpha,
          const int &maxit,
          const double &maxerr){
  
  // init 
  const int p=xedges.n_cols,q=zedges.n_cols,nedges=edgesmat.n_rows;
  mat Sigma = Sig0;
  mat theta = theta0;
  mat Sig0ld(p,p);
  mat th0old(q,p);
  vec mdss(nedges);
  uvec mdord(nedges);
  uvec outl(nedges);
  vec thress(nedges);
  vec vloc(nedges);

  // consts 
  double h   = std::max(floor(alpha*nedges),double((nedges+p+1)/2));
  double err = std::numeric_limits<double>::max();

  // loop till convergence or maximal loop nbr reached
  for(int i=0;i<maxit;i++){
    
    // reset
    Sig0ld = Sigma;
    th0old = theta;
    
    // c steps
    mdss  = mds(xedges,zedges,Sigma,theta);// mahalanobis distances  
    mdord = sort_index(mdss,"ascend");     // order of mdss
    cstep(mdord,xedges,zedges,h,Sigma,theta);    // do a c step 
  
    // compute error and check 
    err = norm(Sigma-Sig0ld,"fro") / (1+norm(Sig0ld,"fro"));
    err = std::max(err,norm(theta-th0old,"fro") / (1+norm(th0old,"fro")));
    
    // error
    if(err < maxerr) break;
  }
  
  // set solutions
  Sig0   = Sigma;
  theta0 = theta;
}

// function to scale for matrix A
void scpp(mat &A, const int &ii, const vec &aa){
  const int pr=A.n_rows,pc=A.n_cols;
  if(ii==1){ vec zz(pc); for(int k=0; k<pr; k++){ zz = A.row(k).t()/double(aa(k)); A.row(k) = zz.t(); }}
  if(ii==2){ vec zz(pr); for(int k=0; k<pc; k++){ zz = A.col(k)/double(aa(k)); A.col(k) = zz; }}
}

// function to align quantiles
void align(const mat &xedges,
           const mat &zedges,
                 mat &Sigma,
           const mat &theta,
           const mat &L,
           const mat &edgesmat,
           const mat &Li){
  
  // init
  const int p=xedges.n_cols;
  const int nedges=edgesmat.n_rows;
  vec mdss(nedges), thress(nedges), vloc(nedges);
  double qcmed = R::qchisq(0.5,p,TRUE,FALSE);  // aligning chi sq median 
  double scf; 
  
  thress = edgethresh(edgesmat,L,Li,1);
  mdss   = mds(xedges,zedges,Sigma,theta);
  vloc   = pow(mdss,2) / thress;
  scf    = median(vloc) / qcmed;
  Sigma  = scf*Sigma ;
}


// function to align quantiles with giso 
//[[Rcpp::export]]
void gisoalign(const mat &xo,
               const mat &zo,
                     mat &Sigma,
               const mat &theta,
               const mat &L){
  
  // init
  const int p    = xo.n_cols;
  const mat Sigi = pinv(Sigma);
  double qcmed   = R::qchisq(0.5,p,TRUE,FALSE);  // aligning chi sq median 
  double scf; 
  mat vv(1,p);
  mat mu    = zo * theta;
  mat gisox      = giso(xo-mu,L);
  const int niso = gisox.n_rows;
  vec mdss(niso);
  
  // compute mahalanobis distances 
  for(int i=0;i<niso;i++){ 
      vv      = gisox.row(i);
      mdss(i) = abs(as_scalar(vv*Sigi*vv.t()));
  }
  
  // align 
  scf    = median(mdss) / qcmed;
  Sigma  = scf*Sigma ;
}




// function to compute solution given multiple initial estimators
//[[Rcpp::export]]
List detmcdBest(const mat &x,
                const mat &z,
                const cube &Sigs0,
                const cube &thetas0,
                const mat &xscales,
                const int &ninit,
                const mat &L,
                const mat &edgesmat,
                const mat &Li,
                const double &alpha,
                const int &maxit,
                const double &maxerr){
  //
  // TODO::
  // implement approximation to Li through reduced svd or such and only calculate once
  
  
  // init 
  const int p=x.n_cols,q=z.n_cols;
  const int nedges=edgesmat.n_rows;
  double qcrew = R::qchisq(0.975,p,TRUE,FALSE); // reweighting chi sq quantile 
  vec mdss(nedges), thress(nedges), scales0(p), obj(ninit), vloc(nedges);
  mat Sigma(p,p), theta(q,p), Scl(p,p), thcl(q,p);
  mat xze(nedges,p), xscel(nedges,p), xedges(nedges,p), zedges(nedges,q);
  cube Sc(p,p,ninit), thc(q,p,ninit);
  uvec outl(nedges);
  uword besind;
  
  // edge differences
  xedges = edgediff(x,edgesmat, L);
  zedges = edgediff(z,edgesmat, L);
  
  for(int i=0;i<ninit;i++){
      // get initial estimators & scale
      Scl      = Sigs0.slice(i);
      thcl     = thetas0.slice(i);
      xscel    = xedges;
      scales0  = xscales.col(i);
      
      // scaling 
      scpp(Scl,1,scales0), scpp(Scl,2,scales0);
      scpp(thcl,2,scales0); scpp(xscel,2,scales0);

      // compute mcd solutions
      detmcd(xscel,zedges,Scl,thcl,L,edgesmat,alpha, maxit,maxerr);
  
      // save
      Sc.slice(i)  = Scl;
      thc.slice(i) = thcl;
      obj(i)       = log(abs(det(Scl)));
  }

  // get solution with smallest objective value 
  besind = obj.index_min();
  Sigma  = Sc.slice(besind);
  theta  = thc.slice(besind);
  
  // get rid of scale 
  scpp(Sigma,1,1/scales0), scpp(Sigma,2,1/scales0);
  scpp(theta,2,1/scales0);
  
  // align quantiles
  align(xedges,zedges,Sigma,theta,L,edgesmat,Li);
  //gisoalign(x,z,Sigma,theta,L);
  
  // reweight
  thress = edgethresh(edgesmat,L,Li,qcrew);
  mdss   = mds(xedges,zedges,Sigma,theta);
  reweightStep(xedges,zedges,Sigma,theta,thress,mdss,outl);
  
  
  // align quantiles
  align(xedges,zedges,Sigma,theta,L,edgesmat,Li);
  //gisoalign(x,z,Sigma,theta,L);
  
  //Rcout << "afterafter: " << Sigma << "\n";


  // return best solution 
  List sol = List::create(Rcpp::Named("Sigma") = Sigma,
                          Rcpp::Named("theta") = theta,
                          Rcpp::Named("outlIndicator") = outl,
                          Rcpp::Named("Xedges_weighted") = xedges,
                          Rcpp::Named("Zedges_weighted") = zedges);
  
  return(sol);
}

