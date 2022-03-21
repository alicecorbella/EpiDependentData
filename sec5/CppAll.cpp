#include <Rcpp.h>
#include <numeric>
#include <unordered_set>
#include <algorithm>
#include <string>
#include <stdio.h>      /* printf */
#include <math.h>       /* lgamma */
using namespace Rcpp;
//[[Rcpp::plugins(cpp11)]]


/*-----------------------------------------------------------
 Function SEEIIRdisCtv solves on intervals dt the system of equation
 of a SEEIIR model with time varying transmission, takes as a iput:
 - begin   start approximation
 - end     start approximation
 - spd     steps per day
 - X0      state value at begin, number of people in each S E I R
 - Theta   vector of parameters, named with values bt, gmm, sgm
 - k       increasing decreasing factor during some intervals
 - Tk      vector of length (end-begin)*spd od 0s and 1s. The 1s 
           identify the intervals during which the factor k is applied
 ------------------------------------------------------------- */

// [[Rcpp::export]]
NumericMatrix SEEIIRdisCtv(double begin, double end, double spd,
                           NumericVector X0, NumericVector Theta, 
                           double k, NumericVector Tk) {
  int sizemat=std::floor((end-begin)*spd)+1;  
  const double dt = 1/spd;
  NumericMatrix mat(sizemat, 7);
  mat(0,0)=begin;
  mat(0,1)=X0[0];
  mat(0,2)=X0[1]/2;
  mat(0,3)=X0[1]/2;
  mat(0,4)=X0[2]/2;
  mat(0,5)=X0[2]/2;
  mat(0,6)=X0[3]; 
  int N = X0[0]+X0[1]+X0[2]+X0[3];
  for (int t=1;t<sizemat;t++){
    mat(t,1) = mat(t-1,1) - dt*Theta["bt"]*(Tk[t-1]*k+1)*mat(t-1,1)*((mat(t-1,4)+mat(t-1,5))/N);
    mat(t,2) = mat(t-1,2) + dt*Theta["bt"]*(Tk[t-1]*k+1)*mat(t-1,1)*((mat(t-1,4)+mat(t-1,5))/N) - dt*Theta["sgm"]*mat(t-1,2);
    mat(t,3) = mat(t-1,3) + dt*Theta["sgm"]*mat(t-1,2) - dt*Theta["sgm"]*mat(t-1,3);
    mat(t,4) = mat(t-1,4) + dt*Theta["sgm"]*mat(t-1,3) - dt*Theta["gmm"]*mat(t-1,4);
    mat(t,5) = mat(t-1,5) + dt*Theta["gmm"]*mat(t-1,4) - dt*Theta["gmm"]*mat(t-1,5);
    mat(t,6) = mat(t-1,6) + dt*Theta["gmm"]*mat(t-1,6);
    mat(t,0) = begin + (t)*dt;
  }
  return(mat);
}

/*-----------------------------------------------------------
 Function LamTseeiirCtv solves on intervals dt the system of equation
 of a SEEIIR time-varying model and retun the number of new infected in 
 each interval of length outdt, it is like applying the function 
 at given values, takes as a ipnut:
 - tm_begin   start transmission approximation
 - tm_end     end transmission approximation
 - tm_spd     steps per days at which we need to make the transmission approximation
 - tm_X0      state value at begin, number of people in each S E I R
 - tm_Theta   vector of parameters, named with values bt, gmm, sgm
 - outdt   length (in days) over which we want to count incidence
 - tm_k       increasing decreasing factor during some intervals
 - tm_Tk      vector of length (end-begin)*spd od 0s and 1s. The 1s 
 identify the intervals during which the factor k is applied
 ------------------------------------------------------------- */
// [[Rcpp::export]]
NumericVector LamTseeiirCtv(double tm_begin, double tm_end, double tm_spd,
                            NumericVector tm_X0, NumericVector tm_Theta, 
                            int outdt , double tm_k, NumericVector tm_Tk){
  int sizevec=std::floor((tm_end-tm_begin)/outdt);
  int thin=std::floor(outdt*tm_spd);
  NumericVector lT(sizevec);
  NumericMatrix x;
  //call the functuon SEIRdisC to approximate the model
  x=SEEIIRdisCtv(tm_begin, tm_end, tm_spd, tm_X0, tm_Theta, tm_k, tm_Tk);
  for (int s=0;s<(sizevec);s++){
    lT[s] = x(s*thin,1) - x((s+1)*thin,1);
  }
  return(lT);
}




/* ----------------------------------------------------------
fdisC is a function that takes as imput the distribution and the 
paramerers of a positive continuous distribution and the deifnition
of the length of the steps and the number of steps and reports 
its discrete version

inputs:
- NAMEDIST : the name of the distribution, one of : 
exponential, gamma, lognormal, weibull
- PARAMS   : a vector of different sizes according to the variable selected
changes length and meaning for each distribution
exponential  -> rate;
gamma        -> shape, rate; so that E(gamma)= shape/rate
lognormal    -> logmu, sigmalog; (sigmalog=sd on the normal scale)
weibull      -> sahape, scale
- SIZESTEP : steps of the t


----------------------------------------------------------*/

// [[Rcpp::export]]
NumericVector fdisC(std::string NAMEDIST, 
                    NumericVector PARAMS,
                    double SIZESTEP, int NSTEPS){
  NumericVector fx(NSTEPS);
  NumericVector px(NSTEPS);
  NumericVector xs(NSTEPS+1);
  NumericVector Fx(NSTEPS+1);
  bool error=false;
  for (int t=0; t<=NSTEPS; t++){
    xs[t] = SIZESTEP*t;
  }
  if(NAMEDIST=="exponential"){
    double rate=PARAMS[0];
    Fx = pexp(xs, rate);
  }else if(NAMEDIST=="gamma"){
    double shape=PARAMS[0];
    double scale=1/PARAMS[1];
    Fx = pgamma(xs, shape, scale);
  }else if(NAMEDIST=="lognormal"){
    double lmu=PARAMS[0];
    double sgl=PARAMS[1];
    Fx = plnorm(xs, lmu, sgl);
  }else if(NAMEDIST=="weibull"){
    double shape=PARAMS[0];
    double scale=PARAMS[1];
    Fx = pweibull(xs, shape, scale);
  }else{
    std::cout<< "NAMEDIST should be one of \n exponential, gamma, lognormal, weibull";
    error=true;
  }
  if(error==true){
    px = NumericVector::create(NSTEPS+1,NAN);
  }else{
    for (int t=0; t<NSTEPS; t++){
      fx[t] = Fx[t+1]-Fx[t];
    }
    px = fx/sum(fx);
  }
  return(px);
}

/* ----------------------------------------------------------
ManyMultinomC is a function that takes as imput a vector of probability, 
and a vecotor of sizes, then for each of the sized the function samples 
a multinomial distribution with those defined probabilities and returns 
it in the form of an iteger matrix with nrows=length of probability 
vecors and ncols=length of the size vector

inputs:
- probs : the vector of probabilities that musr sum to 1
- sizes : vector of the sizes
----------------------------------------------------------*/

// [[Rcpp::export]]
IntegerMatrix ManyMultinomC(NumericVector probs, IntegerVector sizes) {
  int k = probs.size();
  int s = sizes.size();
  IntegerMatrix matans(k,s);
  IntegerVector ans(k);
  for (int i=0;i<s;i++){
    rmultinom(sizes[i], probs.begin(), k, ans.begin());
    matans(_,i)  =ans;
  }
  return(matans);
}

/* ----------------------------------------------------------
conC is a function that takes as imput an iteger matrix 
of sampled multinomial distribution with differetn sizes with 
nrows=length of probability vecors and ncols=length of the size vector
and returns its crossdiagonal sums (i.e. the number of people that 
have  had event 1 in column i and will wait time t before event 2
so eventually who have event 2 at whtat time)
/ / / /
/ / / /
/ / / /

inputs:
- mulmat : the matrix resultig from the multinomial
----------------------------------------------------------*/

// [[Rcpp::export]]
IntegerVector conC(IntegerMatrix mulmat){
  int T = mulmat.cols();
  int S = mulmat.rows();
  IntegerVector out(T);
  for (int t=0;t<T;t++){
    int xt=0;
    for (int s=0; s<(std::min(t+1,S));s++){
      xt += mulmat(s,(t-s));
    }
    out[t] = xt;
  }
  return(out);
}


/* ----------------------------------------------------------
projC puts togetehr the two fucntions ManyMultinomC and 
convC and makes only one function that takes as imput the prob
of the waiting time and gives as output the delayed distribution

inputs:
- probs : the vector of probabilities that must sum to 1
- sizes : vector of the sizes
----------------------------------------------------------*/

// [[Rcpp::export]]
IntegerVector projC(NumericVector probs, IntegerVector sizes){
  int k = probs.size();
  int s = sizes.size();
  IntegerMatrix mulmat(k,s);
  IntegerVector ans(k);
  for (int i=0;i<s;i++){
    rmultinom(sizes[i], probs.begin(), k, ans.begin());
    mulmat(_,i)  =ans;
  }
  int T = mulmat.cols();
  int S = mulmat.rows();
  IntegerVector out(T);
  for (int t=0;t<T;t++){
    int xt=0;
    for (int s=0; s<(std::min(t+1,S));s++){
      xt += mulmat(s,(t-s));
    }
    out[t] = xt;
  }
  return(out);
}


/* ----------------------------------------------------------
projCdet in the deterministic version of projCdet:
takes as imput the prob of the waiting time and the initial averages
and gives as output the delayed distribution

inputs:
- probs : the vector of probabilities that must sum to 1
- sizes : vector of the averages before delay
----------------------------------------------------------*/

// [[Rcpp::export]]
NumericVector projCdet(NumericVector probs, NumericVector sizes){
  int k = probs.size();
  int s = sizes.size();
  NumericMatrix mulmat(k,s);
  NumericVector ans(k);
  for (int i=0;i<s;i++){
    for (int j=0; j<k;j++){
      mulmat(j,i)=sizes(i)*probs(j);
    }
  }
  int T = mulmat.cols();
  int S = mulmat.rows();
  NumericVector out(T);
  for (int t=0;t<T;t++){
    double xt=0;
    for (int s=0; s<(std::min(t+1,S));s++){
      xt += mulmat(s,(t-s));
    }
    out[t] = xt;
  }
  return(out);
}


/* ----------------------------------------------------------
deconC is a function that takes as imput an iteger matrix 
of sampled multinomial distribution with differetn sizes with 
nrows=length of probability vecors and ncols=length of the size vector
and returns its diagonal sums (i.e. the number of people that 
have  had event 1 in column i and have waited time t before event 2
so eventually who had event 1 at whtat time)
\ \ \ \ 
\ \ \ \
\ \ \ \ 
\ \ \ \ 

inputs:
- mulmat : the matrix resultig from the multinomial
----------------------------------------------------------*/

// [[Rcpp::export]]
IntegerVector deconC(IntegerMatrix mulmat){
  int T = mulmat.cols();
  int S = mulmat.rows();
  IntegerVector out(T);
  for (int t=0;t<T;t++){
    int xt=0;
    for (int s=0; s<(std::min(T-t,S));s++){
      xt += mulmat(s,(t+s));
    }
    out[t] = xt;
  }
  return(out);
}


/* ----------------------------------------------------------
backprojC puts together the two fucntions ManyMultinomC and 
convC and makes only one function that takes as imput the prob
of the waiting time and gives as output the 'earlied' distribution

inputs:
- probs : the vector of probabilities that musr sum to 1
- sizes : vector of the sizes
----------------------------------------------------------*/

// [[Rcpp::export]]
IntegerVector backprojC(NumericVector probs, IntegerVector sizes){
  int k = probs.size();
  int s = sizes.size();
  IntegerVector out(s);
  IntegerMatrix matans(k,s);
  IntegerVector ans(k);
  for (int i=0;i<s;i++){
    rmultinom(sizes[i], probs.begin(), k, ans.begin());
    matans(_,i)  =ans;
  }
  out = deconC(matans);
  return(out);
}

// same function as above though applied to all the columns of a matrix
// [[Rcpp::export]]
IntegerMatrix backprojCmat(NumericVector probs, IntegerMatrix sizes){
  IntegerMatrix out(sizes.rows(), sizes.cols());
  int N = sizes.cols();
  for (int n=0;n<N;n++){
    out(_,n) = backprojC(probs, sizes(_,n));
  }
  return(out);
}

/* ----------------------------------------------------------
 dbbCmN density of a beta binomial with numtiple Ns
 inputs:
 - x    : integer value at which to evaluate the density
 - N    : vector of sizes
 - u    : alpha paramtere of the beta binomial
 - v    : beta paramter of the beta binomial
 
 ----------------------------------------------------------*/
// [[Rcpp::export]]
NumericVector dlbbCmN(int x, NumericVector N,
                     double u, double v){
  int S=N.size();
  double inf = std::numeric_limits<double>::infinity();
  NumericVector out(S);
  out = lgamma(N+1)-lgamma(N-x+1)-lgamma(x+1)+lgamma(x+u)+lgamma(N-x+v)+
    lgamma(u+v)-(lgamma(u+N+v)+lgamma(u)+lgamma(v));
  return(out);
}

// [[Rcpp::export]]
double dlbbC(int x, double N,
                      double u, double v){
  double inf = std::numeric_limits<double>::infinity();
  double out;
  out = lgamma(N+1)-lgamma(N-x+1)-lgamma(x+1)+lgamma(x+u)+lgamma(N-x+v)+
    lgamma(u+v)-(lgamma(u+N+v)+lgamma(u)+lgamma(v));
  return(out);
}

/* ----------------------------------------------------------
 cllC conditional log likelihood of the hospitalization given 
 the ICU admissions, fnction per se which is the most expenseive since
 it is thwough simulation
 inputs:
 - xi0      : weekly number of new infections, 
 - del0toH  : vector of relevant delay between infecion and Hospitalization,  
 - delHtoIC : vector of relevant delay between Hospitalization and IC admission,    
 - th0toH   : probability of hospitalization given infection,  
 - thHtoIC  : probability of IC admissions given hospitalization  
 - epsH     : precision of Hospitalization detection, 
 - dH       : vector of Hospitalization detection,  
 - dIC      : vector of IC admission detection, 
 - yH       : data: vector of Hospitalization,  
 - yIC      : data: vector of IC admissions,  
 - NP       : number of simulations.
 ----------------------------------------------------------*/
// [[Rcpp::export]]
double cllC(NumericVector xi0, 
           NumericVector del0toH, NumericVector delHtoIC, 
           double th0toH, double thHtoIC, double epsH, 
           NumericVector dH, NumericVector dIC, 
           IntegerVector yH, IntegerVector yIC, int NP){
  int t_end =xi0.size();
  NumericVector lH =(projCdet(del0toH , xi0))*th0toH;
  NumericVector lIC=(projCdet(delHtoIC, lH))*thHtoIC;
  IntegerMatrix XICgyIC(t_end, NP);
  IntegerMatrix XHtoICgsxHtoIC(t_end, NP);
  IntegerMatrix xHgxHtIC(t_end, NP);
  NumericMatrix lPyHpar(t_end, NP);
  NumericVector lPyH(NP);
  double ll=0;
  for (int t=0;t<t_end;t++){
    XICgyIC(t,_)=Rcpp::rpois(NP, ((1-dIC[t])*lIC[t]))+ yIC[t];
  }  
  for(int n=0;n<NP; n++){
    XHtoICgsxHtoIC(_,n)= backprojC(delHtoIC, XICgyIC(_,n));
  }
  for (int t=0;t<t_end;t++){
    Rcpp::NumericVector xi(NP);
    xi=XHtoICgsxHtoIC(t,_);
    xHgxHtIC(t,_)= Rcpp::rpois(NP, (1-thHtoIC)*lH[t]) + xi;
    Rcpp::NumericVector xb(NP);
    xb=xHgxHtIC(t,_);
    lPyHpar(t,_) = dlbbCmN(yH[t],xb, (dH[t]/(1-dH[t]))*epsH,epsH);
  }
  for(int n=0;n<NP; n++){
    NumericVector x= lPyHpar(_,n);
    lPyH[n] = std::accumulate(x.begin(),x.end(), 0.0);
  }
  ll = log(mean(exp(lPyH)));
  return(ll);
}

// [[Rcpp::export]]
double cllCloop(NumericVector xi0, 
            NumericVector del0toH, NumericVector delHtoIC, 
            double th0toH, double thHtoIC, double epsH, 
            NumericVector dH, NumericVector dIC, 
            IntegerVector yH, IntegerVector yIC, int NP){
  int t_end =xi0.size();
  NumericVector lH =(projCdet(del0toH , xi0))*th0toH;
  NumericVector lIC=(projCdet(delHtoIC, lH))*thHtoIC;
  IntegerMatrix XICgyIC(t_end, NP);
  IntegerMatrix XHtoICgsxHtoIC(t_end, NP);
  IntegerMatrix xHgxHtIC(t_end, NP);
  // NumericMatrix lPyHpar(t_end, NP);
  NumericVector lPyH(NP);
  double ll=0;
  for(int n=0;n<NP; n++){
    lPyH[n]=0;
    for (int t=0;t<t_end;t++){
      XICgyIC(t,n)=R::rpois(((1-dIC[t])*lIC[t]))+ yIC[t];
    }
    XHtoICgsxHtoIC(_,n)= backprojC(delHtoIC, XICgyIC(_,n));
    for (int t=0;t<t_end;t++){
      xHgxHtIC(t,n)= R::rpois((1-thHtoIC)*lH[t]) + XHtoICgsxHtoIC(t,n);
      lPyH[n]+=dlbbC(yH[t], xHgxHtIC(t,n), (dH[t]/(1-dH[t]))*epsH,epsH);
    }
  }   
  ll = log(mean(exp(lPyH)));
  return(ll);
}


// JointLogLikelihood
// [[Rcpp::export]]
double jllC(NumericVector xi0, 
            NumericVector del0toH, NumericVector delHtoIC, 
            double th0toH, double thHtoIC, double epsH, 
            NumericVector dH, NumericVector dIC, 
            IntegerVector yH, IntegerVector yIC, int NP){
  int t_end =xi0.size();
  NumericVector lH =(projCdet(del0toH , xi0))*th0toH;
  NumericVector lIC=(projCdet(delHtoIC, lH))*thHtoIC;
  IntegerMatrix XICgyIC(t_end, NP);
  IntegerMatrix XHtoICgsxHtoIC(t_end, NP);
  IntegerMatrix xHgxHtIC(t_end, NP);
  NumericMatrix lPyHpar(t_end, NP);
  NumericVector lPyH(NP);
  double llic=0;
  double ll=0;
  for (int t=0;t<t_end;t++){
    XICgyIC(t,_)=Rcpp::rpois(NP, ((1-dIC[t])*lIC[t]))+ yIC[t];
    llic+=R::dpois(yIC[t], lIC[t], true);
  }  
  for(int n=0;n<NP; n++){
    XHtoICgsxHtoIC(_,n)= backprojC(delHtoIC, XICgyIC(_,n));
  }
  for (int t=0;t<t_end;t++){
    Rcpp::NumericVector xi(NP);
    xi=XHtoICgsxHtoIC(t,_);
    xHgxHtIC(t,_)= Rcpp::rpois(NP, (1-thHtoIC)*lH[t]) + xi;
    Rcpp::NumericVector xb(NP);
    xb=xHgxHtIC(t,_);
    lPyHpar(t,_) = dlbbCmN(yH[t],xb, (dH[t]/(1-dH[t]))*epsH,epsH);
  }
  for(int n=0;n<NP; n++){
    NumericVector x= lPyHpar(_,n);
    lPyH[n] = std::accumulate(x.begin(),x.end(), 0.0);
  }
  ll = log(mean(exp(lPyH)))+llic;
  return(ll);
}
