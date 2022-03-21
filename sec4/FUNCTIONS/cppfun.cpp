#include <Rcpp.h>
#include <numeric>
#include <unordered_set>
#include <algorithm>
#include <string>
using namespace Rcpp;
//[[Rcpp::plugins(cpp11)]]


/*-----------------------------------------------------------
Function SEIRdisC solves on intervals dt the system of equation
 of a SEIR model, takes as a iput:
 - begin   start approximation
 - end     start approximation
 - spd     steps per day
 - X0      state value at begin, number of people in each S E I R
 - Theta   vector of parameters, named with values bt, gmm, sgm
------------------------------------------------------------- */


// [[Rcpp::export]]
NumericMatrix SEIRdisC(double begin, double end, double spd,
                       NumericVector X0, NumericVector Theta) {
  int sizemat=std::floor((end-begin)*spd)+1;  
  const double dt = 1/spd;
  NumericMatrix mat(sizemat, 5);
  mat(0,0)=begin;
  mat(0,1)=X0[0];
  mat(0,2)=X0[1];
  mat(0,3)=X0[2];
  mat(0,4)=X0[3]; 
  int N = X0[0]+X0[1]+X0[2]+X0[3];
  for (int t=1;t<sizemat;t++){
    mat(t,1) = mat(t-1,1) - dt*Theta["bt"]*mat(t-1,1)*(mat(t-1,3)/N);
    mat(t,2) = mat(t-1,2) + dt*Theta["bt"]*mat(t-1,1)*(mat(t-1,3)/N) - dt*Theta["sgm"]*mat(t-1,2);
    mat(t,3) = mat(t-1,3) + dt*Theta["sgm"]*mat(t-1,2) - dt*Theta["gmm"]*mat(t-1,3);
    mat(t,4) = mat(t-1,4) + dt*Theta["gmm"]*mat(t-1,3);
    mat(t,0) = begin + (t)*dt;
  }
  return(mat);
}



/*-----------------------------------------------------------
 Function LamTseirC solves on intervals dt the system of equation
of a SEIR model and retun the number of new infected in each interval
of length outdt, it is like applying the function at give values
, takes as a ipnut:
- tm_begin   start transmission approximation
- tm_end     end transmission approximation
- tm_spd     steps per days at which we need to make the transmission approximation
- tm_X0      state value at begin, number of people in each S E I R
- tm_Theta   vector of parameters, named with values bt, gmm, sgm
- outdt   length (in days) over which we want to count incidence
 
------------------------------------------------------------- */


// [[Rcpp::export]]
NumericVector LamTseirC(double tm_begin, double tm_end, double tm_spd,
                        NumericVector tm_X0, NumericVector tm_Theta, 
                        int outdt){
  int sizevec=std::floor((tm_end-tm_begin)/outdt);
  int thin=std::floor(outdt*tm_spd);
  NumericVector lT(sizevec);
  NumericMatrix x;
  //call the functuon SEIRdisC to approximate the model
  x=SEIRdisC(tm_begin, tm_end, tm_spd, tm_X0, tm_Theta);
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




/* ----------------------------------------------------------
loglikesevPF: log likelihood of the severity of the process, 
 this is a function of:
 - lamt_tm : number of new infections per time-step (obtained straight from 
             the transmission parameter)
 - outdt   : time steps at which Dinft is inputed: must be the same at which
           the data are inputed (length in number of days)
 - yHO   : number of hospitalisation
 - yIC   : number of IC admissions 
 - NP    : number of particles for the approximation
 
 - parEtoH  : parameters of the time from E to H (vecotr either of size 1 or 2)
 - parHtoIC : parameters of the time from H to IC(vecotr either of size 1 or 2)
 - thEtoH   : prob of going from E to H 
 - thHtoIC  : prob of going from H to IC
 - detH     : prob of detection in H 
 - detIC    : prob of detection in IC
 - NMdistrEtoH : character with the name of the distribution of the time fromm E to H
 - NMdistrHtoIC : character with the name of the distribution of the time from H to IC
----------------------------------------------------------*/

// [[Rcpp::export]]
double loglikesevPF(NumericVector lamt_tm, 
                    IntegerVector yH, IntegerVector yIC, int NP, 
                    double thEtoH, double thHtoIC, 
                    NumericVector delEtoH, NumericVector delHtoIC, 
                    NumericVector detH, NumericVector detIC){
  int t_end = yH.size();
  NumericVector lH =(projCdet(delEtoH , lamt_tm))*thEtoH;                 // average rates 
  NumericVector lIC=(projCdet(delHtoIC, lH))*thHtoIC;
  NumericVector dlIC=lIC*detIC;
  double lPyIC=0;
  for(int i=0; i<t_end;i++){
    lPyIC += R::dpois(yIC[i], dlIC[i], true) ;
  }

  IntegerMatrix XICgyIC(t_end,NP);
  for(int i=0; i<t_end;i++){
    XICgyIC(i,_)= rpois(NP, ((1-detIC[i])*lIC[i]) )+ yIC[i] ;
  }
  IntegerMatrix XHtoICgsxHtoIC(t_end, NP);
  for (int p=0;p<NP; p++){
    XHtoICgsxHtoIC(_,p)=backprojC(delHtoIC, XICgyIC(_,p));
  }
  IntegerMatrix xHgxHtIC(t_end,NP);
  for(int i=0; i<t_end;i++){
    Rcpp::NumericVector xi(NP);
    xi=XHtoICgsxHtoIC(i,_);
    xHgxHtIC(i,_)= rpois(NP, (1-thHtoIC)*lH[i]) +  xi;
  }
  NumericVector lPyHpar(NP,0.0);
  double lPyHsum=0;
  for(int p=0; p<NP; p++){
      for(int i=0; i<t_end;i++){
        lPyHpar[p] += R::dbinom(yH[i], xHgxHtIC(i,p), detH[i], true) ;
    }
    lPyHsum+=exp(lPyHpar[p]);
  }
  
  double lPyH;
  lPyH=log(lPyHsum/NP);
  double lL;
  lL=lPyIC+lPyH;
  return(lL);
}



// [[Rcpp::export]]
double loglikesevIND(NumericVector lamt_tm, 
                     IntegerVector yH, IntegerVector yIC, 
                     double thEtoH, double thHtoIC, 
                     NumericVector delEtoH, NumericVector delHtoIC,  
                     NumericVector detH, NumericVector detIC){
  int t_end = yH.size();
  NumericVector lH =(projCdet(delEtoH , lamt_tm))*thEtoH;                 // average rates 
  NumericVector lIC=(projCdet(delHtoIC, lH))*thHtoIC;
  NumericVector dlIC=lIC*detIC;
  NumericVector dlH=lH*detH;
  double lPyIC=0;
  for(int i=0; i<t_end;i++){
    lPyIC += R::dpois(yIC[i], dlIC[i], true) ;
  }
  double lPyH=0;
  for(int i=0; i<t_end;i++){
    lPyH += R::dpois(yH[i], dlH[i], true) ;
  }
  double lL;
  lL=lPyIC+lPyH;
  return(lL);
}

