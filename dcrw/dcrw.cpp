/* 
 * can we fit a difference correlated random walk model to a GLATOS fish? 
 *                        Cahill March 2023
 */

#include <TMB.hpp>
using namespace density; 

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(y); 
  DATA_VECTOR(dt); 
  
  PARAMETER(logSdlat); 
  PARAMETER(logSdlon); 
  PARAMETER(logSdgamma); 
  PARAMETER_VECTOR(gamma); 
  
  Type sdLat = exp(logSdlat);
  Type sdLon = exp(logSdlon);
  Type sdGamma = exp(logSdgamma);
  
  Type jnll = 0.0;
  
  for(int i = 1; i < gamma.size(); ++i){
    jnll -= dnorm(gamma(i), gamma(i-1), dt(i-1)*sdGamma, true);
  }
  
  matrix<Type> covs(2,2);
  vector<Type> tmp(2);
  
  for(int i = 2; i < y.cols(); ++i){
    tmp = y.col(i) - (y.col(i-1) + Type(gamma(i))*(dt(i-1)/dt(i-2))*(y.col(i-1)-y.col(i-2)));
    covs << dt(i-1)*dt(i-1)*sdLon*sdLon, 0.0, 0.0, dt(i-1)*dt(i-1)*sdLat*sdLat;
    MVNORM_t<Type> nll_dens(covs);
    jnll += nll_dens(tmp);
  }
  ADREPORT(sdLat);
  ADREPORT(sdLon);
  ADREPORT(sdGamma);
  
  return jnll;
}

