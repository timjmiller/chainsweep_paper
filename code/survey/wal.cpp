#include <TMB.hpp>
#include <iostream>

#define see(object) std::cout << #object ":\n" << object << "\n";

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  DATA_VECTOR(wal);
  DATA_VECTOR(lens);
  DATA_VECTOR(pred_lens);
  int n_pred = pred_lens.size();
  
  PARAMETER_VECTOR(lw_pars);
  Type nll = 0.0;

  nll -= dnorm(log(wal), lw_pars(0) + lw_pars(1) * log(lens) - 0.5*exp(2.0*lw_pars(2)), exp(lw_pars(2)), 1).sum();  
	if(n_pred>0)
	{
	  vector<Type> pred_log_wal = lw_pars(0) + lw_pars(1) * log(pred_lens);
	  REPORT(pred_log_wal);
	  ADREPORT(pred_log_wal);
  }
	return(nll);
}

