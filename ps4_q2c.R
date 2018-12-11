# libraries
library(parallel)

# default arguments
args_list = list(
  n_cores=1,
  mc_rep=1e4,
  sigma = 1
)

## get parameters from command line
args = commandArgs(trailingOnly = TRUE)
print(args)

# functions for finding named arguments
args_to_list = function(args){
  ind = grep('=', args)  
  args_list = sapply(args[ind], strsplit,'=')
  names(args_list) = sapply(args_list, function(x) x[1])
  
  args_list = lapply(args_list, function(x) as.numeric(x[2]))
  args_list
}

# get named arguments
args_list_in = args_to_list(args)

# update non default arguments
ignored = c()
for ( arg in names(args_list_in) ) {
  # Check for unknown argument
  if ( is.null(args_list[[arg]]) ) {
    ignored = c(ignored, arg)
  } else{
    # update if known
    args_list[[arg]] = args_list_in[[arg]]
  }
}

source('ps4_q2_funcs.R')

# Parameters: -----------------------------------------------------------------
n= 1000; p = 100
beta = c(rep(0.1, 10),rep(0,p-10) )
dim(beta) = c(p,1)
sigma = args_list$sigma
n_cores = args_list$n_cores
mc_rep = args_list$mc_rep

func4= function(rho, sigma, mc_rep=1e4)
{
  # Sigma: 
  Sigma = rho * beta %*% t(beta)
  diag(Sigma) = 1
  Sigma
  R = chol(Sigma)
  
  # Here is an X for testing: ---------------------------------------------------
  X = matrix( rnorm(n*p), n, p) %*%  R
  # head(X)
  # Monte Carlo simulation function for fixed X:-------- ------------------------
  P = sim_beta(X, beta, sigma =sigma, mc_rep = mc_rep)
  
  # evaluate: -------------------------------------------------------------------
  P_eval = evaluate(P, 1:10)
  # P_eval
  
  
  all0 =
    lapply( c('holm', 'bonferroni', 'BH', 'BY'), function(x){
      evaluate( apply(P, 2, p.adjust, method = x), tp_ind = 1:10)
    })
  all = rbindlist(all0)
  all[ , method := c('holm', 'bonferroni', 'BH', 'BY') ]
  #all
  est = as.data.frame(t(cbind(t(all$fwer), t(all$fdr), t(all$sens), t(all$spec))))
  colnames(est) = 'est'
  est
  se = as.data.frame(t(cbind(t(all$fwer_se), t(all$fdr_se), t(all$sens_se), t(all$spec_se))))
  colnames(se) = 'se'
  se
  metric = as.data.frame(c(rep("fwer",4), rep("fdr",4), rep("sens",4), rep("spec",4)))
  colnames(metric) = 'metric'
  metric
  method = as.data.frame(rep(c("holm","bonferroni","BH","BY"),4))
  colnames(method) = 'method'
  method
  rho1 = as.data.frame(rep(rho, 16))
  colnames(rho1) = 'rho'
  rho1
  sigma1 = as.data.frame(rep(sigma,16))
  sigma1
  colnames(sigma1) = 'sigma'
  sigma1
  cbind(rho1, sigma1, metric, method, est, se)
}

rho = seq(-0.75,0.75,0.25)
ps4_q2_result = mclapply(rho, function(x) func4(x,sigma, mc_rep))
ps4_q2_result = rbindlist(ps4_q2_result)
ps4_q2_result
