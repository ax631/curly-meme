library(doParallel)
library(foreach)
library(parallel)
library(data.table)
library(tidyverse)
source('ps4_q2_funcs.R')

ncores=4
cl = makeCluster(ncores)
registerDoParallel(cl)

f= function(sigma)
{
  library(data.table)
  library(parallel)
  library(tidyverse)
  source('ps4_q2_funcs.R')
  func4= function(rho, sigma)
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
    P = sim_beta(X, beta, sigma =sigma, mc_rep = 100)
    
    # evaluate: -------------------------------------------------------------------
    P_eval = evaluate(P, 1:10)
    # P_eval
    
    
    all0 =
      lapply( c('holm', 'bonferroni', 'BH', 'BY'), function(x){
        evaluate( apply(P, 2, p.adjust, method = x), tp_ind = 1:10)
      })
    all = rbindlist(all0)
    all[ , method := c('holm', 'bonferroni', 'BH', 'BY') ]
    # all
    est = as.data.frame(t(cbind(t(all$fwer), t(all$fdr), t(all$sens), t(all$spec))))
    colnames(est) = 'est'
    se = as.data.frame(t(cbind(t(all$fwer_se), t(all$fdr_se), t(all$sens_se), t(all$spec_se))))
    colnames(se) = 'se'
    metric = as.data.frame(c(rep("fwer",4), rep("fdr",4), rep("sens",4), rep("spec",4)))
    colnames(metric) = 'metric'
    method = as.data.frame(rep(c("holm","bonferroni","BH","BY"),4))
    colnames(method) = 'method'
    rho1 = as.data.frame(rep(rho, 16))
    colnames(rho1) = 'rho'
    sigma1 = as.data.frame(rep(sigma,16))
    colnames(sigma1) = 'sigma'
    cbind(rho1, sigma1, metric, method, est, se)
  }
  
  ps4_q2_result = lapply(rho, function(x) func4(x,sigma))
  ps4_q2_result = rbindlist(ps4_q2_result)
  ps4_q2_result
}


# Parameters: -----------------------------------------------------------------
n= 1000; p = 100
beta = c(rep(0.1, 10),rep(0,p-10) )
dim(beta) = c(p,1)
rho  = seq(-0.75,0.75,0.25)
sigma = c(0.25,0.5,1)
results_q4b = foreach(sigma = c(0.25, 0.5, 1)) %dopar% { f(sigma)
}

rm(list = ls())
# results_q4b = foreach(rho = seq(-0.75,0.75,0.25)) %dopar% {
results_q4b = rbindlist(results_q4b)