#' CARMA_spike_slab_noEM_rcpp Function
#'
#' Perform CARMA analysis using a Spike-and-Slab prior without Expectation-Maximization (EM).
#'
#' @param z Numeric vector representing z-scores.
#' @param ld Numeric matrix representing the linkage disequilibrium (LD) matrix.
#' @param lambda Numeric, regularization parameter controlling the strength of the L1 penalty.
#' @param rho.index Numeric, threshold for computing credible sets in CARMA analysis.
#' @param Max.Model.Dim Numeric, maximum allowed dimension for the causal model.
#' @param all.iter Integer, the total number of iterations to run the CARMA analysis.
#' @param all.inner.iter Integer, the number of inner iterations in each CARMA iteration.
#' @param epsilon.threshold Numeric, threshold for convergence in CARMA iterations.
#' @param num.causal Integer, number of causal variants to be selected in the final model.
#' @param y.var Numeric, variance parameter for the marginal distribution of the response variable.
#' @param tau Numeric, tuning parameter controlling the degree of sparsity in the Spike-and-Slab prior.
#' @param outlier.switch Logical, whether to consider outlier models in the analysis.
#' @param outlier.BF.index Numeric, Bayes Factor threshold for identifying outliers.
#' @param path_to_src Character, the path to the source files.
#' @param sparse_fun Logical, indicating whether to use sparse matrix operations for efficiency.
#'
#' @return A list containing the following components:
#' \item{PIPs}{Vector of Posterior Inclusion Probabilities (PIPs) for each genetic variant.}
#' \item{Credible_set}{Matrix of variants in the credible set based on the given threshold.}
#' \item{Credible_model}{NULL, as this component is not implemented in the function.}
#' \item{Outliers}{Data frame containing information about identified outlier models.}
#' \item{all.C.list}{List containing intermediate results of the CARMA analysis.}
#'
#' @examples
#' # Example usage:
#' result <- CARMA_spike_slab_noEM_rcpp(z = your_genetic_data, ld = your_ld_matrix, path_to_src = "path/to/source/files")
#'
#' @references
#' Provide relevant references if applicable.
#'
#' @seealso
#' Other functions related to CARMA analysis.
#'
#' @keywords
#' CARMA, Spike-and-Slab, Bayesian analysis
#'
#' @export
CARMA_spike_slab_noEM_rcpp <- function(z, ld, lambda = 1, rho.index = 0.99, Max.Model.Dim = 2e+5, 
                                       all.iter = 3, all.inner.iter = 10, epsilon.threshold = 1e-5, num.causal = 10, 
                                       y.var = 1, tau = 0.04, outlier.switch = TRUE, outlier.BF.index = 1/3.2, 
                                       path_to_src, sparse_fun = FALSE) {
  
  t0=Sys.time()
  
  #effect.size.prior='Spike-slab'
  
  library(Rcpp)
  source(paste0(path_to_src,"/CARMA_miscs.R"))
  sourceCpp(paste0(path_to_src,"/ind_Normal_fixed_sigma_marginal.cpp"))
  sourceCpp(paste0(path_to_src,"/outlier_ind_Normal_marginal.cpp"))
  
  if (sparse_fun){
    library(Matrix)
    source(paste0(path_to_src,"/new_MCS_fun_sparse_matrix.R"))
  } else {
    source(paste0(path_to_src,"/new_MCS_fun_vector.R"))
  }
  
  
  p_snp=length(z)
  epsilon.list<-epsilon.threshold*p_snp
  all.epsilon.threshold<-epsilon.threshold*p_snp
  
  print("ordering pizza...")
  all.C.list<-MCS_modified(z,ld,epsilon=epsilon.list,
                           Max.Model.Dim=Max.Model.Dim,lambda = lambda,
                           outlier.switch=outlier.switch,tau=tau,
                           num.causal = num.causal,y.var=y.var,
                           inner.all.iter = all.inner.iter,outlier.BF.index=outlier.BF.index)
  ########Running CARMA######## 
  for(g in 1:all.iter){ 
    
    ac1=all.C.list[[1]][[1]]
    previous.result<-mean(ac1[1:round(length(ac1)/4)])
    
    
    all.C.list<-MCS_modified(z=z,ld,input.conditional.S.list = all.C.list[[4]],
                             Max.Model.Dim=Max.Model.Dim,
                             y.var=y.var,num.causal = num.causal,epsilon=epsilon.list,
                             outlier.switch=outlier.switch,tau=tau,
                             lambda = lambda,
                             inner.all.iter = all.inner.iter,
                             outlier.BF.index=outlier.BF.index)

    ac1=all.C.list[[1]][[1]]
    difference<-abs(previous.result-mean(ac1[1:round(length(ac1)/4)]))
    if(difference<all.epsilon.threshold){
      break
    }
  }
  
  results.list<-list()
  pip=PIP.func(all.C.list[[1]][[1]],all.C.list[[1]][[2]],p_snp)
  credible.set<-credible.set.fun.improved(pip,ld,rho=rho.index)
  credible.model<-NULL
  results.list[[1]]<-pip
  results.list[[2]]<-credible.set
  results.list[[3]]<-credible.model
  results.list[[4]]<-all.C.list[[4]]
  
  results.list[[5]]<-all.C.list
  
  names(results.list)<-c('PIPs','Credible_set','Credible_model','Outliers','all.C.list')
  
  t1=Sys.time()-t0
  print("pizza ordered!")
  print(paste0('FM time is ',round(t1,2)," sec"))
  
  return(results.list)
  
}


