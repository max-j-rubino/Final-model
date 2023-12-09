Replicates = function(sd_beta_ind_v, alpha_v,vr_V) {  
  
  replicates = expand.grid(sd_beta_ind_v, alpha_v,vr_V)
  colnames(replicates) = c('sd_beta_ind_v', 'alpha_v','vr_V')
  
  return(replicates)
} 