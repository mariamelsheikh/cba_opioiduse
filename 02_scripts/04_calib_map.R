
# MAP time depend and no time depend

library(here)
library(lhs) #Latin hypercube sampling

set.seed(123456)
source(here("02_scripts/03a_fun_calib_td.R"))
source(here("02_scripts/03b_fun_calib_ntd.R"))


# Functions

## sample_prior

sample_prior <- function(n, time_depend) { 
  
  if(time_depend == 1) {
    npar <- n_params
    dt_calib <- calib_param_tbl
    vec_names <- v_params_to_be_calib
  } else {
    npar <- n_params_ntd
    dt_calib <- calib_param_tbl_ntd
    vec_names <- v_params_to_be_calib_ntd
  }
  
  draws0 <- randomLHS(n = n, k = npar)
  
  draws <- data.frame(matrix(NA, nrow = n, ncol = npar,
                             dimnames = list(NULL, vec_names)))
  
  for (i in 1:n_params){
    
    draws[, i] <- assign(
      paste0(vec_names[i]),
      qunif(p = draws0[,i],
            min = dt_calib$lb[dt_calib$var_params_cycle %in% vec_names[i]],
            max = dt_calib$ub[dt_calib$var_params_cycle %in% vec_names[i]])
    )
  }
  return(as.matrix(draws))
}


## prior density

l_prior <- function(v_parms_calib, time_depend) {
  
  if(time_depend == 1) {
    dt_calib <- calib_param_tbl
    vec_names <- v_params_to_be_calib
  } else {
    dt_calib <- calib_param_tbl_ntd
    vec_names <- v_params_to_be_calib_ntd
  }
  
  if(is.null(dim(v_parms_calib))) 
    v_parms_calib <- t(v_parms_calib)
  
  lprior <- rep(0, nrow(v_parms_calib))
  
  for (i in 1:n_params){
    
    lprior <- lprior + dunif(x = v_parms_calib[,i],
                             min = dt_calib$lb[dt_calib$var_params_cycle %in% vec_names[i]],
                             max = dt_calib$ub[dt_calib$var_params_cycle %in% vec_names[i]],
                             log = T)
  }
  return(lprior)
}

## likelihood

l_likelihood <- function(v_parms_calib, time_depend) {
  
  if(is.null(dim(v_parms_calib))) v_parms_calib <- t(v_parms_calib)
  
  llik <- rep(0,nrow(v_parms_calib))
  
  for(j in 1:nrow(v_parms_calib)) {
    jj <- tryCatch( {
      #Run model
      
      if(time_depend == 1) {
        res_j <- mod_calib(as.numeric(v_parms_calib[j,]))
      } else {
        res_j <- mod_calib_ntd(as.numeric(v_parms_calib[j,]))
      }
      
      # calculate deaths likelihood
      llik[j] <- llik[j] + sum(dgamma(x = c(calib_target_tbl$target[calib_target_tbl$group == "total_deaths"]),
                                      shape = 1,
                                      rate = 1/res_j[["num_deaths"]],
                                      log = TRUE))
      
      # calculate prev likelihood
      
      llik[j] <- llik[j] + sum(dnorm(x = c(calib_target_tbl$dist_param_1[calib_target_tbl$group == "prev_on_opioidsrx"])/c(calib_target_tbl$dist_param_2[calib_target_tbl$group == "prev_on_opioidsrx"]),
                                     mean = res_j[["wmean_prop_rx_opioid"]],
                                     sd = res_j[["wmean_prop_rx_opioid"]],
                                     log = TRUE))
      
     
      # calculate oat counts likelihood
      
      llik[j] <- llik[j] + sum(dgamma(x = c(calib_target_tbl$target[calib_target_tbl$group == "total_oat"]),
                                      shape = 1,
                                      rate = 1/res_j[["num_oat"]],
                                      log = TRUE))
      
    }, error = function(e) NA)
    if(is.na(jj)) { llik[j] <- -Inf } 
  }
  return(llik)
}

# l_post

l_post <- function(par_vector) {
  return(l_prior(par_vector, time_depend = 1) + l_likelihood(par_vector, time_depend = 1))
}

l_post_notimedep <- function(par_vector) {
  return(l_prior(par_vector, time_depend = 0) + l_likelihood(par_vector, time_depend = 0))
}




# MAP

Sys.time() # "2023-04-18 00:44:51 EDT"
opt_params_notimedepend <- suppressMessages(optim(calib_param_tbl_ntd$basevalue, 
                                                  l_post_notimedep, control = list(fnscale = -1),
                                     method = "Nelder-Mead"))
Sys.time() # "2023-04-18 00:45:36 EDT"

opt_params <- suppressMessages(optim(calib_param_tbl$basevalue, 
                                                  l_post, control = list(fnscale = -1),
                                                  method = "Nelder-Mead"))
Sys.time() # "2023-04-18 00:50:09 EDT"


saveRDS(opt_params, file = here("01_data/map_calib_td.RDS"))
saveRDS(opt_params_notimedepend, file = here("01_data/map_calib_ntd.RDS"))




