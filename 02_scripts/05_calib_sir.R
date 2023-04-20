library(here)
source(here("02_scripts/03a_fun_calib_td.R"))
set.seed(123456)

# Make changes from MAP then run SIR

trans_prob_tbl_new <- cbind.data.frame(rep(v_state_names, each = n_states), 
                                       rep(v_state_names, n_states))

names(trans_prob_tbl_new) <- c("var_from", "var_to")

trans_prob_tbl_new <- trans_prob_tbl_new %>% 
  full_join(., trans_prob_tbl %>% 
              select(var_from, var_to, basevalue), by = c("var_from", "var_to")) %>% 
  replace(is.na(.), 0)


trans_prob_tbl_new$basevalue[trans_prob_tbl_new$var_from == "BN_PN" &
                               trans_prob_tbl_new$var_to == "BPO_MISUSE"] <- 0.00220011
trans_prob_tbl_new$basevalue[trans_prob_tbl_new$var_from == "BPO_CHRONIC" &
                               trans_prob_tbl_new$var_to == "BO_OD_RX"] <- 0.00160398
trans_prob_tbl_new$basevalue[trans_prob_tbl_new$var_from == "BS_DETOX" &
                               trans_prob_tbl_new$var_to == "BI_ILLICIT"] <- 0.22181118
trans_prob_tbl_new$basevalue[trans_prob_tbl_new$var_from == "BS_DETOX" &
                               trans_prob_tbl_new$var_to == "BS_OAT_INI"] <- 0.82000000
trans_prob_tbl_new$basevalue[trans_prob_tbl_new$var_from == "BPO_PALLIATIVE" &
                               trans_prob_tbl_new$var_to == "BO_DEATH"] <- 0.38250000
trans_prob_tbl_new$basevalue[trans_prob_tbl_new$var_from == "BS_OAT_MAINT" &
                               trans_prob_tbl_new$var_to == "BS_OAT_MAINT"] <- 0.71000000


m_P_map <- matrix(trans_prob_tbl_new$basevalue, byrow = T,
                   nrow = n_states, ncol = n_states,
                   dimnames = list(v_state_names, v_state_names))


for (i in 1:n_states){ 
  for (j in 1:n_states){
    m_P_map[i, j] <- ifelse(m_P_map[i, j] == 999, 1 - (rowSums(m_P_map)[i] - 999), m_P_map[i, j])
  }
}

# arr_P_map <- gen_trans_prob_arr(mat_P = m_P_map)


calib_param_tbl$basevalue[calib_param_tbl$var_params == "p__BN_PN__BPO_MISUSE"] <- 0.00220011
calib_param_tbl$lb[calib_param_tbl$var_params == "p__BN_PN__BPO_MISUSE"] <- 0.00220011*0.95
calib_param_tbl$ub[calib_param_tbl$var_params == "p__BN_PN__BPO_MISUSE"] <- 0.00220011*1.05

calib_param_tbl$basevalue[calib_param_tbl$var_params == "p__BPO_CHRONIC__BO_OD_RX"] <- 0.00160398
calib_param_tbl$lb[calib_param_tbl$var_params == "p__BPO_CHRONIC__BO_OD_RX"] <- 0.00160398*0.95
calib_param_tbl$ub[calib_param_tbl$var_params == "p__BPO_CHRONIC__BO_OD_RX"] <- 0.00160398*1.05

calib_param_tbl$basevalue[calib_param_tbl$var_params == "p__BS_DETOX__BI_ILLICIT"] <- 0.22181118
calib_param_tbl$lb[calib_param_tbl$var_params == "p__BS_DETOX__BI_ILLICIT"] <- 0.22181118*0.95
calib_param_tbl$ub[calib_param_tbl$var_params == "p__BS_DETOX__BI_ILLICIT"] <- 0.22181118*1.05

calib_param_tbl$basevalue[calib_param_tbl$var_params == "p__BS_DETOX__BS_OAT_INI"] <- 0.82000000
calib_param_tbl$lb[calib_param_tbl$var_params == "p__BS_DETOX__BS_OAT_INI"] <- 0.82000000*0.95
calib_param_tbl$ub[calib_param_tbl$var_params == "p__BS_DETOX__BS_OAT_INI"] <- 0.82000000*1.05

calib_param_tbl$basevalue[calib_param_tbl$var_params == "p__BPO_PALLIATIVE__BO_DEATH"] <- 0.38250000
calib_param_tbl$lb[calib_param_tbl$var_params == "p__BPO_PALLIATIVE__BO_DEATH"] <- 0.38250000*0.95
calib_param_tbl$ub[calib_param_tbl$var_params == "p__BPO_PALLIATIVE__BO_DEATH"] <- 0.38250000*1.05

calib_param_tbl$basevalue[calib_param_tbl$var_params == "p__BS_OAT_MAINT__BS_OAT_MAINT"] <- 0.71000000
calib_param_tbl$lb[calib_param_tbl$var_params == "p__BS_OAT_MAINT__BS_OAT_MAINT"] <- 0.71000000*0.95
calib_param_tbl$ub[calib_param_tbl$var_params == "p__BS_OAT_MAINT__BS_OAT_MAINT"] <- 0.71000000*1.05

# Functions

## sample_prior

sample_prior <- function(n) { 
  
  draws0 <- randomLHS(n = n, k = n_params)
  
  draws <- data.frame(matrix(NA, nrow = n, ncol = n_params,
                             dimnames = list(NULL, v_params_to_be_calib)))
  
  for (i in 1:n_params){
    
    draws[, i] <- assign(
      paste0(v_params_to_be_calib[i]),
      qunif(p = draws0[,i],
            min = calib_param_tbl$lb[calib_param_tbl$var_params_cycle %in% v_params_to_be_calib[i]],
            max = calib_param_tbl$ub[calib_param_tbl$var_params_cycle %in% v_params_to_be_calib[i]])
    )
  }
  return(as.matrix(draws))
}


## likelihood

l_likelihood <- function(v_parms_calib) {
  
  if(is.null(dim(v_parms_calib))) v_parms_calib <- t(v_parms_calib)
  
  llik <- rep(0,nrow(v_parms_calib))
  
  for(j in 1:nrow(v_parms_calib)) {
    jj <- tryCatch( {
      #Run model

      res_j <- mod_calib(as.numeric(v_parms_calib[j,]))
      
      # calculate deaths likelihood
      llik[j] <- llik[j] + sum(dgamma(x = c(calib_target_tbl$target[calib_target_tbl$group == "total_deaths"]),
                                      shape = 1,
                                      rate = 1/res_j[["num_deaths"]],
                                      log = TRUE))
      
      # calculate prev likelihood
      
      llik[j] <- llik[j] + sum(dnorm(x = c(calib_target_tbl$dist_param_1[calib_target_tbl$group == "prev_on_opioidsrx"])/
                                       c(calib_target_tbl$dist_param_2[calib_target_tbl$group == "prev_on_opioidsrx"]),
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



# SIR


# For each parameter set in samp, calculate the log-likelihood
n_samples <- 10000
uncalib_samples <- sample_prior(n_samples)
v_llik_samps <- rep(NA,n_samples) #Pre-allocate vector

library(progress)
pb <- progress_bar$new(
  format = "  downloading [:bar] :percent :elapsed eta: :eta",
  total = n_samples, clear = FALSE, width= 60)

for(i in 1:n_samples){
  pb$tick()
  Sys.sleep(1 / n_samples)
  
  v_llik_samps[i] <- suppressWarnings(suppressMessages((l_likelihood(uncalib_samples[i, ]))))
  
  if(i/100==round(i/100,0)) { 
    cat('\r',paste(i/n_samples*100,"% done",sep="")) 
  } 
} # ~2hrs

lik_samps <- exp(v_llik_samps - max(v_llik_samps))

m_resamples <- uncalib_samples[sample(length(lik_samps),
                                      size = length(lik_samps), 
                                      prob = lik_samps,
                                      replace=TRUE),]

#Calculate the number of unique parameter sets
SIR_unique_sets <- nrow(unique(m_resamples))
SIR_unique_sets # 6398

saveRDS(m_resamples, file = here("01_data/calib_samples_sir.RDS"))


uncalib_resamples <- uncalib_samples[sample(n_samples,
                                            size = n_samples,
                                            replace=TRUE),]
SIR_unique_sets_uncalib <- nrow(unique(uncalib_resamples))
SIR_unique_sets_uncalib # 6296

saveRDS(uncalib_resamples, file = here("01_data/uncalib_sample_sir.RDS"))

v_yrs_prev_rx <- c(2015:2018)
yrs_prev_rx <- length(v_yrs_prev_rx)
v_yrs_deaths <- c(2016:2020)
yrs_deaths <- length(v_yrs_deaths)
v_yrs_oddeaths <- c(2016:2021)
yrs_oddeaths <- length(v_yrs_oddeaths)
v_yrs_oat <- c(2018:2021)
yrs_oat <- length(v_yrs_oat)

m_prev_rx_use_uncalib <- m_prev_rx_use <- matrix(0, nrow = n_samples, ncol = yrs_prev_rx,
                        dimnames = list(samp = 1:n_samples,
                                        year = v_yrs_prev_rx))

m_deaths_uncalib <- m_deaths <- matrix(0, nrow = n_samples, ncol = yrs_deaths,
                   dimnames = list(samp = 1:n_samples,
                                   year = v_yrs_deaths))

m_od_deaths_uncalib <- m_od_deaths <- matrix(0, nrow = n_samples, ncol = yrs_oddeaths,
                      dimnames = list(samp = 1:n_samples,
                                      year = v_yrs_oddeaths))
m_oat_uncalib <- m_oat <- matrix(0, nrow = n_samples, ncol = yrs_oat,
                dimnames = list(samp = 1:n_samples,
                                year = v_yrs_oat))

pb <- progress_bar$new(
  format = "  downloading [:bar] :percent :elapsed eta: :eta",
  total = nrow(m_resamples), clear = FALSE, width= 60)

for (j in 1:nrow(m_resamples)){
  
  pb$tick()
  Sys.sleep(1 / nrow(m_resamples))
  
  calibmod <- suppressWarnings(suppressMessages(mod_calib(as.numeric(m_resamples[j, ]))))
  calibmod_uncalib <- suppressWarnings(suppressMessages(mod_calib(as.numeric(uncalib_resamples[j, ]))))
  
  
  for (i in 1:yrs_deaths){
    yr <- (v_yrs_deaths[1] - 1) + i
    m_deaths[j, i] <- calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                              year_mon_cycle_tbl$mon == 12] + 1,
                                   "BO_DEATH"] -
      calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
                                              year_mon_cycle_tbl$mon == 12] + 1,
                   "BO_DEATH"]
    
    m_deaths_uncalib[j, i] <- calibmod_uncalib$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                              year_mon_cycle_tbl$mon == 12] + 1,
                                   "BO_DEATH"] -
      calibmod_uncalib$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
                                              year_mon_cycle_tbl$mon == 12] + 1,
                   "BO_DEATH"]
  }
  
  for (i in 1:yrs_oddeaths){
    yr <- (v_yrs_oddeaths[1] - 1) + i 
    m_od_deaths[j, i] <- calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                                 year_mon_cycle_tbl$mon == 12] + 1,
                                      "BO_OD_DEATH"] -
      calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
                                              year_mon_cycle_tbl$mon == 12] + 1,
                   "BO_OD_DEATH"]
    m_od_deaths_uncalib[j, i] <- calibmod_uncalib$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                                 year_mon_cycle_tbl$mon == 12] + 1,
                                      "BO_OD_DEATH"] -
      calibmod_uncalib$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
                                              year_mon_cycle_tbl$mon == 12] + 1,
                   "BO_OD_DEATH"]
  }
  
  # Calculate prev
  
  prop_opioids_rx_uncalib <- prop_opioids_rx <- data.frame(matrix(NA, byrow = T, nrow = 12, ncol = yrs_prev_rx,
                                       list(NULL, paste0("prop_opioids_rx_", str_sub(v_yrs_prev_rx, 3, 4)))))
  
  tot_pop_tbl_uncalib <- tot_pop_tbl <- data.frame(matrix(NA, byrow = T, nrow = 12, ncol = yrs_prev_rx,
                                   list(NULL, paste0("tot_pop_", str_sub(v_yrs_prev_rx, 3, 4)))))
  
  for (i in 1:yrs_prev_rx) {
    
    yr <- (v_yrs_prev_rx[1] - 1) + i
    prop_opioids_rx[, i] <- assign(
      paste0("prop_opioids_rx_", str_sub(yr, 3, 4)),
      rowSums(calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,
                           c("BPO_ACUTE", "BPO_CHRONIC", "BPO_CANCER",
                             "BPO_OTHER", "BPO_PALLIATIVE")]) /
        rowSums(calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
      
    )
    
    tot_pop_tbl[, i] <- assign(
      paste0("tot_pop_", str_sub(yr, 3, 4)),
      rowSums(calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
    )
    
    
    prop_opioids_rx_uncalib[, i] <- assign(
      paste0("prop_opioids_rx_", str_sub(yr, 3, 4)),
      rowSums(calibmod_uncalib$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,
                           c("BPO_ACUTE", "BPO_CHRONIC", "BPO_CANCER",
                             "BPO_OTHER", "BPO_PALLIATIVE")]) /
        rowSums(calibmod_uncalib$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
      
    )
    
    tot_pop_tbl_uncalib[, i] <- assign(
      paste0("tot_pop_", str_sub(yr, 3, 4)),
      rowSums(calibmod_uncalib$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
    )
  }
  
  
  prop_opioids_rx <- prop_opioids_rx %>% 
    pivot_longer(1:4, names_to = "grp", values_to = "target_val") %>% 
    arrange(grp) %>% 
    mutate(year = rep(v_yrs_prev_rx, each = 12),
           year_mon = paste(year, month.abb[1:12], sep = "_")) %>% 
    bind_cols(., tot_pop_tbl %>% 
                pivot_longer(1:4, names_to = "tot_pop", values_to = "tot_pop_val") %>% 
                arrange(tot_pop) %>% dplyr::select(-tot_pop)) %>% 
    group_by(year) %>% 
    summarise(target_val = weighted.mean(target_val, tot_pop_val)) %>% 
    ungroup()
  
  m_prev_rx_use[j, ] <- prop_opioids_rx$target_val
  
  
  prop_opioids_rx_uncalib <- prop_opioids_rx_uncalib %>% 
    pivot_longer(1:4, names_to = "grp", values_to = "target_val") %>% 
    arrange(grp) %>% 
    mutate(year = rep(v_yrs_prev_rx, each = 12),
           year_mon = paste(year, month.abb[1:12], sep = "_")) %>% 
    bind_cols(., tot_pop_tbl_uncalib %>% 
                pivot_longer(1:4, names_to = "tot_pop", values_to = "tot_pop_val") %>% 
                arrange(tot_pop) %>% dplyr::select(-tot_pop)) %>% 
    group_by(year) %>% 
    summarise(target_val = weighted.mean(target_val, tot_pop_val)) %>% 
    ungroup()
  
  m_prev_rx_use_uncalib[j, ] <- prop_opioids_rx_uncalib$target_val
  
  
  # calculate oat
  for (i in 1:yrs_oat){
    yr <- (v_yrs_oat[1] - 1) + i 
    m_oat[j, i] <- sum(calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                               year_mon_cycle_tbl$mon == 6] + 1,
                                    c("BS_OAT_INI", "BS_OAT_MAINT", "BR_OAT_INI", "BR_OAT_MAINT")])
    
    m_oat_uncalib[j, i] <- sum(calibmod_uncalib$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                               year_mon_cycle_tbl$mon == 6] + 1,
                                    c("BS_OAT_INI", "BS_OAT_MAINT", "BR_OAT_INI", "BR_OAT_MAINT")])
    
  }
}

saveRDS(m_prev_rx_use, file = here("01_data/prev_outcome_data_sir.RDS"))
saveRDS(m_od_deaths, file = here("01_data/oddeath_outcome_data_sir.RDS"))
saveRDS(m_deaths, file = here("01_data/death_outcome_data_sir.RDS"))
saveRDS(m_oat, file = here("01_data/oat_outcome_data_sir.RDS"))

saveRDS(m_prev_rx_use_uncalib, file = here("01_data/prev_outcome_data_sir_uncalib.RDS"))
saveRDS(m_od_deaths_uncalib, file = here("01_data/oddeath_outcome_data_sir_uncalib.RDS"))
saveRDS(m_deaths_uncalib, file = here("01_data/death_outcome_data_sir_uncalib.RDS"))
saveRDS(m_oat_uncalib, file = here("01_data/oat_sir_uncalib.RDS"))


# uncalib sample ----------------------------------------------------------


# uncalib_resamples <- uncalib_samples[sample(n_samples,
#                                             size = n_samples,
#                                             replace=TRUE),]
# saveRDS(uncalib_resamples, file = "uncalib_sample_sir.RDS")
# 
# m_prev_rx_use_uncalib <- matrix(0, nrow = n_samples, ncol = 4,
#                                 dimnames = list(samp = 1:n_samples,
#                                                 year = c("2015", "2016", "2017", "2018")))
# 
# m_deaths_uncalib <- matrix(0, nrow = n_samples, ncol = 5,
#                            dimnames = list(samp = 1:n_samples,
#                                            year = c("2016", "2017", "2018", "2019", "2020")))
# 
# 
# m_od_deaths_uncalib <- matrix(0, nrow = n_samples, ncol = 6,
#                               dimnames = list(samp = 1:n_samples,
#                                               year = c("2016", "2017", "2018", "2019",
#                                                        "2020", "2021")))
# 
# m_cancer_uncalib <- matrix(0, nrow = n_samples, ncol = 4,
#                            dimnames = list(samp = 1:n_samples,
#                                            year = c("2015","2016", "2017", "2018")))
# m_oat_uncalib <- matrix(0, nrow = n_samples, ncol = 4,
#                         dimnames = list(samp = 1:n_samples,
#                                         year = c("2018", "2019", "2020","2021")))
# 
# 
# pb <- progress_bar$new(
#   format = "  downloading [:bar] :percent :elapsed eta: :eta",
#   total = n_samples, clear = FALSE, width= 60)
# 
# for (j in 1:nrow(uncalib_samples)){
# 
#   pb$tick()
#   Sys.sleep(1 / nrow(uncalib_samples))
# 
#   calibmod <- suppressWarnings(suppressMessages(mod_calib(as.numeric(uncalib_resamples[j, ]))))
# 
# 
#   for (i in 1:length(2016:2020)){
#     yr <- 2015 + i
#     m_deaths_uncalib[j, i] <- calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
#                                                                       year_mon_cycle_tbl$mon == 12] + 1,
#                                            "BO_DEATH"] -
#       calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
#                                               year_mon_cycle_tbl$mon == 12] + 1,
#                    "BO_DEATH"]
#   }
# 
#   for (i in 1:length(2016:2021)){
#     yr <- 2015 + i
#     m_od_deaths_uncalib[j, i] <- calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
#                                                                          year_mon_cycle_tbl$mon == 12] + 1,
#                                               "BO_OD_DEATH"] -
#       calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
#                                               year_mon_cycle_tbl$mon == 12] + 1,
#                    "BO_OD_DEATH"]
#   }
# 
#   # Calculate prev
# 
#   prop_opioids_rx <- data.frame(matrix(NA, byrow = T, nrow = 12, ncol = length(2015:2018),
#                                        list(NULL, paste0("prop_opioids_rx_", str_sub(2015:2018, 3, 4)))))
# 
#   tot_pop_tbl <- data.frame(matrix(NA, byrow = T, nrow = 12, ncol = length(2015:2018),
#                                    list(NULL, paste0("tot_pop_", str_sub(2015:2018, 3, 4)))))
# 
#   for (i in 1:length(2015:2018)) {
# 
#     yr <- 2014 + i
#     prop_opioids_rx[, i] <- assign(
#       paste0("prop_opioids_rx_", str_sub(yr, 3, 4)),
#       rowSums(calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,
#                            c("BPO_ACUTE", "BPO_CHRONIC", "BPO_CANCER",
#                              "BPO_OTHER", "BPO_PALLIATIVE")]) /
#         rowSums(calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
# 
#     )
# 
#     tot_pop_tbl[, i] <- assign(
#       paste0("tot_pop_", str_sub(yr, 3, 4)),
#       rowSums(calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
#     )
#   }
# 
#   prop_opioids_rx <- prop_opioids_rx %>%
#     pivot_longer(1:4, names_to = "grp", values_to = "target_val") %>%
#     arrange(grp) %>%
#     mutate(year = rep(2015:2018, each = 12),
#            year_mon = paste(year, month.abb[1:12], sep = "_")) %>%
#     bind_cols(., tot_pop_tbl %>%
#                 pivot_longer(1:4, names_to = "tot_pop", values_to = "tot_pop_val") %>%
#                 arrange(tot_pop) %>% dplyr::select(-tot_pop)) %>%
#     group_by(year) %>%
#     summarise(target_val = weighted.mean(target_val, tot_pop_val)) %>%
#     ungroup()
# 
#   m_prev_rx_use_uncalib[j, ] <- prop_opioids_rx$target_val
# 
# 
#   # calculate cancer
# 
#   for (i in 1:length(2015:2018)){
#     yr <- 2014 + i
#     m_cancer_uncalib[j, i] <- sum(calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
#                                                                           year_mon_cycle_tbl$mon == 6] + 1,
#                                                c("BN_CANCER", "BPO_CANCER")])
#   }
# 
#   # calculate oat
#   for (i in 1:length(2018:2021)){
#     yr <- 2017 + i
#     m_oat_uncalib[j, i] <- sum(calibmod$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
#                                                                        year_mon_cycle_tbl$mon == 6] + 1,
#                                             c("BS_OAT_INI", "BS_OAT_MAINT", "BR_OAT_INI", "BR_OAT_MAINT")])
# 
#   }
# 
# }
# 
# saveRDS(m_prev_rx_use_uncalib, file = "prev_outcome_data_sir_uncalib_td.RDS")
# saveRDS(m_od_deaths_uncalib, file = "oddeath_outcome_data_sir_uncalib_td.RDS")
# saveRDS(m_deaths_uncalib, file = "death_outcome_data_sir_uncalib_td.RDS")
# saveRDS(m_cancer_uncalib, file = "cancer_sir_unclaib_td.RDS")
# saveRDS(m_oat_uncalib, file = "oat_sir_uncalib_td.RDS")
# nrow(unique(uncalib_resamples))




