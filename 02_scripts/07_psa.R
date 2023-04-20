# PSA
library(here)
strt <- Sys.time()
source(here("01_fun_data.R"))
sir_sample <- readRDS(here("01_data/calib_samples_sir.RDS"))


ss_t1 <- year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == 2016 & year_mon_cycle_tbl$mon == 1] 
ss_t2 <- year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == 2020 & year_mon_cycle_tbl$mon == 3]
ss_t3 <- year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == 2023 & year_mon_cycle_tbl$mon == 1]

nalx_t1 <- year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == 2017 &
                                      year_mon_cycle_tbl$mon == 1]


gen_arr_nalox <- function(a_P) {
  
  arr_P_nalx <- a_P
  
  for (i in nalx_t1:180){
    m_p_this_cycle_int <- arr_P_nalx[,,i]
    
    m_p_this_cycle_int["BO_OD_ILLICIT", "BO_OD_DEATH"] <- m_p_this_cycle_int["BO_OD_ILLICIT", "BO_OD_DEATH"]*0.95
    m_p_this_cycle_int["BO_OD_ILLICIT", "BI_ILLICIT"] <- 1-(rowSums(m_p_this_cycle_int)["BO_OD_ILLICIT"] - m_p_this_cycle_int["BO_OD_ILLICIT", "BI_ILLICIT"])
    
    m_p_this_cycle_int["BO_OD_RX", "BO_OD_DEATH"] <- m_p_this_cycle_int["BO_OD_RX", "BO_OD_DEATH"]*0.97
    m_p_this_cycle_int["BO_OD_RX", "BPO_MISUSE"] <- 1-(rowSums(m_p_this_cycle_int)["BO_OD_RX"] - m_p_this_cycle_int["BO_OD_RX", "BPO_MISUSE"])
    
    arr_P_nalx[,,i] <- m_p_this_cycle_int
  }
  
  return(arr_P_nalx)
  
}

gen_arr_ss <- function(a_P){
  
  arr_P_ss <- a_P
  
  
  for (i in ss_t1:180){
    
    m_p_this_cycle_int <- arr_P_ss[,,i]
    
    if (i == ss_t1){
      
      m_p_this_cycle_int["BS_OAT_SS", "BN_PN"] <- 0.0046 # 20216
      m_p_this_cycle_int["BS_OAT_SS", "BI_ILLICIT"] <- 0.0438 # 2016
      m_p_this_cycle_int["BS_OAT_SS", "BO_OD_ILLICIT"] <- 0.0021 #2016
      m_p_this_cycle_int["BS_OAT_SS", "BO_DEATH"] <- m_p_this_cycle_int["BPO_CHRONIC", "BO_DEATH"]
      m_p_this_cycle_int["BS_OAT_SS", "BS_OAT_SS"] <- 1 - rowSums(m_p_this_cycle_int)["BS_OAT_SS"]
      
      
      m_p_this_cycle_int["BS_SS", "BN_PN"] <- 0.0046 # 2016
      m_p_this_cycle_int["BS_SS", "BI_ILLICIT"] <- 0.0438 # 2016
      m_p_this_cycle_int["BS_SS", "BO_OD_ILLICIT"] <- 0.0024 # 2016
      m_p_this_cycle_int["BS_SS", "BO_DEATH"] <- m_p_this_cycle_int["BPO_CHRONIC", "BO_DEATH"]
      m_p_this_cycle_int["BS_SS", "BS_SS"] <- 1 - rowSums(m_p_this_cycle_int)["BS_SS"]
      
      arr_P_ss[,,i] <- m_p_this_cycle_int
    } else if (i >= ss_t3) {
      
      m_p_this_cycle_int["BS_OAT_INI", "BS_OAT_SS"] <- 0.003  # starts in 2023
      m_p_this_cycle_int["BS_OAT_INI", "BI_ILLICIT"] <-  1 - (rowSums(m_p_this_cycle_int)["BS_OAT_INI"] - m_p_this_cycle_int["BS_OAT_INI", "BI_ILLICIT"])
      
      
      m_p_this_cycle_int["BS_OAT_MAINT", "BS_OAT_SS"] <- 0.003 # starts in 2023
      m_p_this_cycle_int["BS_OAT_MAINT", "BI_ILLICIT"] <-  1 - (rowSums(m_p_this_cycle_int)["BS_OAT_MAINT"] - m_p_this_cycle_int["BS_OAT_MAINT", "BI_ILLICIT"])
      
      
      m_p_this_cycle_int["BO_OD_ILLICIT", "BS_OAT_SS"] <- 0.005 #2023
      m_p_this_cycle_int["BO_OD_ILLICIT", "BS_SS"] <- 0.005 # 2023
      m_p_this_cycle_int["BO_OD_ILLICIT", "BI_ILLICIT"] <- 1 - (rowSums(m_p_this_cycle_int)["BO_OD_ILLICIT"] - m_p_this_cycle_int["BO_OD_ILLICIT", "BI_ILLICIT"])
      
      
      m_p_this_cycle_int["BO_MOD_BI", "BR_OAT_SS"] <- 0.01 # in 2023
      m_p_this_cycle_int["BO_MOD_BI", "BR_SS"] <- 0.001 # in 2023
      m_p_this_cycle_int["BO_MOD_BI", "BR_ILLICIT"] <- 1 - (rowSums(m_p_this_cycle_int)["BO_MOD_BI"] - m_p_this_cycle_int["BO_MOD_BI", "BR_ILLICIT"]) # in 2023
      
      m_p_this_cycle_int["BR_ILLICIT", "BR_SS"] <- 0.005 #2023
      m_p_this_cycle_int["BR_ILLICIT", "BR_ILLICIT"] <- 1 - (rowSums(m_p_this_cycle_int)["BR_ILLICIT"] - m_p_this_cycle_int["BR_ILLICIT", "BR_ILLICIT"])
      
      
      m_p_this_cycle_int["BR_OAT_SS", "BR_ILLICIT"] <- 0.01 # 2023
      m_p_this_cycle_int["BR_OAT_SS", "BR_OD_ILLICIT"] <- 0.0182 # 2023
      m_p_this_cycle_int["BR_OAT_SS", "BO_DEATH"] <- m_p_this_cycle_int["BR_OAT_MAINT", "BO_DEATH"]
      m_p_this_cycle_int["BR_OAT_SS", "BR_OAT_SS"] <- 1 - rowSums(m_p_this_cycle_int)["BR_OAT_SS"]
      
      
      m_p_this_cycle_int["BR_SS", "BR_ILLICIT"] <- 0.01 # 2023
      m_p_this_cycle_int["BR_SS", "BR_OD_ILLICIT"] <- 0.0182 # 2023
      m_p_this_cycle_int["BR_SS", "BO_DEATH"] <- m_p_this_cycle_int["BR_OAT_MAINT", "BO_DEATH"]
      m_p_this_cycle_int["BR_SS", "BR_SS"] <- 1 - rowSums(m_p_this_cycle_int)["BR_SS"]
      
      
      m_p_this_cycle_int["BR_OAT_MAINT", "BR_OAT_SS"] <- 0.005 #2023
      m_p_this_cycle_int["BR_OAT_MAINT", "BR_OAT_MAINT"] <- 1 - (rowSums(m_p_this_cycle_int)["BR_OAT_MAINT"] - m_p_this_cycle_int["BR_OAT_MAINT", "BR_OAT_MAINT"])
      
      arr_P_ss[,,i] <- m_p_this_cycle_int
    } else {arr_P_ss[,,i] <- m_p_this_cycle_int}
  }
  
  return(arr_P_ss)
}

gen_arr_pg <- function(a_P){
  arr_P_pres_guid <- a_P
  
  for (i in 29:180){
    m_p_this_cycle_int <- arr_P_pres_guid[,,i]
    
    m_p_this_cycle_int["BN_PN", "BPO_ACUTE"] <- m_p_this_cycle_int["BN_PN", "BPO_ACUTE"]*0.57
    m_p_this_cycle_int["BN_PN", "BN_PN"] <- 1-(rowSums(m_p_this_cycle_int)["BN_PN"] - m_p_this_cycle_int["BN_PN", "BN_PN"])
    
    m_p_this_cycle_int["BN_PN", "BPO_ACUTE"] <- m_p_this_cycle_int["BN_PN", "BPO_ACUTE"]*0.57
    m_p_this_cycle_int["BN_PN", "BN_PN"] <- 1-(rowSums(m_p_this_cycle_int)["BN_PN"] - m_p_this_cycle_int["BN_PN", "BN_PN"])
    
    
    m_p_this_cycle_int["BN_CANCER", "BPO_CANCER"] <- m_p_this_cycle_int["BN_CANCER", "BPO_CANCER"]*0.452
    m_p_this_cycle_int["BN_CANCER", "BN_CANCER"] <- 1-(rowSums(m_p_this_cycle_int)["BN_CANCER"] - m_p_this_cycle_int["BN_CANCER", "BN_CANCER"])
    
    m_p_this_cycle_int["BPO_ACUTE", "BPO_CHRONIC"] <- m_p_this_cycle_int["BPO_ACUTE", "BPO_CHRONIC"] * 0.78
    m_p_this_cycle_int["BPO_ACUTE", "BN_PN"] <- 1-(rowSums(m_p_this_cycle_int)["BPO_ACUTE"] - m_p_this_cycle_int["BPO_ACUTE", "BN_PN"])
    
    m_p_this_cycle_int["BPO_CHRONIC", "BN_CHRONIC"] <- m_p_this_cycle_int["BPO_CHRONIC", "BN_CHRONIC"] *1.04
    m_p_this_cycle_int["BPO_CHRONIC", "BPO_CHRONIC"] <- 1-(rowSums(m_p_this_cycle_int)["BPO_CHRONIC"] - m_p_this_cycle_int["BPO_CHRONIC", "BPO_CHRONIC"])
    
    arr_P_pres_guid[,,i] <- m_p_this_cycle_int
  }
  return(arr_P_pres_guid)
}

gen_arr_all <- function(a_P){
  
  arr_P_all_int <- a_P
  
  for (i in nalx_t1:180) {
    
    m_p_this_cycle_int <- arr_P_all_int[,,i]
    
    m_p_this_cycle_int["BO_OD_ILLICIT", "BO_OD_DEATH"] <- m_p_this_cycle_int["BO_OD_ILLICIT", "BO_OD_DEATH"]*0.95
    m_p_this_cycle_int["BO_OD_ILLICIT", "BI_ILLICIT"] <- 1-(rowSums(m_p_this_cycle_int)["BO_OD_ILLICIT"] - m_p_this_cycle_int["BO_OD_ILLICIT", "BI_ILLICIT"])
    
    m_p_this_cycle_int["BO_OD_RX", "BO_OD_DEATH"] <- m_p_this_cycle_int["BO_OD_RX", "BO_OD_DEATH"]*0.97
    m_p_this_cycle_int["BO_OD_RX", "BPO_MISUSE"] <- 1-(rowSums(m_p_this_cycle_int)["BO_OD_RX"] - m_p_this_cycle_int["BO_OD_RX", "BPO_MISUSE"])
    
    if (i >= 29){
      
      m_p_this_cycle_int["BN_PN", "BPO_ACUTE"] <- m_p_this_cycle_int["BN_PN", "BPO_ACUTE"]*0.57
      m_p_this_cycle_int["BN_PN", "BN_PN"] <- 1-(rowSums(m_p_this_cycle_int)["BN_PN"] - m_p_this_cycle_int["BN_PN", "BN_PN"])
      
      m_p_this_cycle_int["BN_PN", "BPO_ACUTE"] <- m_p_this_cycle_int["BN_PN", "BPO_ACUTE"]*0.57
      m_p_this_cycle_int["BN_PN", "BN_PN"] <- 1-(rowSums(m_p_this_cycle_int)["BN_PN"] - m_p_this_cycle_int["BN_PN", "BN_PN"])
      
      m_p_this_cycle_int["BN_CANCER", "BPO_CANCER"] <- m_p_this_cycle_int["BN_CANCER", "BPO_CANCER"]*0.452
      m_p_this_cycle_int["BN_CANCER", "BN_CANCER"] <- 1-(rowSums(m_p_this_cycle_int)["BN_CANCER"] - m_p_this_cycle_int["BN_CANCER", "BN_CANCER"])
      
      m_p_this_cycle_int["BPO_ACUTE", "BPO_CHRONIC"] <- m_p_this_cycle_int["BPO_ACUTE", "BPO_CHRONIC"] * 0.78
      m_p_this_cycle_int["BPO_ACUTE", "BN_PN"] <- 1-(rowSums(m_p_this_cycle_int)["BPO_ACUTE"] - m_p_this_cycle_int["BPO_ACUTE", "BN_PN"])
      
      m_p_this_cycle_int["BPO_CHRONIC", "BN_CHRONIC"] <- m_p_this_cycle_int["BPO_CHRONIC", "BN_CHRONIC"] *1.04
      m_p_this_cycle_int["BPO_CHRONIC", "BPO_CHRONIC"] <- 1-(rowSums(m_p_this_cycle_int)["BPO_CHRONIC"] - m_p_this_cycle_int["BPO_CHRONIC", "BPO_CHRONIC"])
      
      arr_P_all_int[,,i] <- m_p_this_cycle_int
      
    } else { arr_P_all_int[,,i] <- m_p_this_cycle_int}
    
  }
  return(arr_P_all_int)
}


gen_arr_no <- function(m_P_t) {
  arr_P_t <- array(0, dim=c(n_states, n_states, n_cycles),
                   dimnames = list(v_state_names, v_state_names, 1:n_cycles))
  
  p_od_illicit <- m_P_t["BO_OD_ILLICIT", "BI_ILLICIT"] # 1-rest for transitions from OD
  p_illicit_illicit <- m_P_t["BI_ILLICIT", "BI_ILLICIT"]  # 1-rest for transitions from illicit
  
  for (t in 1:n_cycles) {
    
    m_P_this_cycle <- m_P_t
    
    if(t %in% year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year %in% c(2015, 2016, 2017)]){
      
      arr_P_t[,,t] <- m_P_this_cycle
      
    } else if (t %in% year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year %in% c(2018, 2019, 2020)]) {
      
      m_P_this_cycle["BO_OD_ILLICIT","BO_OD_DEATH"] <- m_P_t["BO_OD_ILLICIT", "BO_OD_DEATH"] * (1.0055)^(t-36)
      m_P_this_cycle["BO_OD_ILLICIT", "BI_ILLICIT"] <- 1-(rowSums(m_P_this_cycle)["BO_OD_ILLICIT"] - p_od_illicit)
      
      m_P_this_cycle["BI_ILLICIT", "BO_OD_ILLICIT"] <- m_P_t["BI_ILLICIT", "BO_OD_ILLICIT"] * (1.0055)^(t-36)
      m_P_this_cycle["BI_ILLICIT", "BI_ILLICIT"] <- 1-(rowSums(m_P_this_cycle)["BI_ILLICIT"] - p_illicit_illicit)
      
    } else {
      
      m_P_this_cycle["BO_OD_ILLICIT","BO_OD_DEATH"] <- m_P_t["BO_OD_ILLICIT", "BO_OD_DEATH"] * (1.0055)^(36)
      m_P_this_cycle["BO_OD_ILLICIT", "BI_ILLICIT"] <- 1-(rowSums(m_P_this_cycle)["BO_OD_ILLICIT"] - p_od_illicit)
      
      m_P_this_cycle["BI_ILLICIT", "BO_OD_ILLICIT"] <- m_P_t["BI_ILLICIT", "BO_OD_ILLICIT"] * (1.0055)^(36)
      m_P_this_cycle["BI_ILLICIT", "BI_ILLICIT"] <- 1-(rowSums(m_P_this_cycle)["BI_ILLICIT"] - p_illicit_illicit)
      
    }
    arr_P_t[,,t] <- m_P_this_cycle
  }
  
  return(arr_P_t)
}



n_psa_iter <- 10000

modnames_abr <- c("no","nalox", "ss", "pg", "all")
pref_abr <- c("cost_diff", "cost_diff_per", "death_diff",
              "death_diff_per", "oddeath_diff", "oddeath_diff_per",
              "tot_death", "tot_oddeath")

vect_names_no <- vect_names_nalox <- vect_names_ss <- vect_names_pg <- vect_names_all <- rep(NA, 14)
for (i in 1:14) { 
  vect_names_no[i] <-  paste0("inci_oddeath_", modnames_abr[1], "_", i)
  vect_names_nalox[i] <-  paste0("inci_oddeath_", modnames_abr[2], "_", i)
  vect_names_ss[i] <-  paste0("inci_oddeath_", modnames_abr[3], "_", i)
  vect_names_pg[i] <-  paste0("inci_oddeath_", modnames_abr[4], "_", i)
  vect_names_all[i] <-  paste0("inci_oddeath_", modnames_abr[5], "_", i)
}

colnames_m_psa <- c(paste0(pref_abr[1], "_", modnames_abr[-1]),
                    paste0(pref_abr[2], "_", modnames_abr[-1]),
                    paste0(pref_abr[3], "_", modnames_abr[-1]),
                    paste0(pref_abr[4], "_", modnames_abr[-1]),
                    paste0(pref_abr[5], "_", modnames_abr[-1]),
                    paste0(pref_abr[6], "_", modnames_abr[-1]),
                    paste0(pref_abr[7], "_", modnames_abr),
                    paste0(pref_abr[8], "_", modnames_abr),
                    vect_names_no,
                    vect_names_nalox,
                    vect_names_ss, vect_names_pg,vect_names_all)

m_outcomes_psa <- matrix(0, nrow = n_psa_iter,
                         ncol = length(colnames_m_psa),
                         dimnames = list(1:n_psa_iter,
                                         colnames_m_psa))


sir_sample <- as.data.frame(sir_sample)

names(sir_sample) <- substr(names(sir_sample), start = 4, stop = nchar(names(sir_sample)) - 3)
names(sir_sample) <- sub("__", " ", names(sir_sample))


# sampling parameters - creating matrix of the parameter
trans_prob_tbl <- trans_prob_tbl %>% 
  select(var_from, var_to, basevalue, lb, ub) %>% 
  mutate(vars_params = paste(var_from, var_to))

v_params_psa <- trans_prob_tbl %>% 
  filter(!is.na(lb)) %>% 
  filter(basevalue < 1 & basevalue > 0) %>% 
  filter(!(vars_params %in% names(sir_sample)))

m_params_psa <- matrix(0, nrow = n_psa_iter,
                       ncol = nrow(v_params_psa) + length(names(sir_sample)),
                       dimnames = list(1:n_psa_iter, c(v_params_psa$vars_params, names(sir_sample)))) 


for (var_param in colnames(m_params_psa)) {
  
  if(var_param %in% names(sir_sample)){
    m_params_psa[, var_param] <- sample(sir_sample[, var_param],
                                        size = n_psa_iter, replace = TRUE)
  } else {
    m_params_psa[, var_param] <- runif(n = n_psa_iter,
                                       min = v_params_psa$lb[v_params_psa$vars_params == var_param],
                                       max = v_params_psa$ub[v_params_psa$vars_params == var_param])
  }
}


# PSA

library(progress)

pb <- progress_bar$new(
  format = "  downloading [:bar] :percent :elapsed eta: :eta",
  total = n_psa_iter, clear = FALSE, width= 60)

for (m in 1:n_psa_iter){
  
  pb$tick()
  Sys.sleep(1/n_psa_iter)
  
  
  temp_tbl <- trans_prob_tbl_new %>% 
    mutate(vars_params = paste(var_from, var_to))
  
  vector_test <- m_params_psa[m,]
  
  for (i in 1:nrow(temp_tbl)){
    temp_tbl$basevalue[i] <- ifelse(temp_tbl$vars_params[i] %in% names(vector_test), 
                                    vector_test[names(vector_test) == temp_tbl$vars_params[i]],
                                    temp_tbl$basevalue[i])
  }
  
  m_P_temp <- matrix(temp_tbl$basevalue, byrow = T,
                     nrow = n_states, ncol = n_states,
                     dimnames = list(v_state_names, v_state_names))
  
  for (i in 1:n_states){
    for (j in 1:n_states){
      m_P_temp[i, j] <- ifelse(m_P_temp[i, j] == 999, 1 - (rowSums(m_P_temp)[i] - 999), m_P_temp[i, j])
    }
  }
  
  arr_P_temp <- gen_arr_no(m_P_t = m_P_temp)
  arr_P_pres_guid_1 <- gen_arr_pg(a_P = arr_P_temp)
  arr_P_ss_1 <- gen_arr_ss(a_P = arr_P_temp)
  arr_P_nalx_1 <- gen_arr_nalox(a_P = arr_P_temp)
  arr_P_all_int_1 <- gen_arr_all(a_P = arr_P_ss_1)
  
  modbasecase <- markov_mod(num_cycles = n_cycles,
                            vec_state_names = v_state_names,
                            vec_m_0 = v_m_0,
                            mat_P = NA,
                            array_P = arr_P_temp,
                            vec_cost_states = v_cost_states,
                            disc_fac =0.03, #annual discount factor
                            t_0_cost =0,
                            inc = pop_increase_tbl_new$increase)
  
  mod_pres_guid_1 <- markov_mod(num_cycles = n_cycles,
                                vec_state_names = v_state_names,
                                vec_m_0 = v_m_0,
                                mat_P = NA,
                                array_P = arr_P_pres_guid_1,
                                vec_cost_states = v_cost_states,
                                disc_fac =0.03, #annual discount factor
                                t_0_cost =0,
                                inc = pop_increase_tbl_new$increase)
  
  mod_ss_1 <- markov_mod_ss(num_cycles = n_cycles,
                            vec_state_names = v_state_names,
                            vec_m_0 = v_m_0,
                            mat_P = NA,
                            array_P = arr_P_ss_1,
                            vec_cost_states = v_cost_states,
                            disc_fac =0.03, #annual discount factor
                            t_0_cost = 0,
                            inc = pop_increase_tbl_new$increase)
  
  mod_nalx_1 <- markov_mod(num_cycles = n_cycles,
                           vec_state_names = v_state_names,
                           vec_m_0 = v_m_0,
                           mat_P = NA,
                           array_P = arr_P_nalx_1,
                           vec_cost_states = v_cost_states,
                           disc_fac =0.03, #annual discount factor
                           t_0_cost = 20000000,
                           inc = pop_increase_tbl_new$increase)
  
  mod_all_int_1 <- markov_mod_ss(num_cycles = n_cycles,
                                 vec_state_names = v_state_names,
                                 vec_m_0 = v_m_0,
                                 mat_P = NA,
                                 array_P = arr_P_all_int_1,
                                 vec_cost_states = v_cost_states,
                                 disc_fac =0.03, #annual discount factor
                                 t_0_cost = 20000000,
                                 inc = pop_increase_tbl_new$increase)
  
  cost_diff_temp <- c(mod_nalx_1$total_net_present_cost - modbasecase$total_net_present_cost,
                      mod_ss_1$total_net_present_cost - modbasecase$total_net_present_cost,
                      mod_pres_guid_1$total_net_present_cost - modbasecase$total_net_present_cost,
                      mod_all_int_1$total_net_present_cost - modbasecase$total_net_present_cost)
  
  cost_diff_per_temp <- c(round((cost_diff_temp/c(modbasecase$total_net_present_cost))*100,2))
  
  death_diff_temp <- c(mod_nalx_1$total_deaths - modbasecase$total_deaths,
                       mod_ss_1$total_deaths - modbasecase$total_deaths,
                       mod_pres_guid_1$total_deaths - modbasecase$total_deaths,
                       mod_all_int_1$total_deaths - modbasecase$total_deaths)
  
  death_diff_per_temp <- c(round((death_diff_temp/c(modbasecase$total_deaths))*100 ,2))
  
  oddeath_diff_temp <- c(mod_nalx_1$total_od_deaths - modbasecase$total_od_deaths,
                         mod_ss_1$total_od_deaths - modbasecase$total_od_deaths,
                         mod_pres_guid_1$total_od_deaths - modbasecase$total_od_deaths,
                         mod_all_int_1$total_od_deaths - modbasecase$total_od_deaths)
  
  oddeath_diff_per_temp <- c(round((oddeath_diff_temp/c(modbasecase$total_od_deaths))*100,2))
  
  total_death_temp <- c(modbasecase$total_deaths,
                        mod_nalx_1$total_deaths,
                        mod_ss_1$total_deaths,
                        mod_pres_guid_1$total_deaths,
                        mod_all_int_1$total_deaths)
  
  total_oddeath_temp <- c(modbasecase$total_od_deaths,
                          mod_nalx_1$total_od_deaths,
                          mod_ss_1$total_od_deaths,
                          mod_pres_guid_1$total_od_deaths,
                          mod_all_int_1$total_od_deaths)
  
  m_outcomes_psa[m, ] <- c(cost_diff_temp,
                           cost_diff_per_temp,
                           death_diff_temp,
                           death_diff_per_temp,
                           oddeath_diff_temp,
                           oddeath_diff_per_temp,
                           total_death_temp,
                           total_oddeath_temp,
                           modbasecase$inci_oddeaths,
                           mod_nalx_1$inci_oddeaths,
                           mod_ss_1$inci_oddeaths,
                           mod_pres_guid_1$inci_oddeaths,
                           mod_all_int_1$inci_oddeaths)
}

saveRDS(m_outcomes_psa, file = here("01_data/m_outcomes_psa.RDS"))
