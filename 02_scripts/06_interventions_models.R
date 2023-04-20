source(here("02_scripts/01_fun_data.R"))


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


# mod_basecase ---------------------------------------------------------

# Make changes from MAP 

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

arr_P_map <- gen_trans_prob_arr(mat_P = m_P_map)

mod_basecase <- markov_mod(num_cycles = n_cycles,
                             vec_state_names = v_state_names,
                             vec_m_0 = v_m_0,
                             mat_P = NA,
                             array_P = arr_P_map,
                             vec_cost_states = v_cost_states,
                             disc_fac = 0.03, #annual discount factor
                             t_0_cost = 0,
                             inc = pop_increase_tbl_new$increase)


arr_P_map_ntd <- array(0, dim=c(n_states, n_states, n_cycles),
                       dimnames = list(v_state_names, v_state_names, 1:n_cycles))

for (i in 1:180){
    
  arr_P_map_ntd[,,i] <- m_P_map
}
    
  
mod_basecase_2 <- markov_mod(num_cycles = n_cycles,
                             vec_state_names = v_state_names,
                             vec_m_0 = v_m_0,
                             mat_P = NA,
                             array_P = arr_P_map_ntd,
                             vec_cost_states = v_cost_states,
                             disc_fac = 0.03, #annual discount factor
                             t_0_cost = 0,
                             inc = pop_increase_tbl_new$increase)

# mod_pres_guid_1 ---------------------------------------------------------

arr_P_pres_guid_1 <- gen_arr_pg(arr_P_map)  

mod_pres_guid_1 <- markov_mod(num_cycles = n_cycles,
                              vec_state_names = v_state_names,
                              vec_m_0 = v_m_0,
                              mat_P = NA,
                              array_P = arr_P_pres_guid_1,
                              vec_cost_states = v_cost_states,
                              disc_fac = 0.03, #annual discount factor
                              t_0_cost = 0,
                              inc = pop_increase_tbl_new$increase)

# mod_pres_guid_2 ---------------------------------------------------------

arr_P_pres_guid_2 <- gen_arr_pg(arr_P_map_ntd)
mod_pres_guid_2 <- markov_mod(num_cycles = n_cycles,
                              vec_state_names = v_state_names,
                              vec_m_0 = v_m_0,
                              mat_P = NA,
                              array_P = arr_P_pres_guid_2,
                              vec_cost_states = v_cost_states,
                              disc_fac = 0.03, #annual discount factor
                              t_0_cost = 0,
                              inc = pop_increase_tbl_new$increase)

# mod_ss_1 ----------------------------------------------------------------
arr_P_ss_1 <- gen_arr_ss(arr_P_map)

mod_ss_1 <- markov_mod_ss(num_cycles = n_cycles,
                          vec_state_names = v_state_names,
                          vec_m_0 = v_m_0,
                          mat_P = NA,
                          array_P = arr_P_ss_1,
                          vec_cost_states = v_cost_states,
                          disc_fac = 0.03, #annual discount factor
                          t_0_cost = 0,
                          inc = pop_increase_tbl_new$increase)

# mod_ss_2 ----------------------------------------------------------------

arr_P_ss_2 <- gen_arr_ss(arr_P_map_ntd)

mod_ss_2 <- markov_mod_ss(num_cycles = n_cycles,
                          vec_state_names = v_state_names,
                          vec_m_0 = v_m_0,
                          mat_P = NA,
                          array_P = arr_P_ss_2,
                          vec_cost_states = v_cost_states,
                          disc_fac = 0.03, #annual discount factor
                          t_0_cost = 0,
                          inc = pop_increase_tbl_new$increase)


# mod_nalx_1 --------------------------------------------------------------

arr_P_nalx_1 <- gen_arr_nalox(arr_P_map)
mod_nalx_1 <- markov_mod(num_cycles = n_cycles,
                         vec_state_names = v_state_names,
                         vec_m_0 = v_m_0,
                         mat_P = NA,
                         array_P = arr_P_nalx_1,
                         vec_cost_states = v_cost_states,
                         disc_fac = 0.03, #annual discount factor
                         t_0_cost = 20000000,
                         inc = pop_increase_tbl_new$increase)

# mod_nalx_2 --------------------------------------------------------------


arr_P_nalx_2 <- gen_arr_nalox(arr_P_map_ntd)

mod_nalx_2 <- markov_mod(num_cycles = n_cycles,
                         vec_state_names = v_state_names,
                         vec_m_0 = v_m_0,
                         mat_P = NA,
                         array_P = arr_P_nalx_2,
                         vec_cost_states = v_cost_states,
                         disc_fac = 0.03, #annual discount factor
                         t_0_cost = 20000000,
                         inc = pop_increase_tbl_new$increase)



# mod_all_int_1 -----------------------------------------------------------

arr_P_all_int_1 <- gen_arr_all(arr_P_map)
mod_all_int_1 <- markov_mod_ss(num_cycles = n_cycles,
                               vec_state_names = v_state_names,
                               vec_m_0 = v_m_0,
                               mat_P = NA,
                               array_P = arr_P_all_int_1,
                               vec_cost_states = v_cost_states,
                               disc_fac = 0.03, #annual discount factor
                               t_0_cost = 20000000,
                               inc = pop_increase_tbl_new$increase)

# mod_all_int_2 -----------------------------------------------------------

arr_P_all_int_2 <- gen_arr_all(arr_P_map_ntd)
mod_all_int_2 <- markov_mod_ss(num_cycles = n_cycles,
                               vec_state_names = v_state_names,
                               vec_m_0 = v_m_0,
                               mat_P = NA,
                               array_P = arr_P_all_int_2,
                               vec_cost_states = v_cost_states,
                               disc_fac = 0.03, #annual discount factor
                               t_0_cost = 20000000,
                               inc = pop_increase_tbl_new$increase)

