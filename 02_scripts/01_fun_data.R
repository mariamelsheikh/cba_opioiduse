# Packages ----------------------------------------------------------------
library(here)
library(tidyverse)
library(readxl)
library(lubridate)
library(expm) #for matrix operations


# Load Data --------------------------------------------------------------------

ini_states_tbl_org <- read_excel(here('01_data/markov_model_parameters_preprior.xlsx'), sheet = "Initial population parameters")
trans_prob_tbl <- read_excel(here('01_data/markov_model_parameters_preprior.xlsx'), sheet = "transition probabilities ra")
costs_tbl <- read_excel(here('01_data/markov_model_parameters_preprior.xlsx'), sheet = "Cost")
pop_increase_tbl <- read_excel(here('01_data/markov_model_parameters_preprior.xlsx'), sheet = "population per cycle")


# Parameters --------------------------------------------------------------

start_date <- ymd("2015-01-01")
end_date <- ymd("2030-01-01")

cycle_length <- 1/12
n_cycles <- (interval(start_date, end_date) / years(1)) / cycle_length # 180 cycles

v_state_names <- ini_states_tbl_org$`Variable `
n_states <- length(v_state_names)


# Initial population ------------------------------------------------------

v_m_0 <- ini_states_tbl_org$basevalue
names(v_m_0) <- v_state_names


# Costs -------------------------------------------------------------------

# creating a vector for costs
v_cost_states <- costs_tbl$Cost
names(v_cost_states) <- costs_tbl$`Variable `


# population increase -----------------------------------------------------

pop_increase_tbl_new <- cbind.data.frame(year = rep(pop_increase_tbl$Year,
                                                    each = 12),
                                         cycle = 1:n_cycles,
                                         increase = rep(pop_increase_tbl$increase_per_month,
                                                        each = 12))


# create a transition probability matrix -------------------------------------------

trans_prob_tbl_new <- cbind.data.frame(rep(v_state_names, each = n_states), 
                                         rep(v_state_names, n_states))

names(trans_prob_tbl_new) <- c("var_from", "var_to")
names(trans_prob_tbl) <- c("var_from", "var_to",
                 "basevalue", "lb", "ub", "range")

trans_prob_tbl_new <- trans_prob_tbl_new %>% 
    full_join(., trans_prob_tbl %>% 
                select(var_from, var_to, basevalue), by = c("var_from", "var_to")) %>% 
    replace(is.na(.), 0)

m_P <- matrix(trans_prob_tbl_new$basevalue, byrow = T,
                nrow = n_states, ncol = n_states,
                dimnames = list(v_state_names, v_state_names))

for (i in 1:n_states){ 
  
    for (j in 1:n_states){
      m_P[i, j] <- ifelse(m_P[i, j] == 999, 1 - (rowSums(m_P)[i] - 999), m_P[i, j])
    }
  }
  

# rowSums(m_P) # all rows add up to 1
# sum(m_P < 0) # 0 negative numbers


# transition probability array --------------------------------------------

year_mon_cycle_tbl <- cbind.data.frame(
  cycle = 1:n_cycles,
  year = rep(seq(from = 2015, to = 2029), each = 12),
  mon = rep(1:12, 15)
)

# this function generates array with time dependent probabilities for od-illicit -> od_death and illicit -> od_illicit
# it can take on either a matrix or an array (if we want to make changes to an already existing array)

gen_trans_prob_arr <- function(mat_P){
  
  arrP <- array(0, dim=c(n_states, n_states, n_cycles),
                 dimnames = list(v_state_names, v_state_names, 1:n_cycles))
  
  for (t in 1:n_cycles) {
    
    m_P_this_cycle <- mat_P
    
    if(t %in% year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year %in% c(2015, 2016, 2017)]){
      
      arrP[,,t] <- m_P_this_cycle
      
    } else if (t %in% year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year %in% c(2018, 2019, 2020)]) {
      
      m_P_this_cycle["BO_OD_ILLICIT","BO_OD_DEATH"] <- mat_P["BO_OD_ILLICIT", "BO_OD_DEATH"] * (1.005)^(t-36)
      m_P_this_cycle["BO_OD_ILLICIT", "BI_ILLICIT"] <- 1-(rowSums(m_P_this_cycle)["BO_OD_ILLICIT"] - mat_P["BO_OD_ILLICIT", "BI_ILLICIT"])
      
      m_P_this_cycle["BI_ILLICIT", "BO_OD_ILLICIT"] <- mat_P["BI_ILLICIT", "BO_OD_ILLICIT"] * (1.005)^(t-36)
      m_P_this_cycle["BI_ILLICIT", "BI_ILLICIT"] <- 1-(rowSums(m_P_this_cycle)["BI_ILLICIT"] - mat_P["BI_ILLICIT", "BI_ILLICIT"])
      
    } else {
      
      m_P_this_cycle["BO_OD_ILLICIT","BO_OD_DEATH"] <- mat_P["BO_OD_ILLICIT", "BO_OD_DEATH"] * (1.005)^(36)
      m_P_this_cycle["BO_OD_ILLICIT", "BI_ILLICIT"] <- 1-(rowSums(m_P_this_cycle)["BO_OD_ILLICIT"] - mat_P["BO_OD_ILLICIT", "BI_ILLICIT"])
      
      m_P_this_cycle["BI_ILLICIT", "BO_OD_ILLICIT"] <- mat_P["BI_ILLICIT", "BO_OD_ILLICIT"] * (1.005)^(36)
      m_P_this_cycle["BI_ILLICIT", "BI_ILLICIT"] <- 1-(rowSums(m_P_this_cycle)["BI_ILLICIT"] - mat_P["BI_ILLICIT", "BI_ILLICIT"])
      
    }
    
    arrP[,,t] <- m_P_this_cycle
  }
  
  return(arrP)
}

arr_P <- gen_trans_prob_arr(mat_P = m_P)


# Model function ----------------------------------------------------------

markov_mod <- function(num_cycles = n_cycles,
                       vec_state_names = v_state_names,
                       vec_m_0 = v_m_0, 
                       mat_P = m_P, 
                       array_P = NA,
                       vec_cost_states = v_cost_states,
                       disc_fac =0.03, #annual discount factor
                       t_0_cost =0,
                       inc = pop_increase_tbl_new$increase){
  
  
  if(max(is.na(mat_P)) + max(is.na(array_P)) != 1){
    stop("must provide m_P or arr_P but not both")
  }
  
  
  n_states <- length(vec_state_names)
  
  
  if(max(!is.na(mat_P))){
    if(min(rowSums(mat_P)) - 1 > 1e-5 | max(rowSums(mat_P)) - 1 > 1e-5) { 
      stop("rows of transition matrix must all sum to one")
    }
  } else{
    
    for(t in 1:num_cycles){
      row_sums = rowSums(array_P[,,t])
      if(min(row_sums) - 1 > 1e-5 | max(row_sums) - 1 > 1e-5){
        stop("rows of the transition matrices in the transition 
            array do not all sum to 1 for all cycles")
      }
    }
  }
  
  m_M <- matrix(0, nrow = num_cycles + 1, ncol = n_states,
                dimnames = list(0:num_cycles, vec_state_names))
  
  m_cost_outcome_by_cycle_by_state <- matrix(0, nrow = num_cycles, ncol = n_states)
  
  m_M[1, ] <- vec_m_0
  
  for (t in 1:num_cycles){
    
    if(max(!is.na(array_P))){
      mat_P <- array_P[,,t]
    }
    
    m_M[t+1, ] <- m_M[t, ] %*% mat_P
    m_M[t+1, "BN_PN"] <- m_M[t+1, "BN_PN"] + inc[t]

  }
  
  m_extra_cycle <- m_M[181, ] %*% mat_P
  total_od_deaths_181 <- m_extra_cycle[,"BO_OD_DEATH"]
  extra_od_deaths <- total_od_deaths_181 - m_M[181,"BO_OD_DEATH"]
  total_mod_bi_181 <- m_extra_cycle[,"BO_MOD_BI"]
  extra_mod_bi <- total_mod_bi_181 - m_M[181,"BO_MOD_BI"]
  total_sev_bi_181 <- m_extra_cycle[,"BO_SEVERE_BI"]
  extra_sev_bi <- total_mod_bi_181 - m_M[181,"BO_SEVERE_BI"]
  extra_odseq_costs <- (extra_mod_bi * v_cost_states["BO_MOD_BI"]) +
    (extra_sev_bi * v_cost_states["BO_SEVERE_BI"])
  
  v_cost <- m_M[2:181, ] %*% vec_cost_states
  
  v_discountweight <- 1/((1 + disc_fac*cycle_length))^(0:num_cycles)
  
  
  total_net_present_cost <- t(v_cost) %*% (v_discountweight[-1]) + extra_odseq_costs + t_0_cost
  
  total_deaths <- m_M[181, "BO_DEATH"]
  total_od_deaths <- m_M[181, "BO_OD_DEATH"] + extra_od_deaths
  
  inci_oddeaths <- rep(NA, 14)
  
  for (i in 1:14){
    yr <- 2015 + i 
    inci_oddeaths[i] <- m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                       year_mon_cycle_tbl$mon == 12] + 1, "BO_OD_DEATH"] -
      m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
                                     year_mon_cycle_tbl$mon == 12] + 1,"BO_OD_DEATH"]}
  
  
  return(list(m_M = m_M,
              total_deaths = total_deaths,
              total_od_deaths = total_od_deaths,
              inci_oddeaths = inci_oddeaths,
              total_net_present_cost = total_net_present_cost,
              extra_od_deaths = extra_od_deaths,
              extra_mod_bi = extra_mod_bi,
              extra_sev_bi = extra_sev_bi
  ))
}

# Running basecase model - time dependent - pre-calibration ---------------------------------

mod_basecase <- markov_mod(num_cycles = n_cycles,
                           vec_state_names = v_state_names,
                           vec_m_0 = v_m_0, 
                           mat_P = NA, 
                           array_P = arr_P,
                           vec_cost_states = v_cost_states,
                           disc_fac =0.03, #annual discount factor
                           t_0_cost =0,
                           inc = pop_increase_tbl_new$increase)

# Running basecase model - no time dependent - pre-calibration ---------------------------------

mod_basecase_notimedep <- markov_mod(num_cycles = n_cycles,
                                     vec_state_names = v_state_names,
                                     vec_m_0 = v_m_0,
                                     mat_P = m_P,
                                     array_P = NA,
                                     vec_cost_states = v_cost_states,
                                     disc_fac =0.03, #annual discount factor
                                     t_0_cost =0,
                                     inc = pop_increase_tbl_new$increase)

# function for markov model for safer supply ---------------------------------

markov_mod_ss <- function(num_cycles = n_cycles,
                          vec_state_names = v_state_names,
                          vec_m_0 = v_m_0,
                          mat_P = NA, 
                          array_P = arr_P,
                          vec_cost_states = v_cost_states,
                          disc_fac =0.03, #annual discount factor
                          t_0_cost =0,
                          inc = pop_increase_tbl_new$increase) {
  
  
  if(max(is.na(mat_P)) + max(is.na(array_P)) != 1){
    stop("must provide m_P or arr_P but not both")
  }
  
  
  n_states <- length(vec_state_names)
  
  
  if(max(!is.na(mat_P))){
    if(min(rowSums(mat_P)) - 1 > 1e-5 | max(rowSums(mat_P)) - 1 > 1e-5) { 
      stop("rows of transition matrix must all sum to one")
    }
  } else{
    
    for(t in 1:num_cycles){
      row_sums = rowSums(array_P[,,t])
      if(min(row_sums) - 1 > 1e-5 | max(row_sums) - 1 > 1e-5){
        stop("rows of the transition matrices in the transition 
            array do not all sum to 1 for all cycles")
      }
    }
  }
  
  m_M <- matrix(0, nrow = num_cycles + 1, ncol = n_states,
                dimnames = list(0:num_cycles, vec_state_names))
  
  m_cost_outcome_by_cycle_by_state <- matrix(0, nrow = num_cycles, ncol = n_states)
  
  m_M[1, ] <- vec_m_0
  
  for (t in 1:n_cycles){
    
    mat_P <- array_P[,,t]
    
    if(t == ss_t1){
      m_M[t+1, ] <- m_M[t, ] %*% mat_P
      m_M[t+1, "BN_PN"] <- m_M[t+1, "BN_PN"] + pop_increase_tbl_new$increase[t]
      m_M[t+1, "BS_OAT_MAINT"] <- m_M[t+1, "BS_OAT_MAINT"] - 1000
      m_M[t+1, "BS_OAT_SS"] <- m_M[t+1, "BS_OAT_SS"] + 1000
      
    } else if(t == ss_t2){
      m_M[t+1, ] <- m_M[t, ] %*% mat_P
      m_M[t+1, "BN_PN"] <- m_M[t+1, "BN_PN"] + pop_increase_tbl_new$increase[t]
      m_M[t+1, "BS_OAT_MAINT"] <- m_M[t+1, "BS_OAT_MAINT"] - 10000
      m_M[t+1, "BS_OAT_SS"] <- m_M[t+1, "BS_OAT_SS"] + 10000
      
    } else{
      
      m_M[t+1, ] <- m_M[t, ] %*% mat_P
      m_M[t+1, "BN_PN"] <- m_M[t+1, "BN_PN"] + pop_increase_tbl_new$increase[t]
    }
  }
  
  m_extra_cycle <- m_M[181, ] %*% mat_P
  total_od_deaths_181 <- m_extra_cycle[,"BO_OD_DEATH"]
  extra_od_deaths <- total_od_deaths_181 - m_M[181,"BO_OD_DEATH"]
  total_mod_bi_181 <- m_extra_cycle[,"BO_MOD_BI"]
  extra_mod_bi <- total_mod_bi_181 - m_M[181,"BO_MOD_BI"]
  total_sev_bi_181 <- m_extra_cycle[,"BO_SEVERE_BI"]
  extra_sev_bi <- total_mod_bi_181 - m_M[181,"BO_SEVERE_BI"]
  extra_odseq_costs <- (extra_mod_bi * v_cost_states["BO_MOD_BI"]) +
    (extra_sev_bi * v_cost_states["BO_SEVERE_BI"])
  
  v_cost <- m_M[2:181, ] %*% vec_cost_states
  
  v_discountweight <- 1/((1 + disc_fac*cycle_length))^(0:num_cycles)
  
  total_net_present_cost <- t(v_cost) %*% (v_discountweight[-1]) + extra_odseq_costs + t_0_cost
  
  total_deaths <- m_M[181, "BO_DEATH"]
  total_od_deaths <- m_M[181, "BO_OD_DEATH"] + extra_od_deaths
  
  inci_oddeaths <- rep(NA, 14)
  
  for (i in 1:14){
    yr <- 2015 + i 
    inci_oddeaths[i] <- m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                       year_mon_cycle_tbl$mon == 12] + 1, "BO_OD_DEATH"] -
      m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
                                     year_mon_cycle_tbl$mon == 12] + 1,"BO_OD_DEATH"]}
  
  
  return(list(m_M = m_M,
              total_deaths = total_deaths,
              total_od_deaths = total_od_deaths,
              inci_oddeaths = inci_oddeaths,
              total_net_present_cost = total_net_present_cost,
              extra_od_deaths = extra_od_deaths,
              extra_mod_bi = extra_mod_bi,
              extra_sev_bi = extra_sev_bi
  ))
}


