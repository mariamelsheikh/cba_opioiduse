# "One-way sensitivity analysis - pre-calibration for status quo model - time dependent"


# Packages + source scripts + data ----------------------------------------
library(here)
library(plotly)
theme_set(theme_minimal())
source(here("02_scripts/ontario/01_fun_data_on.R"))


# Initial population ------------------------------------------------------

##### creating the empty matrix to be populated by the OWS results
m_owsa_ini_pop <- matrix(0, 
                         nrow = nrow(ini_states_tbl_on),
                         ncol = 8,
                         dimnames = list(
                           rows = paste0("ini_pop_",v_state_names),
                           cols = c("input_low", "input_high",
                                    "costs_low", "costs_high",
                                    "deaths_low", "deaths_high",
                                    "od_deaths_low", "od_deaths_high")
                         ))

##### loop for OWS results 
for (state in v_state_names){
  
  t_ini_pop_temp <- ini_states_tbl_on %>% 
    mutate(basevalue = ifelse(`Variable ` == state, lb, basevalue))
  
  mod <- markov_mod(num_cycles = n_cycles,
                    vec_state_names = v_state_names,
                    vec_m_0 = t_ini_pop_temp$basevalue,
                    mat_P = m_P,
                    array_P = NA,
                    vec_cost_states = v_cost_states,
                    disc_fac = 0.03, t_0_cost = 0,
                    inc = pop_increase_tbl_new$increase)
  
  deaths_lb <- mod$m_M[(n_cycles + 1), "BO_DEATH"]
  od_deaths_lb <- mod$m_M[(n_cycles + 1), "BO_OD_DEATH"] + mod$extra_od_deaths
  costs_lb <- mod$total_net_present_cost
  
  t_ini_pop_temp <- ini_states_tbl_on %>% 
    mutate(basevalue = ifelse(`Variable ` == state, ub, basevalue))
  
  mod <- markov_mod(num_cycles = n_cycles,
                    vec_state_names = v_state_names,
                    vec_m_0 = t_ini_pop_temp$basevalue,
                    mat_P = m_P,
                    array_P = NA,
                    vec_cost_states = v_cost_states,
                    disc_fac = 0.03, t_0_cost = 0,
                    inc = pop_increase_tbl_new$increase)
  
  deaths_ub <- mod$m_M[(n_cycles + 1), "BO_DEATH"]
  od_deaths_ub <- mod$m_M[(n_cycles + 1), "BO_OD_DEATH"] + mod$extra_od_deaths
  costs_ub <- mod$total_net_present_cost
  
  m_owsa_ini_pop[paste0("ini_pop_",state), ] <- c(ini_states_tbl_on$lb[ini_states_tbl_on$`Variable ` == state],
                                                  ini_states_tbl_on$ub[ini_states_tbl_on$`Variable ` == state],
                                                  costs_lb, costs_ub, 
                                                  deaths_lb, deaths_ub,
                                                  od_deaths_lb, od_deaths_ub)
}


t_owsa_ini_pop_ntd <- as_tibble(m_owsa_ini_pop) %>% 
  mutate(costs_range = abs(costs_high - costs_low),
         deaths_range = abs(deaths_high - deaths_low),
         od_deaths_range = abs(od_deaths_high - od_deaths_low),
         name = rownames(m_owsa_ini_pop))


# Transition Probabilities ------------------------------------------------

trans_prob_tbl_ows <- cbind.data.frame(rep(v_state_names, each = n_states), 
                                       rep(v_state_names, n_states))

names(trans_prob_tbl_ows) <- c("var_from", "var_to")
names(trans_prob_tbl) <- c("var_from", "var_to",
                           "basevalue", "lb", "ub", "range")

trans_prob_tbl_ows <- trans_prob_tbl_ows %>% 
  full_join(., trans_prob_tbl %>% 
              select(var_from, var_to, basevalue, lb, ub), by = c("var_from", "var_to")) %>% 
  replace(is.na(.), 0)

v_var1 <- unique(trans_prob_tbl_ows$var_from[!(trans_prob_tbl_ows$var_from %in% c("BO_OD_DEATH", "BO_DEATH"))])
v_var2 <- trans_prob_tbl_ows$var_to[trans_prob_tbl_ows$basevalue == 999]


dsa_fun <- function(var1, var2) {
  
  m_temp <- matrix(0,
                   nrow = nrow(trans_prob_tbl_ows[trans_prob_tbl_ows$var_from == var1, ]),
                   ncol = 8 + 4,
                   dimnames = list(
                     rows = paste0("p_", trans_prob_tbl_ows$var_to[trans_prob_tbl_ows$var_from == var1]),
                     cols = c("input_bc", "input_low", "input_high",
                              "costs_low", "costs_high", "cost_bc",
                              "deaths_low", "deaths_high", "deaths_bc",
                              "od_deaths_low", "od_deaths_high", "od_deaths_bc")
                   ))
  
  for(state_to in unique(trans_prob_tbl_ows$var_to)){
    
    t_p_temp_full <- trans_prob_tbl_ows %>% 
      mutate(basevalue = ifelse(var_to == state_to & 
                                  var_from == var1 &
                                  var_to != var2, lb, basevalue))
    
    m_P_temp <- matrix(t_p_temp_full$basevalue, byrow = T,
                       nrow = n_states, ncol = n_states,
                       dimnames = list(v_state_names, v_state_names))
    
    for (i in 1:n_states){
      for (j in 1:n_states){
        m_P_temp[i, j] <- ifelse(m_P_temp[i, j] == 999, 1 - (rowSums(m_P_temp)[i] - 999), m_P_temp[i, j])
      }
    }
    
    mod <- markov_mod(array_P = NA, mat_P = m_P_temp)
    
    deaths_lb <- mod$m_M[(n_cycles + 1), "BO_DEATH"]
    od_deaths_lb <- mod$m_M[(n_cycles + 1), "BO_OD_DEATH"] + mod$extra_od_deaths
    costs_lb <- mod$total_net_present_cost
    
    t_p_temp_full <- trans_prob_tbl_ows %>% 
      mutate(basevalue = ifelse(var_to == state_to & 
                                  var_from == var1 &
                                  var_to != var2, ub, basevalue))
    
    m_P_temp <- matrix(t_p_temp_full$basevalue, byrow = T,
                       nrow = n_states, ncol = n_states,
                       dimnames = list(v_state_names, v_state_names))
    
    
    for (i in 1:n_states){
      for (j in 1:n_states){
        m_P_temp[i, j] <- ifelse(m_P_temp[i, j] == 999, 1 - (rowSums(m_P_temp)[i] - 999), m_P_temp[i, j])
      }
    }
    
    mod <- markov_mod(array_P = NA, mat_P = m_P_temp)
    
    
    deaths_ub <- mod$m_M[(n_cycles + 1), "BO_DEATH"]
    od_deaths_ub <- mod$m_M[(n_cycles + 1), "BO_OD_DEATH"] + mod$extra_od_deaths
    costs_ub <- mod$total_net_present_cost
    
    
    t_p_temp_full <- trans_prob_tbl_ows
    
    m_P_temp <- matrix(t_p_temp_full$basevalue, byrow = T,
                       nrow = n_states, ncol = n_states,
                       dimnames = list(v_state_names, v_state_names))
    
    
    for (i in 1:n_states){
      for (j in 1:n_states){
        m_P_temp[i, j] <- ifelse(m_P_temp[i, j] == 999, 1 - (rowSums(m_P_temp)[i] - 999), m_P_temp[i, j])
      }
    }
    
    
    mod <- markov_mod(array_P = NA, mat_P = m_P_temp)
    
    
    deaths_bc_ows <- mod$m_M[(n_cycles + 1), "BO_DEATH"]
    od_deaths_bc_ows <- mod$m_M[(n_cycles + 1), "BO_OD_DEATH"] + mod$extra_od_deaths
    costs_bc_ows <- mod$total_net_present_cost
    
    
    
    
    m_temp[paste0("p_", state_to), ] <- c(trans_prob_tbl_ows$basevalue[trans_prob_tbl_ows$var_from == var1 &
                                                                         trans_prob_tbl_ows$var_to == state_to],
                                          trans_prob_tbl_ows$lb[trans_prob_tbl_ows$var_from == var1 &
                                                                  trans_prob_tbl_ows$var_to == state_to],
                                          trans_prob_tbl_ows$ub[trans_prob_tbl_ows$var_from == var1 &
                                                                  trans_prob_tbl_ows$var_to == state_to],
                                          costs_lb, costs_ub,costs_bc_ows,
                                          deaths_lb, deaths_ub,deaths_bc_ows,
                                          od_deaths_lb, od_deaths_ub, od_deaths_bc_ows)
  }
  
  return(as_tibble(m_temp) %>% 
           mutate(group = var1,
                  cost_range = abs(costs_high - costs_low),
                  deaths_range = abs(deaths_high - deaths_low),
                  od_deaths_range = abs(od_deaths_high - od_deaths_low),
                  name = rownames(m_temp)))
}

ows_tbl_prob_ntd <- data.frame()

for (i in 1:length(v_var1)){
  ows_tbl_prob_ntd <- rbind.data.frame(ows_tbl_prob_ntd, dsa_fun(v_var1[i], v_var2[i]))
}

write.csv(t_owsa_ini_pop_ntd, file = "01_data/ontario/ows_tbl_ini_pop_ntd_on.csv")
write.csv(ows_tbl_prob_ntd, file = "01_data/ontario/ows_tbl_prob_ntd_on.csv")




