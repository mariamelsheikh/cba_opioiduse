# Here im only looking at transition probabilities that resulted in more than 2.5% change in number of deaths and number of od deaths
#  and im not calibrating all the other time dependent probabilities just the one at the beginning
library(here)
source(here("02_scripts/01_fun_data.R"))
ows_tbl_prob <- read.csv(here("01_data/ows_tbl_prob.csv"))
calib_target_tbl <- read_excel(here("01_data/markov_model_parameters_preprior.xlsx"), sheet = "calibration_targets")


#Calibration target
target  <- calib_target_tbl$target
names(target) <- paste0(calib_target_tbl$group, "_", calib_target_tbl$year)

# deaths_bc <- mod_basecase$m_M[181, "BO_DEATH"]
# od_deaths_bc <- mod_basecase$m_M[181, "BO_OD_DEATH"] + mod_basecase$extra_od_deaths 

# subset parameters form ows_tbl_prob that resulted in more than 2.5% difference in deaths and od deaths
ows_tbl_prob_subset5 <- ows_tbl_prob %>%
  filter((abs(deaths_bc - deaths_low) / deaths_bc) > 0.025 |
           (abs(deaths_high - deaths_bc) / deaths_bc) > 0.025 |
           (abs(od_deaths_bc - od_deaths_low) / od_deaths_bc) > 0.025 |
           (abs(od_deaths_high - od_deaths_bc) / od_deaths_bc) > 0.025) %>% 
  mutate(var_from = group,
         var_to = substr(name, 3, nchar(name)),
         var_params = paste0("p__", var_from, "__", var_to)) %>% 
  dplyr::select(var_params) %>% 
  distinct()


##### Creating the table containing parameters to be calibrated, and
#####  ensuring that if lb and ub weren't specified from literature to have them +/- 5%

#creating table with the time-dependent probabilities

calib_param_tbl_timdep <- data.frame()
for (i in 1:180){
  calib_param_tbl_timdep <- rbind.data.frame(calib_param_tbl_timdep,
                                             data.frame(cycle = i,
                                                           p__BI_ILLICIT__BO_OD_ILLICIT = arr_P["BI_ILLICIT", "BO_OD_ILLICIT", i],
                                                        p__BO_OD_ILLICIT__BO_OD_DEATH = arr_P["BO_OD_ILLICIT", "BO_OD_DEATH", i]))
}

calib_param_tbl_timdep <- calib_param_tbl_timdep %>% 
  distinct(p__BI_ILLICIT__BO_OD_ILLICIT, 
           p__BO_OD_ILLICIT__BO_OD_DEATH,
           .keep_all = T) %>% 
  pivot_longer(2:3, names_to = "var_params", values_to = "basevalue") %>% 
  mutate(var_params_cycle = paste(var_params, cycle, sep = "__")) %>% 
  select(-cycle) %>% 
  mutate(lb = basevalue * 0.95,
         ub = basevalue * 1.05)


calib_param_tbl_org <- as.data.frame(arr_P[,,1]) %>% 
  pivot_longer(1:31, names_to = "var_to", values_to = "basevalue") %>% 
  mutate(var_from = rep(rownames(m_P), each = n_states)) %>% 
  filter(basevalue != 0) %>% 
  left_join(., trans_prob_tbl %>% 
              select(var_to, var_from,
                     basevalue, lb, ub, range),
            by = c("var_to", "var_from", "basevalue")) %>%
  filter(!is.na(lb)) %>% 
  mutate(lb = ifelse(range %in% c("no"), basevalue*0.95, lb),
         ub = ifelse(range %in% c("no", "informed_l"), basevalue * 1.05, ub)) %>% 
  mutate(var_params = paste0("p__", var_from,"__",var_to),
         var_params_cycle = paste(var_params, 1, sep = "__")) %>%
  select(var_params,var_params_cycle, basevalue, lb, ub)

# Final table with only parameters to be calibrated: 
# those form ows_tbl_prob_subset5, all transition probabilities to death/OD_death
calib_param_tbl <- ows_tbl_prob_subset5 %>% 
  bind_rows(., calib_param_tbl_org %>% 
              select(var_params) %>% 
              filter(grepl("DEATH", var_params))) %>% 
  left_join(., calib_param_tbl_org, by = "var_params") %>% 
  arrange(var_params) %>% 
  filter(!is.na(var_params_cycle)) %>% 
  filter(basevalue != 1) %>%  # removes absorbing states
  distinct()
# %>% 
#   bind_rows(., calib_param_tbl_timdep)

v_params_to_be_calib <- calib_param_tbl$var_params_cycle
n_params <- length(v_params_to_be_calib)

# v_parms_calib <- calib_param_tbl$basevalue[calib_param_tbl$var_params_cycle == v_params_to_be_calib]
# model function that takes in vector of parameters and runs markov model with calibration targets

mod_calib <- function(v_parms_calib) {
  
  v_parms_calib <- suppressWarnings({
    data.frame(v_parms_calib) %>% 
      mutate(params = v_params_to_be_calib) %>% 
      rename(value = v_parms_calib) %>% 
      separate(params, into = c("temp", "var_from",
                                "var_to", "cycle"),
               sep = "__" ) %>%
      dplyr::select(-temp)
  })
  
  v_parms_calib$cycle <- as.numeric(v_parms_calib$cycle)
  temp_trans_prob_tbl <- trans_prob_tbl_new %>% 
      left_join(., v_parms_calib %>% select(-cycle),
              by = c("var_from", "var_to")) %>% 
    mutate(basevalue = ifelse(!is.na(value), value, basevalue)) %>% 
    select(-value)
  
  m_P_temp <- matrix(temp_trans_prob_tbl$basevalue, byrow = T,
                nrow = n_states, ncol = n_states,
                dimnames = list(v_state_names, v_state_names))
  
  
  for (i in 1:n_states){ 
    for (j in 1:n_states){
      m_P_temp[i, j] <- ifelse(m_P_temp[i, j] == 999, 1 - (rowSums(m_P_temp)[i] - 999), m_P_temp[i, j])
    }
  }
  
  
  arr_P_calib <- gen_trans_prob_arr(mat_P = m_P_temp)
  
  
 m_M <- markov_mod(num_cycles = n_cycles,
                      vec_state_names = v_state_names,
                      vec_m_0 = v_m_0, 
                      mat_P = NA, 
                      array_P = arr_P_calib,
                      vec_cost_states = v_cost_states,
                      disc_fac =0.03, #annual discount factor
                      t_0_cost =0,
                      inc = pop_increase_tbl_new$increase)$m_M
  
  #### calculate targets
  
  # number of deaths (excluding OD-deaths)
  yr1_deaths <- 2016
  yrs_deaths <- length(2016:2020)
  num_deaths <- rep(NA, yrs_deaths)
  
  for (i in 1:yrs_deaths){
    yr <- (yr1_deaths - 1) + i
    num_deaths[i] <- m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                    year_mon_cycle_tbl$mon == 12] + 1, "BO_DEATH"] -
      m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
                                     year_mon_cycle_tbl$mon == 12] + 1, "BO_DEATH"]
  }
  
  # number of OD-deaths
  yr1_oddeaths <- 2016
  yrs_oddeaths <- length(2016:2021)
  num_od_deaths <- rep(NA, yrs_oddeaths)
  
  for (i in 1:yrs_oddeaths){
    yr <- (yr1_oddeaths - 1) + i 
    num_od_deaths[i] <- m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                       year_mon_cycle_tbl$mon == 12] + 1, "BO_OD_DEATH"] -
      m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
                                     year_mon_cycle_tbl$mon == 12] + 1, "BO_OD_DEATH"]
    
  }
  
  
  # prevalence of people on prescribed opioids
  yr1_prev_rx <- 2015
  yrlast_prev_rx <- 2018
  yrs_prev_rx <- length(2015:2018)
  prop_opioids_rx <- data.frame(matrix(NA, byrow = T, nrow = 12, ncol = yrs_prev_rx,
                                       list(NULL, paste0("prop_opioids_rx_", str_sub(yr1_prev_rx:yrlast_prev_rx, 3, 4)))))
  
  tot_pop_tbl <- data.frame(matrix(NA, byrow = T, nrow = 12, ncol = yrs_prev_rx,
                                   list(NULL, paste0("tot_pop_", str_sub(yr1_prev_rx:yrlast_prev_rx, 3, 4)))))
  
  for (i in 1:yrs_prev_rx) {
    
    yr <- (yr1_prev_rx - 1) + i
    prop_opioids_rx[, i] <- assign(
      paste0("prop_opioids_rx_", str_sub(yr, 3, 4)),
      rowSums(m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,
                  c("BPO_ACUTE", "BPO_CHRONIC", "BPO_CANCER",
                    "BPO_OTHER", "BPO_PALLIATIVE")]) /
        rowSums(m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
      
    )
    
    tot_pop_tbl[, i] <- assign(
      paste0("tot_pop_", str_sub(yr, 3, 4)),
      rowSums(m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
    )
  }
  
  prop_opioids_rx <- prop_opioids_rx %>% 
    pivot_longer(1:4, names_to = "grp", values_to = "target_val") %>% 
    arrange(grp) %>% 
    mutate(year = rep(yr1_prev_rx:yrlast_prev_rx, each = 12),
           year_mon = paste(year, month.abb[1:12], sep = "_")) %>% 
    bind_cols(., tot_pop_tbl %>% 
                pivot_longer(1:4, names_to = "tot_pop", values_to = "tot_pop_val") %>% 
                arrange(tot_pop) %>% dplyr::select(-tot_pop)) %>% 
    group_by(year) %>% 
    summarise(target_val = weighted.mean(target_val, tot_pop_val)) %>% 
    ungroup()
  
  return(list(num_deaths = num_deaths,
              num_od_deaths = num_od_deaths,
              wmean_prop_rx_opioid = prop_opioids_rx$target_val,
              m_M = m_M)
  )
}
