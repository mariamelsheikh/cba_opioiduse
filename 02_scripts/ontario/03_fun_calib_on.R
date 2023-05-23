# Here im only looking at transition probabilities that resulted in more than 2.5% change in number of deaths and number of od deaths
#  and im not calibrating all the other time dependent probabilities just the one at the beginning
library(here)
source(here("02_scripts/ontario/01_fun_data_on.R"))
ows_tbl_prob <- read.csv(here("01_data/ontario/ows_tbl_prob_on.csv"))
calib_target_mon_tbl <- read_excel(here("01_data/ontario/on_markov_model_parameters_preprior.xlsx"), sheet = "mon_calib_targets")
calib_target_yrly_tbl <- read_excel(here("01_data/ontario/on_markov_model_parameters_preprior.xlsx"), sheet = "yrly_calib_targets")


#Calibration target
calib_target_tbl <- calib_target_mon_tbl %>% 
  pivot_longer(2:4, names_to = "group", values_to = "value") %>% 
  filter(!is.na(value)) %>% 
  separate(1, c("year", "mon", "day"), sep = "-") %>% 
  select(-day) %>% 
  mutate(mon = as.numeric(mon),
         year = as.numeric(year)) %>% 
  bind_rows(., calib_target_yrly_tbl %>% 
              pivot_longer(2:5, names_to = "group",
                           values_to = "value") %>% 
              mutate(year = as.numeric(Year)) %>% 
              mutate(mon = NA) %>% select(-Year)) 


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
for (i in 1:n_cycles){
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
         ub = ifelse(range %in% c("no", "informed_l"),
                     basevalue * 1.05, ub)) %>% 
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

calc_targets_on <- function(m_M_fun, m_P_fun){
  #### calculate targets
  
  # number of deaths (excluding OD-deaths)
  yr1_deaths <- 2010
  yrs_deaths <- length(2010:2021)
  num_deaths <- rep(NA, yrs_deaths)
  
  num_deaths[1] <- m_M_fun[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == 2010 &
                                                      year_mon_cycle_tbl$mon == 12] + 1, "BO_DEATH"]
  
  for (i in 2:yrs_deaths){
    yr <- (yr1_deaths - 1) + i
    num_deaths[i] <- m_M_fun[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                    year_mon_cycle_tbl$mon == 12] + 1, "BO_DEATH"] -
      m_M_fun[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == (yr - 1) &
                                     year_mon_cycle_tbl$mon == 12] + 1, "BO_DEATH"]
  }
  
  tot_deaths <- cbind.data.frame(year = calib_target_yrly_tbl$Year[!is.na(calib_target_yrly_tbl$total_overall_deaths)],
                                 tot_deaths = num_deaths)
  
  # number of OD-deaths
  yrs_mons_oddeath <- nrow(calib_target_mon_tbl %>%
                             select(`Month, Year`, total_od_deaths) %>% 
                             filter(!is.na(total_od_deaths)))
  
  inc_oddeath <- calib_target_mon_tbl %>%
    select(`Month, Year`, total_od_deaths) %>% 
    filter(!is.na(total_od_deaths)) %>% 
    mutate(new_od_deaths = NA,
           tot_pop_val = NA) %>% 
    separate(`Month, Year`, into = c("year", "mon", "day"), sep = "-") %>%
    mutate(year = as.numeric(year),
           mon = as.numeric(mon)) %>% 
    select(-c(day, total_od_deaths)) %>% 
    left_join(., year_mon_cycle_tbl, by = c("year", "mon"))
  
  for (i in 1:nrow(inc_oddeath)) {
    
    this_cycle <- inc_oddeath[i, "cycle"] + 1
    
    
    inc_oddeath[i, "new_od_deaths"] <- m_M_fun[as.numeric(this_cycle) + 1, "BO_OD_DEATH"] -
      m_M_fun[as.numeric(this_cycle), "BO_OD_DEATH"]
    
    inc_oddeath[i, "tot_pop_val", i] <- sum(m_M_fun[as.numeric(this_cycle), ])
    
  }
  
  inc_oddeath <- inc_oddeath %>% 
    mutate(mon_inc = new_od_deaths / tot_pop_val)
  
  yrly_inc_oddeath <- inc_oddeath %>% 
    filter(!(year %in% c(2022))) %>% 
    group_by(year) %>% 
    summarise(weighted_inc = weighted.mean(mon_inc, tot_pop_val)) %>% 
    ungroup()
  
  
  # prevalence of people on prescribed opioids
  yrs_mons_prev_rx <- nrow(calib_target_mon_tbl %>%
                             select(`Month, Year`, total_rx_opioids) %>% 
                             filter(!is.na(total_rx_opioids)))
  
  
  prev_opioids_rx <- calib_target_mon_tbl %>%
    select(`Month, Year`, total_rx_opioids) %>% 
    filter(!is.na(total_rx_opioids)) %>% 
    mutate(total_rx_opioids = NA,
           tot_pop_val = NA) %>% 
    separate(`Month, Year`, into = c("year", "mon", "day"), sep = "-") %>%
    mutate(year = as.numeric(year),
           mon = as.numeric(mon)) %>% 
    select(-day) %>% 
    left_join(., year_mon_cycle_tbl, by = c("year", "mon"))
  
  
  for (i in 1:nrow(prev_opioids_rx)) {
    
    this_cycle <- prev_opioids_rx[i, "cycle"] + 1
    
    
    prev_opioids_rx[i, "total_rx_opioids"] <- sum(m_M_fun[as.numeric(this_cycle),
                                                      c("BPO_ACUTE", "BPO_CHRONIC",
                                                        "BPO_CANCER","BPO_PALLIATIVE")])
    
    
    prev_opioids_rx[i, "tot_pop_val", i] <- sum(m_M_fun[as.numeric(this_cycle), ])
    
  }
  
  prev_opioids_rx <- prev_opioids_rx %>% 
    mutate(mon_prev = total_rx_opioids / tot_pop_val)
  
  yrly_prev_opioids_rx <- prev_opioids_rx %>% 
    filter(!(year %in% c(2012, 2022))) %>% 
    group_by(year) %>% 
    summarise(weighted_prev = weighted.mean(mon_prev, tot_pop_val)) %>% 
    ungroup()
  
  
  # prevalence of OAT numbers
  yrs_mons_oat <- nrow(calib_target_mon_tbl %>%
                         select(`Month, Year`, total_oat) %>% 
                         filter(!is.na(total_oat)))
  
  
  prev_oat <- calib_target_mon_tbl %>%
    select(`Month, Year`, total_oat) %>% 
    filter(!is.na(total_oat)) %>% 
    mutate(total_oat = NA,
           tot_pop_val = NA) %>% 
    separate(`Month, Year`, into = c("year", "mon", "day"), sep = "-") %>%
    mutate(year = as.numeric(year),
           mon = as.numeric(mon)) %>% 
    select(-day) %>% 
    left_join(., year_mon_cycle_tbl, by = c("year", "mon"))
  
  
  for (i in 1:nrow(prev_oat)) {
    
    this_cycle <- prev_oat[i, "cycle"] + 1
    
    
    prev_oat[i, "total_oat"] <- sum(m_M_fun[as.numeric(this_cycle),
                                        c("BS_OAT_INI", "BS_OAT_MAINT",
                                          "BS_OAT_SS","BS_SS", "BR_OAT_INI",
                                          "BR_OAT_MAINT", "BR_OAT_SS", "BR_SS")])
    
    
    prev_oat[i, "tot_pop_val", i] <- sum(m_M_fun[as.numeric(this_cycle), ])
    
  }
  
  prev_oat <- prev_oat %>% 
    mutate(mon_prev = total_oat / tot_pop_val)
  
  yrly_prev_oat <- prev_oat %>% 
    filter(!(year %in% c(2013, 2022))) %>% 
    group_by(year) %>% 
    summarise(weighted_prev = weighted.mean(mon_prev, tot_pop_val)) %>% 
    ungroup()
  
  
  # incidence of cancer - yrly
  yr1_inc_cancer <- 2010
  yrlast_inc_cancer <- 2022
  yrs_inc_cancer <- length(2010:2022)
  inc_cancer <- data.frame(matrix(NA, byrow = T, nrow = 12,
                                  ncol = yrs_inc_cancer,
                                  list(NULL, paste0("new_cancer_",
                                                    str_sub(yr1_inc_cancer:yrlast_inc_cancer, 3, 4)))))
  
  tot_pop_tbl <- data.frame(matrix(NA, byrow = T, nrow = 12, ncol = yrs_inc_cancer,
                                   list(NULL, paste0("tot_pop_",
                                                     str_sub(yr1_inc_cancer:yrlast_inc_cancer, 3, 4)))))
  
  for (i in 1:yrs_inc_cancer) {
    
    yr <- (yr1_inc_cancer - 1) + i
    inc_cancer[, i] <- assign(
      paste0("new_cancer_", str_sub(yr, 3, 4)),
      ((m_M_fun[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,
            "BN_PN"] * m_P_fun["BN_PN", "BN_CANCER"]) +
         (m_M_fun[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,
              "BN_CHRONIC"] * m_P_fun["BN_CHRONIC", "BN_CANCER"]) +
         (m_M_fun[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,
                  "BI_ILLICIT"] * m_P_fun["BI_ILLICIT", "BN_CANCER"]))
    )
    
    tot_pop_tbl[, i] <- assign(
      paste0("tot_pop_", str_sub(yr, 3, 4)),
      rowSums(m_M_fun[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
    )
  }
  
  inc_cancer <- inc_cancer %>%
    pivot_longer(1:yrs_inc_cancer, names_to = "grp",
                 values_to = "new_cases") %>%
    arrange(grp) %>%
    mutate(year = rep(yr1_inc_cancer:yrlast_inc_cancer, each = 12),
           year_mon = paste(year, month.abb[1:12], sep = "_")) %>%
    bind_cols(., tot_pop_tbl %>%
                pivot_longer(1:yrs_inc_cancer, names_to = "tot_pop",
                             values_to = "tot_pop_val") %>%
                arrange(tot_pop) %>% dplyr::select(-tot_pop)) %>%
    mutate(mon_incidence = new_cases / tot_pop_val,
           temp_mon_num = rep(c(1:12), yrs_inc_cancer)) %>%
    left_join(., year_mon_cycle_tbl, by = c("year",
                                            "temp_mon_num" = "mon"))
  
  inc_cancer_yrly <- inc_cancer %>% 
    mutate(jul_pop_val = ifelse(grepl("Jul", year_mon), tot_pop_val, NA)) %>% 
    group_by(year) %>% 
    summarise(new_cases = sum(new_cases),
              weighted_inc = weighted.mean(mon_incidence, tot_pop_val),
              jul_tot_pop_val =  sum(jul_pop_val, na.rm = T),
              jul_yrly_inc = new_cases/jul_tot_pop_val) %>% 
    ungroup() ##### not sure how to calculate the annual total population value, 
  # im not going to use incidence i will use total number of new cases
  
  # cancer prevalence 
  yr1_prev_cancer <- 2015
  yrlast_prev_cancer <- 2022
  yrs_prev_cancer <- length(2015:2022)
  prev_cancer <- data.frame(matrix(NA, byrow = T, nrow = 12, ncol = yrs_prev_cancer,
                                   list(NULL, paste0("tot_cancer_cases_",
                                                     str_sub(yr1_prev_cancer:yrlast_prev_cancer, 3, 4))))) %>% 
    pivot_longer(1:yrs_prev_cancer, names_to = "grp", values_to = "tot_cancer_cases") %>% 
    arrange(grp) %>% 
    mutate(year = rep(yr1_prev_cancer:yrlast_prev_cancer, each = 12),
           year_mon = paste(year, month.abb[1:12], sep = "_"),
           mon = rep(1:12, yrs_prev_cancer)) %>% 
    bind_cols(., data.frame(matrix(NA, byrow = T, nrow = 12, ncol = yrs_prev_cancer,
                                   list(NULL, paste0("tot_pop_",
                                                     str_sub(yr1_prev_cancer:yrlast_prev_cancer, 3, 4))))) %>% 
                pivot_longer(1:yrs_prev_cancer, names_to = "tot_pop", values_to = "tot_pop_val") %>% 
                arrange(tot_pop) %>% dplyr::select(-tot_pop)) %>% 
    left_join(., year_mon_cycle_tbl, by = c("year", "mon"))
  
  
  y0 <- sum(m_M_fun[prev_cancer$cycle[1] + 1, c("BN_CANCER", "BPO_CANCER")]) + 
    (0.7*(sum(m_M_fun[prev_cancer$cycle[1] + 1, c("BN_PALLIATIVE", "BPO_PALLIATIVE")])))
  
  
  prev_cancer[1, "tot_cancer_cases"] <- inc_cancer$new_cases[inc_cancer$cycle == prev_cancer$cycle[1]] + 
    y0 - (
      m_M_fun[prev_cancer$cycle[1] + 1, "BN_PN"]*m_P_fun["BN_CANCER", "BO_DEATH"] + 
        m_M_fun[prev_cancer$cycle[1] + 1, "BN_CANCER"]*m_P_fun["BN_CANCER", "BN_PN"] +
        m_M_fun[prev_cancer$cycle[1] + 1, "BPO_CANCER"]*m_P_fun["BPO_CANCER", "BO_DEATH"]
    )
  
  prev_cancer[1, "tot_pop_val"] <- sum(m_M_fun[1,])
  
  for (i in 2:nrow(prev_cancer)) {
    
    this_cycle <- prev_cancer$cycle[i]
    
    prev_cancer[i, "tot_cancer_cases"] <- inc_cancer$new_cases[inc_cancer$cycle == this_cycle] + 
      prev_cancer[this_cycle, "tot_cancer_cases"] - (
        m_M_fun[this_cycle + 1, "BN_PN"]*m_P_fun["BN_CANCER", "BO_DEATH"] + 
          m_M_fun[this_cycle + 1, "BN_CANCER"]*m_P_fun["BN_CANCER", "BN_PN"] +
          m_M_fun[this_cycle + 1, "BPO_CANCER"]*m_P_fun["BPO_CANCER", "BO_DEATH"]
      )
    
    prev_cancer[i, "tot_pop_val"] <- sum(m_M_fun[this_cycle + 1,])
  }
  
  prev_cancer <- prev_cancer %>% 
    mutate(mon_prev = tot_cancer_cases / tot_pop_val)
  
  prev_cancer_yrly <- prev_cancer %>% 
    group_by(year) %>% 
    summarise(weighted_prev = weighted.mean(mon_prev, tot_pop_val)) %>% 
    ungroup()
  
  return(list(tot_deaths = tot_deaths,
              inc_oddeath = inc_oddeath,
              prev_opioids_rx = prev_opioids_rx,
              prev_oat = prev_oat,
              inc_cancer = inc_cancer_yrly,
              prev_cancer = prev_cancer)
  )
}



mod_calib_on <- function(v_parms_calib) {
  
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
  
  calculated_targets <- calc_targets_on(m_M_fun = m_M,
                                        m_P_fun = m_P_temp)
  
  
  # # number of deaths (excluding OD-deaths)
  # yr1_deaths <- 2010
  # yrs_deaths <- length(2010:2021)
  # num_deaths <- rep(NA, yrs_deaths)
  # 
  # for (i in 1:yrs_deaths){
  #   yr <- (yr1_deaths - 1) + i
  #   num_deaths[i] <- m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
  #                                                   year_mon_cycle_tbl$mon == 12] + 1, "BO_DEATH"] -
  #     m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
  #                                    year_mon_cycle_tbl$mon == 12] + 1, "BO_DEATH"]
  # }
  # 
  # # number of OD-deaths
  # yrs_mons_oddeath <- nrow(calib_target_mon_tbl %>%
  #                            select(`Month, Year`, total_od_deaths) %>% 
  #                            filter(!is.na(total_od_deaths)))
  # 
  # inc_oddeath <- calib_target_mon_tbl %>%
  #   select(`Month, Year`, total_od_deaths) %>% 
  #   filter(!is.na(total_od_deaths)) %>% 
  #   mutate(new_od_deaths = NA,
  #          tot_pop_val = NA) %>% 
  #   separate(`Month, Year`, into = c("year", "mon", "day"), sep = "-") %>%
  #   mutate(year = as.numeric(year),
  #          mon = as.numeric(mon)) %>% 
  #   select(-c(day, total_od_deaths)) %>% 
  #   left_join(., year_mon_cycle_tbl, by = c("year", "mon"))
  # 
  # for (i in 1:nrow(inc_oddeath)) {
  #   
  #   this_cycle <- inc_oddeath[i, "cycle"] + 1
  #   
  #   
  #   inc_oddeath[i, "new_od_deaths"] <- m_M[as.numeric(this_cycle) + 1, "BO_OD_DEATH"] -
  #     m_M[as.numeric(this_cycle), "BO_OD_DEATH"]
  # 
  #   inc_oddeath[i, "tot_pop_val", i] <- sum(m_M[as.numeric(this_cycle), ])
  #   
  # }
  # 
  # inc_oddeath <- inc_oddeath %>% 
  #   mutate(mon_inc = new_od_deaths / tot_pop_val)
  # 
  # yrly_inc_oddeath <- inc_oddeath %>% 
  #   filter(!(year %in% c(2022))) %>% 
  #   group_by(year) %>% 
  #   summarise(weighted_inc = weighted.mean(mon_inc, tot_pop_val)) %>% 
  #   ungroup()
  # 
  # 
  # # prevalence of people on prescribed opioids
  # yrs_mons_prev_rx <- nrow(calib_target_mon_tbl %>%
  #                       select(`Month, Year`, total_rx_opioids) %>% 
  #                       filter(!is.na(total_rx_opioids)))
  # 
  # 
  # prev_opioids_rx <- calib_target_mon_tbl %>%
  #   select(`Month, Year`, total_rx_opioids) %>% 
  #   filter(!is.na(total_rx_opioids)) %>% 
  #   mutate(total_rx_opioids = NA,
  #          tot_pop_val = NA) %>% 
  #   separate(`Month, Year`, into = c("year", "mon", "day"), sep = "-") %>%
  #   mutate(year = as.numeric(year),
  #          mon = as.numeric(mon)) %>% 
  #   select(-day) %>% 
  #   left_join(., year_mon_cycle_tbl, by = c("year", "mon"))
  # 
  # 
  # for (i in 1:nrow(prev_opioids_rx)) {
  #   
  #   this_cycle <- prev_opioids_rx[i, "cycle"] + 1
  #   
  # 
  #   prev_opioids_rx[i, "total_rx_opioids"] <- sum(m_M[as.numeric(this_cycle),
  #                                                         c("BPO_ACUTE", "BPO_CHRONIC",
  #                                                         "BPO_CANCER","BPO_PALLIATIVE")])
  #   
  #   
  #   prev_opioids_rx[i, "tot_pop_val", i] <- sum(m_M[as.numeric(this_cycle), ])
  #     
  # }
  # 
  # prev_opioids_rx <- prev_opioids_rx %>% 
  #   mutate(mon_prev = total_rx_opioids / tot_pop_val)
  # 
  # yrly_prev_opioids_rx <- prev_opioids_rx %>% 
  #   filter(!(year %in% c(2012, 2022))) %>% 
  #   group_by(year) %>% 
  #   summarise(weighted_prev = weighted.mean(mon_prev, tot_pop_val)) %>% 
  #   ungroup()
  # 
  # 
  # # prevalence of OAT numbers
  # yrs_mons_oat <- nrow(calib_target_mon_tbl %>%
  #                            select(`Month, Year`, total_oat) %>% 
  #                            filter(!is.na(total_oat)))
  # 
  # 
  # prev_oat <- calib_target_mon_tbl %>%
  #   select(`Month, Year`, total_oat) %>% 
  #   filter(!is.na(total_oat)) %>% 
  #   mutate(total_oat = NA,
  #          tot_pop_val = NA) %>% 
  #   separate(`Month, Year`, into = c("year", "mon", "day"), sep = "-") %>%
  #   mutate(year = as.numeric(year),
  #          mon = as.numeric(mon)) %>% 
  #   select(-day) %>% 
  #   left_join(., year_mon_cycle_tbl, by = c("year", "mon"))
  # 
  # 
  # for (i in 1:nrow(prev_oat)) {
  #   
  #   this_cycle <- prev_oat[i, "cycle"] + 1
  #   
  #   
  #   prev_oat[i, "total_oat"] <- sum(m_M[as.numeric(this_cycle),
  #                                       c("BS_OAT_INI", "BS_OAT_MAINT",
  #                                         "BS_OAT_SS","BS_SS", "BR_OAT_INI",
  #                                         "BR_OAT_MAINT", "BR_OAT_SS", "BR_SS")])
  #   
  #   
  #   prev_oat[i, "tot_pop_val", i] <- sum(m_M[as.numeric(this_cycle), ])
  #   
  # }
  # 
  # prev_oat <- prev_oat %>% 
  #   mutate(mon_prev = total_oat / tot_pop_val)
  # 
  # yrly_prev_oat <- prev_oat %>% 
  #   filter(!(year %in% c(2013, 2022))) %>% 
  #   group_by(year) %>% 
  #   summarise(weighted_prev = weighted.mean(mon_prev, tot_pop_val)) %>% 
  #   ungroup()
  # 
  # 
  # # incidence of cancer - yrly
  # yr1_inc_cancer <- 2010
  # yrlast_inc_cancer <- 2022
  # yrs_inc_cancer <- length(2010:2022)
  # inc_cancer <- data.frame(matrix(NA, byrow = T, nrow = 12,
  #                                 ncol = yrs_inc_cancer,
  #                               list(NULL, paste0("new_cancer_",
  #                                                 str_sub(yr1_inc_cancer:yrlast_inc_cancer, 3, 4)))))
  # 
  # tot_pop_tbl <- data.frame(matrix(NA, byrow = T, nrow = 12, ncol = yrs_inc_cancer,
  #                                  list(NULL, paste0("tot_pop_",
  #                                                    str_sub(yr1_inc_cancer:yrlast_inc_cancer, 3, 4)))))
  # 
  # for (i in 1:yrs_inc_cancer) {
  # 
  #   yr <- (yr1_inc_cancer - 1) + i
  #   inc_cancer[, i] <- assign(
  #     paste0("new_cancer_", str_sub(yr, 3, 4)),
  #     ((m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,
  #           "BN_PN"] * m_P_temp["BN_PN", "BN_CANCER"]) +
  #       (m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,
  #            "BN_CHRONIC"] * m_P_temp["BN_CHRONIC", "BN_CANCER"]))
  #   )
  # 
  #   tot_pop_tbl[, i] <- assign(
  #     paste0("tot_pop_", str_sub(yr, 3, 4)),
  #     rowSums(m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
  #   )
  # }
  # 
  # inc_cancer <- inc_cancer %>%
  #   pivot_longer(1:yrs_inc_cancer, names_to = "grp",
  #                values_to = "new_cases") %>%
  #   arrange(grp) %>%
  #   mutate(year = rep(yr1_inc_cancer:yrlast_inc_cancer, each = 12),
  #          year_mon = paste(year, month.abb[1:12], sep = "_")) %>%
  #   bind_cols(., tot_pop_tbl %>%
  #               pivot_longer(1:yrs_inc_cancer, names_to = "tot_pop",
  #                            values_to = "tot_pop_val") %>%
  #               arrange(tot_pop) %>% dplyr::select(-tot_pop)) %>%
  #   mutate(mon_incidence = new_cases / tot_pop_val,
  #          temp_mon_num = rep(c(1:12), yrs_inc_cancer)) %>%
  #   left_join(., year_mon_cycle_tbl, by = c("year",
  #                                           "temp_mon_num" = "mon"))
  # 
  # inc_cancer_yrly <- inc_cancer %>% 
  #   mutate(jul_pop_val = ifelse(grepl("Jul", year_mon), tot_pop_val, NA)) %>% 
  #   group_by(year) %>% 
  #   summarise(new_cases = sum(new_cases),
  #             weighted_inc = weighted.mean(mon_incidence, tot_pop_val),
  #             jul_tot_pop_val =  sum(jul_pop_val, na.rm = T),
  #             jul_yrly_inc = new_cases/jul_tot_pop_val) %>% 
  #   ungroup() ##### not sure how to calculate the annual total population value, 
  # # im not going to use incidence i will use total number of new cases
  # 
  # # cancer prevalence 
  # yr1_prev_cancer <- 2015
  # yrlast_prev_cancer <- 2022
  # yrs_prev_cancer <- length(2015:2022)
  # prev_cancer <- data.frame(matrix(NA, byrow = T, nrow = 12, ncol = yrs_prev_cancer,
  #                                 list(NULL, paste0("tot_cancer_cases_",
  #                                                   str_sub(yr1_prev_cancer:yrlast_prev_cancer, 3, 4))))) %>% 
  #   pivot_longer(1:yrs_prev_cancer, names_to = "grp", values_to = "tot_cancer_cases") %>% 
  #   arrange(grp) %>% 
  #   mutate(year = rep(yr1_prev_cancer:yrlast_prev_cancer, each = 12),
  #          year_mon = paste(year, month.abb[1:12], sep = "_"),
  #          mon = rep(1:12, yrs_prev_cancer)) %>% 
  #   bind_cols(., data.frame(matrix(NA, byrow = T, nrow = 12, ncol = yrs_prev_cancer,
  #                                  list(NULL, paste0("tot_pop_",
  #                                                    str_sub(yr1_prev_cancer:yrlast_prev_cancer, 3, 4))))) %>% 
  #               pivot_longer(1:yrs_prev_cancer, names_to = "tot_pop", values_to = "tot_pop_val") %>% 
  #               arrange(tot_pop) %>% dplyr::select(-tot_pop)) %>% 
  #   left_join(., year_mon_cycle_tbl, by = c("year", "mon"))
  #   
  # 
  # y0 <- sum(m_M[prev_cancer$cycle[1] + 1, c("BN_CANCER", "BPO_CANCER")]) + 
  #   (0.7*(sum(m_M[prev_cancer$cycle[1] + 1, c("BN_PALLIATIVE", "BPO_PALLIATIVE")])))
  # 
  # 
  # prev_cancer[1, "tot_cancer_cases"] <- inc_cancer$new_cases[inc_cancer$cycle == prev_cancer$cycle[1]] + 
  #   y0 - (
  #   m_M[prev_cancer$cycle[1] + 1, "BN_PN"]*m_P_temp["BN_CANCER", "BO_DEATH"] + 
  #     m_M[prev_cancer$cycle[1] + 1, "BN_CANCER"]*m_P_temp["BN_CANCER", "BN_PN"] +
  #     m_M[prev_cancer$cycle[1] + 1, "BPO_CANCER"]*m_P_temp["BPO_CANCER", "BO_DEATH"]
  # )
  # 
  # prev_cancer[1, "tot_pop_val"] <- sum(m_M[1,])
  # 
  # for (i in 2:nrow(prev_cancer)) {
  #   
  #   this_cycle <- prev_cancer$cycle[i]
  #   
  #   prev_cancer[i, "tot_cancer_cases"] <- inc_cancer$new_cases[inc_cancer$cycle == this_cycle] + 
  #     prev_cancer[this_cycle, "tot_cancer_cases"] - (
  #     m_M[this_cycle + 1, "BN_PN"]*m_P_temp["BN_CANCER", "BO_DEATH"] + 
  #       m_M[this_cycle + 1, "BN_CANCER"]*m_P_temp["BN_CANCER", "BN_PN"] +
  #       m_M[this_cycle + 1, "BPO_CANCER"]*m_P_temp["BPO_CANCER", "BO_DEATH"]
  #   )
  #   
  #   prev_cancer[i, "tot_pop_val"] <- sum(m_M[this_cycle + 1,])
  # }
  # 
  # prev_cancer <- prev_cancer %>% 
  #   mutate(mon_prev = tot_cancer_cases / tot_pop_val)
  # 
  # prev_cancer_yrly <- prev_cancer %>% 
  #   group_by(year) %>% 
  #   summarise(weighted_prev = weighted.mean(mon_prev, tot_pop_val)) %>% 
  #   ungroup()
  
  return(list(tot_deaths = calculated_targets$tot_deaths,
              inc_oddeath = calculated_targets$inc_oddeath,
              prev_opioids_rx = calculated_targets$prev_opioids_rx,
              prev_oat = calculated_targets$prev_oat,
              inc_cancer = calculated_targets$inc_cancer,
              prev_cancer = calculated_targets$prev_cancer,
              m_M = m_M)
  )
}
