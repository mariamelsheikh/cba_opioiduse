library(here)

source(here("02_scripts/03a_fun_calib_td.R"))
theme_set(theme_minimal())

opt_params <- readRDS(here("01_data/map_calib_td.RDS"))

mod_opt_map <- mod_calib(as.numeric(opt_params$par))


# prevalence of rx opioids use --------------------------------------------


v_yrlast_prev_rx <- c(2015:2018)
yrs_prev_rx <- length(v_yrlast_prev_rx)

prop_opioids_rx_calib <- prop_opioids_rx_uncalib <- data.frame(matrix(NA, byrow = T,
                                                                      nrow = 12,
                                                                      ncol = yrs_prev_rx,
                                                                      list(NULL,
                                                                           paste0("prop_opioids_rx_",
                                                                                  str_sub(v_yrlast_prev_rx, 3, 4)))))
tot_pop_calib_tbl <- tot_pop_uncalib_tbl <- data.frame(matrix(NA, byrow = T,
                                                              nrow = 12,
                                                              ncol = yrs_prev_rx,
                                                              list(NULL,
                                                                   paste0("tot_pop_",
                                                                          str_sub(v_yrlast_prev_rx, 3, 4)))))

for (i in 1:yrs_prev_rx) {
  yr <- (v_yrlast_prev_rx[1] - 1) + i
  prop_opioids_rx_uncalib[, i] <- assign(
    paste0("prop_opioids_rx_", str_sub(yr, 3, 4)),
    rowSums(mod_basecase$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,
                             c("BPO_ACUTE", "BPO_CHRONIC", "BPO_CANCER",
                               "BPO_OTHER", "BPO_PALLIATIVE")]) /
      rowSums(mod_basecase$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
    
  )
  
  tot_pop_uncalib_tbl[, i] <- assign(
    paste0("tot_pop_", str_sub(yr, 3, 4)),
    rowSums(mod_basecase$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
  )
  
  prop_opioids_rx_calib[, i] <- assign(
    paste0("prop_opioids_rx_", str_sub(yr, 3, 4)),
    rowSums(mod_opt_map$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,
                            c("BPO_ACUTE", "BPO_CHRONIC", "BPO_CANCER",
                              "BPO_OTHER", "BPO_PALLIATIVE")]) /
      rowSums(mod_opt_map$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
    
  )
  
  tot_pop_calib_tbl[, i] <- assign(
    paste0("tot_pop_", str_sub(yr, 3, 4)),
    rowSums(mod_opt_map$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr] + 1,])
  )
}

prop_opioids_rx_uncalib_tbl <- prop_opioids_rx_uncalib %>% 
  pivot_longer(1:4, names_to = "grp", values_to = "target_val") %>% 
  arrange(grp) %>% 
  mutate(year = rep(2015:2018,
                    each = 12),
         year_mon = paste(year, month.abb[1:12], sep = "_")) %>% 
  bind_cols(., tot_pop_uncalib_tbl %>% 
              pivot_longer(1:4, names_to = "tot_pop",
                           values_to = "tot_pop_val") %>% 
              arrange(tot_pop) %>% select(-tot_pop))


prop_opioids_rx_calib_tbl <- prop_opioids_rx_calib %>% 
  pivot_longer(1:4, names_to = "grp", values_to = "target_val") %>% 
  arrange(grp) %>% 
  mutate(year = rep(2015:2018,
                    each = 12),
         year_mon = paste(year, month.abb[1:12], sep = "_")) %>% 
  bind_cols(., tot_pop_calib_tbl %>% 
              pivot_longer(1:4, names_to = "tot_pop",
                           values_to = "tot_pop_val") %>% 
              arrange(tot_pop) %>% select(-tot_pop))


prop_opioids_rx_target_tbl <- calib_target_tbl %>% 
  filter(group == "prev_on_opioidsrx") %>% 
  select(year, target, group) %>% 
  rename(grp = group,
         target_val = target) %>% 
  mutate(year_mon = paste(year, month.abb[6], sep = "_"))

prop_opioids_rx_uncalib_wei_mean <- prop_opioids_rx_uncalib_tbl %>% 
  group_by(year) %>% 
  summarise(target_val = weighted.mean(target_val, tot_pop_val)) %>% 
  ungroup() %>% 
  mutate(grp = "Uncalibrated Model") %>% 
  bind_rows(., prop_opioids_rx_target_tbl %>% 
              select(year, target_val) %>% 
              mutate(grp = "target")) %>% 
  bind_rows(., prop_opioids_rx_calib_tbl %>% 
              group_by(year) %>% 
              summarise(target_val = weighted.mean(target_val, tot_pop_val)) %>% 
              ungroup() %>% 
              mutate(grp = "Calibrated Model"))

prop_opioids_rx_target_tbl$year_mon <- factor(prop_opioids_rx_target_tbl$year_mon,
                                              levels = paste(rep(2015:2018,
                                                                 each = 12),
                                                             month.abb[1:12],
                                                             sep = "_"),
                                              labels = paste(rep(2015:2018,
                                                                 each = 12),
                                                             month.abb[1:12],
                                                             sep = "_"))
prop_opioids_rx_uncalib_tbl$year_mon <- factor(prop_opioids_rx_uncalib_tbl$year_mon,
                                               levels = paste(rep(2015:2018,
                                                                  each = 12),
                                                              month.abb[1:12],
                                                              sep = "_"),
                                               labels = paste(rep(2015:2018,
                                                                  each = 12),
                                                              month.abb[1:12],
                                                              sep = "_"))
prop_opioids_rx_calib_tbl$year_mon <- factor(prop_opioids_rx_calib_tbl$year_mon,
                                             levels = paste(rep(2015:2018,
                                                                each = 12),
                                                            month.abb[1:12],
                                                            sep = "_"),
                                             labels = paste(rep(2015:2018,
                                                                each = 12),
                                                            month.abb[1:12],
                                                            sep = "_"))

prop_opioids_rx_tbl <- prop_opioids_rx_uncalib_tbl %>% 
  mutate(group = "Uncalibrated Model") %>% 
  bind_rows(.,prop_opioids_rx_calib_tbl %>% 
              mutate(group = "Calibrated Model"))

wprop_opioids_rx_tbl <- prop_opioids_rx_target_tbl %>% 
  mutate(grp = "Target") %>% 
  bind_rows(., prop_opioids_rx_uncalib_tbl %>% 
              group_by(year) %>% 
              summarise(target_val = weighted.mean(target_val, tot_pop_val)) %>% 
              ungroup() %>% 
              mutate(year_mon = prop_opioids_rx_target_tbl$year_mon,
                     grp = "Uncalibrated Model")) %>% 
  bind_rows(., prop_opioids_rx_calib_tbl %>% 
              group_by(year) %>% 
              summarise(target_val = weighted.mean(target_val, tot_pop_val)) %>% 
              ungroup() %>% 
              mutate(year_mon = prop_opioids_rx_target_tbl$year_mon,
                     grp = "Calibrated Model"))

prop_opioids_rx_tbl$group <- factor(prop_opioids_rx_tbl$group,
                                    labels = c("Calibrated Model", "Uncalibrated Model"),
                                    levels = c("Calibrated Model", "Uncalibrated Model"))

wprop_opioids_rx_tbl$grp <- factor(wprop_opioids_rx_tbl$grp,
                                   labels = c("Target", "Calibrated Model", "Uncalibrated Model"),
                                   levels = c("Target", "Calibrated Model", "Uncalibrated Model"))

# Deaths ---------------------------------------------------------------

# number of deaths

num_deaths_calib <- num_deaths_uncalib <- rep(NA, length(2016:2020))
for (i in 1:length(2016:2020)){
  yr <- 2015 + i
  num_deaths_uncalib[i] <- mod_basecase$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                                       year_mon_cycle_tbl$mon == 12] + 1,
                                            "BO_DEATH"] -
    mod_basecase$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
                                                year_mon_cycle_tbl$mon == 12] + 1,
                     "BO_DEATH"]
  
  num_deaths_calib[i] <- mod_opt_map$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                                    year_mon_cycle_tbl$mon == 12] + 1,
                                         "BO_DEATH"] -
    mod_opt_map$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
                                               year_mon_cycle_tbl$mon == 12] + 1,
                    "BO_DEATH"]
}

# number of OD-deaths

num_od_deaths_calib <- num_od_deaths_uncalib <- rep(NA, length(2016:2021))
for (i in 1:length(2016:2021)){
  yr <- 2015 + i 
  num_od_deaths_uncalib[i] <- mod_basecase$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                                          year_mon_cycle_tbl$mon == 12] + 1,
                                               "BO_OD_DEATH"] -
    mod_basecase$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
                                                year_mon_cycle_tbl$mon == 12] + 1,
                     "BO_OD_DEATH"]
  
  num_od_deaths_calib[i] <- mod_opt_map$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                                       year_mon_cycle_tbl$mon == 12] + 1,
                                            "BO_OD_DEATH"] -
    mod_opt_map$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr - 1 &
                                               year_mon_cycle_tbl$mon == 12] + 1,
                    "BO_OD_DEATH"]
}

deaths_target_uncalib_tbl <- calib_target_tbl %>% 
  select(year, target, group) %>% 
  rename(target_val = target,
         target = group) %>% 
  filter(target %in% c("total_deaths", "total_od_deaths")) %>% 
  mutate(target = ifelse(target == "total_deaths", "Total deaths",
                         ifelse(target == "total_od_deaths",
                                "Total OD-deaths", NA)),
         group = "Target") %>% 
  bind_rows(., data.frame(target_val = num_deaths_uncalib) %>% 
              mutate(target = "Total deaths") %>% 
              bind_rows(., data.frame(target_val = num_od_deaths_uncalib) %>%
                          mutate(target = "Total OD-deaths")) %>% 
              mutate(year = c(2016:2020, 2016:2021),
                     group = "Uncalibrated Model")) %>% 
  bind_rows(., data.frame(target_val = num_deaths_calib) %>% 
              mutate(target = "Total deaths") %>% 
              bind_rows(., data.frame(target_val = num_od_deaths_calib) %>%
                          mutate(target = "Total OD-deaths")) %>% 
              mutate(year = c(2016:2020, 2016:2021),
                     group = "Calibrated Model"))

deaths_target_uncalib_tbl$group <- factor(deaths_target_uncalib_tbl$group,
                                          labels = c("Target", "Calibrated Model", "Uncalibrated Model"),
                                          levels = c("Target", "Calibrated Model", "Uncalibrated Model"))


# OAT ---------------------------------------------------------------------

num_oat_calib <- num_oat_uncalib <- rep(NA, length(2018:2021))
for (i in 1:length(2018:2021)){
  yr <- 2017 + i 
  num_oat_uncalib[i] <- sum(mod_basecase$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                                        year_mon_cycle_tbl$mon == 6] + 1,
                                             c("BS_OAT_INI", "BS_OAT_MAINT", "BR_OAT_INI", "BR_OAT_MAINT")])
  
  num_oat_calib[i] <- sum(mod_opt_map$m_M[year_mon_cycle_tbl$cycle[year_mon_cycle_tbl$year == yr &
                                                                     year_mon_cycle_tbl$mon == 6] + 1,
                                          c("BS_OAT_INI", "BS_OAT_MAINT", "BR_OAT_INI", "BR_OAT_MAINT")])
}


oat_target_uncalib_tbl <- calib_target_tbl %>% 
  select(year, target, group) %>% 
  rename(target_val = target,
         target = group) %>% 
  filter(target %in% c("total_oat")) %>% 
  mutate(target = ifelse(target == "total_oat",
                         "Total OAT", NA),
         group = "Target") %>% 
  bind_rows(., data.frame(target_val = num_oat_uncalib) %>%
              mutate(target = "Total OAT",
                     year = c(2018:2021),
                     group = "Uncalibrated Model")) %>% 
  bind_rows(.,  data.frame(target_val = num_oat_calib) %>%
              mutate(target = "Total OAT",
                     year = c(2018:2021),
                     group = "Calibrated Model"))

oat_target_uncalib_tbl$group <- factor(oat_target_uncalib_tbl$group,
                                       labels = c("Target", "Calibrated Model", "Uncalibrated Model"),
                                       levels = c("Target", "Calibrated Model", "Uncalibrated Model"))

