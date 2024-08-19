# Loading required libraries 
library(tidyverse); library(rstan);library(loo)

# Loading in and processing data 
max_age <- 75
df <- readRDS("data/processed_serological_data_ageAggregated.rds")

## Setting up the data for the models
foi_mapping_df <- data.frame(year_range = c(paste0(min(df$first_year_exposed), "-1974"), "1975-1984", "1985-1994", "1995-2004", 2005:2021))
foi_mapping_df$foi_number <- seq(1:length(foi_mapping_df$year_range))
foi_mapping_df$rep <- c(1974 - min(df$first_year_exposed) + 1, 10, 10, 10, rep(1, length(2005:2021)))
N <- nrow(df)
A <- length(unique(df$age))
years <- length(min(df$first_year_exposed):max(df$final_year_exposed))
num_unique_foi <- max(foi_mapping_df$foi_number)
foi_year_index <- rep(foi_mapping_df$foi_number, times = foi_mapping_df$rep)
X <- as.matrix(df[c("first_year_exposed_stand", "final_year_exposed_stand", "age")])
alpha <- 0.02
num_age_group_FOIs <- length(unique(df$age_group))
number_ages <- length(unique(df$age))
age_group_index <- rep(1:num_age_group_FOIs, times = c(6, 18-5, max_age - 18))
y <- df$seropositive
tot <- df$total
sens <- 0.892
spec <- 0.988
zika_cross <- 0.1
zika_year <- 2016 - min(df$first_year_exposed) + 1
alpha_mean <- 0.02
alpha_sd <- 0.01
iter <- 2500
chains <- 3

fresh_run <- FALSE

## Model 1: Time-Varying FOI only

### Setting up and running the model
if (fresh_run) {
  model1 <- rstan::stan_model("scripts/angolaFOI_Model1_TimeVaryingFOI.stan")
  model1_stan_data <- list(
    N = N,
    A = A,
    years = years,
    num_unique_foi = num_unique_foi,
    foi_year_index = foi_year_index,
    X = X,
    y = y,
    tot = tot, 
    sens = sens,
    spec = spec,
    zika_year = zika_year,
    zika_cross = zika_cross)
  model1_fit <- rstan::sampling(model1, data = model1_stan_data, iter = iter, chains = chains, cores = 3)
  saveRDS(model1_fit, "outputs/angolaFOI_Model1_TimeVaryingFOI_Outputs.rds")
} else {
  model1_fit <- readRDS("outputs/angolaFOI_Model1_TimeVaryingFOI_Outputs.rds")
}
model1_params <- rstan::extract(model1_fit)
model1_outputs <- data.frame(
  point_estimate = apply(model1_params$lambda, 2, mean),
  lower = apply(model1_params$lambda, 2, quantile, 0.05),
  upper = apply(model1_params$lambda, 2, quantile, 0.95),
  year = foi_mapping_df$year_range,
  zika_cross = zika_cross)

### Calculating LOO-CV
log_lik_1 <- extract_log_lik(model1_fit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_1), cores = 2) 
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)

### Posterior predictive checks
y_draw_model1 <- data.frame(
  point_estimate = apply(model1_params$y_draw, 2, mean),
  lower = apply(model1_params$y_draw, 2, quantile, 0.05),
  upper = apply(model1_params$y_draw, 2, quantile, 0.95),
  year = df$final_year_exposed,
  age = df$age,
  N = df$total) %>% 
  mutate(
    point_estimate = point_estimate / N,
    lower = lower / N,
    upper = upper / N)

model1_fit_plot <- ggplot(y_draw_model1) +
  geom_ribbon(aes(x = age, ymax = upper, ymin = lower), fill = "grey20", alpha = 0.4) +
  geom_point(data = df %>% 
               rename(year = final_year_exposed), aes(x = age, y = IgG_positive/IgG_tested),
             col = "black", pch = 21, fill = "#95C9C2", size = 2) +
  geom_line(aes(x = age, y = point_estimate), col = "grey20", linewidth = 1) +
  facet_wrap(~year) +
  labs(x = "Age (years)", y = "Seroprevalence")

## Model 2: Time-Varying and Age-Specific FOI Modifiers only

### Setting up and running the model
if (fresh_run) {
  model2 <- rstan::stan_model("scripts/angolaFOI_Model2_TimeVaryingFOI_AgeVariableFOI.stan")
  model2_stan_data <- list(
    N = N,
    A = A,
    years = years,
    num_unique_foi = num_unique_foi,
    foi_year_index = foi_year_index,
    X = X,
    num_age_group_FOIs = num_age_group_FOIs,
    number_ages = number_ages,
    age_group_index = age_group_index,
    y = y,
    tot = tot, 
    sens = sens,
    spec = spec,
    zika_year = zika_year,
    zika_cross = zika_cross)
  model2_fit <- rstan::sampling(model2, data = model2_stan_data, iter = iter, chains = chains, cores = 3)
  saveRDS(model2_fit, "outputs/angolaFOI_Model2_TimeVaryingFOI_AgeVariableFOI_Outputs.rds")
} else {
  model2_fit <- readRDS("outputs/angolaFOI_Model2_TimeVaryingFOI_AgeVariableFOI_Outputs.rds")
}
model2_params <- rstan::extract(model2_fit)
model2_outputs <- data.frame(
  point_estimate = apply(model2_params$lambda, 2, mean),
  lower = apply(model2_params$lambda, 2, quantile, 0.05),
  upper = apply(model2_params$lambda, 2, quantile, 0.95),
  year = foi_mapping_df$year_range,
  zika_cross = zika_cross)
apply(model2_params$Mu_age_raw, 2, mean)

### Calculating LOO-CV
log_lik_2 <- extract_log_lik(model2_fit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_2), cores = 2) 
loo_2 <- loo(log_lik_2, r_eff = r_eff, cores = 2)
loo_compare(loo_1, loo_2)

### Posterior predictive checks
y_draw_model2 <- data.frame(
  point_estimate = apply(model2_params$y_draw, 2, mean),
  lower = apply(model2_params$y_draw, 2, quantile, 0.05),
  upper = apply(model2_params$y_draw, 2, quantile, 0.95),
  year = df$final_year_exposed,
  age = df$age,
  N = df$total) %>% 
  mutate(
    point_estimate = point_estimate / N,
    lower = lower / N,
    upper = upper / N)

model2_fit_plot <- ggplot(y_draw_model2) +
  geom_ribbon(aes(x = age, ymax = upper, ymin = lower), fill = "grey20", alpha = 0.4) +
  geom_point(data = df %>% 
               rename(year = final_year_exposed), aes(x = age, y = IgG_positive/IgG_tested),
             col = "black", pch = 21, fill = "#95C9C2", size = 2) +
  geom_line(aes(x = age, y = point_estimate), col = "grey20", linewidth = 1) +
  facet_wrap(~year) +
  labs(x = "Age (years)", y = "Seroprevalence")

## Model 3: Time-Varying FOI and Seroreversion 

### Setting up and running the model
if (fresh_run) {
  model3 <- rstan::stan_model("scripts/angolaFOI_Model3_TimeVaryingFOI_Seroreversion.stan")
  model3_stan_data <- list(
    N = N,
    A = A,
    years = years,
    num_unique_foi = num_unique_foi,
    foi_year_index = foi_year_index,
    X = X,
    y = y,
    tot = tot, 
    sens = sens,
    spec = spec,
    zika_year = zika_year,
    zika_cross = zika_cross,
    alpha_mean = alpha_mean,
    alpha_sd = alpha_sd)
  # model3_fit <- rstan::sampling(model3, data = model3_stan_data, iter = iter, chains = chains, cores = 3)
  model3_fit <- rstan::sampling(model3, data = model3_stan_data, iter = 2500, chains = 1, cores = 1)
  saveRDS(model3_fit, "outputs/angolaFOI_Model3_TimeVaryingFOI_Seroreversion_Outputs.rds")
} else {
  model3_fit <- readRDS("outputs/angolaFOI_Model3_TimeVaryingFOI_Seroreversion_Outputs.rds")
}
model3_params <- rstan::extract(model3_fit)
model3_outputs <- data.frame(
  point_estimate = apply(model3_params$lambda, 2, mean),
  lower = apply(model3_params$lambda, 2, quantile, 0.05),
  upper = apply(model3_params$lambda, 2, quantile, 0.95),
  year = foi_mapping_df$year_range,
  zika_cross = zika_cross)
mean(model3_params$alpha)
median(model3_params$alpha)
quantile(model3_params$alpha, prob = c(0.025, 0.975))

### Calculating LOO-CV
log_lik_3 <- extract_log_lik(model3_fit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_3), cores = 2) 
loo_3 <- loo(log_lik_3, r_eff = r_eff, cores = 2)
loo_compare(loo_1, loo_2)
loo_compare(loo_1, loo_3)
loo_compare(loo_2, loo_3)

### Posterior predictive checks
y_draw_model3 <- data.frame(
  point_estimate = apply(model3_params$y_draw, 2, mean),
  lower = apply(model3_params$y_draw, 2, quantile, 0.05),
  upper = apply(model3_params$y_draw, 2, quantile, 0.95),
  year = df$final_year_exposed,
  age = df$age,
  N = df$total) %>% 
  mutate(
    point_estimate = point_estimate / N,
    lower = lower / N,
    upper = upper / N)

model3_fit_plot <- ggplot(y_draw_model3) +
  geom_ribbon(aes(x = age, ymax = upper, ymin = lower), fill = "grey20", alpha = 0.4) +
  geom_point(data = df %>% 
               rename(year = final_year_exposed), aes(x = age, y = IgG_positive/IgG_tested),
             col = "black", pch = 21, fill = "#95C9C2", size = 2) +
  geom_line(aes(x = age, y = point_estimate), col = "grey20", linewidth = 1) +
  facet_wrap(~year) +
  labs(x = "Age (years)", y = "Seroprevalence")

## Model 4: Time-Varying FOI and Seroreversion and Age-Specific effects

### Setting up and running the model
model4 <- rstan::stan_model("models/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI.stan")
if (fresh_run) {
  model4_stan_data <- list(
    N = N,
    A = A,
    years = years,
    num_unique_foi = num_unique_foi,
    foi_year_index = foi_year_index,
    X = X,
    num_age_group_FOIs = num_age_group_FOIs,
    number_ages = number_ages,
    age_group_index = age_group_index,
    y = y,
    tot = tot, 
    sens = sens,
    spec = spec,
    zika_year = zika_year,
    zika_cross = zika_cross,
    alpha_mean = alpha_mean,
    alpha_sd = alpha_sd)
  model4_fit <- rstan::sampling(model4, data = model4_stan_data, iter = iter, chains = chains, cores = 3)
  saveRDS(model4_fit, "outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_Outputs.rds")
} else {
  model4_fit <- readRDS("outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_Outputs.rds")
}

model4_params <- rstan::extract(model4_fit)
model4_outputs <- data.frame(
  point_estimate = apply(model4_params$lambda, 2, mean),
  lower = apply(model4_params$lambda, 2, quantile, 0.05),
  upper = apply(model4_params$lambda, 2, quantile, 0.95),
  year = foi_mapping_df$year_range,
  zika_cross = zika_cross)
# plot(2005:2021, model4_outputs$point_estimate[5:21], type = "l", col = "black")
# lines(2005:2021, model1_outputs$point_estimate[5:21], col = "blue")
# lines(2005:2021, model2_outputs$point_estimate[5:21], col = "red")
# lines(2005:2021, model3_outputs$point_estimate[5:21], col = "orange")
apply(model4_params$Mu_age_raw, 2, mean)

### Calculating LOO-CV
log_lik_4 <- extract_log_lik(model4_fit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_4), cores = 2) 
loo_4 <- loo(log_lik_4, r_eff = r_eff, cores = 2)
print(loo_4)
loo_compare(loo_4, loo_1) # time-varying FOI only
loo_compare(loo_4, loo_2) # time-varying FOI and age-specific modifiers
loo_compare(loo_4, loo_3) # time-varying FOI and seroreversion
loo_compare(loo_3, loo_1) 

### Posterior predictive checks
y_draw_model4 <- data.frame(
  point_estimate = apply(model4_params$y_draw, 2, mean),
  lower = apply(model4_params$y_draw, 2, quantile, 0.05),
  upper = apply(model4_params$y_draw, 2, quantile, 0.95),
  year = df$final_year_exposed,
  age = df$age,
  N = df$total) %>% 
  mutate(
    point_estimate = point_estimate / N,
    lower = lower / N,
    upper = upper / N)

model4_fit_plot <- ggplot(y_draw_model4) +
  geom_ribbon(aes(x = age, ymax = upper, ymin = lower), fill = "grey20", alpha = 0.4) +
  geom_point(data = df %>% 
               rename(year = final_year_exposed), aes(x = age, y = IgG_positive/IgG_tested),
             col = "black", pch = 21, fill = "#95C9C2", size = 2) +
  geom_line(aes(x = age, y = point_estimate), col = "grey20", linewidth = 1) +
  facet_wrap(~year) +
  labs(x = "Age (years)", y = "Seroprevalence")


cowplot::plot_grid(model1_fit, model2_fit, model3_fit, model4_fit,
                   nrow = 2, ncol = 2)

# Final figure generation

## No Zika cross reactivity
if (fresh_run) {
  model4_stan_data_noCross <- model4_stan_data
  model4_stan_data_noCross$zika_cross <- 0.0
  model4_fit_noCross <- rstan::sampling(model4, data = model4_stan_data_noCross, 
                                        iter = iter, chains = chains, cores = 3)
  saveRDS(model4_fit_noCross, 
          "outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_NoCross_Outputs.rds")
} else {
  model4_fit_noCross <- readRDS("outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_NoCross_Outputs.rds")
}
no_cross_params <- rstan::extract(model4_fit_noCross)
no_cross_plt_foi <- data.frame(
  point_estimate = apply(no_cross_params$lambda, 2, mean),
  lower = apply(no_cross_params$lambda, 2, quantile, 0.05),
  upper = apply(no_cross_params$lambda, 2, quantile, 0.95),
  year = foi_mapping_df$year_range,
  zika_cross = 0)

apply(no_cross_params$Mu_age_raw, 2, mean)

no_cross_plt_foi %>% 
  filter(year %in% paste0(2005:2021)) %>%
  filter(zika_cross == 0) %>%
  summarise(index = which(point_estimate == min(point_estimate)),
            year = year[index],
            median = point_estimate[index],
            lower = lower[index],
            upper = upper[index])

no_cross_plt_foi %>% 
  filter(year %in% paste0(2005:2021)) %>%
  filter(zika_cross == 0) %>%
  summarise(index = which(point_estimate == max(point_estimate)),
            year = year[index],
            median = point_estimate[index],
            lower = lower[index],
            upper = upper[index])

## Low Zika cross reactivity
if (fresh_run) {
  model4_stan_data_lowCross <- model4_stan_data
  model4_stan_data_lowCross$zika_cross <- 0.1
  model4_fit_lowCross <- rstan::sampling(model4, data = model4_stan_data_lowCross, 
                                         iter = iter, chains = chains, cores = 3)
  saveRDS(model4_fit_lowCross, 
          "outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_LowCross_Outputs.rds")
} else {
  model4_fit_lowCross <- readRDS("outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_LowCross_Outputs.rds")
}
low_cross_params <- rstan::extract(model4_fit_lowCross)
low_cross_plt_foi <- data.frame(
  point_estimate = apply(low_cross_params$lambda, 2, mean),
  lower = apply(low_cross_params$lambda, 2, quantile, 0.05),
  upper = apply(low_cross_params$lambda, 2, quantile, 0.95),
  year =  foi_mapping_df$year_range,
  zika_cross = 0.1)

## High Zika cross reactivity
if (fresh_run) {
  model4_stan_data_midCross <- model4_stan_data
  model4_stan_data_midCross$zika_cross <- 0.3
  model4_fit_midCross <- rstan::sampling(model4, data = model4_stan_data_midCross, 
                                         iter = iter, chains = chains, cores = 3)
  saveRDS(model4_fit_midCross, 
          "outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_MidCross_Outputs.rds")
} else{
  model4_fit_midCross <- readRDS("outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_MidCross_Outputs.rds")
}
mid_cross_params <- rstan::extract(model4_fit_midCross)
mid_cross_plt_foi <- data.frame(
  point_estimate = apply(mid_cross_params$lambda, 2, mean),
  lower = apply(mid_cross_params$lambda, 2, quantile, 0.05),
  upper = apply(mid_cross_params$lambda, 2, quantile, 0.95),
  year =  foi_mapping_df$year_range,
  zika_cross = 0.3)

## Summarising and plotting the results
overall <- rbind(no_cross_plt_foi, low_cross_plt_foi, mid_cross_plt_foi) %>%
  filter(year %in% paste0(2005:2021)) %>%
  mutate(year = as.numeric(year)) %>%
  filter(zika_cross == 0 | (zika_cross == 0.1 & year >= 2016) | (zika_cross == 0.3 & year >= 2016)) %>%
  mutate(plot = case_when(zika_cross == 0.0 & year < 2016 ~ "aBefore Zika",
                          zika_cross == 0.3 & year >= 2016 ~ "cAfter Zika - High Cross-Reactivity",
                          zika_cross == 0.1 & year >= 2016 ~ "bAfter Zika - Low Cross-Reactivity",
                          zika_cross == 0.0 & year >= 2016 ~ "cAfter Zika - No Cross-Reactivity"))
extra <- overall[overall$zika_cross == 0.0 & overall$year == 2016, ]
extra$plot <- "aBefore Zika"
overall2 <- rbind(overall, extra)

seroprev <- ggplot(overall2) +
  geom_line(aes(x = year, y = point_estimate, col = plot)) +
  geom_ribbon(aes(x = year, ymax = upper, ymin = lower, fill = plot), alpha = 0.2) +
  labs(x = "Year", y = "Yearly force of infection") +
  theme_bw() +
  scale_fill_manual(values = c("#948D9B", "#9B598E","#D2B1CC", "#7D4271"),
                    labels = c("Pre-Zika Period", "Low Cross-Reactivity",
                               "High Cross-Reactivity", "No Cross-Reactivity"),
                    name = "Degree of Zika Cross-Reactivity") +
  scale_colour_manual(values = c("#948D9B",  "#9B598E", "#D2B1CC","#7D4271"),
                      labels = c("Pre-Zika Period", "Low Cross-Reactivity",
                                 "High Cross-Reactivity", "No Cross-Reactivity"),
                      name = "Degree of Zika Cross-Reactivity") +
  coord_cartesian(xlim = c(2005, 2021)) +
  theme(legend.position = c(0.25, 0.65))

overall2 %>%
  filter(year >= 2005) %>%
  filter(zika_cross == 0) %>%
  summarise(index = which(point_estimate == max(point_estimate)),
            year = year[index],
            median = point_estimate[index],
            lower = lower[index],
            upper = upper[index])
overall2 %>%
  filter(year >= 2005) %>%
  filter(zika_cross == 0) %>%
  summarise(index = which(point_estimate == min(point_estimate)),
            year = year[index],
            median = point_estimate[index],
            lower = lower[index],
            upper = upper[index])

## Separate location FOI calculations

# Loading in and processing data

## Model fitting for Luanda
luanda_df <- readRDS(file = "data/processedLuanda_seroData.rds")
luanda_foi_mapping_df <- data.frame(year_range = c(paste0(min(luanda_df$first_year_exposed), "-1974"), "1975-1984", "1985-1994", "1995-2004", 2005:2021))
luanda_foi_mapping_df$foi_number <- seq(1:length(luanda_foi_mapping_df$year_range))
luanda_foi_mapping_df$rep <- c(1974 - min(luanda_df$first_year_exposed) + 1, 10, 10, 10, rep(1, length(2005:2021)))

luanda_stan_data <- list(
  N = nrow(luanda_df),
  A = length(unique(luanda_df$age)),
  years = length(min(luanda_df$first_year_exposed):max(luanda_df$final_year_exposed)),
  num_unique_foi = max(luanda_foi_mapping_df$foi_number),
  foi_year_index = rep(luanda_foi_mapping_df$foi_number, times = luanda_foi_mapping_df$rep),
  X = as.matrix(luanda_df[c("first_year_exposed_stand", "final_year_exposed_stand", "age")]),
  num_age_group_FOIs = length(unique(luanda_df$age_group)),
  number_ages = length(unique(luanda_df$age)),
  age_group_index = rep(1:num_age_group_FOIs, times = c(6, 18-5, max_age - 18)),
  y = luanda_df$seropositive,
  tot = luanda_df$total, 
  sens = sens,
  spec = spec,
  zika_year = 2016 - min(luanda_df$first_year_exposed) + 1,
  zika_cross = zika_cross,
  alpha_mean = alpha_mean,
  alpha_sd = alpha_sd)

if (fresh_run) {
  fit_luanda <- rstan::sampling(model4, data = luanda_stan_data, iter = 3000, chains = 1)
  saveRDS(fit_luanda, 
          "outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_LowCross_LuandaONLY_Outputs.rds")
  
} else {
  fit_luanda <- readRDS("outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_LowCross_LuandaONLY_Outputs.rds")
}
luanda_params <- rstan::extract(fit_luanda)

## Model fitting for Viana
viana_df <- readRDS("data/processedViana_seroData.rds")
viana_foi_mapping_df <- data.frame(year_range = c(paste0(min(viana_df$first_year_exposed), "-1974"), "1975-1984", "1985-1994", "1995-2004", 2005:2021))
viana_foi_mapping_df$foi_number <- seq(1:length(viana_foi_mapping_df$year_range))
viana_foi_mapping_df$rep <- c(1974 - min(viana_df$first_year_exposed) + 1, 10, 10, 10, rep(1, length(2005:2021)))

viana_stan_data <- list(
  N = nrow(viana_df),
  A = length(unique(viana_df$age)),
  years = length(min(viana_df$first_year_exposed):max(viana_df$final_year_exposed)),
  num_unique_foi = max(viana_foi_mapping_df$foi_number),
  foi_year_index = rep(viana_foi_mapping_df$foi_number, times = viana_foi_mapping_df$rep),
  X = as.matrix(viana_df[c("first_year_exposed_stand", "final_year_exposed_stand", "age")]),
  num_age_group_FOIs = length(unique(viana_df$age_group)),
  number_ages = length(unique(viana_df$age)),
  age_group_index = rep(1:num_age_group_FOIs, times = c(6, 18-5, max_age - 18)),
  y = viana_df$seropositive,
  tot = viana_df$total, 
  sens = sens,
  spec = spec,
  zika_year = 2016 - min(viana_df$first_year_exposed) + 1,
  zika_cross = zika_cross,
  alpha_mean = alpha_mean,
  alpha_sd = alpha_sd)

if (fresh_run) {
  fit_viana <- rstan::sampling(model4, data = viana_stan_data, iter = 3000, chains = 1)
  saveRDS(fit_viana, 
          "outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_LowCross_VianaONLY_Outputs.rds")
  
} else {
  fit_viana <- readRDS("outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_LowCross_VianaONLY_Outputs.rds")
}
viana_params <- rstan::extract(fit_viana)

## Model fitting for Zango
zango_df <- readRDS("data/processedZango_seroData.rds")
zango_foi_mapping_df <- data.frame(year_range = c(paste0(min(zango_df$first_year_exposed), "-1974"), "1975-1984", "1985-1994", "1995-2004", 2005:2021))
zango_foi_mapping_df$foi_number <- seq(1:length(zango_foi_mapping_df$year_range))
zango_foi_mapping_df$rep <- c(1974 - min(zango_df$first_year_exposed) + 1, 10, 10, 10, rep(1, length(2005:2021)))

zango_stan_data <- list(
  N = nrow(zango_df),
  A = length(min(zango_df$age):max(zango_df$age)),
  years = length(min(zango_df$first_year_exposed):max(zango_df$final_year_exposed)),
  num_unique_foi = max(zango_foi_mapping_df$foi_number),
  foi_year_index = rep(zango_foi_mapping_df$foi_number, times = zango_foi_mapping_df$rep),
  X = as.matrix(zango_df[c("first_year_exposed_stand", "final_year_exposed_stand", "age")]),
  num_age_group_FOIs = length(unique(zango_df$age_group)),
  number_ages = length(min(zango_df$age):max(zango_df$age)) + 1,
  age_group_index = rep(1:num_age_group_FOIs, times = c(6, 18-5, max_age - 18)),
  y = zango_df$seropositive,
  tot = zango_df$total, 
  sens = sens,
  spec = spec,
  zika_year = zika_year,
  zika_cross = zika_cross,
  alpha_mean = alpha_mean,
  alpha_sd = alpha_sd)

if (fresh_run) {
  fit_zango <- rstan::sampling(model4, data = zango_stan_data, iter = 3000, chains = 1)
  saveRDS(fit_zango, 
          "outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_LowCross_ZangoONLY_Outputs.rds")
  
} else {
  fit_zango <- readRDS("outputs/angolaFOI_Model4_TimeVaryingFOI_Seroreversion_AgeVariableFOI_LowCross_ZangoONLY_Outputs.rds")
}
zango_params <- rstan::extract(fit_zango)

plt_foi_luanda <- data.frame(
  location = "Luanda",
  point_estimate = apply(luanda_params$lambda, 2, mean),
  lower = apply(luanda_params$lambda, 2, quantile, 0.05),
  upper = apply(luanda_params$lambda, 2, quantile, 0.95),
  year = luanda_foi_mapping_df$year_range)
plt_foi_viana <- data.frame(
  location = "Viana",
  point_estimate = apply(viana_params$lambda, 2, mean),
  lower = apply(viana_params$lambda, 2, quantile, 0.05),
  upper = apply(viana_params$lambda, 2, quantile, 0.95),
  year = viana_foi_mapping_df$year_range)
plt_foi_zango <- data.frame(
  location = "Zango",
  point_estimate = apply(zango_params$lambda, 2, mean),
  lower = apply(zango_params$lambda, 2, quantile, 0.05),
  upper = apply(zango_params$lambda, 2, quantile, 0.95),
  year = zango_foi_mapping_df$year_range)
plt_foi_overall <- rbind(plt_foi_luanda, plt_foi_viana, plt_foi_zango) %>%
  filter(year %in% paste0(2005:2021)) %>%
  mutate(year = as.numeric(year))

mean(luanda_params$alpha)
mean(viana_params$alpha)
mean(zango_params$alpha)

location_specific_plot <- ggplot(plt_foi_overall) +
  geom_line(aes(x = year, y = point_estimate, col = location)) +
  geom_ribbon(aes(x = year, ymax = upper, ymin = lower, fill = location), alpha = 0.2) +
  labs(x = "Year", y = "Yearly force of infection") +
  theme_bw() +
  facet_grid(.~location, scales = "free_y") +
  coord_cartesian(xlim = c(2005, 2021), ylim = c(0, 0.45)) +
  theme(strip.background = element_rect(fill = "white"),
        legend.position = "none")

x <- cowplot::plot_grid(seroprev, location_specific_plot, nrow  = 2, labels = c("A", "B"))
ggsave(plot = x,
       filename = "updated_sero_figure.pdf",
       width = 7.5,
       height = 5.75)
