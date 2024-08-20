functions {
  real prob_infected_time_model_single_age(
    int first_year_exposed,
    int final_year_exposed,
    vector foi_vector,
    vector age_modifier,
    real seroreversion_rate) { 
      real prob = 0.0;
      int age = final_year_exposed - first_year_exposed;
      for(j in 1:age) {
        real foi = foi_vector[first_year_exposed + j] * age_modifier[j];
        real lambda_over_both = foi / (foi + seroreversion_rate);
        real e_lower = exp(-(foi + seroreversion_rate));
        prob = lambda_over_both + e_lower * (prob - lambda_over_both);
      }
     return prob;
    }
}

data {
  int<lower=0> N;                // number of data points (years * age groups sampled in each year)
  int<lower=1> A;                // total number of ages in the dataset

  int<lower=0> years;                   // length of yearly foi vector - i.e. number of years covered by the data
  int<lower=0> num_unique_foi;          // number of unique fois (decadal until 2005, and then yearly onwards).
  int<lower=0> foi_year_index[years];  // mapping each year to a particular FOI
  int X[N,3];                           // first year of exposure, and final year of exposure, and age
  int<lower=0> num_age_group_FOIs;     // number of age-groups to estimate FOI modifier for
  int<lower=0> number_ages;            // number of ages in the dataset     
  int<lower=0> age_group_index[number_ages];     // mapping each age-group to a particular FOI modifier

  int y[N];                      // number testing positive in each age-year group
  int tot[N];                    // total number tested in each age-year group
  real<lower=0,upper=1> sens;    // diagnostic sensitivity
  real<lower=0,upper=1> spec;    // diagnostic specificity
  
  int zika_year;                 // year after which we need to account for zika cross-reactivity
  real zika_cross;               // zika cross-reactivity - proportion of dengue seropositive samples that are actually zika seropositive
  real<lower=0,upper=1> alpha_mean;   // seroreversion parameter
  real<lower=0,upper=1> alpha_sd;   // seroreversion parameter
}

parameters {
  real<lower=0.0001, upper = 1> lambda[num_unique_foi]; // number of yearly FOIs that can be estimated from the data
  real<lower=0> Mu_age_raw[num_age_group_FOIs - 1];         // scaler that determines age-specific fraction of yearly FOI
  real<lower=0, upper=1> alpha;
}

transformed parameters {
  real <lower = 0, upper = 1> true_seropositive[N];      // proportion truly seropositive for each age-year group
  real <lower = 0, upper = 1> test_seropositive[N];      // proportion testing seropositive for each age-year group
  real <lower = 0, upper = 1> test_seropositive_zika_adjusted[N]; // proportion testing seropositive for each age-year group post-zika

  vector[years] modelled_lambda;
  for (i in 1:years) {
    modelled_lambda[i] = lambda[foi_year_index[i]];
  }
  real  Mu_age[num_age_group_FOIs] = append_array(Mu_age_raw, {1.0});
  vector[number_ages] modelled_age_modifier;
  for (i in 1:number_ages) {
    modelled_age_modifier[i] = Mu_age[age_group_index[i]];
  }
  for(i in 1:N) {
    true_seropositive[i] = prob_infected_time_model_single_age(
      X[i, 1], 
      X[i, 2],
      modelled_lambda,
      modelled_age_modifier,
      alpha);
    test_seropositive[i] = ((1 - zika_cross) * sens * true_seropositive[i]) + ((1 - spec) * (1 - true_seropositive[i]));
    test_seropositive_zika_adjusted[i] = (sens * (true_seropositive[i] / (1 - zika_cross))) + ((1 - spec) * (1 - (true_seropositive[i] / (1 - zika_cross))));
  }
}

model {
  // Priors
  lambda ~ normal(0, 5); 
  alpha ~ normal(alpha_mean, alpha_sd);
  for(k in 1:num_age_group_FOIs) {
    Mu_age[k] ~ normal(1, 1);
  }
  for(j in 1:N) {
    if (X[j, 2] >= zika_year) {
      y[j] ~ binomial(tot[j], test_seropositive_zika_adjusted[j]);
    } else {
      y[j] ~ binomial(tot[j], test_seropositive[j]);
    }
  }
}

generated quantities {
 // posterior predictive check 
 real y_draw[N];
 real log_lik[N];
 for(q in 1:N) {
    if (X[q, 2] >= zika_year) {
       y_draw[q] = binomial_rng(tot[q], test_seropositive_zika_adjusted[q]);
    } else {
       y_draw[q] = binomial_rng(tot[q], test_seropositive[q]);
    }
  }
  // log likelihood for model comparison in loo package
  for (n in 1:N) {
    if (X[n, 2] >= zika_year) {
      log_lik[n] = binomial_lpmf(y[n] | tot[n], test_seropositive_zika_adjusted[n]);
    } else {
      log_lik[n] = binomial_lpmf(y[n] | tot[n], test_seropositive[n]);
    }
  }
}

