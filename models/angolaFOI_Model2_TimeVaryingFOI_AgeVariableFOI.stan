
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
}

parameters {
  real<lower=0.0001, upper = 1> lambda[num_unique_foi]; // number of yearly FOIs that can be estimated from the data
  real<lower=0> Mu_age_raw[num_age_group_FOIs - 1];         // scaler that determines age-specific fraction of yearly FOI
}

transformed parameters {
  real <lower = 0, upper = 1> true_seropositive[N];      // proportion truly seropositive for each age-year group
  real <lower = 0, upper = 1> test_seropositive[N];      // proportion testing seropositive for each age-year group
  real <lower = 0, upper = 1> test_seropositive_zika_adjusted[N]; // proportion testing seropositive for each age-year group post-zika

  // real<lower=0> modelled_lambda[years];   // mapping each year to a particular FOI
  // for (i in 1:years) {
  //   modelled_lambda[i] = lambda[foi_year_index[i]];
  // }
  real<lower=0> modelled_lambda[years] = lambda[foi_year_index];
  real  Mu_age[num_age_group_FOIs] = append_array(Mu_age_raw, {1.0});
  real<lower=0> modelled_age_modifier[number_ages] = Mu_age[age_group_index];
  // real<lower=0> modelled_age_modifier[number_ages];
  // for (i in 1:number_ages) {
  //   modelled_age_modifier[i] = Mu_age[age_group_index[i]];
  // }
  
  for(i in 1:N) {
    true_seropositive[i] = 1 - exp(-sum(to_vector(modelled_lambda[X[i,1]:X[i,2]]) .* to_vector(modelled_age_modifier[1:X[i,3]+1])));
    test_seropositive[i] = (sens * true_seropositive[i]) + ((1 - spec) * (1 - true_seropositive[i]));
    test_seropositive_zika_adjusted[i] = (sens * (true_seropositive[i] / (1 - zika_cross))) + ((1 - spec) * (1 - (true_seropositive[i] / (1 - zika_cross))));
  }
  
}

model {
  // Priors
  lambda ~ normal(0, 0.5); 
  for(k in 1:num_age_group_FOIs) {
    Mu_age[k] ~ normal(1, 0.5);
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
/// posterior predictive check 
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

