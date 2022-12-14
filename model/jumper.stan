data{
  int yrs; //number of years to make estimates
  int years[yrs];//year indexes for years to make estimates
  int a_yrs;//number of years of observed august abundance estimates
  int a_years[a_yrs];  //year indexes for observed august abundance estimates
  int s_yrs;//number of years of observed september abundance estimates
  int s_years[s_yrs];//year indexes for observed september abundance estimates
  int obs_yrs;//number of years with both observed august and september abundance estimates
  int obs_years[obs_yrs];//indexes for years with both observed august and september abundance 
  vector<lower=0>[a_yrs] N_aug_obs;//observed august abundance point estimate (median of pdf)
  vector<lower=0>[a_yrs] N_aug_sd_obs;//observed august abundance estimate lognormal sd
  vector<lower=0>[s_yrs] N_sept_obs;//observed september abundance estimate (median of pdf)
  vector<lower=0>[s_yrs] N_sept_sd_obs;//observed september abundance estimate lognormal sd
  int tr[obs_yrs];//observed trap count between august and september estimates
}
parameters{
  vector[yrs-1] z_eps_sept_obs;//scaled process error in september abundance
  real<lower=0> N_sept_sd;//september abundance process error sd
  real<lower=0> N_sept_1;//september abundance in first year
  vector[yrs] z_eps_tr;//scaled random effect for proportion trapped between august and september
  vector[yrs-1] z_eps_aug;//scaled random effect for august abundance as proportion of september
  real<lower=0> p_tr_sd;//process error sd for proportion trapped between august and september
  real<lower=0> p_aug_sd;//process error sd for august abundance as proportion of september
  real mu_tr;//logit mean august abundance as proportion of september
  real p_aug_1;//logit mean proportion trapped between august and september
}
transformed parameters{
  vector<lower=0>[yrs] N_sept;//state estimate of september abundance
  vector<lower=0>[yrs] N_aug;//state estimate of august abundnace
  vector<lower=0>[yrs] inc;//difference between august and september state estimate of abundance
  vector<lower=0,upper=1>[yrs] p_aug;//annual august proportion of september abundance
  vector<lower=0,upper=1>[yrs] p_tr;//annual proportion trapped between august and september
  N_sept[1] = N_sept_1;
  N_sept[2:yrs] = exp(log(N_sept_1) + cumulative_sum(z_eps_sept_obs * N_sept_sd));
  p_aug[1] = p_aug_1;
  p_aug[2:yrs] = inv_logit(logit(p_aug_1) + cumulative_sum(z_eps_aug * p_aug_sd));
  p_tr = inv_logit(mu_tr + z_eps_tr * p_tr_sd);
  N_aug = N_sept .* p_aug;
  inc = N_sept - N_aug;
}
model{
  //liklihood for binomial jumper ratio	
  for(i in 1:obs_yrs){
    tr[i] ~ poisson(p_tr[obs_years[i]] * inc[obs_years[i]]);
  }
  //likelihood for observed august abundance
  for(i in 1 : a_yrs){
    N_aug_obs[i] ~ lognormal(log(N_aug[a_years[i]]), N_aug_sd_obs[i]);
  }
  //likelihood for observed september abundance
  for(i in 1 : s_yrs){
    N_sept_obs[i]  ~ lognormal(log(N_sept[s_years[i]]), N_sept_sd_obs[i]);
  }
  //Priors
  z_eps_sept_obs ~ std_normal();
  N_sept_sd ~ std_normal();
  N_sept_1 ~ lognormal(0,6);
  z_eps_aug ~ std_normal();
  z_eps_tr ~ std_normal();
  p_tr_sd ~ cauchy(0,2.5);
  p_aug_sd ~ std_normal();
  mu_tr ~ normal(0,3);
  p_aug_1 ~ beta(1,1);
}
generated quantities{
  real hier_z_eps_tr;
  real hier_p_tr;
  hier_z_eps_tr = normal_rng(0,1);
  hier_p_tr = inv_logit(mu_tr + hier_z_eps_tr * p_tr_sd);
}