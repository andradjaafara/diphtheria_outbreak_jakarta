functions {

  // this is the model functions
  // includes x_r and x_i variables as per https://jrmihalj.github.io/estimating-transmission-by-fitting-mechanistic-models-in-Stan/

  real[] diphModel(real t,
                   real[] ymod,
                   real[] params,
                   real[] x_r,
                   int[] x_i){
  
    real dydt[12];
    real t1=344;
	real t2=359;
	real t3=372;
	
	if (t > t3 && t <= 9999){
		dydt[1] = -params[1] * params[19] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[1]/params[23] + params[17] * params[14] * ymod[3] - params[10];
		dydt[2] = params[1] * params[19] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[1]/params[23] - params[11] * ymod[2] - params[9] * ymod[2];
		dydt[3] = params[16] * params[11] * ymod[2] - params[6] * params[12] * ymod[3] - (1 - params[6]) * params[14] * ymod[3];
		dydt[4] = (1 - params[16]) * params[11] * ymod[2] + params[11] * ymod[9] - params[14] * ymod[4];
		dydt[5] = (1 - params[6] - params[17]) * params[14] * ymod[3] + params[14] * ymod[4] - params[14] * ymod[5];
		dydt[6] = (1 - params[20]) * params[14] * ymod[5] - params[15] * ymod[6];
		dydt[7] = params[6] * params[12] * ymod[3] + params[20] * params[14] * ymod[5] + params[9] * (ymod[2] + ymod[9]);
		dydt[8] = params[10] - params[1] * params[19] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[8]/params[23];
		dydt[9] = params[1] * params[19] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[8]/params[23] - params[11] * ymod[9] - params[9] * ymod[9];
		dydt[10] = params[16] * params[11] * ymod[2];
		dydt[11] = (1 - params[16]) * params[11] * ymod[2] + (1 - params[6] - params[17]) * params[14] * ymod[3] + params[14] * ymod[4] + (1 - params[20]) * params[14] * ymod[5];
		dydt[12] = ymod[1]/(ymod[1]+ymod[8]+ymod[7])*params[4]*params[1]*params[19]/(params[11]+params[9]) + ymod[8]/(ymod[1]+ymod[8]+ymod[7])*params[4]*params[1]*params[19]/(params[11]+params[9]) + params[16]*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[9]))*params[1]*params[19]/(params[6]*params[12]+(1-params[6])*params[14]) + ((1-params[16])*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[9]))+(ymod[8]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[9])))*(params[2]*params[1]*params[19]/(params[14])) + (params[16]*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[9]))*((1-params[6]-params[17])*params[14]/(params[6]*params[12]+(1-params[6])*params[14]))+((1-params[16])*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[9]))+(ymod[8]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[9])))*(params[14]/(params[14])))*(params[3]*params[1]*params[19]/(params[14])) + (params[16]*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[9]))*((1-params[6]-params[17])*params[14]/(params[6]*params[12]+(1-params[6])*params[14]))+((1-params[16])*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[9]))+(ymod[8]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[9])))*(params[14]/(params[14])))*(params[3]*params[1]*params[19]/(params[15]))*(1-params[20]);
	} else if (t > t2 && t <= t3){
		dydt[1] = -params[1] * params[7] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[1]/params[23] + params[17] * params[14] * ymod[3] - params[10];
		dydt[2] = params[1] * params[7] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[1]/params[23] - params[11] * ymod[2] - params[8] * ymod[2];
		dydt[3] = params[16] * params[11] * ymod[2] - params[6] * params[13] * ymod[3] - (1 - params[6]) * params[14] * ymod[3];
		dydt[4] = (1 - params[16]) * params[11] * ymod[2] + params[11] * ymod[9] - params[14] * ymod[4];
		dydt[5] = (1 - params[6] - params[17]) * params[14] * ymod[3] + params[14] * ymod[4] - params[14] * ymod[5];
		dydt[6] = (1 - params[20]) * params[14] * ymod[5] - params[15] * ymod[6];
		dydt[7] = params[6] * params[13] * ymod[3] + params[20] * params[14] * ymod[5] + params[8] * (ymod[2] + ymod[9]);
		dydt[8] = params[10] - params[1] * params[7] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[8]/params[23];
		dydt[9] = params[1] * params[7] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[8]/params[23] - params[11] * ymod[9] - params[8] * ymod[9];
		dydt[10] = params[16] * params[11] * ymod[2];
		dydt[11] = (1 - params[16]) * params[11] * ymod[2] + (1 - params[6] - params[17]) * params[14] * ymod[3] + params[14] * ymod[4] + (1 - params[20]) * params[14] * ymod[5];
		dydt[12] = ymod[1]/(ymod[1]+ymod[8]+ymod[7])*params[4]*params[1]*params[7]/(params[11]+params[8]) + ymod[8]/(ymod[1]+ymod[8]+ymod[7])*params[4]*params[1]*params[7]/(params[11]+params[8]) + params[16]*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8]))*params[1]*params[7]/(params[6]*params[13]+(1-params[6])*params[14]) + ((1-params[16])*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8]))+(ymod[8]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8])))*(params[2]*params[1]*params[7]/(params[14])) + (params[16]*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8]))*((1-params[6]-params[17])*params[14]/(params[6]*params[13]+(1-params[6])*params[14]))+((1-params[16])*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8]))+(ymod[8]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8])))*(params[14]/(params[14])))*(params[3]*params[1]*params[7]/(params[14])) + (params[16]*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8]))*((1-params[6]-params[17])*params[14]/(params[6]*params[13]+(1-params[6])*params[14]))+((1-params[16])*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8]))+(ymod[8]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8])))*(params[14]/(params[14])))*(params[3]*params[1]*params[7]/(params[15]))*(1-params[20]);
	} else if (t > t1 && t <= t2){
		dydt[1] = -params[1] * params[7] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[1]/params[23] + params[17] * params[14] * ymod[3] - params[18];
		dydt[2] = params[1] * params[7] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[1]/params[23] - params[11] * ymod[2] - params[8] * ymod[2];
		dydt[3] = params[16] * params[11] * ymod[2] - params[6] * params[13] * ymod[3] - (1 - params[6]) * params[14] * ymod[3];
		dydt[4] = (1 - params[16]) * params[11] * ymod[2] + params[11] * ymod[9] - params[14] * ymod[4];
		dydt[5] = (1 - params[6] - params[17]) * params[14] * ymod[3] + params[14] * ymod[4] - params[14] * ymod[5];
		dydt[6] = (1 - params[20]) * params[14] * ymod[5] - params[15] * ymod[6];
		dydt[7] = params[6] * params[13] * ymod[3] + params[20] * params[14] * ymod[5] + params[8] * (ymod[2] + ymod[9]);
		dydt[8] = params[18] - params[1] * params[7] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[8]/params[23];
		dydt[9] = params[1] * params[7] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[8]/params[23] - params[11] * ymod[9] - params[8] * ymod[9];
		dydt[10] = params[16] * params[11] * ymod[2];
		dydt[11] = (1 - params[16]) * params[11] * ymod[2] + (1 - params[6] - params[17]) * params[14] * ymod[3] + params[14] * ymod[4] + (1 - params[20]) * params[14] * ymod[5];
		dydt[12] = ymod[1]/(ymod[1]+ymod[8]+ymod[7])*params[4]*params[1]*params[7]/(params[11]+params[8]) + ymod[8]/(ymod[1]+ymod[8]+ymod[7])*params[4]*params[1]*params[7]/(params[11]+params[8]) + params[16]*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8]))*params[1]*params[7]/(params[6]*params[13]+(1-params[6])*params[14]) + ((1-params[16])*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8]))+(ymod[8]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8])))*(params[2]*params[1]*params[7]/(params[14])) + (params[16]*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8]))*((1-params[6]-params[17])*params[14]/(params[6]*params[13]+(1-params[6])*params[14]))+((1-params[16])*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8]))+(ymod[8]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8])))*(params[14]/(params[14])))*(params[3]*params[1]*params[7]/(params[14])) + (params[16]*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8]))*((1-params[6]-params[17])*params[14]/(params[6]*params[13]+(1-params[6])*params[14]))+((1-params[16])*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8]))+(ymod[8]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[8])))*(params[14]/(params[14])))*(params[3]*params[1]*params[7]/(params[15]))*(1-params[20]);
	} else{
		dydt[1] = -params[1] * params[19] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[1]/params[23] + params[17] * params[14] * ymod[3] - params[18];
		dydt[2] = params[1] * params[19] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[1]/params[23] - params[11] * ymod[2] - params[21] * ymod[2];
		dydt[3] = params[16] * params[11] * ymod[2] - params[5] * params[12] * ymod[3] - (1 - params[5]) * params[14] * ymod[3];
		dydt[4] = (1 - params[16]) * params[11] * ymod[2] + params[11] * ymod[9] - params[14] * ymod[4] - params[22] * ymod[4];
		dydt[5] = (1 - params[5] - params[17]) * params[14] * ymod[3] + params[14] * ymod[4] - params[14] * ymod[5];
		dydt[6] = (1 - params[20]) * params[14] * ymod[5] - params[15] * ymod[6];
		dydt[7] = params[5] * params[12] * ymod[3] + params[20] * params[14] * ymod[5] + params[21] * (ymod[2] + ymod[9]) + params[22] * ymod[4];
		dydt[8] = params[18] - params[1] * params[19] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[8]/params[23];
		dydt[9] = params[1] * params[19] * (ymod[3] + params[2] * ymod[4] + params[3] * (ymod[5] + ymod[6]) + params[4] * (ymod[2] + ymod[9])) * ymod[8]/params[23] - params[11] * ymod[9] - params[21] * ymod[9];
		dydt[10] = params[16] * params[11] * ymod[2];
		dydt[11] = (1 - params[16]) * params[11] * ymod[2] + (1 - params[5] - params[17]) * params[14] * ymod[3] + params[14] * ymod[4] + (1 - params[20]) * params[14] * ymod[5];
		dydt[12] = ymod[1]/(ymod[1]+ymod[8]+ymod[7])*params[4]*params[1]*params[19]/(params[11]+params[21]) + ymod[8]/(ymod[1]+ymod[8]+ymod[7])*params[4]*params[1]*params[19]/(params[11]+params[21]) + params[16]*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[21]))*params[1]*params[19]/(params[5]*params[12]+(1-params[5])*params[14]) + ((1-params[16])*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[21]))+(ymod[8]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[21])))*(params[2]*params[1]*params[19]/(params[14]+params[22])) + (params[16]*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[21]))*((1-params[5]-params[17])*params[14]/(params[5]*params[12]+(1-params[5])*params[14]))+((1-params[16])*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[21]))+(ymod[8]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[21])))*(params[14]/(params[14]+params[22])))*(params[3]*params[1]*params[19]/(params[14])) + (params[16]*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[21]))*((1-params[5]-params[17])*params[14]/(params[5]*params[12]+(1-params[5])*params[14]))+((1-params[16])*(ymod[1]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[21]))+(ymod[8]/(ymod[1]+ymod[8]+ymod[7]))*(params[11]/(params[11]+params[21])))*(params[14]/(params[14]+params[22])))*(params[3]*params[1]*params[19]/(params[15]))*(1-params[20]);
	}
    return dydt;
  }
}

data {
  int<lower = 1> n_obs; // number of days sampled
  int<lower = 1> n_params; // number of model parameters
  int<lower = 1> n_difeq; // number of differential equations in the system
  int<lower = 1> n_fake; // to generate "predicted"/"unsampled" data
  int<lower = 1> n_samples; // number of data samples
  int<lower = 1> N; // total number of population
  int y[n_obs]; // the poisson distributed data
  real t0; // initial time point (zero)
  real ts[n_obs]; // time points that were sampled
  real fake_ts[n_fake]; // time points for "predicted"/"unsampled" data
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real<lower = 0, upper = 1> beta; // infection rate
  real<lower = 0, upper = 1> theta1; // infectiousness of primary carriers relative to symptomatic
  real<lower = 0, upper = 1> theta2; // infectiousness of secondary/chronic carriers relative to symptomatic
  real<lower = 0, upper = 1> theta3; // infectiousness of exposed/prodormal infected relative to symptomatic
  real<lower = 0, upper = 1> p1; // reporting rate at the beginning of the outbreak
  real<lower = 0, upper = 1> p2; // reporting rate after declaration
  real<lower = 0, upper = 1> kappa; // proportion of the reduced contact after declaration + during school break
  real<lower = 0, upper = 1> sigmaE_1; // removal of exposed/prodormal as the effect of contact tracing during peak
  real<lower = 0, upper = 1> sigmaE_2; // removal of exposed/prodormal as the effect of contact tracing during decline
  real<lower = 0, upper = 1> S0; // initial proportion of susceptible hosts
}

transformed parameters {
  real y0[n_difeq]; // initial conditions for all model states
  real inci[n_samples];
  real R0;
  vector[126] inci2;
  vector[126] inci3;
  vector[44] inci4;
  vector[82] inci5;
  vector[126] inci6;
  vector[n_fake] inci1;
  real inci_samp[n_samples];
  real y_hat[n_fake, n_difeq]; // output from the ODE solver
  real<lower = 0> params[n_params];

  params[1] = beta;
  params[2] = theta1;
  params[3] = theta2;
  params[4] = theta3;
  params[5] = p1;
  params[6] = p2;
  params[7] = kappa;
  params[8] = sigmaE_1;
  params[9] = sigmaE_2;
  params[10] = 40000 * S0 * 1.974; //rho
  params[11] = 1.0/3.0; //tau
  params[12] = 1/(3.88+2); //gammaI_1
  params[13] = 1/(1.12+2); //gammaI_2
  params[14] = 1.0/18.0; //gammaC1
  params[15] = 1.0/18.0; //gammaC3
  params[16] = 0.70; //delta
  params[17] = 0.05; //eta
  params[18] = 0; //rho1
  params[19] = 1; //kappa1
  params[20] = 0.95; //epsilon
  params[21] = 0; //sigmaE
  params[22] = 0; //sigmaC
  params[23] = N; //total population
  
  R0 = theta3*beta/params[11] + params[16]*beta/(p1*params[12]+(1-p1)*params[14]) + (1-params[16])*theta1*beta/params[14] + (params[16]*(1-p1-params[17])*params[14]/(p1*params[12]+(1-p1)*params[14])+(1-params[16]))*theta2*beta/params[14] + (params[16]*(1-p1-params[17])*params[14]/(p1*params[12]+(1-p1)*params[14])+(1-params[16]))*theta2*beta/params[15]*(1-params[20]);
  
  y0[1] = N*S0-1;
  y0[2] = 0;
  y0[3] = 1;
  y0[4] = 0;
  y0[5] = 0;
  y0[6] = 0;
  y0[7] = 0;
  y0[8] = N-N*S0-1;
  y0[9] = 0;
  y0[10] = 0;
  y0[11] = 0;
  y0[12] = 1;
  
  y_hat = integrate_ode_rk45(diphModel, y0, t0, fake_ts, params, x_r, x_i);
  
  inci1 = to_vector(y_hat[,10]);
  inci2 = append_row(0.0,inci1[1:125]);
  inci3 = inci1 - inci2;
  inci4 = inci3[1:44]*p1;
  inci5 = inci3[45:126]*p2;
  inci6 = append_row(inci4,inci5);
  
  for (j in 1:n_samples){
	inci[j] = inci6[7*(j-1)+1]+inci6[7*(j-1)+2]+inci6[7*(j-1)+3]+inci6[7*(j-1)+4]+inci6[7*(j-1)+5]+inci6[7*(j-1)+6]+inci6[7*(j-1)+7];
  }
  
  //for(m in 1:n_samples)
  //  inci_samp[m] = inci[m];
  inci_samp = inci;
}

model {
  R0 ~ gamma(37,14);
  S0 ~ beta(21.66,143.12);
  theta1 ~ beta(15.0/2.0,46.0/2.0);
  theta2 ~ beta(1,1);
  theta3 ~ beta(1,1);
  kappa ~ beta(1,1);
  p1 ~ beta(15,15);
  p2 ~ beta(27,6);
  sigmaE_1 ~ beta(1,1);
  sigmaE_2 ~ beta(1,1);
  
  y ~ poisson(inci_samp);
  
  //for (t in 1:n_samples)
  //	y[t] ~ poisson(inci_samp[t]);
}

generated quantities {
  // generate log_likelihood;
  vector[n_obs] log_lik;
  for (t in 1:n_samples) {
    log_lik[t] = poisson_lpmf(y[t] | inci_samp[t]);
  }
}
