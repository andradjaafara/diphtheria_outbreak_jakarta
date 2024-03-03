### DIPHTHERIA MCMC FUNCTIONS ###
# SET OF FUNCTIONS FOR JAKARTA DIPHTHERIA OUTBREAK MODELLING BASED ON: http://sbfnk.github.io/mfiidd/introduction.html

### CREATE FITMODEL FUNCTION FOR DIPHTHERIA OUTBREAK MODELLING

## create a simple deterministic model model with constant population size

model_name <- "Modelling diphtheria outbreak in Jakarta"
model_state.names <- c("S","E","I","C1","C2","C3","R","V","Ev","In","As","Rt")
model_theta.names <- c("beta","S0","theta1","theta2","theta3","p1","p2","kappa","sigmaE_1","sigmaE_2")

model_simulateDeterministic <- function(theta,init.state,times) {
  
  model_ode <- function(time, state, parameters) {
    
    ## parameters
    beta = parameters[["beta"]]
    S0 = parameters[["S0"]]
    theta1 = parameters[["theta1"]]
    theta2 = parameters[["theta2"]]
    theta3 = parameters[["theta3"]]
    p1 = parameters[["p1"]]
    p2 = parameters[["p2"]]
    kappa = parameters[["kappa"]]
    sigmaE_1 = parameters[["sigmaE_1"]]
    sigmaE_2 = parameters[["sigmaE_2"]]
    rho = 40000 * S0 * 1.974 # 40000 roughly average per day vaccinated; times susceptible times 1.974531 ratio S in children to overall S
    tau = 1/3
    gammaI_1 = 1/(3.88+2) # added two days for clearance process from medication
    gammaI_2 = 1/(1.12+2) # added two days for clearance process from medication
    gammaC1 = 1/18
    gammaC3 = 1/18
    delta = 0.7 # alter between 5% and 70%
    eta = 0.05
    rho1 = 0
    kappa1 = 1
    epsilon = 0.95 # changed from previously 80%
    sigmaE = 0
    sigmaC = 0
    
    ## states
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    C1 <- state[["C1"]]
    C2 <- state[["C2"]]
    C3 <- state[["C3"]]
    R <- state[["R"]]
    V <- state[["V"]]
    Ev <- state[["Ev"]]
    In <- state[["In"]]
    As <- state[["As"]]
    Rt <- state[["Rt"]]
    
    N <- S + E + I + C1 + C2 + C3 + R + V
    
    if (time > 372 & time <= 9999) {
      dS  = -beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho
      dE  = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_2 * E
      dI  = delta * tau * E - p2 * gammaI_1 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_1 * I + epsilon * gammaC1 * C2 + sigmaE_2 * (E + Ev)
      dV  = rho - beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_2 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
    } else if (time > 359 & time <= 372) {
      dS  = -beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho
      dE  = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_1 * E
      dI  = delta * tau * E - p2 * gammaI_2 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_2 * I + epsilon * gammaC1 * C2 + sigmaE_1 * (E + Ev)
      dV  = rho - beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_1 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + V/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*beta*kappa/(p2*gammaI_2+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(theta1*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC3))*(1-epsilon)
    } else if (time > 344 & time <= 359) {
      dS  = -beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho1
      dE  = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_1 * E
      dI  = delta * tau * E - p2 * gammaI_2 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_2 * I + epsilon * gammaC1 * C2 + sigmaE_1 * (E + Ev)
      dV  = rho1 - beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_1 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + V/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*beta*kappa/(p2*gammaI_2+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(theta1*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC3))*(1-epsilon)
    } else {
      dS  = -beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho1
      dE  = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE * E
      dI  = delta * tau * E - p1 * gammaI_1 * I - (1 - p1) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1 - sigmaC * C1
      dC2 = (1 - p1 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p1 * gammaI_1 * I + epsilon * gammaC1 * C2 + sigmaE * (E + Ev) + sigmaC * C1
      dV  = rho1 - beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p1 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE) + V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE) + delta*(S/(S+V+R))*(tau/(tau+sigmaE))*beta*kappa1/(p1*gammaI_1+(1-p1)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(theta1*beta*kappa1/(gammaC1+sigmaC)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE))*((1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(gammaC1/(gammaC1+sigmaC)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE))*((1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(gammaC1/(gammaC1+sigmaC)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
    }
    
    return(list(c(dS, dE, dI, dC1, dC2, dC3, dR, dV, dEv, dIn, dAs, dRt)))
  }
  
  trajectory <- data.frame(ode(y = init.state,
                               times = times,
                               func = model_ode,
                               parms = theta,
                               method = "lsoda"))
  
  return(trajectory)
}

## function to compute log-prior
model_prior <- function(theta, log = FALSE) {
  
  beta = theta[["beta"]]
  S0 = theta[["S0"]]
  theta1 = theta[["theta1"]]
  theta2 = theta[["theta2"]]
  theta3 = theta[["theta3"]]
  #rho = theta[["rho"]]
  p1 = theta[["p1"]]
  p2 = theta[["p2"]]
  kappa = theta[["kappa"]]
  sigmaE_1 = theta[["sigmaE_1"]]
  sigmaE_2 = theta[["sigmaE_2"]]
  tau = 1/3
  gammaI_1 = 1/(3.88+2)
  gammaI_2 = 1/(1.12+2)
  gammaC1 = 1/18
  gammaC3 = 1/18
  delta = 0.7 # alter between 5% and 70%
  eta = 0.05
  rho1 = 0
  kappa1 = 1
  epsilon = 0.95
  sigmaE = 0
  sigmaC = 0
  
  R0 = theta3*beta/tau + delta*beta/(p1*gammaI_1+(1-p1)*gammaC1) + (1-delta)*theta1*beta/gammaC1 +
    (delta*(1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1)+(1-delta))*theta2*beta/gammaC1 +
    (delta*(1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1)+(1-delta))*theta2*beta/gammaC3*(1-epsilon)
  
  ## normal prior on R0: N[6,1]
  log.prior.R0 <- dnorm(R0, mean = 4.5, sd = 3, log = TRUE)
  #log.prior.R0 <- dpois(R0, lambda = 4.5, log = TRUE)
  log.prior.S0 <- dbeta(theta[["S0"]], shape1 = 0.662, shape2 = 3.25516, log = TRUE)
  log.prior.theta1 <- dbeta(theta[["theta1"]], shape1=3, shape2=21, log = TRUE)
  log.prior.theta2 <- dbeta(theta[["theta2"]], shape1=1, shape2=1, log = TRUE)
  log.prior.theta3 <- dbeta(theta[["theta3"]], shape1=1, shape2=1, log = TRUE)
  #log.prior.rho <- dbeta(theta[["rho"]], shape1=1, shape2=1, log = TRUE)
  log.prior.kappa <- dbeta(theta[["kappa"]], shape1=1, shape2=1, log = TRUE)
  log.prior.p1 <- dbeta(theta[["p1"]], shape1=15, shape2=15, log = TRUE)
  log.prior.p2 <- dbeta(theta[["p2"]], shape1=27, shape2=6, log = TRUE)
  log.prior.sigmaE_1 <- dbeta(theta[["sigmaE_1"]], shape1=1, shape2=1, log = TRUE)
  log.prior.sigmaE_2 <- dbeta(theta[["sigmaE_2"]], shape1=1, shape2=1, log = TRUE)
  
  #log.sum <- log.prior.R0 + log.prior.S0 + log.prior.theta1 + log.prior.theta2 + log.prior.theta3 + log.prior.p1 + log.prior.p2 +
  #  log.prior.rho + log.prior.kappa + log.prior.sigmaE_1 + log.prior.sigmaE_2
  
  log.sum <- log.prior.R0 + log.prior.S0 + log.prior.theta1 + log.prior.theta2 + log.prior.theta3 + log.prior.p1 + log.prior.p2 +
    log.prior.kappa + log.prior.sigmaE_1 + log.prior.sigmaE_2
  
  return(ifelse(log, log.sum, exp(log.sum)))
}

## function to compute the likelihood of one data point
model_pointLike <- function(data.point, model.point, theta, log = FALSE){
  
  ## the prevalence is observed through a Poisson process
  return(dpois(x = data.point[["obs"]],
               lambda = model.point[["I"]],
               log = log))
}

## function to generate observation from a model simulation
model_genObsPoint <- function(model.point, theta){
  
  ## the prevalence is observed through a Poisson process
  obs.point <- rpois(n = 1, lambda = model.point[["I"]])
  
  return(c(obs = obs.point))
}

## create deterministic model fitmodel
diphmodel_more_params_final <- fitR::fitmodel(name = model_name,
                                              state.names = model_state.names,
                                              theta.names = model_theta.names,
                                              simulate = model_simulateDeterministic,
                                              dprior = model_prior,
                                              rPointObs = model_genObsPoint,
                                              dPointObs = model_pointLike)

### CREATE FUNCTIONS FOR TO CALCULATE LOG POSTERIOR MCMC
## calculate log-likelihood from model estimates compared to data; modified from original to accomodate weekly data
dTrajObs_mod <- function (fitmodel, theta, init.state, data, par2, log = FALSE) 
{
  
  tmin = 300
  tmax = max(times_of_observed)
  
  if (tmin>(min(times_of_observed))){
    cat("\n")
    cat("**************************************************************************\n")
    cat("**************************************************************************\n")
    cat("The time-of-introduction is _after_ the first cases appeared!",t0,min(times_of_observed)-time_binning,"\n")
    cat("**************************************************************************\n")
    cat("**************************************************************************\n")
    cat("\n")
  }
  tmin = min(t0,min(times_of_observed))
  
  times = seq(tmin,tmax)
  
  if ("S0" %in% names(theta)) {
    npop = 10400000
    init.state["I"]  = 1
    init.state["S"]  = npop * theta["S0"] - init.state["I"]
    init.state["E"]  = 0
    init.state["C1"] = 0
    init.state["C2"] = 0
    init.state["C3"] = 0
    init.state["R"]  = 0     
    init.state["V"]  = (npop - init.state["S"])
    init.state["Ev"] = 0
  }
  traj <- fitmodel$simulate(theta, init.state, times)
  
  data.diff = diff(traj$In) 
  
  diff.t = 344 - (par2[1] + 1)
  diff.t2 = diff.t + 1
  data.diff[1:diff.t] = data.diff[1:diff.t]*theta[["p1"]]
  data.diff[diff.t2:length(data.diff)] = data.diff[diff.t2:length(data.diff)]*theta[["p2"]]
  
  diff.t3 = 302 - (par2[1] + 1) + 1
  
  data.diff2 = data.diff[diff.t3:length(data.diff)]
  diff.matrix = matrix(data.diff2,ncol = 7,byrow = T)
  incidence_predicted = rowSums(diff.matrix)
  traj2 = data.frame(time= 1:length(incidence_predicted), I = incidence_predicted)
  
  dens <- 0
  for (i in 1:nrow(data)) {
    data.point <- unlist(data[i, ])
    model.point <- unlist(traj2[i, ])
    dens <- dens + fitmodel$dPointObs(data.point = data.point, 
                                      model.point = model.point, theta = theta, log = TRUE)
  }
  return(ifelse(log, dens, exp(dens)))
}

## calculate log-likelihood posterior including prior log-likelihood and data
my_dLogPosterior_mod <- function(fitmodel, theta, init.state, data, par2) {
  
  # calculate the `fitmodel` log-prior for parameter `theta`
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  
  # calculate the `fitmodel` log-likelihood for parameter `theta` and
  # initial state `init.state`, with respect to the data set `data`.
  log.likelihood <- dTrajObs_mod(fitmodel, theta, init.state, data, par2, log = TRUE)
  
  # return the logged posterior probability
  log.posterior <- log.prior + log.likelihood
  return(log.posterior)
}

## diphtheria log-likelihood posterior calculation; par2 input can be modified c(t0, p, px)
my_dLogPosterior_diphtheria_more_params_final <- function(theta, par2=c(300,0.5,0.8)) {
  
  npop = 10400000
  I_0  = 1
  S_0  = npop * 0.15 - I_0
  E_0  = 0
  C1_0 = 0
  C2_0 = 0
  C3_0 = 0
  R_0  = 0      
  V_0  = npop - S_0 - I_0
  Ev_0 = 0
  
  init = c(S=S_0, E=E_0, I=I_0, C1=C1_0, C2=C2_0, C3=C3_0, R=R_0, V=V_0, Ev=Ev_0, In=0, As=0, Rt=0)
  
  return(my_dLogPosterior_mod(fitmodel = diphmodel_more_params_final,
                              theta = theta,
                              init.state = init,
                              data = diphdata,
                              par2=par2))
  
}

process_model_parameters_deterministic <- function(param,param.post,times_of_observed,seed){
  set.seed(seed)
  sample.index.CI <- sample(1:3400,5000,replace=TRUE)
  param.samples <- param.post[sample.index.CI,1:11]
  init.state = c(S=0, E=0, I=0, C1=0, C2=0, C3=0, R=0, V=0, Ev=0, In=0, As=0, Rt=0)
  t0 = 300
  tmin = 300
  tmax = max(times_of_observed)
  
  xtick = c(319,335,349,366,380,397,411,425)
  xlabel = c("15-Nov","01-Dec","15-Dec","01-Jan","15-Jan","01-Feb","15-Feb","01-Mar")
  
  if (tmin>(min(times_of_observed))){
    cat("\n")
    cat("**************************************************************************\n")
    cat("**************************************************************************\n")
    cat("The time-of-introduction is _after_ the first cases appeared!",t0,min(times_of_observed)-time_binning,"\n")
    cat("**************************************************************************\n")
    cat("**************************************************************************\n")
    cat("\n")
  }
  tmin = min(t0,min(times_of_observed))
  
  times.output = seq(tmin,tmax)
  
  if ("S0" %in% names(param)) {
    npop = 10400000
    init.state["I"]  = 1
    init.state["S"]  = npop * param["S0"] - init.state["I"]
    init.state["E"]  = 0
    init.state["C1"] = 0
    init.state["C2"] = 0
    init.state["C3"] = 0
    init.state["R"]  = 0     
    init.state["V"]  = (npop - init.state["S"] - 1)
    init.state["Ev"] = 0
  }
  traj <- diphmodel_more_params_final$simulate(param, init.state, times.output)
  
  data.diff = diff(traj$In)[-(1:2)]
  R_eff = diff(traj$Rt)
  R_eff = R_eff[10:length(R_eff)]
  
  diff.t = 344 - (t0)
  diff.t2 = diff.t + 1
  data.diff[1:diff.t] = data.diff[1:diff.t]*as.numeric(param[6])
  data.diff[diff.t2:length(data.diff)] = data.diff[diff.t2:length(data.diff)]*as.numeric(param[7])
  
  diff.t3 = 1
  
  data.diff2 = data.diff[diff.t3:length(data.diff)]
  diff.matrix = matrix(data.diff2,ncol = 7,byrow = T)
  incidence_predicted = rowSums(diff.matrix)
  traj2 = data.frame(time= 1:length(incidence_predicted), I = incidence_predicted)
  
  CI.samples <- matrix(nrow=5000,ncol=length(incidence_predicted))
  CI.samples.total <- vector()
  CI.samples.all <- matrix(nrow=5000,ncol=length(incidence_predicted))
  CI.samples.total.all <- vector()
  Reff.samples <- matrix(nrow=5000,ncol=119)
  
  for (i in 1:5000){
    init.state = c(S=0, E=0, I=0, C1=0, C2=0, C3=0, R=0, V=0, Ev=0, In=0, As=0, Rt=0)
    param.CI <- param.samples[i,]
    # param.CI["rho"] = 0 # without ORI
    # param.CI[9] = 0; param.CI[10] = 0; param.CI[11] = 0; param.CI[12] = 0 # without contact tracing
    # param.CI["rho"] = (40000 * param.CI[["S0"]] * 1.974)*1.5 # ORI longer and wider
    if ("S0" %in% names(param.CI)) {
      npop = 10400000
      init.state["I"]  = 1
      init.state["S"]  = npop * param.CI[["S0"]] - init.state["I"]
      init.state["E"]  = 0
      init.state["C1"] = 0
      init.state["C2"] = 0
      init.state["C3"] = 0
      init.state["R"]  = 0     
      init.state["V"]  = (npop - init.state["S"] - 1)
      init.state["Ev"] = 0
    }
    traj <- diphmodel_more_params_final$simulate(param.CI, init.state, times.output)
    
    data.diff = diff(traj$In)[-(1:2)]
    R_eff = diff(traj$Rt)
    R_eff = R_eff[10:length(R_eff)]
    data.diff3 = data.diff
    
    diff.t = 344 - (t0)
    diff.t2 = diff.t + 1
    data.diff[1:diff.t] = data.diff[1:diff.t]*param.CI[[6]]
    data.diff[diff.t2:length(data.diff)] = data.diff[diff.t2:length(data.diff)]*param.CI[[7]]
    
    diff.t3 = 1
    
    data.diff2 = data.diff[diff.t3:length(data.diff)]
    data.diff4 = data.diff3[diff.t3:length(data.diff3)]
    diff.matrix = matrix(data.diff2,ncol = 7,byrow = T)
    diff.matrix2 = matrix(data.diff4,ncol = 7,byrow = T)
    
    incidence_predicted.CI = rowSums(diff.matrix)
    incidence_predicted.CI.all = rowSums(diff.matrix2)
    
    CI.samples[i,] <- incidence_predicted.CI
    CI.samples.total[i] <- sum(incidence_predicted.CI)
    CI.samples.all[i,] <- incidence_predicted.CI.all
    CI.samples.total.all[i] <- sum(incidence_predicted.CI.all)
    Reff.samples[i,] <- R_eff
  }
  
  simulation.summary = data.frame(x1=numeric(),
                                  x2=numeric(),
                                  x3=numeric(),
                                  x4=numeric(),
                                  x5=numeric(),
                                  x6=numeric(),
                                  x7=numeric())
  
  for (i in 1:length(incidence_predicted.CI)){
    dat = CI.samples[,i]
    median.dat = median(dat)
    mean.dat = mean(dat)
    q.25.dat = quantile(dat,0.25)
    q.75.dat = quantile(dat,0.75)
    q.2.5.dat = quantile(dat,0.025)
    q.97.5.dat = quantile(dat,0.975)
    summary = c(i,median.dat,mean.dat,q.25.dat,q.75.dat,q.2.5.dat,q.97.5.dat)
    simulation.summary = rbind(simulation.summary,summary)
  }
  
  simulation.summary.all = data.frame(x1=numeric(),
                                      x2=numeric(),
                                      x3=numeric(),
                                      x4=numeric(),
                                      x5=numeric(),
                                      x6=numeric(),
                                      x7=numeric())
  
  for (i in 1:length(incidence_predicted.CI.all)){
    dat = CI.samples.all[,i]
    median.dat = median(dat)
    mean.dat = mean(dat)
    q.25.dat = quantile(dat,0.25)
    q.75.dat = quantile(dat,0.75)
    q.2.5.dat = quantile(dat,0.025)
    q.97.5.dat = quantile(dat,0.975)
    summary.all = c(i,median.dat,mean.dat,q.25.dat,q.75.dat,q.2.5.dat,q.97.5.dat)
    simulation.summary.all = rbind(simulation.summary.all,summary.all)
  }
  
  
  Reff.summary = data.frame(x1=numeric(),
                            x2=numeric(),
                            x3=numeric(),
                            x4=numeric(),
                            x5=numeric(),
                            x6=numeric(),
                            x7=numeric())
  
  for (i in 1:119){
    dat = Reff.samples[,i]
    median.dat = median(dat)
    mean.dat = mean(dat)
    q.25.dat = quantile(dat,0.25)
    q.75.dat = quantile(dat,0.75)
    q.2.5.dat = quantile(dat,0.025)
    q.97.5.dat = quantile(dat,0.975)
    summary = c(i,median.dat,mean.dat,q.25.dat,q.75.dat,q.2.5.dat,q.97.5.dat)
    Reff.summary = rbind(Reff.summary,summary)
  }
  
  total.summary = data.frame(x1=numeric(),
                             x2=numeric(),
                             x3=numeric(),
                             x4=numeric(),
                             x5=numeric(),
                             x6=numeric())
  
  median.dat = median(CI.samples.total)
  mean.dat = mean(CI.samples.total)
  q.25.dat = quantile(CI.samples.total,0.25)
  q.75.dat = quantile(CI.samples.total,0.75)
  q.2.5.dat = quantile(CI.samples.total,0.025)
  q.97.5.dat = quantile(CI.samples.total,0.975)
  summary = c(median.dat,mean.dat,q.25.dat,q.75.dat,q.2.5.dat,q.97.5.dat)
  total.summary = rbind(total.summary,summary)
  
  total.summary.all = data.frame(x1=numeric(),
                                 x2=numeric(),
                                 x3=numeric(),
                                 x4=numeric(),
                                 x5=numeric(),
                                 x6=numeric())
  
  median.dat = median(CI.samples.total.all)
  mean.dat = mean(CI.samples.total.all)
  q.25.dat = quantile(CI.samples.total.all,0.25)
  q.75.dat = quantile(CI.samples.total.all,0.75)
  q.2.5.dat = quantile(CI.samples.total.all,0.025)
  q.97.5.dat = quantile(CI.samples.total.all,0.975)
  summary.all = c(median.dat,mean.dat,q.25.dat,q.75.dat,q.2.5.dat,q.97.5.dat)
  total.summary.all = rbind(total.summary.all,summary.all)
  
  colnames(simulation.summary) = c("time","median","mean","q25","q75","q2.5","q97.5")
  simulation.summary$time <- times_of_observed
  colnames(simulation.summary.all) = c("time","median","mean","q25","q75","q2.5","q97.5")
  simulation.summary.all$time <- times_of_observed
  colnames(Reff.summary) = c("time","median","mean","q25","q75","q2.5","q97.5")
  Reff.summary$time <- times.output[11:length(times.output)]
  colnames(total.summary) = c("median","mean","q25","q75","q2.5","q97.5")
  colnames(total.summary.all) = c("median","mean","q25","q75","q2.5","q97.5")
  
  return(list(simulation.summary,
              total.summary,
              simulation.summary.all,
              total.summary.all,
              Reff.summary))
  
  # write.csv(simulation.summary,simulation.summary.name,row.names = F)
  # write.csv(total.summary,total.summary.name,row.names = F)
  # write.csv(simulation.summary.all,simulation.summary.all.name,row.names = F)
  # write.csv(total.summary.all,total.summary.all.name,row.names = F)
  # write.csv(Reff.summary,Reff.summary.name,row.names = F)
}

process_model_deterministic_plot <- function(data,
                                             times_of_observed,
                                             incidence_observed){
  simulation.summary <- data[[1]]
  total.summary <- data[[2]]
  simulation.summary.all <- data[[3]]
  total.summary.all <- data[[4]]
  Reff.summary <- data[[5]]
  
  xtick = c(319,335,349,366,380,397,411,425)
  xlabel = c("15-Nov","01-Dec","15-Dec","01-Jan","15-Jan","01-Feb","15-Feb","01-Mar")
  
  tmin = 300
  tmax = max(times_of_observed)
  times.output = seq(tmin,tmax)
  
  par(mar=c(3, 3, 2, 3))
  par(mgp=c(5,1,0))
  ymax=80
  plot(times_of_observed
       ,incidence_observed
       ,cex.main=0.8
       ,ylim=c(0,1.05*ymax)
       # ,xlab="Time"
       # ,ylab="Incidence"
       ,cex=2,pch=16,xaxt="n")
  title(ylab="Incidence", line=2, cex.lab=1)
  title(xlab="Time", line=2, cex.lab=1)
  polygon(c(times_of_observed, rev(times_of_observed)), c(simulation.summary$q97.5, rev(simulation.summary$q2.5)),
          col = "rosybrown1", border = NA)
  polygon(c(times_of_observed, rev(times_of_observed)), c(simulation.summary$q75, rev(simulation.summary$q25)),
          col = "indianred1", border = NA)
  points(times_of_observed
         ,incidence_observed
         ,cex.main=0.8
         ,ylim=c(0,1.2*ymax)
         ,xlab="Time"
         ,ylab="Incidence"
         ,cex=2,pch=16)
  lines(times_of_observed
        ,simulation.summary$median
        ,col="red4"
        ,lwd=3)
  lines(times_of_observed
        ,simulation.summary$mean
        ,col="red4"
        ,lwd=3,lty=2)
  abline(v=344.5,lwd=3,lty=2,col="blue")
  text(346.5,1,"Start of catch up vaccinations program (ORI)",cex=0.75,adj=0)
  
  abline(v=365,lwd=3,lty=2,col="green")
  text(367,1.05*ymax,"Start of new year 2018",cex=0.75,adj=0)
  
  color_transparent <- adjustcolor('forestgreen', alpha.f = 0.3)
  polygon(c(357,372,372,357),c(-5,-5,100005,100005),col=color_transparent,border=NA)
  
  #abline(v=333,lwd=3,lty=2,col="orange")
  #text(335,90,"Declaration of diphtheria outbreak in Jakarta/Indonesia",cex=0.75,adj=0)
  
  legend(390,1.05*ymax,c("50% CI","95% CI","School Holiday"),fill=c("indianred1","rosybrown1",color_transparent),
         border=c("indianred1","rosybrown1",color_transparent),bty="n",cex=1)
  
  axis(1, at=xtick, labels=xlabel)
  
  par(new=TRUE)
  plot(times.output[11:length(times.output)],Reff.summary$median,lwd=4,xlab="",ylab="",ylim=c(0,4),
       axes=FALSE,type="l",col="orange",lty=2)
  mtext("Effective reproduction number (Rt)",side=4,col="black",line=1.7)
  mtext("Incidence",side=2,col="black",line=1.7)
  axis(4,ylim=c(0,4),col="black",col.axis="black",las=1)
  abline(h=1,col="orange",lwd=3)  
  
  p_output <- recordPlot()
  return(p_output)
}

#### From here are different models for different intervention scenarios:
# model_baseline
# model_without_ORI
# model_without_CT
# model_delayed_ORI
# model_more_ORI

# Baseline model

model_baseline <- function(theta,init.state,times) {
  
  model_ode <- function(time, state, parameters) {
    
    ## parameters
    beta = parameters[["beta"]]
    S0 = parameters[["S0"]]
    theta1 = parameters[["theta1"]]
    theta2 = parameters[["theta2"]]
    theta3 = parameters[["theta3"]]
    p1 = parameters[["p1"]]
    p2 = parameters[["p2"]]
    kappa = parameters[["kappa"]]
    sigmaE_1 = parameters[["sigmaE_1"]]
    sigmaE_2 = parameters[["sigmaE_2"]]
    rho = 40000 * S0 * 1.974 # 40000 roughly average per day vaccinated; times susceptible times 1.974531 ratio S in children to overall S
    tau = 1/3
    gammaI_1 = 1/(3.88+2) # added two days for clearance process from medication
    gammaI_2 = 1/(1.12+2) # added two days for clearance process from medication
    gammaC1 = 1/18
    gammaC3 = 1/18
    delta = 0.7 # alter between 5% and 70%
    eta = 0.05
    rho1 = 0
    kappa1 = 1
    epsilon = 0.95 # changed from previously 80%
    sigmaE = 0
    sigmaC = 0
    
    ## states
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    C1 <- state[["C1"]]
    C2 <- state[["C2"]]
    C3 <- state[["C3"]]
    R <- state[["R"]]
    V <- state[["V"]]
    Ev <- state[["Ev"]]
    In <- state[["In"]]
    As <- state[["As"]]
    Rt <- state[["Rt"]]
    
    N <- S + E + I + C1 + C2 + C3 + R + V
    
    if (time > 372 & time <= 9999) {
      dS  = -beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho
      dE  = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_2 * E
      dI  = delta * tau * E - p2 * gammaI_1 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_1 * I + epsilon * gammaC1 * C2 + sigmaE_2 * (E + Ev)
      dV  = rho - beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_2 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
    } else if (time > 359 & time <= 372) {
      dS  = -beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho
      dE  = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_1 * E
      dI  = delta * tau * E - p2 * gammaI_2 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_2 * I + epsilon * gammaC1 * C2 + sigmaE_1 * (E + Ev)
      dV  = rho - beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_1 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + V/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*beta*kappa/(p2*gammaI_2+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(theta1*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC3))*(1-epsilon)
    } else if (time > 344 & time <= 359) {
      dS  = -beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho1
      dE  = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_1 * E
      dI  = delta * tau * E - p2 * gammaI_2 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_2 * I + epsilon * gammaC1 * C2 + sigmaE_1 * (E + Ev)
      dV  = rho1 - beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_1 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + V/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*beta*kappa/(p2*gammaI_2+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(theta1*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC3))*(1-epsilon)
    } else {
      dS  = -beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho1
      dE  = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE * E
      dI  = delta * tau * E - p1 * gammaI_1 * I - (1 - p1) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1 - sigmaC * C1
      dC2 = (1 - p1 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p1 * gammaI_1 * I + epsilon * gammaC1 * C2 + sigmaE * (E + Ev) + sigmaC * C1
      dV  = rho1 - beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p1 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE) + V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE) + delta*(S/(S+V+R))*(tau/(tau+sigmaE))*beta*kappa1/(p1*gammaI_1+(1-p1)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(theta1*beta*kappa1/(gammaC1+sigmaC)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE))*((1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(gammaC1/(gammaC1+sigmaC)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE))*((1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(gammaC1/(gammaC1+sigmaC)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
    }
    
    return(list(c(dS, dE, dI, dC1, dC2, dC3, dR, dV, dEv, dIn, dAs, dRt)))
  }
  
  trajectory <- data.frame(ode(y = init.state,
                               times = times,
                               func = model_ode,
                               parms = theta,
                               method = "lsoda"))
  
  return(trajectory)
}

# Scenario 1: without ORI

model_without_ORI <- function(theta,init.state,times) {
  
  model_ode <- function(time, state, parameters) {
    
    ## parameters
    beta = parameters[["beta"]]
    S0 = parameters[["S0"]]
    theta1 = parameters[["theta1"]]
    theta2 = parameters[["theta2"]]
    theta3 = parameters[["theta3"]]
    p1 = parameters[["p1"]]
    p2 = parameters[["p2"]]
    kappa = parameters[["kappa"]]
    sigmaE_1 = parameters[["sigmaE_1"]]
    sigmaE_2 = parameters[["sigmaE_2"]]
    rho = 0
    tau = 1/3
    gammaI_1 = 1/(3.88+2) # added two days for clearance process from medication
    gammaI_2 = 1/(1.12+2) # added two days for clearance process from medication
    gammaC1 = 1/18
    gammaC3 = 1/18
    delta = 0.7 # alter between 5% and 70%
    eta = 0.05
    rho1 = 0
    kappa1 = 1
    epsilon = 0.95 # changed from previously 80%
    sigmaE = 0
    sigmaC = 0
    
    ## states
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    C1 <- state[["C1"]]
    C2 <- state[["C2"]]
    C3 <- state[["C3"]]
    R <- state[["R"]]
    V <- state[["V"]]
    Ev <- state[["Ev"]]
    In <- state[["In"]]
    As <- state[["As"]]
    Rt <- state[["Rt"]]
    
    N <- S + E + I + C1 + C2 + C3 + R + V
    
    if (time > 372 & time <= 9999) {
      dS  = -beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho
      dE  = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_2 * E
      dI  = delta * tau * E - p2 * gammaI_1 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_1 * I + epsilon * gammaC1 * C2 + sigmaE_2 * (E + Ev)
      dV  = rho - beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_2 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
    } else if (time > 359 & time <= 372) {
      dS  = -beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho
      dE  = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_1 * E
      dI  = delta * tau * E - p2 * gammaI_2 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_2 * I + epsilon * gammaC1 * C2 + sigmaE_1 * (E + Ev)
      dV  = rho - beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_1 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + V/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*beta*kappa/(p2*gammaI_2+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(theta1*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC3))*(1-epsilon)
    } else if (time > 344 & time <= 359) {
      dS  = -beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho1
      dE  = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_1 * E
      dI  = delta * tau * E - p2 * gammaI_2 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_2 * I + epsilon * gammaC1 * C2 + sigmaE_1 * (E + Ev)
      dV  = rho1 - beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_1 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + V/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*beta*kappa/(p2*gammaI_2+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(theta1*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC3))*(1-epsilon)
    } else {
      dS  = -beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho1
      dE  = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE * E
      dI  = delta * tau * E - p1 * gammaI_1 * I - (1 - p1) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1 - sigmaC * C1
      dC2 = (1 - p1 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p1 * gammaI_1 * I + epsilon * gammaC1 * C2 + sigmaE * (E + Ev) + sigmaC * C1
      dV  = rho1 - beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p1 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE) + V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE) + delta*(S/(S+V+R))*(tau/(tau+sigmaE))*beta*kappa1/(p1*gammaI_1+(1-p1)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(theta1*beta*kappa1/(gammaC1+sigmaC)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE))*((1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(gammaC1/(gammaC1+sigmaC)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE))*((1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(gammaC1/(gammaC1+sigmaC)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
    }
    
    return(list(c(dS, dE, dI, dC1, dC2, dC3, dR, dV, dEv, dIn, dAs, dRt)))
  }
  
  trajectory <- data.frame(ode(y = init.state,
                               times = times,
                               func = model_ode,
                               parms = theta,
                               method = "lsoda"))
  
  return(trajectory)
}

# Scenario 2: without contact tracing and treatment of contacts

model_without_CT <- function(theta,init.state,times) {
  
  model_ode <- function(time, state, parameters) {
    
    ## parameters
    beta = parameters[["beta"]]
    S0 = parameters[["S0"]]
    theta1 = parameters[["theta1"]]
    theta2 = parameters[["theta2"]]
    theta3 = parameters[["theta3"]]
    p1 = parameters[["p1"]]
    p2 = parameters[["p2"]]
    kappa = parameters[["kappa"]]
    sigmaE_1 = 0
    sigmaE_2 = 0
    rho = 40000 * S0 * 1.974 # 40000 roughly average per day vaccinated; times susceptible times 1.974531 ratio S in children to overall S
    tau = 1/3
    gammaI_1 = 1/(3.88+2) # added two days for clearance process from medication
    gammaI_2 = 1/(1.12+2) # added two days for clearance process from medication
    gammaC1 = 1/18
    gammaC3 = 1/18
    delta = 0.7 # alter between 5% and 70%
    eta = 0.05
    rho1 = 0
    kappa1 = 1
    epsilon = 0.95 # changed from previously 80%
    sigmaE = 0
    sigmaC = 0
    
    ## states
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    C1 <- state[["C1"]]
    C2 <- state[["C2"]]
    C3 <- state[["C3"]]
    R <- state[["R"]]
    V <- state[["V"]]
    Ev <- state[["Ev"]]
    In <- state[["In"]]
    As <- state[["As"]]
    Rt <- state[["Rt"]]
    
    N <- S + E + I + C1 + C2 + C3 + R + V
    
    if (time > 372 & time <= 9999) {
      dS  = -beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho
      dE  = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_2 * E
      dI  = delta * tau * E - p2 * gammaI_1 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_1 * I + epsilon * gammaC1 * C2 + sigmaE_2 * (E + Ev)
      dV  = rho - beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_2 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
    } else if (time > 359 & time <= 372) {
      dS  = -beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho
      dE  = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_1 * E
      dI  = delta * tau * E - p2 * gammaI_2 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_2 * I + epsilon * gammaC1 * C2 + sigmaE_1 * (E + Ev)
      dV  = rho - beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_1 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + V/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*beta*kappa/(p2*gammaI_2+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(theta1*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC3))*(1-epsilon)
    } else if (time > 344 & time <= 359) {
      dS  = -beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho1
      dE  = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_1 * E
      dI  = delta * tau * E - p2 * gammaI_2 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_2 * I + epsilon * gammaC1 * C2 + sigmaE_1 * (E + Ev)
      dV  = rho1 - beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_1 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + V/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*beta*kappa/(p2*gammaI_2+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(theta1*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC3))*(1-epsilon)
    } else {
      dS  = -beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho1
      dE  = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE * E
      dI  = delta * tau * E - p1 * gammaI_1 * I - (1 - p1) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1 - sigmaC * C1
      dC2 = (1 - p1 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p1 * gammaI_1 * I + epsilon * gammaC1 * C2 + sigmaE * (E + Ev) + sigmaC * C1
      dV  = rho1 - beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p1 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE) + V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE) + delta*(S/(S+V+R))*(tau/(tau+sigmaE))*beta*kappa1/(p1*gammaI_1+(1-p1)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(theta1*beta*kappa1/(gammaC1+sigmaC)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE))*((1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(gammaC1/(gammaC1+sigmaC)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE))*((1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(gammaC1/(gammaC1+sigmaC)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
    }
    
    return(list(c(dS, dE, dI, dC1, dC2, dC3, dR, dV, dEv, dIn, dAs, dRt)))
  }
  
  trajectory <- data.frame(ode(y = init.state,
                               times = times,
                               func = model_ode,
                               parms = theta,
                               method = "lsoda"))
  
  return(trajectory)
}

# Scenario 3: delayed vaccinations 2 weeks

model_delayed_ORI <- function(theta,init.state,times) {
  
  model_ode <- function(time, state, parameters) {
    
    ## parameters
    beta = parameters[["beta"]]
    S0 = parameters[["S0"]]
    theta1 = parameters[["theta1"]]
    theta2 = parameters[["theta2"]]
    theta3 = parameters[["theta3"]]
    p1 = parameters[["p1"]]
    p2 = parameters[["p2"]]
    kappa = parameters[["kappa"]]
    sigmaE_1 = parameters[["sigmaE_1"]]
    sigmaE_2 = parameters[["sigmaE_2"]]
    rho = 40000 * S0 * 1.974 # 40000 roughly average per day vaccinated; times susceptible times 1.974531 ratio S in children to overall S
    tau = 1/3
    gammaI_1 = 1/(3.88+2) # added two days for clearance process from medication
    gammaI_2 = 1/(1.12+2) # added two days for clearance process from medication
    gammaC1 = 1/18
    gammaC3 = 1/18
    delta = 0.7 # alter between 5% and 70%
    eta = 0.05
    rho1 = 0
    kappa1 = 1
    epsilon = 0.95 # changed from previously 80%
    sigmaE = 0
    sigmaC = 0
    
    ## states
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    C1 <- state[["C1"]]
    C2 <- state[["C2"]]
    C3 <- state[["C3"]]
    R <- state[["R"]]
    V <- state[["V"]]
    Ev <- state[["Ev"]]
    In <- state[["In"]]
    As <- state[["As"]]
    Rt <- state[["Rt"]]
    
    N <- S + E + I + C1 + C2 + C3 + R + V
    
    if (time > 372 & time <= 9999) {
      dS  = -beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho
      dE  = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_2 * E
      dI  = delta * tau * E - p2 * gammaI_1 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_1 * I + epsilon * gammaC1 * C2 + sigmaE_2 * (E + Ev)
      dV  = rho - beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_2 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
    } else if (time > 359 & time <= 372) {
      dS  = -beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho1
      dE  = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_1 * E
      dI  = delta * tau * E - p2 * gammaI_2 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_2 * I + epsilon * gammaC1 * C2 + sigmaE_1 * (E + Ev)
      dV  = rho1 - beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_1 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + V/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*beta*kappa/(p2*gammaI_2+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(theta1*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC3))*(1-epsilon)
    } else if (time > 344 & time <= 359) {
      dS  = -beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho1
      dE  = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_1 * E
      dI  = delta * tau * E - p2 * gammaI_2 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_2 * I + epsilon * gammaC1 * C2 + sigmaE_1 * (E + Ev)
      dV  = rho1 - beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_1 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + V/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*beta*kappa/(p2*gammaI_2+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(theta1*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC3))*(1-epsilon)
    } else {
      dS  = -beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho1
      dE  = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE * E
      dI  = delta * tau * E - p1 * gammaI_1 * I - (1 - p1) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1 - sigmaC * C1
      dC2 = (1 - p1 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p1 * gammaI_1 * I + epsilon * gammaC1 * C2 + sigmaE * (E + Ev) + sigmaC * C1
      dV  = rho1 - beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p1 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE) + V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE) + delta*(S/(S+V+R))*(tau/(tau+sigmaE))*beta*kappa1/(p1*gammaI_1+(1-p1)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(theta1*beta*kappa1/(gammaC1+sigmaC)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE))*((1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(gammaC1/(gammaC1+sigmaC)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE))*((1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(gammaC1/(gammaC1+sigmaC)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
    }
    
    return(list(c(dS, dE, dI, dC1, dC2, dC3, dR, dV, dEv, dIn, dAs, dRt)))
  }
  
  trajectory <- data.frame(ode(y = init.state,
                               times = times,
                               func = model_ode,
                               parms = theta,
                               method = "lsoda"))
  
  return(trajectory)
}

# Scenario 4: more ORI (twice coverage + longer)

model_more_ORI <- function(theta,init.state,times) {
  
  model_ode <- function(time, state, parameters) {
    
    ## parameters
    beta = parameters[["beta"]]
    S0 = parameters[["S0"]]
    theta1 = parameters[["theta1"]]
    theta2 = parameters[["theta2"]]
    theta3 = parameters[["theta3"]]
    p1 = parameters[["p1"]]
    p2 = parameters[["p2"]]
    kappa = parameters[["kappa"]]
    sigmaE_1 = parameters[["sigmaE_1"]]
    sigmaE_2 = parameters[["sigmaE_2"]]
    rho = 2 * 40000 * S0 * 1.974 # 80000 roughly average per day vaccinated (twice the normal); times susceptible times 1.974531 ratio S in children to overall S
    tau = 1/3
    gammaI_1 = 1/(3.88+2) # added two days for clearance process from medication
    gammaI_2 = 1/(1.12+2) # added two days for clearance process from medication
    gammaC1 = 1/18
    gammaC3 = 1/18
    delta = 0.7 # alter between 5% and 70%
    eta = 0.05
    rho1 = 0
    kappa1 = 1
    epsilon = 0.95 # changed from previously 80%
    sigmaE = 0
    sigmaC = 0
    
    ## states
    S <- state[["S"]]
    E <- state[["E"]]
    I <- state[["I"]]
    C1 <- state[["C1"]]
    C2 <- state[["C2"]]
    C3 <- state[["C3"]]
    R <- state[["R"]]
    V <- state[["V"]]
    Ev <- state[["Ev"]]
    In <- state[["In"]]
    As <- state[["As"]]
    Rt <- state[["Rt"]]
    
    N <- S + E + I + C1 + C2 + C3 + R + V
    
    if (time > 372 & time <= 9999) {
      dS  = -beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho
      dE  = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_2 * E
      dI  = delta * tau * E - p2 * gammaI_1 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_1 * I + epsilon * gammaC1 * C2 + sigmaE_2 * (E + Ev)
      dV  = rho - beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_2 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE_2) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*beta*kappa1/(p2*gammaI_1+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(theta1*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_2))*((1-p2-eta)*gammaC1/(p2*gammaI_1+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_2))+(V/(S+V+R))*(tau/(tau+sigmaE_2)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
    } else if (time > 359 & time <= 372) {
      dS  = -beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho
      dE  = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_1 * E
      dI  = delta * tau * E - p2 * gammaI_2 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_2 * I + epsilon * gammaC1 * C2 + sigmaE_1 * (E + Ev)
      dV  = rho - beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_1 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + V/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*beta*kappa/(p2*gammaI_2+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(theta1*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC3))*(1-epsilon)
    } else if (time > 344 & time <= 359) {
      dS  = -beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho1
      dE  = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE_1 * E
      dI  = delta * tau * E - p2 * gammaI_2 * I - (1 - p2) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1
      dC2 = (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p2 * gammaI_2 * I + epsilon * gammaC1 * C2 + sigmaE_1 * (E + Ev)
      dV  = rho1 - beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE_1 * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p2 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + V/(S+V+R)*theta3*beta*kappa/(tau+sigmaE_1) + delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*beta*kappa/(p2*gammaI_2+(1-p2)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(theta1*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE_1))*((1-p2-eta)*gammaC1/(p2*gammaI_2+(1-p2)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE_1))+(V/(S+V+R))*(tau/(tau+sigmaE_1)))*(gammaC1/(gammaC1)))*(theta2*beta*kappa/(gammaC3))*(1-epsilon)
    } else {
      dS  = -beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N + eta * gammaC1 * I - rho1
      dE  = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * S/N - tau * E - sigmaE * E
      dI  = delta * tau * E - p1 * gammaI_1 * I - (1 - p1) * gammaC1 * I
      dC1 = (1 - delta) * tau * E + tau * Ev - gammaC1 * C1 - sigmaC * C1
      dC2 = (1 - p1 - eta) * gammaC1 * I + gammaC1 * C1 - gammaC1 * C2
      dC3 = (1 - epsilon) * gammaC1 * C2 - gammaC3 * C3
      dR  = p1 * gammaI_1 * I + epsilon * gammaC1 * C2 + sigmaE * (E + Ev) + sigmaC * C1
      dV  = rho1 - beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N
      dEv = beta * kappa1 * (I + theta1 * C1 + theta2 * (C2 + C3) + theta3 * (E + Ev)) * V/N - tau * Ev - sigmaE * Ev
      dIn = delta * tau * E
      dAs = (1 - delta) * tau * E + (1 - p1 - eta) * gammaC1 * I + gammaC1 * C1 + (1 - epsilon) * gammaC1 * C2
      dRt  = S/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE) + V/(S+V+R)*theta3*beta*kappa1/(tau+sigmaE) + delta*(S/(S+V+R))*(tau/(tau+sigmaE))*beta*kappa1/(p1*gammaI_1+(1-p1)*gammaC1) + ((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(theta1*beta*kappa1/(gammaC1+sigmaC)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE))*((1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(gammaC1/(gammaC1+sigmaC)))*(theta2*beta*kappa1/(gammaC1)) + (delta*(S/(S+V+R))*(tau/(tau+sigmaE))*((1-p1-eta)*gammaC1/(p1*gammaI_1+(1-p1)*gammaC1))+((1-delta)*(S/(S+V+R))*(tau/(tau+sigmaE))+(V/(S+V+R))*(tau/(tau+sigmaE)))*(gammaC1/(gammaC1+sigmaC)))*(theta2*beta*kappa1/(gammaC3))*(1-epsilon)
    }
    
    return(list(c(dS, dE, dI, dC1, dC2, dC3, dR, dV, dEv, dIn, dAs, dRt)))
  }
  
  trajectory <- data.frame(ode(y = init.state,
                               times = times,
                               func = model_ode,
                               parms = theta,
                               method = "lsoda"))
  
  return(trajectory)
}

#### Below is the code to run simulations estimating impact of different interventions
calculate_diff <- function(posterior_samples,sim_model,seed.set){
  set.seed(seed.set)
  sample.index.CI <- sample(1:3400,5000,replace=TRUE)
  # param.post <- read.csv("output/trace-both-100k-full-final-mod-6-phase4.csv", header = T)
  param.post <- posterior_samples
  
  param.samples <- param.post[sample.index.CI,1:11]
  
  CI.samples <- matrix(nrow=5000,ncol=length(incidence_observed))
  CI.samples.total <- vector()
  CI.samples.all <- matrix(nrow=5000,ncol=length(incidence_observed))
  CI.samples.total.all <- vector()
  Reff.samples <- matrix(nrow=5000,ncol=119)
  
  CI.samples.baseline <- matrix(nrow=5000,ncol=length(incidence_observed))
  CI.samples.total.baseline <- vector()
  CI.samples.all.baseline <- matrix(nrow=5000,ncol=length(incidence_observed))
  CI.samples.total.all.baseline <- vector()
  Reff.samples.baseline <- matrix(nrow=5000,ncol=119)
  
  tmin = 300
  tmax = max(times_of_observed)
  t0 <- 300
  times.output = seq(tmin,tmax)
  
  for (i in 1:5000){
    init.state = c(S=0, E=0, I=0, C1=0, C2=0, C3=0, R=0, V=0, Ev=0, In=0, As=0, Rt=0)
    param.CI <- param.samples[i,]
    if ("S0" %in% names(param.CI)) {
      npop = 10400000
      init.state["I"]  = 1
      init.state["S"]  = npop * param.CI[["S0"]] - init.state["I"]
      init.state["E"]  = 0
      init.state["C1"] = 0
      init.state["C2"] = 0
      init.state["C3"] = 0
      init.state["R"]  = 0     
      init.state["V"]  = (npop - init.state["S"])
      init.state["Ev"] = 0
    }
    traj <- sim_model(param.CI, init.state, times.output)
    
    data.diff = diff(traj$In)[-(1:2)]
    R_eff = diff(traj$Rt)
    R_eff = R_eff[10:length(R_eff)]
    data.diff3 = data.diff
    
    diff.t = 344 - (t0)
    diff.t2 = diff.t + 1
    data.diff[1:diff.t] = data.diff[1:diff.t]*param.CI[[6]]
    data.diff[diff.t2:length(data.diff)] = data.diff[diff.t2:length(data.diff)]*param.CI[[7]]
    
    diff.t3 = 1
    
    data.diff2 = data.diff[diff.t3:length(data.diff)]
    data.diff4 = data.diff3[diff.t3:length(data.diff3)]
    diff.matrix = matrix(data.diff2,ncol = 7,byrow = T)
    diff.matrix2 = matrix(data.diff4,ncol = 7,byrow = T)
    
    incidence_predicted.CI = rowSums(diff.matrix)
    incidence_predicted.CI.all = rowSums(diff.matrix2)
    
    CI.samples[i,] <- incidence_predicted.CI
    CI.samples.total[i] <- sum(incidence_predicted.CI)
    CI.samples.all[i,] <- incidence_predicted.CI.all
    CI.samples.total.all[i] <- sum(incidence_predicted.CI.all)
    Reff.samples[i,] <- R_eff
    
    init.state = c(S=0, E=0, I=0, C1=0, C2=0, C3=0, R=0, V=0, Ev=0, In=0, As=0, Rt=0)
    param.CI <- param.samples[i,]
    if ("S0" %in% names(param.CI)) {
      npop = 10400000
      init.state["I"]  = 1
      init.state["S"]  = npop * param.CI[["S0"]] - init.state["I"]
      init.state["E"]  = 0
      init.state["C1"] = 0
      init.state["C2"] = 0
      init.state["C3"] = 0
      init.state["R"]  = 0     
      init.state["V"]  = (npop - init.state["S"])
      init.state["Ev"] = 0
    }
    traj <- model_baseline(param.CI, init.state, times.output)
    
    data.diff = diff(traj$In)[-(1:2)]
    R_eff = diff(traj$Rt)
    R_eff = R_eff[10:length(R_eff)]
    data.diff3 = data.diff
    
    diff.t = 344 - (t0)
    diff.t2 = diff.t + 1
    data.diff[1:diff.t] = data.diff[1:diff.t]*param.CI[[6]]
    data.diff[diff.t2:length(data.diff)] = data.diff[diff.t2:length(data.diff)]*param.CI[[7]]
    
    diff.t3 = 1
    
    data.diff2 = data.diff[diff.t3:length(data.diff)]
    data.diff4 = data.diff3[diff.t3:length(data.diff3)]
    diff.matrix = matrix(data.diff2,ncol = 7,byrow = T)
    diff.matrix2 = matrix(data.diff4,ncol = 7,byrow = T)
    
    incidence_predicted.CI = rowSums(diff.matrix)
    incidence_predicted.CI.all = rowSums(diff.matrix2)
    
    CI.samples.baseline[i,] <- incidence_predicted.CI
    CI.samples.total.baseline[i] <- sum(incidence_predicted.CI)
    CI.samples.all.baseline[i,] <- incidence_predicted.CI.all
    CI.samples.total.all.baseline[i] <- sum(incidence_predicted.CI.all)
    Reff.samples.baseline[i,] <- R_eff
  }
  
  diff.intv.all <- CI.samples.total.all - CI.samples.total.all.baseline
  diff.intv <- CI.samples.total - CI.samples.total.baseline
  
  diff.intv.all.summary = data.frame(x1=numeric(),
                                     x2=numeric(),
                                     x3=numeric(),
                                     x4=numeric(),
                                     x5=numeric(),
                                     x6=numeric())
  
  median.dat = median(diff.intv.all)
  mean.dat = mean(diff.intv.all)
  q.25.dat = quantile(diff.intv.all,0.25)
  q.75.dat = quantile(diff.intv.all,0.75)
  q.2.5.dat = quantile(diff.intv.all,0.025)
  q.97.5.dat = quantile(diff.intv.all,0.975)
  summary.all = c(median.dat,mean.dat,q.25.dat,q.75.dat,q.2.5.dat,q.97.5.dat)
  diff.intv.all.summary = rbind(diff.intv.all.summary,summary.all)
  
  diff.intv.summary = data.frame(x1=numeric(),
                                 x2=numeric(),
                                 x3=numeric(),
                                 x4=numeric(),
                                 x5=numeric(),
                                 x6=numeric())
  
  median.dat = median(diff.intv)
  mean.dat = mean(diff.intv)
  q.25.dat = quantile(diff.intv,0.25)
  q.75.dat = quantile(diff.intv,0.75)
  q.2.5.dat = quantile(diff.intv,0.025)
  q.97.5.dat = quantile(diff.intv,0.975)
  summary = c(median.dat,mean.dat,q.25.dat,q.75.dat,q.2.5.dat,q.97.5.dat)
  diff.intv.summary = rbind(diff.intv.summary,summary)
  
  colnames(diff.intv.all.summary) = c("median","mean","q25","q75","q2.5","q97.5")
  colnames(diff.intv.summary) = c("median","mean","q25","q75","q2.5","q97.5")
  
  return(list(diff.intv.all.summary=diff.intv.all.summary,
              diff.intv.summary=diff.intv.summary,
              diff.intv.all=diff.intv.all,
              diff.intv=diff.intv))
  
  # write.csv(diff.intv.all.summary,difference.all,row.names = F)
  # write.csv(diff.intv.summary,difference,row.names = F)
  # write.csv(diff.intv.all,difference.long.all,row.names = F)
  # write.csv(diff.intv,difference.long,row.names = F)
}
