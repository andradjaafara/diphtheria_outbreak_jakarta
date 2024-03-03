data {
  int<lower=0> N; // Number of age groups
  int<lower=0> P; // Number of phases
  matrix[N, N] X[P]; // Contact matrix
  vector[N] y[P]; // Incidence vector for each phase
  int z[P, N]; // Counts data matrix for each phase
  vector<lower=0>[N] pop; // Population size for each age group
}

parameters {  
  vector<lower=0, upper=1>[N] s; // Susceptibility for each age group
}

transformed parameters {
  vector[N] lambda[P]; // Force of infections
  simplex[N] theta[P]; // Probabilities for each category for each phase
  
  for (p in 1:P){
    for (n in 1:N){
      lambda[p][n] = (s[n] * sum(X[p][n] .* to_row_vector(y[p])));
	}
	theta[p] = lambda[p] / sum(lambda[p]);
  }
}

model {
  // Model fitting for each phase
  for (p in 1:P) {
    z[p] ~ multinomial(theta[p]);
  }
  
  s ~ beta(1, 1); // Prior for susceptibility
}

generated quantities {
  // Expected count for each category in each phase
  int expected_z[P, N];
  
  for (p in 1:P)
    expected_z[p] = multinomial_rng(theta[p], sum(z[p]));
}
