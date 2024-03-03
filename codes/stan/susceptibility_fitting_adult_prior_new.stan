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
  vector[N] c[P];      // Estimated incidence
  simplex[N] theta[P]; // Probabilities for each category for each phase
  
  for (p in 1:P){
    for (n in 1:N){
      vector[N] y_over_pop = y[p] ./ pop;  // Element-wise division
      lambda[p][n] = sum(X[p][n] .* to_row_vector(y_over_pop));
      c[p][n] = s[n] * lambda[p][n] * pop[n]; // Calculating estimated incidence
    }
    theta[p] = c[p] ./ sum(c[p]); // Calculating theta using the proportion of c
  }
}

model {
  // Model fitting for each phase
  for (p in 1:P) {
    z[p] ~ multinomial(theta[p]);
  }
  
  s[1] ~ beta(30, 90); // Prior for the first element of susceptibility (0-4 years old)
  s[5] ~ beta(50,50); // Prior for the sixth element of susceptibility (20-59 years old)
  s[6] ~ beta(50,50); // Prior for the sixth element of susceptibility (60+ years old)
  for (n in 2:4) {
    s[n] ~ beta(1, 1); // Prior for the remaining elements
  }
}

generated quantities {
  // Expected count for each category in each phase
  int expected_z[P, N];
  
  for (p in 1:P)
    expected_z[p] = multinomial_rng(theta[p], sum(z[p]));
}
