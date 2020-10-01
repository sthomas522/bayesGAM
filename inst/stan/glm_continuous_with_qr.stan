
data {
 // Define variables in data
 // Number of observations (an integer)
 int<lower=0> N;
 // Number of fixed parameters
 int<lower=0> p;
 // Number of eps param
 int<lower=0> r;
 // Variables
 vector[N] y; // gaussian
 // fixed effects design matrix
 matrix[N, p] X;
 // indicator whether to use QR decomposition
 int<lower=0,upper=1> qr; // 0 = no, 1 = yes
 // family number:  1=gaussian, 2=binomial, 3=poisson
 int<lower=1, upper=3> famnum;
 // link number
 int linknum;
  // offset
 vector[N] offset;

 ////////////////////////////////////////////////////////
 // hyperpriors
 ////////////////////////////////////////////////////////

 // epsilon
 int epsnum[r+1];
 int eps_max_params;
 matrix[r, eps_max_params] eps_param;

 // beta
 int betanum[p+1];
 int beta_max_params;
 matrix[p, beta_max_params] beta_param;
}

transformed data {
  // qr for X
  matrix[N, p] Q_x;
  matrix[p, p] R_x;
  matrix[p, p] R_x_inverse;

  // thin and scale the QR decomposition
  Q_x = qr_Q(X)[, 1:p] * sqrt(N - 1);
  R_x = qr_R(X)[1:p, ] / sqrt(N - 1);
  R_x_inverse = inverse(R_x);

}

parameters {
 // Define parameters to estimate
 vector[p] theta_b;
 vector<lower=0>[r] eps;
}

model {
 // Prior part of Bayesian inference

   for (k1 in 1:p) {
     if (betanum[k1] == 1) {
       theta_b[k1] ~ normal(beta_param[k1, 1], beta_param[k1, 2]);
     } else if (betanum[k1] == 2) {
       theta_b[k1] ~ student_t(beta_param[k1, 1], beta_param[k1, 2], beta_param[k1, 3]);
     }
   }

 // Gaussian
 if (famnum == 1) {
   for (k2 in 1:r) {
      if (epsnum[k2] == 1) {
       eps[k2] ~ normal(eps_param[k2, 1], eps_param[k2, 2]);
     } else if (epsnum[k2] == 2) {
       eps[k2] ~ student_t(eps_param[k2, 1], eps_param[k2, 2], eps_param[k2, 3]);
     }
   }


    // Identity
    if (linknum == 1) {
      if (qr == 1) {
        y ~ normal(Q_x*theta_b + offset, eps[1]);
      } else {
        y ~ normal(X*theta_b + offset, eps[1]);
      }
    // Log
    } else if (linknum == 2) {
      if (qr == 1) {
        y ~ normal(exp(Q_x*theta_b + offset), eps[1]);
      } else {
        y ~ normal(exp(X*theta_b + offset), eps[1]);
      }
    // Inverse
    } else if (linknum == 3) {
      if (qr == 1) {
        y ~ normal(inv(Q_x*theta_b + offset), eps[1]);
      } else {
        y ~ normal(inv(X*theta_b + offset), eps[1]);
      }
    }

  }
}

generated quantities {
  vector[p] beta;
  vector[N] log_lik;
  if (qr == 1) {
    beta = R_x_inverse * theta_b; // coefficients on x
  } else {
    beta = theta_b;
  }
  
    
    // extract log lik
    for (n in 1:N) {
       if (linknum == 1) {
         log_lik[n] = normal_lpdf(y[n] | X[n, ]*beta, eps);
       }
       // log link
       else if (linknum == 2) {
         log_lik[n] = normal_lpdf(exp(y[n]) | X[n, ]*beta, eps);
       }
       // inverse link
       else if (linknum == 3) {
         log_lik[n] = normal_lpdf(inv(y[n]) | X[n, ]*beta, eps);
       }
    }
}


