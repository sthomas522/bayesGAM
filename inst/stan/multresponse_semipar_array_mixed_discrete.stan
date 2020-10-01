// discrete mixed no random intercept
data {
 // Define variables in data
 // Number of observations (an integer)
 int<lower=0> N;
 // Number of fixed parameters
 int<lower=0> p;
 // Number of lambda param
 int<lower=0> q;
 // Number of eps param
 int<lower=0> r;
 // Number of random effect parameters or knots
 int<lower=0> nk;
 int<lower=0> ny;
 // matrix[N, ny] y;
 int y[N, ny];
 // Variables
 // fixed effects design matrix
 matrix[N, p] X;
 // max col of Z
 int max_col;
 // number of columns for each Z matrix
 int zvars[q+1];
 // random effects design matrix
 matrix[N, nk] Z;
 // random effects random intercept matrix
 // matrix[N, max_col] Zarray[q];

 // indicator whether to use QR decomposition
 int<lower=0,upper=1> qr; // 0 = no, 1 = yes
 // indicator whether to split QR decomposition across multiple matrices
 int<lower=0,upper=1> qrsplit;

 // family number:  1=gaussian, 2=binomial, 3=poisson
 int<lower=1, upper=3> famnum;
 // link number
 int linknum;
 // offset - TODO currently unused
 // vector[N] offset;


  // TODO:  feedback
  // matrix ny x q:  1's or 0's to include/exclude individual Z

   ////////////////////////////////////////////////////////
 // hyperpriors
 ////////////////////////////////////////////////////////

 // beta
 int betanum[p*r+1];
 int beta_max_params;
 matrix[p*r, beta_max_params] beta_param;

  // lambda
 int lambdanum[q*r+1];
 int lambda_max_params;
 matrix[q*r, lambda_max_params] lambda_param;

}

transformed data {
  // qr for X
  matrix[N, p] Q_x;
  matrix[p, p] R_x;
  matrix[p, p] R_x_inverse;

  // qr for Z knots matrices
  matrix[N, nk] Q_z;
  matrix[nk, nk] R_z;
  matrix[nk, nk] R_z_inverse;

  // qr for Z with array of matrices
  matrix[N, max_col] qz[q];
  matrix[max_col, max_col] rz[q];
  matrix[max_col, max_col] rzinv[q];

  // thin and scale the QR decomposition
  Q_x = qr_Q(X)[, 1:p] * sqrt(N - 1);
  R_x = qr_R(X)[1:p, ] / sqrt(N - 1);
  R_x_inverse = inverse(R_x);

  // thin and scale the QR decomposition
  Q_z = qr_Q(Z)[, 1:nk] * sqrt(N - 1);
  R_z = qr_R(Z)[1:nk, ] / sqrt(N - 1);
  R_z_inverse = inverse(R_z);

}

parameters {
 // Define parameters to estimate
 vector[p] theta_b[ny];

 // u = lambda * tau
 // random intercept
 vector[nk] tau[ny];

 // residual sd
 vector<lower=0>[q] lambda[ny];
 // vector<lower=0>[r] eps;

}

transformed parameters {
  vector[nk] theta_u[ny];
  vector[N] yhat[ny];

  // TODO: fix this loop
  for (jj in 1:ny) {
    for (l4 in 1:1) {
      int i = 1;
      for (j4 in 1:q) {
        for (k4 in 1:zvars[j4]) {
          theta_u[jj][i] = tau[jj][i] * lambda[jj][j4];
          i = i + 1;
        }
      }
    }
  }

  // adding restriction design matrix for feedback loop
  for (jj in 1:ny) {
    yhat[jj] = Q_x*theta_b[jj] + Q_z*theta_u[jj];
    // yhat[jj] = Q_x*theta_b[jj] + Q_z * diag_matrix(feedback[jj]) *theta_u[jj] ;
    // yhat[jj] = Q_x*theta_b[jj] + diag_post_multiply(Q_z, feedback[jj]) *theta_u[jj] ;
  }

}

model {
 // Prior part of Bayesian inference

 // fixed effects prior
 // for (k1 in 1:p) {
 //   for (jj in 1:ny) {
 //     theta_b[jj] ~ normal(0, 1e6);
 //   }
//  }

  // nested loop for multvariate response
  for (j1 in 1:r) {
   for (k1 in 1:p) {
     if (betanum[k1*j1] == 1) {
       theta_b[j1, k1] ~ normal(beta_param[k1*j1, 1], beta_param[k1*j1, 2]);
     } else if (betanum[k1] == 2) {
       theta_b[j1, k1] ~ student_t(beta_param[k1*j1, 1], beta_param[k1*j1, 2], beta_param[k1*j1, 3]);
     }
   }
  }
    for (j1 in 1:r) {
     for (k1 in 1:q) {
       if (lambdanum[k1*j1] == 1) {
         lambda[j1, k1] ~ normal(lambda_param[k1*j1, 1], lambda_param[k1*j1, 2]);
       } else if (lambdanum[k1*j1] == 2) {
         lambda[j1, k1] ~ student_t(lambda_param[k1*j1, 1], lambda_param[k1*j1, 2], lambda_param[k1*j1, 3]);
       }
     }
    }

  for (jj in 1:ny)
    tau[jj] ~ normal(0, 1);

  if (famnum == 2) {
   for (jj in 1:ny) {
      // logit link
       if (linknum == 4) {
         y[, jj] ~ bernoulli_logit(yhat[jj]);
       }
        // probit link
       else if (linknum == 5) {
        y[, jj] ~ bernoulli(Phi(yhat[jj]));
       }
       // cauchit link
       else if (linknum == 6) {
        y[, jj] ~ bernoulli(atan(yhat[jj]) / pi() + 0.5);
       }
       // log link
       else if (linknum == 2) {
        y[, jj] ~ bernoulli(exp(yhat[jj]));
       }
       // cloglog link
       else if (linknum == 7) {
        y[, jj] ~ bernoulli(inv_cloglog(yhat[jj]));
       }
    }
   } else if (famnum == 3) {
     for (jj in 1:ny) {
         // log link
         if (linknum == 2) {
           y[, jj] ~ poisson_log(yhat[jj]);
         // identity link
         } else if (linknum == 1) {
           y[, jj] ~ poisson(yhat[jj]);
         // sqrt link
         } else if (linknum == 8) {
           y[, jj] ~ poisson(square(yhat[jj]));
         }
     }
   }

}
generated quantities {

  vector[nk] u[ny];
  vector[p] beta[ny];
  vector[N] log_lik[ny];

  if (qr == 1) {
    for (jj in 1:ny) {
      beta[jj] = R_x_inverse * theta_b[jj];
    }
    for (jj in 1:ny) {
      u[jj] = R_z_inverse * theta_u[jj];
    }
  } else
    {
      for (jj in 1:ny)
      {
        beta[jj] = theta_b[jj];
        u[jj] = theta_u[jj];
      }
    }
    
  // extract log_lik
  for (jj in 1:ny) {
   for (n in 1:N) {
     // Binomial
     if (famnum == 2) {
       // logit link
       if (linknum == 4) {
          log_lik[jj][n] = bernoulli_logit_lpmf(y[n] | X[n, ]*beta[jj] + Z[n,]*u[jj]);
        // probit link
       } else if (linknum == 5) {
          log_lik[jj][n] = bernoulli_lpmf( y[n] | Phi(X[n, ]*beta[jj] + Z[n,]*u[jj] ));
        // cauchit link
       } else if (linknum == 6) {
         log_lik[jj][n] = bernoulli_lpmf( y[n] | atan(X[n, ]*beta[jj] + Z[n,]*u[jj])/pi() + 0.5);
        // log link
       } else if (linknum == 2) {
         log_lik[jj][n] = bernoulli_lpmf( y[n] | exp(X[n, ]*beta[jj] + Z[n,]*u[jj]));
        // cloglog link
       } else if (linknum == 7) {
         log_lik[jj][n] = bernoulli_lpmf( y[n] | inv_cloglog(X[n, ]*beta[jj] + Z[n,]*u[jj]));
       }
    
    
      // Poisson
     } else if (famnum == 3) {
       // log link
       if (linknum == 2) {
         log_lik[jj][n] = poisson_log_lpmf(y[n] | X[n, ]*beta[jj] + Z[n,]*u[jj]);
        // identity link
       } else if (linknum == 1) {
         log_lik[jj][n] = poisson_lpmf(y[n] | X[n, ]*beta[jj] + Z[n,]*u[jj]);
        // sqrt link
       } else if (linknum == 8) {
         log_lik[jj][n] = poisson_lpmf(y[n] | square(X[n, ]*beta[jj] + Z[n,]*u[jj]));
       }
     }
   }
  }

  

}

