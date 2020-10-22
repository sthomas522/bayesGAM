// discrete no random effects
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
 // indicator whether to use QR decomposition
 int<lower=0,upper=1> qr; // 0 = no, 1 = yes
 // indicator whether to split QR decomposition across multiple matrices
 // int<lower=0,upper=1> qrsplit;

 // family number:  1=gaussian, 2=binomial, 3=poisson
 int<lower=1, upper=3> famnum;
 // link number
 int linknum;
 // offset - TODO currently unused
 // vector[N] offset;

 ////////////////////////////////////////////////////////
 // hyperpriors
 ////////////////////////////////////////////////////////

 // beta
 int betanum[p*r+1];
 int beta_max_params;
 matrix[p*r, beta_max_params] beta_param;
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
 vector[p] theta_b[ny];

 // u = lambda * tau
 // random intercept
 // vector[nk] tau[ny];

 // residual sd
 // vector<lower=0>[q] lambda[ny];
 // vector<lower=0>[r] eps;
}

transformed parameters {
  // vector[nk] theta_u[ny];
  // vector[nk] u[ny];
  vector[p] beta[ny];
  
  if (qr == 1) {
    for (jj in 1:ny) {
      beta[jj] = R_x_inverse * theta_b[jj];
    }
  } else
    {
      for (jj in 1:ny)
      {
        beta[jj] = theta_b[jj];
      }
    }
    
}

model {
 // Prior part of Bayesian inference

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

  if (famnum == 2) {
    vector[N] yhat[ny];
  
    for (jj in 1:ny) {
      yhat[jj] = Q_x*theta_b[jj];
    }

   // prior on random effects variance
   // random intercepts
   for (jj in 1:ny) {
     // lambda[jj] ~ student_t(35, 0, 1e0);
     // tau[jj] ~ normal(0, 1);

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
    vector[N] yhat[ny];
  
    for (jj in 1:ny) {
      yhat[jj] = Q_x*theta_b[jj];
    }

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
  vector[N] log_lik[ny];
  
  // extract log_lik
  for (jj in 1:ny) {
   for (n in 1:N) {
     // Binomial
     if (famnum == 2) {
       // logit link
       if (linknum == 4) {
          log_lik[jj][n] = bernoulli_logit_lpmf(y[n] | X[n, ]*beta[jj]);
        // probit link
       } else if (linknum == 5) {
          log_lik[jj][n] = bernoulli_lpmf( y[n] | Phi(X[n, ]*beta[jj] ));
        // cauchit link
       } else if (linknum == 6) {
         log_lik[jj][n] = bernoulli_lpmf( y[n] | atan(X[n, ]*beta[jj])/pi() + 0.5);
        // log link
       } else if (linknum == 2) {
         log_lik[jj][n] = bernoulli_lpmf( y[n] | exp(X[n, ]*beta[jj]));
        // cloglog link
       } else if (linknum == 7) {
         log_lik[jj][n] = bernoulli_lpmf( y[n] | inv_cloglog(X[n, ]*beta[jj]));
       }
    
    
      // Poisson
     } else if (famnum == 3) {
       // log link
       if (linknum == 2) {
         log_lik[jj][n] = poisson_log_lpmf(y[n] | X[n, ]*beta[jj]);
        // identity link
       } else if (linknum == 1) {
         log_lik[jj][n] = poisson_lpmf(y[n] | X[n, ]*beta[jj]);
        // sqrt link
       } else if (linknum == 8) {
         log_lik[jj][n] = poisson_lpmf(y[n] | square(X[n, ]*beta[jj]));
       }
     }
   }
  }

}

