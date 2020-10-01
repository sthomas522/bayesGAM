// continuous no random effects
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
 matrix[N, ny] y;
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

 // epsilon
 int epsnum[r+1];
 int eps_max_params;
 matrix[r, eps_max_params] eps_param;

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
 vector[nk] tau[ny];

 // residual sd
 vector<lower=0>[q] lambda[ny];
 vector<lower=0>[r] eps;

}

transformed parameters {
  vector[nk] theta_u[ny];
  vector[N] yhat[ny];

  for (jj in 1:ny) {
    yhat[jj] = Q_x*theta_b[jj];
  }

}

model {
 // Prior part of Bayesian inference

 // fixed effects prior
 // for (k1 in 1:p) {
//    for (jj in 1:ny) {
//      theta_b[jj] ~ normal(0, 1e6);
//    }
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

  if (famnum == 1) {
    for (k2 in 1:r) {
      if (epsnum[k2] == 1) {
       eps[k2] ~ normal(eps_param[k2, 1], eps_param[k2, 2]);
     } else if (epsnum[k2] == 2) {
       eps[k2] ~ student_t(eps_param[k2, 1], eps_param[k2, 2], eps_param[k2, 3]);
     }
   }


   // prior on random effects variance
   // random intercepts
   for (jj in 1:ny) {
     lambda[jj] ~ student_t(35, 0, 1e0);
     tau[jj] ~ normal(0, 1);
     // eps[jj][r] ~ student_t(1, 0, 25);
    //  y[, jj] ~ normal(yhat[jj], eps[jj]);

      // identity link
     if (linknum == 1) {
       y[, jj] ~ normal(yhat[jj], eps[jj]);
     }
     // log link
     else if (linknum == 2) {
       y[, jj] ~ normal(exp(yhat[jj]), eps[jj]);
     }
     // inverse link
     else if (linknum == 3) {
       y[, jj] ~ normal(inv(yhat[jj]), eps[jj]);
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
  } else
    {
      for (jj in 1:ny)
      {
        beta[jj] = theta_b[jj];
      }
    }
    
  // extract log lik
  for (jj in 1:ny) {
    for (n in 1:N) {
         if (linknum == 1) {
           log_lik[jj][n] = normal_lpdf(y[n] | X[n, ]*beta[jj], eps[jj]);
         }
         // log link
         else if (linknum == 2) {
           log_lik[jj][n] = normal_lpdf(exp(y[n]) | X[n, ]*beta[jj], eps[jj]);
         }
         // inverse link
         else if (linknum == 3) {
           log_lik[jj][n] = normal_lpdf(inv(y[n]) | X[n, ]*beta[jj], eps[jj]);
         }
    }
  }

}

