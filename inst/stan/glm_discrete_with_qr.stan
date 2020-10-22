
data {
 // Define variables in data
 // Number of observations (an integer)
 int<lower=0> N;
 // Number of fixed parameters
 int<lower=0> p;
 // Variables
 int y[N]; // binomial, poisson
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
}

transformed parameters {
  vector[p] beta;

  if (qr == 1) {
    beta = R_x_inverse * theta_b; // coefficients on x
  } else {
    beta = theta_b;
  }
}

model {
   for (k1 in 1:p) {
     if (betanum[k1] == 1) {
       theta_b[k1] ~ normal(beta_param[k1, 1], beta_param[k1, 2]);
     } else if (betanum[k1] == 2) {
       theta_b[k1] ~ student_t(beta_param[k1, 1], beta_param[k1, 2], beta_param[k1, 3]);
     }
   }

 // Binomial
 if (famnum == 2) {
   // logit link
   if (linknum == 4) {
      if (qr == 1) {
        y ~ bernoulli_logit(Q_x * theta_b + offset);
      } else {
        y ~ bernoulli_logit(X * theta_b + offset);
      }
    // probit link
   } else if (linknum == 5) {
     if (qr == 1) {
        y ~ bernoulli(Phi(Q_x * theta_b + offset));
      } else {
        y ~ bernoulli(Phi(X * theta_b + offset));
      }
    // cauchit link
   } else if (linknum == 6) {
     if (qr == 1) {
        y ~ bernoulli(atan(Q_x * theta_b + offset) / pi() + 0.5);
      } else {
        y ~ bernoulli(atan(X * theta_b + offset) / pi() + 0.5);
      }
    // log link
   } else if (linknum == 2) {
     if (qr == 1) {
        y ~ bernoulli(exp(Q_x * theta_b + offset));
      } else {
        y ~ bernoulli(exp(X * theta_b + offset));
      }
    // cloglog link
   } else if (linknum == 7) {
     if (qr == 1) {
        y ~ bernoulli(inv_cloglog(Q_x * theta_b + offset));
      } else {
        y ~ bernoulli(inv_cloglog(X * theta_b + offset));
      }
   }


  // Poisson
 } else if (famnum == 3) {
   // log link
   if (linknum == 2) {
      if (qr == 1) {
        y ~ poisson_log(Q_x * theta_b + offset);
      } else {
        y ~ poisson_log(X * theta_b + offset);
      }
    // identity link
   } else if (linknum == 1) {
     if (qr == 1) {
        y ~ poisson(Q_x * theta_b + offset);
      } else {
        y ~ poisson(X * theta_b + offset);
      }
    // sqrt link
   } else if (linknum == 8) {
     if (qr == 1) {
        y ~ poisson(square(Q_x * theta_b + offset));
      } else {
        y ~ poisson(square(X * theta_b + offset));
      }
   }
 }
}

generated quantities {
  vector[N] log_lik;

  // extract log_lik
   for (n in 1:N) {
     // Binomial
     if (famnum == 2) {
       // logit link
       if (linknum == 4) {
          log_lik[n] = bernoulli_logit_lpmf(y[n] | X[n, ]*beta);
        // probit link
       } else if (linknum == 5) {
          log_lik[n] = bernoulli_lpmf( y[n] | Phi(X[n, ]*beta));
        // cauchit link
       } else if (linknum == 6) {
         log_lik[n] = bernoulli_lpmf( y[n] | atan(X[n, ]*beta)/pi() + 0.5);
        // log link
       } else if (linknum == 2) {
         log_lik[n] = bernoulli_lpmf( y[n] | exp(X[n, ]*beta));
        // cloglog link
       } else if (linknum == 7) {
         log_lik[n] = bernoulli_lpmf( y[n] | inv_cloglog(X[n, ]*beta));
       }
    
    
      // Poisson
     } else if (famnum == 3) {
       // log link
       if (linknum == 2) {
         log_lik[n] = poisson_log_lpmf(y[n] | X[n, ]*beta);
        // identity link
       } else if (linknum == 1) {
         log_lik[n] = poisson_lpmf(y[n] | X[n, ]*beta);
        // sqrt link
       } else if (linknum == 8) {
         log_lik[n] = poisson_lpmf(y[n] | square(X[n, ]*beta));
       }
     }
   }
   
}


