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
  // Number of knots i.e. random effects
 int<lower=0> nk;
 // number of columns for each Z matrix
 int zvars[q+1];
 // Variables
 vector[N] y;
 // fixed effects design matrix
 matrix[N, p] X;
 // random effects design matrix
 matrix[N, nk] Z;
 // max col of Z
 int max_col;
 // array of Z
 // real Zarray[N, max_col, q];
 // matrix[N, max_col] Zarray[q];

 // indicator whether to use QR decomposition
 int<lower=0,upper=1> qr; // 0 = no, 1 = yes
 // indicator whether to split QR decomposition across multiple matrices
 int<lower=0,upper=1> qrsplit;
// family number:  1=gaussian, 2=binomial, 3=poisson
 int<lower=1, upper=3> famnum;
 // link number
 int linknum;
 // offset
 vector[N] offset;

 ////////////////////////////////////////////////////////
 // hyperpriors
 ////////////////////////////////////////////////////////

 // lambda
 int lambdanum[q+1];
 int lambda_max_params;
 matrix[q, lambda_max_params] lambda_param;

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

  // qr for Z
  matrix[N, nk] Q_z;
  matrix[nk, nk] R_z;
  matrix[nk, nk] R_z_inverse;

  // qr for Z with array of matrices
  matrix[N, max_col] qz[q];
  matrix[max_col, max_col] rz[q];
  matrix[max_col, max_col] rzinv[q];

  // matrix[max_col, max_col] R_z_inverse_temp;
  // matrix[N, max_col] Ztemp;

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
 vector[p] theta_b;
 vector[nk] tau;
 vector<lower=0>[q] lambda;
 vector<lower=0>[r] eps;
}

transformed parameters {
  vector[nk] theta_u;

// TODO: fix this loop
  for (l4 in 1:1) {
    int i = 1;
    for (j4 in 1:q) {
      for (k4 in 1:zvars[j4]) {
        theta_u[i] = tau[i] * lambda[j4];
        i = i + 1;
      }
    }
  }

}

model {
 // Prior part of Bayesian inference
 for (k3 in 1:q) {
   if (lambdanum[k3] == 1) {
     lambda[k3] ~ normal(lambda_param[k3, 1], lambda_param[k3, 2]);
   } else if (lambdanum[k3] == 2) {
     lambda[k3] ~ student_t(lambda_param[k3, 1], lambda_param[k3, 2], lambda_param[k3, 3]);
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


   for (k1 in 1:p) {
     if (betanum[k1] == 1) {
       theta_b[k1] ~ normal(beta_param[k1, 1], beta_param[k1, 2]);
     } else if (betanum[k1] == 2) {
       theta_b[k1] ~ student_t(beta_param[k1, 1], beta_param[k1, 2], beta_param[k1, 3]);
     }
   }

   for (k2 in 1:nk)
     tau[k2] ~ normal(0, 1);

    // identity link
    if (linknum == 1) {
      if (qr == 1) {
        y ~ normal(Q_x*theta_b + Q_z*theta_u + offset, eps[1]);

      } else {
        y ~ normal(X*theta_b + Z*theta_u + offset, eps[1]);
      }
    // log link
    } else if (linknum == 2) {
      if (qr == 1) {
        y ~ normal(exp(Q_x*theta_b + Q_z*theta_u + offset), eps[1]);

      } else {
        y ~ normal(exp(X*theta_b + Z*theta_u + offset), eps[1]);
      }
    // inverse link
    } else if (linknum == 3) {
      if (qr == 1) {
        y ~ normal(inv(Q_x*theta_b + Q_z*theta_u + offset), eps[1]);
      } else {
        y ~ normal(inv(X*theta_b + Z*theta_u + offset), eps[1]);
      }
    }

 }
}

generated quantities {
  vector[p] beta;
  vector[nk] u;
  vector[N] log_lik;

  if (qr == 1) {
    beta = R_x_inverse * theta_b; // coefficients on x
    u = R_z_inverse * theta_u;
  } else {
    beta = theta_b;
    u = theta_u;
  }
  
    // extract log lik
    for (n in 1:N) {
       if (linknum == 1) {
         log_lik[n] = normal_lpdf(y[n] | X[n, ]*beta + Z[n,]*u, eps);
       }
       // log link
       else if (linknum == 2) {
         log_lik[n] = normal_lpdf(exp(y[n]) | X[n, ]*beta + Z[n,]*u, eps);
       }
       // inverse link
       else if (linknum == 3) {
         log_lik[n] = normal_lpdf(inv(y[n]) | X[n, ]*beta + Z[n,]*u, eps);
       }
    }
}

