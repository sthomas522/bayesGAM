// Continuous multresponse
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
 // max col of Z
 int max_col;
 // number of columns for each Z matrix
 int zvars[q+1];
 // number of columns for Z intercept
 int<lower=0> nrandint;
 // number of columns for Z nonparametric
 int<lower=0> nnp;
 // random effects design matrix
 matrix[N, nnp] Znp;
 // random effects random intercept matrix
 matrix[N, nrandint] Zint;
 // matrix[N, max_col] Zarray[q];
 matrix[N, nk] Z;

 // indicator whether to use QR decomposition
 int<lower=0,upper=1> qr; // 0 = no, 1 = yes
 // indicator whether to split QR decomposition across multiple matrices
 int<lower=0,upper=1> qrsplit;
 // indicator of multivariate independence
 int<lower=0,upper=1> mvindep;
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

  // lambda for nonparametric
 int lambdanum[q*r+1];
 int lambda_max_params;
 matrix[q*r, lambda_max_params] lambda_param;

 // number of off-diagonal
 int a_num_offdiagonal;
 int anum[a_num_offdiagonal+1];
 int a_max_params;
 matrix[a_num_offdiagonal, a_max_params] a_param;
}

transformed data {
  // qr for X
  matrix[N, p] Q_x;
  matrix[p, p] R_x;
  matrix[p, p] R_x_inverse;

  // qr for Z knots matrices
  matrix[N, nnp] Q_z;
  matrix[nnp, nnp] R_z;
  matrix[nnp, nnp] R_z_inverse;

  // qr for Z with array of matrices
  // matrix[N, max_col] qz[q];
  // matrix[max_col, max_col] rz[q];
  // matrix[max_col, max_col] rzinv[q];

  // thin and scale the QR decomposition X
  Q_x = qr_Q(X)[, 1:p] * sqrt(N - 1);
  R_x = qr_R(X)[1:p, ] / sqrt(N - 1);
  R_x_inverse = inverse(R_x);

  // thin and scale the QR decomposition
  if (nnp > 0) {
    Q_z = qr_Q(Znp)[, 1:nnp] * sqrt(N - 1);
    R_z = qr_R(Znp)[1:nnp, ] / sqrt(N - 1);
    R_z_inverse = inverse(R_z);    
  } 


}

parameters {
 // Define parameters to estimate
 vector[p] theta_b[ny];

 // TODO: reverse indices
 matrix[nrandint, ny] trans_u_random;
 vector<lower=0>[ny] lambda_random;

 // nonparametric, if any
 vector[nnp] tau[ny];

 // residual sd
 vector<lower=0>[q-1] lambda[ny];
 vector<lower=0>[r] eps;

 // a parameters
 vector[mvindep ? 0 : a_num_offdiagonal] a;
}

transformed parameters {
  vector[nnp] theta_u[ny];
  vector[N] yhat[ny];

  /////////////////////////////////////////////////////////////
  // random intercept
  // create diagonal covariance matrix for now
  matrix[ny, ny] sigma_u_random;

  // local block
  for (ll2 in 1:1)
  {
    matrix[ny, ny] L;
    matrix[ny, ny] Dhalf;

    // assign LDLT decomposition
    Dhalf = diag_matrix(lambda_random);
    L = diag_matrix(rep_vector(1.0, ny));
    for (ll in 1:1) {
      int iter = 1;
       for (ii in 1:ny) {
        for (jj in 1:ny) {
          if (jj > ii) {
            if (mvindep == 1) {
              L[jj, ii] = 0;
            } else {
              L[jj, ii] = a[iter];  
            }
            iter = iter + 1;
          }
        }
      }
    }
    

  // sigma_u_random = L * Dhalf * Dhalf * (L');
  sigma_u_random = tcrossprod(L * Dhalf);
  
  }
  

  /////////////////////////////////////////////////////////////
  // alternate approach
  if (q >= 2) {
    for (l4 in 1:ny) {
      int i = 1;
      for (j4 in 2:q) {
        for (k4 in 1:zvars[j4]) {
          theta_u[l4][i] = tau[l4][i] * lambda[l4][(j4-1)];
          i = i + 1;
        }
      }
    }    
  }

  // multivariate response
  for (jj in 1:ny) {
   yhat[jj] = Q_x*theta_b[jj] + Zint*col(trans_u_random, jj);
  }
  
   // add if nonparametric terms present
  if (q >= 2) {
    for (jj in 1:ny) {
      yhat[jj] = yhat[jj] + Q_z[jj]*theta_u[jj]; 
    }    
  }

}

model {
 // Prior part of Bayesian inference

 // off diagonal w prior input
 if (mvindep == 0) {
   for (jj in 1:a_num_offdiagonal) {
     // a[jj] ~ normal(0, 1e6);
    if (anum[jj] == 1) {
       a[jj] ~ normal(a_param[jj, 1], a_param[jj, 2]);
     } else if (anum[jj] == 2) {
       a[jj] ~ student_t(a_param[jj, 1], a_param[jj, 2], a_param[jj, 3]);
     }
   }
 }

 for (jj in 1:nrandint) {
    trans_u_random[jj] ~ multi_normal(rep_vector(0.0, ny), sigma_u_random);
 }

 for (j1 in 1:r) {
   for (k1 in 1:p) {
     if (betanum[k1*j1] == 1) {
       theta_b[j1, k1] ~ normal(beta_param[k1*j1, 1], beta_param[k1*j1, 2]);
     } else if (betanum[k1] == 2) {
       theta_b[j1, k1] ~ student_t(beta_param[k1*j1, 1], beta_param[k1*j1, 2], beta_param[k1*j1, 3]);
     }
   }
  }

  // nested loop for multvariate response
  for (j1 in 1:r) {
    if (lambdanum[1] == 1) {
      lambda_random[j1] ~ normal(lambda_param[j1, 1], lambda_param[j1, 2]);
    } else if (lambdanum[1] == 2) {
      lambda_random[j1] ~ student_t(lambda_param[j1, 1], lambda_param[j1, 2], lambda_param[j1, 3]);
    }
    // TODO restore nonparametric

    if (q >= 2) {
      for (k1 in 1:(q-1)) {
       if (lambdanum[k1*j1] == 1) {
         lambda[j1, k1] ~ normal(lambda_param[k1*j1, 1], lambda_param[k1*j1, 2]);
       } else if (lambdanum[k1*j1] == 2) {
         lambda[j1, k1] ~ student_t(lambda_param[k1*j1, 1], lambda_param[k1*j1, 2], lambda_param[k1*j1, 3]);
       }
     }

    }
  }
  
 for (jj in 1:ny)
   tau[jj] ~ normal(0, 1);

  // nested loop for multvariate response
  if (famnum == 1) {

     for (k2 in 1:r) {
        if (epsnum[k2] == 1) {
         eps[k2] ~ normal(eps_param[k2, 1], eps_param[k2, 2]);
       } else if (epsnum[k2] == 2) {
         eps[k2] ~ student_t(eps_param[k2, 1], eps_param[k2, 2], eps_param[k2, 3]);
       }
     }

       // identity link
      for (jj in 1:ny) {
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
  vector[p] beta[ny];
  vector[ny] dhalf_inv;
  matrix[ny, ny] sigma_u_correlation;
  vector[nnp] nonpar[ny];
  vector[nrandint+nnp] u[ny];
  vector[N] log_lik[ny];

  if (qr == 1) {
      for (jj in 1:ny) {
        beta[jj] = R_x_inverse * theta_b[jj];
        if (q >= 2) {
           nonpar[jj] = R_z_inverse * theta_u[jj];
        }

      }
    } else {
      for (jj in 1:ny) {
        beta[jj] = theta_b[jj];
       if (q >= 2) {
         nonpar[jj] = theta_u[jj];
       }
      }
  }


  // correlation matrix
  dhalf_inv = diagonal(sigma_u_random); 
  for (jj in 1:ny) {
    dhalf_inv[jj] = 1 / sqrt(dhalf_inv[jj]);
  }

  sigma_u_correlation = quad_form_diag(sigma_u_random, dhalf_inv);

  for (jj in 1:ny) {
    for (kk in 1:nrandint) {
     u[jj][kk] = trans_u_random[kk][jj];
    }
    if (q >= 2) {
      for (ll in 1:nnp) {
        u[jj][nrandint+ll] = nonpar[jj][ll]; 
      }      
    }

  }
  
  // extract log lik
  for (jj in 1:ny) {
    for (n in 1:N) {
         if (linknum == 1) {
           log_lik[jj][n] = normal_lpdf(y[n] | X[n, ]*beta[jj] + Z[n,]*u[jj], eps[jj]);
         }
         // log link
         else if (linknum == 2) {
           log_lik[jj][n] = normal_lpdf(exp(y[n]) | X[n, ]*beta[jj] + Z[n,]*u[jj], eps[jj]);
         }
         // inverse link
         else if (linknum == 3) {
           log_lik[jj][n] = normal_lpdf(inv(y[n]) | X[n, ]*beta[jj] + Z[n,]*u[jj], eps[jj]);
         }
    }
  }
  
  
}

