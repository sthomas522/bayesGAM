test_that("test_each_stan_prog", {
  # cover all stan progs not covered by other test file
  
  # create sample data
  set.seed(651)
  X <- matrix(rnorm(100*5), ncol=5)
  y1 <- sample(x=c(0,1), size=100, replace=TRUE, prob=c(0.5, 0.5))
  y2 <- sample(x=c(0,1), size=100, replace=TRUE, prob=c(0.5, 0.5))
  idnum <- rep(1:20, each=5)
  dat <- data.frame(X = X, 
                    y1 = y1,
                    y2 = y2, 
                    y1cont = rnorm(100), 
                    y2cont = rnorm(100), 
                    idnum=factor(idnum))
  
  # glm_discrete_mixed_with_qr
  # logistic regression
  f1 <- bayesGAM(y1 ~ X.1 + X.2 + X.3 + X.4 + X.5, 
                 data=dat, family=binomial, 
                 random = ~idnum, 
                 chains = 1, iter = 500)
  cf1 <- coef(f1)
  expect_equal(length(cf1), 27)
  expect_true(sum(abs(cf1)) > 0)
  
  # multresponse_semipar_array
  f2 <- bayesGAM(cbind(y1cont, y2cont) ~ 
                   X.1 + X.2 + X.3 + X.4 + X.5, 
                 data=dat, family=gaussian, 
                 chains=1, iter=500)
  cf2 <- coef(f2)
  expect_equal(length(cf2), 14)
  expect_true(sum(abs(cf2)) > 0)
  
  # multresponse_semipar_array_discrete
  f3 <- bayesGAM(cbind(y1, y2) ~ 
                   X.1 + X.2 + X.3 + X.4 + X.5, 
                 data=dat, family=binomial, 
                 chains=1, iter=500)
  cf3 <- coef(f3)
  expect_equal(length(cf3), 12)
  expect_true(sum(abs(cf3)) > 0)
  
  # multresponse_semipar_array_mixed
  f4 <- bayesGAM(cbind(y1cont, y2cont) ~ 
                   np(X.1, X.2), 
                 data=dat, family=gaussian, 
                 chains=1, iter=500)
  cf4 <- coef(f4)
  expect_equal(length(cf4), 62)
  expect_true(sum(abs(cf4)) > 0)
  
  # multresponse_semipar_array_mixed_discrete
  f5 <- bayesGAM(cbind(y1, y2) ~ X.3 + 
                   np(X.1, X.2), 
                 data=dat, family=binomial, 
                 chains=1, iter=500)
  cf5 <- coef(f5)
  expect_equal(length(cf5), 60)
  expect_true(sum(abs(cf5)) > 0)
  
  # multresponse_semipar_array_mixed_randomint_discrete
  f6 <- bayesGAM(cbind(y1, y2) ~ X.3 + 
                   np(X.1, X.2), 
                 random = ~idnum, 
                 data=dat, family=binomial, 
                 chains=1, iter=500)
  cf6 <- coef(f6)
  expect_equal(length(cf6), 113)
  expect_true(sum(abs(cf6)) > 0)
  
  # random intercept only
  f7 <- bayesGAM(y1cont ~ X.1 + X.2 + X.3, 
                 random = ~idnum, 
                 family=gaussian, data=dat, 
                 chains=1, iter=500)
  cf7 <- coef(f7)
  expect_equal(length(cf7), 26)
  expect_true(sum(abs(cf7)) > 0)
  
  # multresponse_semipar_array_mixed
  f8 <- bayesGAM(cbind(y1, y2) ~ X.3 + 
                   np(X.1, X.2), 
                 data=dat, family=gaussian, 
                 chains=1, iter=500)
  cf8 <- coef(f8)
  expect_equal(length(cf8), 64)
  expect_true(sum(abs(cf8)) > 0)
  
})