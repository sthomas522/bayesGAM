test_that("hmc testing", {
  require(stats); require(graphics); require(SemiPar)
  # Linear regression example
  # linear: identity
  f1 <- mvgamHMC(weight ~ np(height), data = women, iter=500,
                   family = gaussian(link="identity"), 
                 chains=1, seed=521)
  # cf1 <- as.vector(round(coef(f1), 5))
  # chk1 <- c( 1.29652, 0.24450, -80.02591, 3.30562, 
  #            0.02730, 0.02801 , -0.03589, 0.15149)
  # expect_equal(cf1, chk1)
  cf1 <- coef(f1)
  expect_equal(length(cf1), 8)
  expect_true(all(cf1 != 0))
  
  # linear: log
  f2 <- mvgamHMC(weight ~ np(height), data = women, iter=500,
                 family = gaussian(link="log"), 
                 chains=1, beta=st(c(35, 0, 4)), seed=521)
  # cf2 <- as.vector(round(coef(f2), 5))
  # chk2 <- c(0.00494 , 0.25505 , 3.28320, 0.02503, 
  #           0.00005, 0.00011,-0.00041 ,  0.00044)
  # expect_equal(cf2, chk2) 
  cf2 <- coef(f2)
  expect_equal(length(cf2), 8)
  expect_true(all(cf2 != 0))
  
  
  # linear: sqrt
  f3 <- mvgamHMC(weight ~ np(height), data = women, iter=500,
                 family = gaussian(link="sqrt"), 
                 chains=1, seed=521,
                 beta = normal(c(0, 10)))
  # cf3 <- as.vector(round(coef(f3), 5))
  # chk3 <- c(0.71834, 0.73243, 12.94780, -0.13236, 
  #           -0.00408, 0.01758, -0.00070, -0.00021)
  # expect_equal(cf3, chk3)
  cf3 <- coef(f3)
  expect_equal(length(cf3), 8)
  expect_true(all(cf3 != 0))
  
  # binomial models
  set.seed(651)
  X <- matrix(rnorm(100*5), ncol=5)
  y <- sample(x=c(0,1), size=100, replace=TRUE, prob=c(0.5, 0.5))
  
  bdat <- data.frame(y, X)
  
  # binomial: logit
  f4 <- mvgamHMC(y ~ . - 1, data=bdat, chains=1, iter=500,
                 family=binomial(link="logit"),
                    spcontrol = list(qr = TRUE), 
                 seed=521)

  # cf4 <- as.vector(round(coef(f4), 5))
  # chk4 <- c(-0.00923, -0.30050, -0.15579 ,0.07052, 0.11117)
  # expect_equal(cf4, chk4)
  cf4 <- coef(f4)
  expect_equal(length(cf4), 5)
  expect_true(all(cf4 != 0))

  # binomial: logit2 qr=FALSE
  f5 <- mvgamHMC(y ~ . - 1, data=bdat, chains=1, iter=500,
                 family=binomial(link="logit"),
                 spcontrol = list(qr = FALSE), 
                 seed=521)

  # cf5 <- as.vector(round(coef(f5), 5))
  # chk5 <- c(-0.01116 , -0.30224, -0.16636 ,0.06900, 0.11216)
  # expect_equal(cf5, chk5)
  cf5 <- coef(f5)
  expect_equal(length(cf5), 5)
  expect_true(all(cf5 != 0))

  # binomial: probit
  f6 <- mvgamHMC(y ~ . - 1, data=bdat, chains=1, iter=500,
                 family=binomial(link="probit"),
                 spcontrol = list(qr = TRUE), 
                 seed=521)

  # cf6 <- as.vector(round(coef(f6), 5))
  # chk6 <- c(-0.00460, -0.18241,  -0.09817, 0.03969, 0.06609)
  # expect_equal(cf6, chk6)
  cf6 <- coef(f6)
  expect_equal(length(cf6), 5)
  expect_true(all(cf6 != 0))

  # binomial: cauchit
  f7 <- mvgamHMC(y ~ . - 1, data=bdat, chains=1, iter=500,
                 family=binomial(link="cauchit"),
                 spcontrol = list(qr = TRUE), 
                 seed=521)

  # cf7 <- as.vector(round(coef(f7), 5))
  # chk7 <- c(-0.02866, -0.27692, -0.16637,0.06580,0.15567)
  # expect_equal(cf7, chk7)
  cf7 <- coef(f7)
  expect_equal(length(cf7), 5)
  expect_true(all(cf7 != 0))
  

  # binomial: logit3 with . in formula
  dat<- data.frame(abs(X),
                   y=y)
  f8 <- mvgamHMC(y ~ ., data=dat, chains=1, iter=500,
                 family=binomial(link="logit"),
                 spcontrol = list(qr = TRUE), 
                 seed=521)

  # cf8 <- as.vector(round(coef(f8), 5))
  # chk8 <- c(-0.00421 ,-0.10923  ,  0.12493, -0.10314,-0.01567, 0.18157 )
  # expect_equal(cf8, chk8)
  cf8 <- coef(f8)
  expect_equal(length(cf8), 6)
  expect_true(all(cf8 != 0))

  # binomial: log TODO

  # binomial: cloglog
  data("PimaIndiansDiabetes2", package = "mlbench")
  PimaIndiansDiabetes2$y <- ifelse(PimaIndiansDiabetes2$diabetes=="pos",1,0)
  f9 <- mvgamHMC(y ~ pregnant+pedigree+age, chains=1, 
                 data=PimaIndiansDiabetes2, iter=500,
                 family=binomial(link="cloglog"),
                 spcontrol = list(qr = TRUE), seed=521)

  # cf9 <- as.vector(round(coef(f9), 5))
  # chk9 <- c(-2.17828, 0.07697 ,0.79679 ,0.01780)
  # expect_equal(cf9, chk9)
  cf9 <- coef(f9)
  expect_equal(length(cf9), 4)
  expect_true(all(cf9 != 0))

  # Poisson: log
  ## Dobson (1990) Page 93: Randomized Controlled Trial :
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  pdata <- data.frame(counts, outcome, treatment)

  f10 <- mvgamHMC(counts ~ outcome + treatment, 
                  data=pdata, chains=1, iter=500,
                  family = poisson(link="log"),
                       spcontrol = list(qr = TRUE), 
                  seed=521)

  # cf10 <- as.vector(round(coef(f10), 5))
  # chk10 <- c(3.03935, -0.45839,-0.29909,0.00179,-0.00786)
  # expect_equal(cf10, chk10)
  cf10 <- coef(f10)
  expect_equal(length(cf10), 5)
  expect_true(all(cf10 != 0))


  # Poisson: identity
  # counts2 <- sample(0:4, size=9, replace=TRUE)
  # pdata$counts2 <- counts2
  # f11 <- mvgamHMC(counts2 ~ outcome + treatment, data=pdata,
  #                 family = poisson(link="identity"),
  #                 spcontrol = list(qr = TRUE),
  #                 cores=1, chains=1, seed=521)
  # 
  # cf11 <- as.vector(round(coef(f11), 5))
  # chk11 <- c(1.64220, 0.27777 ,1.53792 , 0.07323,-0.54652)
  # expect_equal(cf11, chk11)


  # Poisson: sqrt
  f12 <- mvgamHMC(counts ~ outcome + treatment, 
                  data=pdata, chains=1, iter=500,
                  family = poisson(link="sqrt"),
                  spcontrol = list(qr = TRUE), 
                  seed=521)

  # cf12 <- as.vector(round(coef(f12), 5))
  # chk12 <- c(-1.32924, 0.12008 , 3.86230 ,2.03010 ,-0.16244)
  # expect_equal(cf12, chk12)
  cf12 <- coef(f12)
  expect_equal(length(cf12), 5)
  expect_true(all(cf12 != 0))
  
  
  # bivariate smoothing
  data(scallop)
  set.seed(982)
  f13 <- mvgamHMC(log(tot.catch+1) ~ np(longitude, latitude),
                 data=scallop, cores=1, chains=1, iter=500,
                 seed=521)
  
  # cf13 <- as.vector(round(sum(abs(coef(f13))), 3))
  # chk13 <- round(244.773, 3)
  # expect_equal(cf13, chk13)
  cf13 <- coef(f13)
  expect_equal(length(cf13), 42)
  expect_true(all(cf13 != 0))
  
  
  # # multivariate response
  # f14 <- mvgamHMC(cbind(logvar1, logvar2, logvar3) ~ np(x1, x2) +
  #                  x3 + x4, data = maskdata,
  #                 random = ~ factor(idnum), cores=1, 
  #                 chains=1, iter=100, 
  #                 seed=521)
  # # cf14 <- as.vector(round(sum(abs(coef(f14))), 3))
  # # chk14 <- round(565.785, 3)
  # # expect_equal(cf14, chk14)
  # cf14 <- coef(f14)
  # expect_equal(length(cf14), 198)
  # expect_true(all(cf14 != 0))
  
  # autoregressive
  f15 <- mvgamHMC(lh ~ L(lh), family=gaussian, cores=1, chains=1, 
                  seed=521)
  cf15 <- coef(f15)
  expect_equal(length(cf15), 3)
  expect_true(all(cf15 != 0))
  
  # posterior predict
  set.seed(981)
  pp <- posterior_predict(f13, draws=50)
  ppdim <- dim(pp@pp[[1]])
  expect_equal(ppdim, c(50, 148))
  expect_equal(sum(is.na(pp@pp)), 0)
  
  # chksum <- round(sum(unlist(pp@pp, sum)), 2)
  # expect_equal(chksum, 25741.55)
  
  # predict
  set.seed(432)
  f <- mvgamHMC(weight ~ np(height), data = women,
                family = gaussian, iter=1000)
  newheights <- with(women, rnorm(10, mean = mean(height)), sd=sd(height))
  women2 <- data.frame(height=newheights)
  
  pred <- predict(f, women2, draws=100)
  
  ppdim <- dim(pred@pp[[1]])
  expect_equal(ppdim, c(100, 10))
  expect_equal(sum(is.na(pred@pp)), 0)
  
  # chksum <- round(sum(unlist(pred@pp)), 1)
  # expect_equal(chksum, 135534.9)
  
  # random intercept
  dat <- data.frame(id = rep(1:5, each=10),
                    y1 = rnorm(50),
                    y2 = rexp(50),
                    x = runif(50))
  
  f16 <- mvgamHMC(cbind(y1, y2) ~ np(x), random = ~factor(id), data=dat, 
                chains=1, iter=500)
  
  expect_equal(length(coef(f16)), 55)
})



