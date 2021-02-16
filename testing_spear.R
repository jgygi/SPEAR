library(Rcpp)
sourceCpp("/Users/lg689/Dropbox/COVID19/shared_folder/SPEARcomplete/src/_spear.cpp")
source("/Users/lg689/Dropbox/COVID19/shared_folder/SPEARcomplete/R/helpers.R")
source("/Users/lg689/Dropbox/COVID19/shared_folder/SPEARcomplete/R/spear.R")


seed = 2021
N = 100; Ntest = 5000;P = 100; D = 5;  num_factors = 5; print_out = 100; num_specific = 2; 

max_iter = 1000;  thres_elbo = 0.1; thres_count = 5; pi = .2; 
warm_up  = 50; weights = seq(0, 2, length.out = 20); weights = weights[c(length(weights):1)]

data = dataGen_gassian(N = N, Ntest = Ntest, P = P, D = D, seed = seed, num_factors = num_factors,
                       c = sqrt(1 *log(P*D)/N), pi = pi, eta = 1, num_specific =num_specific, Ymodel = "XY",pi_reg =N/(P*D*log(P * D) * log(P * D)))

X = data$data.tr$X; Y = data$data.tr$Y; 
pattern_samples = data$data.tr$pattern_samples
pattern_features = data$data.tr$pattern_features
functional_path = data$data.tr$functional_path;
ws =weights; nclasses = data$data.tr$nclasses

a0 = 1e-2; b0 = 1e-2; a1 = sqrt(nrow(X)); b1 = sqrt(nrow(X));
thres_elbo = 0.01; thres_count = 5; thres_factor = 1e-8;print_out = 10;
a2= sqrt(nrow(X)); b2 = sqrt(nrow(X)); inits_post_mu = NULL;seed = 1; inits_type = "pca";
family = 0

Xobs = apply(X, c(1,2), function(z) ifelse(is.na(z), 0,1))
Yobs = apply(Y, c(1,2), function(z) ifelse(is.na(z), 0,1))
Z = X

tmp = spear(X = X, Xobs = Xobs, Y = Y,  Yobs = Yobs, Z = Z, ws = c(1, 0.5, 0),
            family = family,  nclasses = nclasses, functional_path = functional_path,
            pattern_samples = pattern_samples, pattern_features = pattern_features,
            num_factors = num_factors,  print_out = print_out, warm_up=warm_up,
            max_iter = max_iter, seed = seed, thres_elbo = thres_elbo,
                     thres_count = thres_count)

meanFactors = matrix(0, ncol = num_factors, nrow = nrow(Y))
for(k in 1:length(pattern_samples)){
  ii = pattern_samples[[k]]
  jj = pattern_features[[k]]
  post_betas = tmp$post_betas[,,k,3]
  meanFactors[ii,] = Z[ii,jj]%*%post_betas[jj,]
}
plot(Y~meanFactors[,1])

require(parallel)
foldid = sample(1:5, N, replace = T)
tmp_cv <- cv.spear(X = X, Xobs = Xobs, Y = Y,  Yobs = Yobs, Z = Z, ws =weights,
                   family = 0,  nclasses = nclasses, functional_path = functional_path, foldid = foldid,
                   pattern_samples = pattern_samples, pattern_features = pattern_features,
                   num_factors = num_factors,  print_out = print_out, warm_up=warm_up,
                   max_iter = max_iter, seed = seed, thres_elbo = thres_elbo,
                   thres_count = thres_count, crossYonly = F, numCores = NULL)

tmp_eval <- cv.evaluation(fitted.obj = tmp_cv, X = X, Y = Y, Z = Z, family = 0, nclasses = nclasses, 
                          pattern_samples = pattern_samples, 
                          pattern_features = pattern_features, nlambda = 20)

tmp_eval$cvm

F0 =cbind(X%*%tmp_cv$results$post_betas[,1,1,20],X%*%tmp_cv$factors_coefs[,1,1,,20])
cor(F0)

require(glmnet)
lasso.fit <-cv.glmnet(X, Y, foldid = foldid)



seed = 2021
N = 200; Ntest = 5000;P = 100; D = 4;  num_factors = 5; print_out = 100; num_specific = 2; family = 2

max_iter = 1000;  thres_elbo = 0.1; thres_count = 5; pi = .2; 
warm_up  = 50; weights = seq(0, 2, length.out = 20); weights = weights[c(length(weights):1)]

data = dataGen_ordinal(N = N, Ntest =Ntest, P =P, levels = 4,
                       D = D, seed = 2020, num_factors = 5, c = sqrt(0 *log(P*D)/N), 
                       pi = 0.2, eta = 5, num_specific =D-2, Ymodel = "XY",
                       pi_reg = 0.05)

X = data$data.tr$X; Y = data$data.tr$Y; 
pattern_samples = data$data.tr$pattern_samples
pattern_features = data$data.tr$pattern_features
functional_path = data$data.tr$functional_path;
ws =weights; nclasses = data$data.tr$nclasses

a0 = 1e-2; b0 = 1e-2; a1 = sqrt(nrow(X)); b1 = sqrt(nrow(X));
thres_elbo = 0.1; thres_count = 5; thres_factor = 1e-8;print_out = 10;
a2= sqrt(nrow(X)); b2 = sqrt(nrow(X)); inits_post_mu = NULL;seed = 1; inits_type = "pca";

Xobs = apply(X, c(1,2), function(z) ifelse(is.na(z), 0,1))
Yobs = apply(Y, c(1,2), function(z) ifelse(is.na(z), 0,1))
Z = X

sourceCpp("/Users/lg689/Dropbox/COVID19/shared_folder/SPEARcomplete/src/_spear.cpp")
source("/Users/lg689/Dropbox/COVID19/shared_folder/SPEARcomplete/R/helpers.R")
source("/Users/lg689/Dropbox/COVID19/shared_folder/SPEARcomplete/R/spear.R")

family = 0
tmp = spear(X = X, Xobs = Xobs, Y = Y,  Yobs = Yobs, Z = Z, ws = c(0),
            family = family,  functional_path = functional_path, nclasses = nclasses,
            pattern_samples = pattern_samples, pattern_features = pattern_features,
            num_factors = num_factors,  print_out = print_out, warm_up=warm_up,
            max_iter = max_iter, seed = seed, thres_elbo = thres_elbo,
            thres_count = thres_count)
tmp$interceptsY
meanFactors = matrix(0, ncol = num_factors, nrow = nrow(Y))
for(k in 1:length(pattern_samples)){
  ii = pattern_samples[[k]]
  jj = pattern_features[[k]]
  post_betas = tmp$post_betas[,,k,1]
  meanFactors[ii,] = Z[ii,jj]%*%post_betas[jj,]
}
par(mfrow = c(1,2))
boxplot(meanFactors[,1]~Y)
plot(meanFactors[,1]~data$data.tr$truth)

cor(meanFactors[,1], data$data.tr$truth)
require(parallel)

foldid = sample(1:5, N, replace = T)
tmp_cv <- cv.spear(X = X, Xobs = Xobs, Y = Y,  Yobs = Yobs, Z = Z, ws = weights,
                   family = 0,  nclasses = nclasses, functional_path = functional_path, foldid = foldid,
                   pattern_samples = pattern_samples, pattern_features = pattern_features,
                   num_factors = num_factors,  print_out = print_out, warm_up=warm_up,
                   max_iter = max_iter, seed = seed, thres_elbo = thres_elbo,
                   thres_count = thres_count, crossYonly = F, numCores = NULL)


tmp_eval <- cv.evaluation(fitted.obj = tmp_cv, X = X, Y = Y, Z = Z, family = 0, nclasses = nclasses, 
                          pattern_samples = pattern_samples, 
                          pattern_features = pattern_features, nlambda = 100)

tmp_eval$cvm
yhat.tr = rep(0, N)
yhat.te = rep(0, Ntest)
for(k in 1:num_patterns){
  ii = pattern_samples[[k]]
  yhat.tr[ii] = Z%*%tmp_eval$reg_coefs[,k,1,which.min(tmp_eval$cvm)]+tmp_eval$intercepts[[1]][which.min(tmp_eval$cvm),1]
}
yhat.te = data$data.te$X%*%tmp_eval$reg_coefs[,k,1,which.min(tmp_eval$cvm)] + tmp_eval$intercepts[[1]][which.min(tmp_eval$cvm),1]

cor(yhat.te, data$data.te$Y, method = "spearman")

lasso.fit = cv.glmnet(X, Y, foldid = foldid)

yhat.lasso.tr = predict(lasso.fit, data$data.tr$X, s = "lambda.min")

yhat.lasso.te = predict(lasso.fit, data$data.te$X, s = "lambda.min")

par(mfrow = c(2,2))
boxplot(yhat.tr~data$data.tr$Y)
boxplot(yhat.lasso.tr~data$data.tr$Y)
boxplot(yhat.te~data$data.te$Y)
boxplot(yhat.lasso.te~data$data.te$Y)



boxplot(yhat.lasso.te~data$data.te$Y)
boxplot(yhat.te~data$data.te$Y)
