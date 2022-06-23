library(readxl)
library(splines)
library(haven)
library(rjags)
library(mcmcplots)
library(runjags)
library(ggplot2)
library(coda)
library(R2OpenBUGS)
library(dplyr)

df <- read_excel("df.xlsx")

df$Time_centered <- scale(df$Time, scale = FALSE)

bspline <-function(X., XL., XR., NDX., BDEG.){
  dx <- (XR. - XL.)/NDX.
  knots <- seq(XL. - BDEG.*dx, XR. + BDEG.*dx, by=dx)
  B <- spline.des(knots, X., BDEG.+1, 0*X.)$design
  B
}

MM.basis <- function (x, xl, xr, ndx, bdeg, pord, decom = 2) {
  B = bspline(x,xl,xr,ndx,bdeg)
  m = ncol(B)
  n = nrow(B)
  D = diff(diag(m), differences=pord)
  P.svd = svd(crossprod(D))
  U = (P.svd$u)[,1:(m-pord)] # eigenvectors
  d = (P.svd$d)[1:(m-pord)] # eigenvalues
  Delta = diag(1/sqrt(d))
  Z = B%*%U%*%Delta
  X = NULL
  for(i in 0:(pord-1)){
    X = cbind(X,x^i)
  }
  list(X = X, Z = Z, d = d, B = B)
}

basis_lmm <- MM.basis(df$Time_centered, min(df$Time_centered)-0.5, max(df$Time_centered)+0.5, 39, 3, 2)
Z <- test$Z
X <- test$X


Nsubjects <- length(levels(as.factor(df$ID)))
num_data <- nrow(df)

first_timepoints <- c(1, 603, 1236, 1828, 2459, 3065, 3606, 4244, 4882, 5518,
                      6076, 6695, 7326, 7954, 8570, 9209, 9843, 10457, 11086, 11578,
                      12119, 12721, 13296, 13929, 14551, 15179, 15810, 16422, 17485)
not_first_timepoints <- 1:num_data
not_first_timepoints <- not_first_timepoints[! not_first_timepoints %in% first_timepoints]


jagsdata_model <- list(response_pos = df$log_positives,
                             response_neg = df$log_negatives,
                             X = X,
                             Z = Z,
                             n = num_data,
                             num.knots = 40,
                             nsubjects = Nsubjects,
                             subject = df$ID,
                             first_timepoints = first_timepoints,
                             not_first_timepoints = not_first_timepoints)
#JAGS Model
model_jags <- function(){
  for (k in 1:n)
  {
    response_pos[k]~dnorm(m_pos[k],taueps_pos[subject[k]])
    response_neg[k]~dnorm(m_neg[k],taueps_neg[subject[k]])
    
    std_resid_neg[k] <- (response_neg[k] - m_neg[k])/sigmaeps_neg[subject[k]]
    std_resid_pos[k] <- (response_pos[k] - m_pos[k])/sigmaeps_pos[subject[k]]
  }
  
  #Likelihood of the model for the first timepoint of each country
  for (k in first_timepoints)
  {
    m_pos[k]<-f_pos[k]+fi_pos[k] 
    m_neg[k]<-f_neg[k]+fi_neg[k] 
    
    f_pos[k]<-b[1]*Z[k,1]+b[2]*Z[k,2]+b[3]*Z[k,3]+
      b[4]*Z[k,4]+b[5]*Z[k,5]+b[6]*Z[k,6]+
      b[7]*Z[k,7]+b[8]*Z[k,8]+b[9]*Z[k,9]+
      b[10]*Z[k,10]+b[11]*Z[k,11]+b[12]*Z[k,12]+
      b[13]*Z[k,13]+b[14]*Z[k,14]+b[15]*Z[k,15]+
      b[16]*Z[k,16]+b[17]*Z[k,17]+b[18]*Z[k,18]+
      b[19]*Z[k,19]+b[20]*Z[k,20]+b[21]*Z[k,21]+
      b[22]*Z[k,22]+b[23]*Z[k,23]+b[24]*Z[k,24]+
      b[25]*Z[k,25]+b[26]*Z[k,26]+b[27]*Z[k,27]+
      b[28]*Z[k,28]+b[29]*Z[k,29]+b[30]*Z[k,30]+
      b[31]*Z[k,31]+b[32]*Z[k,32]+b[33]*Z[k,33]+
      b[34]*Z[k,34]+b[35]*Z[k,35]+b[36]*Z[k,36]+
      b[37]*Z[k,37]+b[38]*Z[k,38]+b[39]*Z[k,39]+
      b[40]*Z[k,40]
    
    fi_pos[k]<-(delta[subject[k],1])*X[k,1]+(delta[subject[k],2])*X[k,2]+
      d[subject[k],1]*Z[k,1]+d[subject[k],2]*Z[k,2]+d[subject[k],3]*Z[k,3]+
      d[subject[k],4]*Z[k,4]+d[subject[k],5]*Z[k,5]+d[subject[k],6]*Z[k,6]+
      d[subject[k],7]*Z[k,7]+d[subject[k],8]*Z[k,8]+d[subject[k],9]*Z[k,9]+
      d[subject[k],10]*Z[k,10]+d[subject[k],11]*Z[k,11]+d[subject[k],12]*Z[k,12]+
      d[subject[k],13]*Z[k,13]+d[subject[k],14]*Z[k,14]+d[subject[k],15]*Z[k,15]+
      d[subject[k],16]*Z[k,16]+d[subject[k],17]*Z[k,17]+d[subject[k],18]*Z[k,18]+
      d[subject[k],19]*Z[k,19]+d[subject[k],20]*Z[k,20]+d[subject[k],21]*Z[k,21]+
      d[subject[k],22]*Z[k,22]+d[subject[k],23]*Z[k,23]+d[subject[k],24]*Z[k,24]+
      d[subject[k],25]*Z[k,25]+d[subject[k],26]*Z[k,26]+d[subject[k],27]*Z[k,27]+
      d[subject[k],28]*Z[k,28]+d[subject[k],29]*Z[k,29]+d[subject[k],30]*Z[k,30]+
      d[subject[k],31]*Z[k,31]+d[subject[k],32]*Z[k,32]+d[subject[k],33]*Z[k,33]+
      d[subject[k],34]*Z[k,34]+d[subject[k],35]*Z[k,35]+d[subject[k],36]*Z[k,36]+
      d[subject[k],37]*Z[k,37]+d[subject[k],38]*Z[k,38]+d[subject[k],39]*Z[k,39]+
      d[subject[k],40]*Z[k,40]
    
    f_neg[k]<-a[1]*Z[k,1]+a[2]*Z[k,2]+a[3]*Z[k,3]+
      a[4]*Z[k,4]+a[5]*Z[k,5]+a[6]*Z[k,6]+
      a[7]*Z[k,7]+a[8]*Z[k,8]+a[9]*Z[k,9]+
      a[10]*Z[k,10]+a[11]*Z[k,11]+a[12]*Z[k,12]+
      a[13]*Z[k,13]+a[14]*Z[k,14]+a[15]*Z[k,15]+
      a[16]*Z[k,16]+a[17]*Z[k,17]+a[18]*Z[k,18]+
      a[19]*Z[k,19]+a[20]*Z[k,20]+a[21]*Z[k,21]+
      a[22]*Z[k,22]+a[23]*Z[k,23]+a[24]*Z[k,24]+
      a[25]*Z[k,25]+a[26]*Z[k,26]+a[27]*Z[k,27]+
      a[28]*Z[k,28]+a[29]*Z[k,29]+a[30]*Z[k,30]+
      a[31]*Z[k,31]+a[32]*Z[k,32]+a[33]*Z[k,33]+
      a[34]*Z[k,34]+a[35]*Z[k,35]+a[36]*Z[k,36]+
      a[37]*Z[k,37]+a[38]*Z[k,38]+a[39]*Z[k,39]+
      a[40]*Z[k,40]
    
    fi_neg[k]<-(delta[subject[k],3])*X[k,1]+(delta[subject[k],4])*X[k,2]+
      c[subject[k],1]*Z[k,1]+c[subject[k],2]*Z[k,2]+c[subject[k],3]*Z[k,3]+
      c[subject[k],4]*Z[k,4]+c[subject[k],5]*Z[k,5]+c[subject[k],6]*Z[k,6]+
      c[subject[k],7]*Z[k,7]+c[subject[k],8]*Z[k,8]+c[subject[k],9]*Z[k,9]+
      c[subject[k],10]*Z[k,10]+c[subject[k],11]*Z[k,11]+c[subject[k],12]*Z[k,12]+
      c[subject[k],13]*Z[k,13]+c[subject[k],14]*Z[k,14]+c[subject[k],15]*Z[k,15]+
      c[subject[k],16]*Z[k,16]+c[subject[k],17]*Z[k,17]+c[subject[k],18]*Z[k,18]+
      c[subject[k],19]*Z[k,19]+c[subject[k],20]*Z[k,20]+c[subject[k],21]*Z[k,21]+
      c[subject[k],22]*Z[k,22]+c[subject[k],23]*Z[k,23]+c[subject[k],24]*Z[k,24]+
      c[subject[k],25]*Z[k,25]+c[subject[k],26]*Z[k,26]+c[subject[k],27]*Z[k,27]+
      c[subject[k],28]*Z[k,28]+c[subject[k],29]*Z[k,29]+c[subject[k],30]*Z[k,30]+
      c[subject[k],31]*Z[k,31]+c[subject[k],32]*Z[k,32]+c[subject[k],33]*Z[k,33]+
      c[subject[k],34]*Z[k,34]+c[subject[k],35]*Z[k,35]+c[subject[k],36]*Z[k,36]+
      c[subject[k],37]*Z[k,37]+c[subject[k],38]*Z[k,38]+c[subject[k],39]*Z[k,39]+
      c[subject[k],40]*Z[k,40]
    
  }
  
  #Likelihood of the model for the other timepoints for each country
  for (k in not_first_timepoints)
  {
    m_pos[k]<-f_pos[k]+fi_pos[k] + rho.time1[subject[k]] * (response_pos[k-1] - f_pos[k-1] - fi_pos[k-1])
    m_neg[k]<-f_neg[k]+fi_neg[k] + rho.time2[subject[k]] * (response_neg[k-1] - f_neg[k-1] - fi_neg[k-1]) 
    
    f_pos[k]<-b[1]*Z[k,1]+b[2]*Z[k,2]+b[3]*Z[k,3]+
      b[4]*Z[k,4]+b[5]*Z[k,5]+b[6]*Z[k,6]+
      b[7]*Z[k,7]+b[8]*Z[k,8]+b[9]*Z[k,9]+
      b[10]*Z[k,10]+b[11]*Z[k,11]+b[12]*Z[k,12]+
      b[13]*Z[k,13]+b[14]*Z[k,14]+b[15]*Z[k,15]+
      b[16]*Z[k,16]+b[17]*Z[k,17]+b[18]*Z[k,18]+
      b[19]*Z[k,19]+b[20]*Z[k,20]+b[21]*Z[k,21]+
      b[22]*Z[k,22]+b[23]*Z[k,23]+b[24]*Z[k,24]+
      b[25]*Z[k,25]+b[26]*Z[k,26]+b[27]*Z[k,27]+
      b[28]*Z[k,28]+b[29]*Z[k,29]+b[30]*Z[k,30]+
      b[31]*Z[k,31]+b[32]*Z[k,32]+b[33]*Z[k,33]+
      b[34]*Z[k,34]+b[35]*Z[k,35]+b[36]*Z[k,36]+
      b[37]*Z[k,37]+b[38]*Z[k,38]+b[39]*Z[k,39]+
      b[40]*Z[k,40]
    
    fi_pos[k]<-(delta[subject[k],1])*X[k,1]+(delta[subject[k],2])*X[k,2]+
      d[subject[k],1]*Z[k,1]+d[subject[k],2]*Z[k,2]+d[subject[k],3]*Z[k,3]+
      d[subject[k],4]*Z[k,4]+d[subject[k],5]*Z[k,5]+d[subject[k],6]*Z[k,6]+
      d[subject[k],7]*Z[k,7]+d[subject[k],8]*Z[k,8]+d[subject[k],9]*Z[k,9]+
      d[subject[k],10]*Z[k,10]+d[subject[k],11]*Z[k,11]+d[subject[k],12]*Z[k,12]+
      d[subject[k],13]*Z[k,13]+d[subject[k],14]*Z[k,14]+d[subject[k],15]*Z[k,15]+
      d[subject[k],16]*Z[k,16]+d[subject[k],17]*Z[k,17]+d[subject[k],18]*Z[k,18]+
      d[subject[k],19]*Z[k,19]+d[subject[k],20]*Z[k,20]+d[subject[k],21]*Z[k,21]+
      d[subject[k],22]*Z[k,22]+d[subject[k],23]*Z[k,23]+d[subject[k],24]*Z[k,24]+
      d[subject[k],25]*Z[k,25]+d[subject[k],26]*Z[k,26]+d[subject[k],27]*Z[k,27]+
      d[subject[k],28]*Z[k,28]+d[subject[k],29]*Z[k,29]+d[subject[k],30]*Z[k,30]+
      d[subject[k],31]*Z[k,31]+d[subject[k],32]*Z[k,32]+d[subject[k],33]*Z[k,33]+
      d[subject[k],34]*Z[k,34]+d[subject[k],35]*Z[k,35]+d[subject[k],36]*Z[k,36]+
      d[subject[k],37]*Z[k,37]+d[subject[k],38]*Z[k,38]+d[subject[k],39]*Z[k,39]+
      d[subject[k],40]*Z[k,40]
    
    f_neg[k]<-a[1]*Z[k,1]+a[2]*Z[k,2]+a[3]*Z[k,3]+
      a[4]*Z[k,4]+a[5]*Z[k,5]+a[6]*Z[k,6]+
      a[7]*Z[k,7]+a[8]*Z[k,8]+a[9]*Z[k,9]+
      a[10]*Z[k,10]+a[11]*Z[k,11]+a[12]*Z[k,12]+
      a[13]*Z[k,13]+a[14]*Z[k,14]+a[15]*Z[k,15]+
      a[16]*Z[k,16]+a[17]*Z[k,17]+a[18]*Z[k,18]+
      a[19]*Z[k,19]+a[20]*Z[k,20]+a[21]*Z[k,21]+
      a[22]*Z[k,22]+a[23]*Z[k,23]+a[24]*Z[k,24]+
      a[25]*Z[k,25]+a[26]*Z[k,26]+a[27]*Z[k,27]+
      a[28]*Z[k,28]+a[29]*Z[k,29]+a[30]*Z[k,30]+
      a[31]*Z[k,31]+a[32]*Z[k,32]+a[33]*Z[k,33]+
      a[34]*Z[k,34]+a[35]*Z[k,35]+a[36]*Z[k,36]+
      a[37]*Z[k,37]+a[38]*Z[k,38]+a[39]*Z[k,39]+
      a[40]*Z[k,40]
    
    fi_neg[k]<-(delta[subject[k],3])*X[k,1]+(delta[subject[k],4])*X[k,2]+
      c[subject[k],1]*Z[k,1]+c[subject[k],2]*Z[k,2]+c[subject[k],3]*Z[k,3]+
      c[subject[k],4]*Z[k,4]+c[subject[k],5]*Z[k,5]+c[subject[k],6]*Z[k,6]+
      c[subject[k],7]*Z[k,7]+c[subject[k],8]*Z[k,8]+c[subject[k],9]*Z[k,9]+
      c[subject[k],10]*Z[k,10]+c[subject[k],11]*Z[k,11]+c[subject[k],12]*Z[k,12]+
      c[subject[k],13]*Z[k,13]+c[subject[k],14]*Z[k,14]+c[subject[k],15]*Z[k,15]+
      c[subject[k],16]*Z[k,16]+c[subject[k],17]*Z[k,17]+c[subject[k],18]*Z[k,18]+
      c[subject[k],19]*Z[k,19]+c[subject[k],20]*Z[k,20]+c[subject[k],21]*Z[k,21]+
      c[subject[k],22]*Z[k,22]+c[subject[k],23]*Z[k,23]+c[subject[k],24]*Z[k,24]+
      c[subject[k],25]*Z[k,25]+c[subject[k],26]*Z[k,26]+c[subject[k],27]*Z[k,27]+
      c[subject[k],28]*Z[k,28]+c[subject[k],29]*Z[k,29]+c[subject[k],30]*Z[k,30]+
      c[subject[k],31]*Z[k,31]+c[subject[k],32]*Z[k,32]+c[subject[k],33]*Z[k,33]+
      c[subject[k],34]*Z[k,34]+c[subject[k],35]*Z[k,35]+c[subject[k],36]*Z[k,36]+
      c[subject[k],37]*Z[k,37]+c[subject[k],38]*Z[k,38]+c[subject[k],39]*Z[k,39]+
      c[subject[k],40]*Z[k,40]
    
  }
  
  
  #Prior for the random parameters of the overall curve
  for (k in 1:num.knots){b[k]~dnorm(0,taub)}
  for (k in 1:num.knots){a[k]~dnorm(0,taua)}
  
  #Prior for the random parameters for the individual deviations
  #from the overall curve
  for (i in 1:nsubjects)
  {for (k in 1:num.knots){d[i,k]~dnorm(0,taud)}}
  for (i in 1:nsubjects)
  {for (k in 1:num.knots){c[i,k]~dnorm(0,tauc)}}
  
  #Prior for monomial parameters of the overall curve
  for (l in 1:4){beta[l]~dnorm(0,0.0001)}
  
  #Prior for monomial parameters of curves describing the individual
  #deviations from the group curve
  # for (i in 1:nsubjects)
  # {for (j in 1:4){delta[i,j]~dnorm(0,taudelta[j])}}
  for (j in 1:nsubjects){
    delta[j,1:4]~dmnorm(beta[1:4], inv_D[1:4, 1:4])
  }
  
  for (k in 1:4) {
    for (l in 1:4) {
      corr[k,l] <- sigma2delta[k,l] / (sigmadelta[k] * sigmadelta[l])
    }
  }
  
  sigma2delta[1:4,1:4]  <- inverse(inv_D[1:4,1:4])
  
  for (i in 1:4) {sigmadelta[i] <- sqrt(sigma2delta[i,i])}
  
  inv_D ~ dwish(4*priorR_D[, ], 5)
  for (l in 1:4) {
    priorR_D[l, l] ~ dgamma(0.5, 0.001)
  }
  
  priorR_D[1, 2] <- 0
  priorR_D[1, 3] <- 0
  priorR_D[1, 4] <- 0
  priorR_D[2, 1] <- 0
  priorR_D[2, 3] <- 0
  priorR_D[2, 4] <- 0
  priorR_D[3, 1] <- 0
  priorR_D[3, 2] <- 0
  priorR_D[3, 4] <- 0
  priorR_D[4, 1] <- 0
  priorR_D[4, 2] <- 0
  priorR_D[4, 3] <- 0
  
  
  
  #Priors of precision parameters
  for (k in 1:nsubjects) {
    taueps_pos[k] <- pow(sigmaeps_pos[k], -2)
    taueps_neg[k] <- pow(sigmaeps_neg[k], -2)
    vareps_pos[k] <- pow(sigmaeps_pos[k], 2)
    vareps_neg[k] <- pow(sigmaeps_neg[k], 2)
    log(sigmaeps_pos[k]) <- u.pos[k]
    log(sigmaeps_neg[k]) <- u.neg[k]
  }
  
  for (k in 1:nsubjects)
  {
    u.pos[k] ~ dnorm(logsigma_respos, tauvarpos)
    u.neg[k] ~ dnorm(logsigma_resneg, tauvarneg)
  }
  
  logsigma_respos ~ dnorm(0, 0.1)
  logsigma_resneg ~ dnorm(0, 0.1)
  
  tauvarpos ~ dgamma(0.001,0.001) 
  tauvarneg ~ dgamma(0.001,0.001) 
  
  varvarpos <- 1/tauvarpos
  varvarneg <- 1/tauvarneg
  
  #prior autocorrelation parameters  
  for (k in 1:nsubjects) {
    rho.time1[k] <- exp(- exp(- lrho.time1[k]))
    rho.time2[k] <- exp(- exp(- lrho.time2[k]))
    
    lrho.time1[k] <- u.rho1[k]
    lrho.time2[k] <- u.rho2[k]
    
  }
  
  for (k in 1:nsubjects)
  {
    u.rho1[k] ~ dnorm(loglog.rho.time1, tauvarrho1)
    u.rho2[k] ~ dnorm(loglog.rho.time2, tauvarrho2)
  }
  
  loglog.rho.time1 ~ dnorm(0, .1)
  loglog.rho.time2 ~ dnorm(0, .1)
  
  tauvarrho1 ~ dgamma(0.001,0.001) 
  tauvarrho2 ~ dgamma(0.001,0.001)
  
  varvarrho1 <- 1/tauvarrho1  
  varvarrho2 <- 1/tauvarrho2 
  
  
  sigma2_b <- 1/taub
  taub ~ dgamma(0.001,0.001)
  
  sigma2_d <- 1/taud
  taud ~ dgamma(0.001,0.001) 
  
  sigma2_a <- 1/taua
  taua ~ dgamma(0.001,0.001) 
  
  sigma2_c <- 1/tauc
  tauc ~ dgamma(0.001,0.001)   
  
  
}

params_model <- c("sigma2_b", "sigma2_d", "sigma2_a", "sigma2_c", 
                        "beta", 
                        "rho.time1", "rho.time2", "loglog.rho.time1", "loglog.rho.time2",
                        "tauvarrho1", "tauvarrho2", "tauvarpos", "tauvarneg",
                        "varvarrho1", "varvarrho2", "varvarpos", "varvarneg",
                        "sigma2delta",
                        "logsigma_respos", "logsigma_resneg", "vareps_pos", "vareps_neg",
                        "corr",
                        "m_pos", "m_neg")

write.model(model_jags, "model_penalized_splines_AR1.bug")

results_model <- run.jags(model="model_penalized_splines_AR1.bug", monitor=params_model,
                          data=jagsdata_model, n.chains=4, method="parallel", jags.refresh = 30,
                          burnin = 100000, adapt = 1000, sample = 50000, thin = 1,
                          modules = c("glm on"), summarise = TRUE, plots = FALSE)


model_summary <- add.summary(results_model)
model_summary

model_mcmc <- as.mcmc.list(results_model)
summary(model_mcmc)
traplot(model_mcmc, parms=c("tauvarrho1", "tauvarrho2", "tauvarpos", "tauvarneg"))
traplot(model_mcmc, parms=c("varvarrho1", "varvarrho2", "varvarpos", "varvarneg"))
traplot(model_mcmc, parms=c("loglog.rho.time1", "loglog.rho.time2"))
traplot(model_mcmc, parms=c("sigma2_b", "sigma2_d", "sigma2_a", "sigma2_c"))
traplot(model_mcmc, parms=c("beta"))

#obtain estimated parameters (variances and autocorrelation for every country)
s <- summary(model_mcmc)
s
m <- s$statistics[,"Mean"]
m
nrow(as.data.frame(m))

rho.time_pos <- as.data.frame(m)[9:38,]
rho.time_neg <- as.data.frame(m)[39:68,]

var_inn_pos <- as.data.frame(m)[89:118,]
var_inn_neg <- as.data.frame(m)[119:148,]

df$rho_pos <- 0
df$rho_neg <- 0

df$var_pos <- 0
df$var_neg <- 0


for(i in 1:nrow(df)) {
  df$rho_pos[i] <- rho.time_pos[df$ID[i]]
  df$rho_neg[i] <- rho.time_neg[df$ID[i]]
  
  df$var_pos[i] <- var_inn_pos[df$ID[i]]
  df$var_neg[i] <- var_inn_neg[df$ID[i]]
}


df$eps_pos <- df$var_pos / (1 - df$rho_pos^2)
df$eps_neg <- df$var_neg / (1 - df$rho_neg^2)

#scatterplot of rho and var for each country
country_specific_effects <- df %>% select(Location, ID, rho_pos, rho_neg, 
                                           var_pos, var_neg, eps_pos, eps_neg)
country_specific_effects <- countrie_specific_effects %>% distinct(Location, ID, rho_pos, rho_neg, 
                                                                   var_pos, var_neg, eps_pos, eps_neg)

ggplot(country_specific_effects, aes(x=var_neg, y=rho_neg)) +
  geom_point() + 
  geom_text(label=country_specific_effects$Location)

ggplot(country_specific_effects, aes(x=var_pos, y=rho_pos)) +
  geom_point() + 
  geom_text(label=country_specific_effects$Location)

country_specific_effects$Location <- factor(country_specific_effects$Location, levels = country_specific_effects$Location[order(country_specific_effects$rho_neg)])

country_specific_effects %>%
  mutate(name = fct_reorder(Location, rho_neg)) %>% ggplot( aes(x=Location, y=rho_neg)) +
  geom_point() +
  coord_flip() +
  xlab("") +
  theme_bw() +
  ylab("Autocorrelation parameter value for number of negative tests")

country_specific_effects$Location <- factor(country_specific_effects$Location, levels = country_specific_effects$Location[order(country_specific_effects$rho_pos)])

country_specific_effects %>%
  mutate(name = fct_reorder(Location, rho_neg)) %>% ggplot( aes(x=Location, y=rho_pos)) +
  geom_point() +
  coord_flip() +
  xlab("") +
  theme_bw() +
  ylab("Autocorrelation parameter value for number of positive tests")

country_specific_effects$Location <- factor(country_specific_effects$Location, levels = country_specific_effects$Location[order(country_specific_effects$var_neg)])

country_specific_effects %>%
  mutate(name = fct_reorder(Location, rho_neg)) %>% ggplot( aes(x=Location, y=var_neg)) +
  geom_point() +
  coord_flip() +
  xlab("") +
  theme_bw() +
  ylab("Variance parameter value for number of negative tests")

country_specific_effects$Location <- factor(country_specific_effects$Location, levels = country_specific_effects$Location[order(country_specific_effects$var_pos)])

country_specific_effects %>%
  mutate(name = fct_reorder(Location, rho_neg)) %>% ggplot( aes(x=Location, y=var_pos)) +
  geom_point() +
  coord_flip() +
  xlab("") +
  theme_bw() +
  ylab("Variance parameter value for number of positive tests")


df$mu_pos <- as.data.frame(m)[169:18255,]
df$mu_neg <- as.data.frame(m)[18256:36342,]


df$resid_neg <- df$log_negatives - df$mu_neg
df$resid_pos <- df$log_positives - df$mu_pos

df$std_resid_neg <- df$resid_neg / (sqrt(df$eps_neg))
df$std_resid_pos <- df$resid_pos / (sqrt(df$eps_pos))

#observed log-scale
colors <- c("Positive tests" = "#387a96", "Negative tests" = "#b32550")
ggplot(df, aes(x = Time, group = Location)) +
  geom_line(aes(y = log_positives, color = "Positive tests")) +
  geom_line(aes(y = log_negatives, color = "Negative tests")) +
  labs(x = "Year", y = "(%)", color = "Legend") +
  scale_color_manual(values = colors) +
  scale_x_continuous("Day") +
  scale_y_continuous("Observed number of tests (log-scale)", expand = c(0, 0), breaks = c(0,5,10)) +
  coord_cartesian(ylim = c(-0.5, 13)) +
  facet_wrap( ~ factor(Location), ncol = 5) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey96"),
    axis.line.x = element_line(colour = "black"),
    axis.title = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm")
  ) 

#observed versus predicted
ggplot(df, aes(x = Time, group = Location)) +
  geom_line(aes(y = log_positives, color = "Positive tests"),
            linetype = "dashed") +
  geom_line(aes(y = mu_pos, color = "Positive tests")) +
  geom_line(aes(y = log_negatives, color = "Negative tests"),
            linetype = "dashed") +
  geom_line(aes(y = mu_neg, color = "Negative tests")) +
  labs(x = "Year", y = "(%)", color = "Legend") +
  scale_color_manual(values = colors) +
  scale_x_continuous("Day") +
  scale_y_continuous("Number of tests (log-scale)", expand = c(0, 0)) +
  coord_cartesian(ylim = c(-0.5, 13)) +
  facet_wrap( ~ factor(Location), ncol = 5) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey96"),
    axis.line.x = element_line(colour = "black"),
    axis.title = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm")
  ) +
  ggtitle("Observed (dashed line) and estimated trajectories (solid line)")

#residuals
ggplot(df, aes(x = Time, group = Location)) +
  geom_point(aes(y = std_resid_pos, color = "Positive tests"), size=0.5) +
  geom_point(aes(y = std_resid_neg, color = "Negative tests"), size=0.5) +
  labs(x = "Year", y = "(%)", color = "Legend") +
  scale_color_manual(values = colors) +
  scale_x_continuous("Day") +
  scale_y_continuous("Residuals", expand = c(0, 0)) +
  coord_cartesian(ylim = c(-1, 1)) +
  facet_wrap( ~ factor(Location), ncol = 5) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey96"),
    axis.line.x = element_line(colour = "black"),
    axis.title = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm")
  ) +
  ggtitle("Residuals")

ggplot(df, aes(x = mu_pos, group = Location)) +
  geom_point(aes(y = std_resid_pos, color = "Positive tests"), size=0.5) +
  labs(x = "Year", y = "(%)", color = "Legend") +
  scale_color_manual(values = colors) +
  scale_x_continuous("Fitted values") +
  scale_y_continuous("Residuals", expand = c(0, 0)) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_wrap( ~ factor(Location), ncol = 5) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey96"),
    axis.line.x = element_line(colour = "black"),
    axis.title = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.position = "none",
    legend.key.size = unit(0.3, "cm")
  ) 

ggplot(df, aes(x = mu_neg, group = Location)) +
  geom_point(aes(y = std_resid_neg, color = "Negative tests"), size=0.5) +
  labs(x = "Year", y = "(%)", color = "Legend") +
  scale_color_manual(values = colors) +
  scale_x_continuous("Fitted values") +
  scale_y_continuous("Residuals", expand = c(0, 0)) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_wrap( ~ factor(Location), ncol = 5) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey96"),
    axis.line.x = element_line(colour = "black"),
    axis.title = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.position = "None",
    legend.key.size = unit(0.3, "cm")
  ) 

# plot acf with ggplot
ggacf <- function(series, main) {
  #significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(series)))  
  a<-acf(series, plot=F, lag.max = 25)
  a.2<-with(a, data.frame(lag, acf))
  g<- ggplot(a.2[-1,], aes(x=lag,y=acf)) + 
    geom_bar(stat = "identity", position = "identity", width = 0.1) + xlab('Lag') + ylab('ACF') +
    ggtitle(main) +
    ylim(-1, 1)
  #geom_hline(yintercept=c(significance_level,-significance_level), lty=3);
  
  # fix scale for integer lags
  if (all(a.2$lag%%1 == 0)) {
    g<- g + scale_x_discrete(limits = seq(1, max(a.2$lag)));
  }
  return(g);
}

library(dplyr)
group1 <- df %>% filter(ID == 1)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 2)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 3)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 4)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 5)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 6)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 7)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 8)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 9)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 10)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 11)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 12)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 13)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 14)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 15)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 16)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 17)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 18)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 19)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 20)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 21)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 22)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 23)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 24)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 25)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 26)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 27)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 28)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 29)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

group1 <- df %>% filter(ID == 30)
ggacf(group1$std_resid_neg_eps, main=paste(group1$Location[1]))
ggacf(group1$std_resid_pos_eps, main=paste(group1$Location[1]))

