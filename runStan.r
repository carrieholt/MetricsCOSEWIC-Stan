#==============================================================================
# Code to run probability of declines over 3 generations using Stan for COSEWIC
# Applied to example stocks of Sockeye Salmon
# Created by: Carrie Holt
# Date last modified: 29 Oct 2021
#==============================================================================


# Package dependencies
library(dplyr)
library(tidyr)
library(rstan)
library(devtools)
library(ggplot2)
install_github("SOLV-Code/MetricsCOSEWIC", dependencies = TRUE,
               build_vignettes = FALSE)
library(MetricsCOSEWIC)

#==============================================================================
# Data
#==============================================================================

# Example data for 22 sockeye stocks are included in MetricsCOSEWIC package
# Quick look at spawners data over the most recent 13 years (3 gens + 1year),
# 2005-2017

SR_Sample %>% filter(Year>(2017-13)) %>%
  ggplot(aes(Year, Spn)) + geom_line(aes(colour=factor(Stock)))


#read in data
data.in <- SR_Sample %>% select(Stock,Year,Spn) %>% rename(DU=Stock,Abd = Spn)
head(data.in)
gen <- 4
yrs.do <- (3 * gen) +1
last.year <- 2017

# Choose a stock to look at: "Stock18", "Stock2" and "Stock3" are examples 
# where posteriors were sensitive to choice of prior (i.e., posteriors from 
# previous implementation of RStanarm with informative priors differed from 
# those when we used the diffuse priors described below)
stk <- "Stock18"
data.in <- data.in %>% filter(DU==stk) %>% mutate(logAbd=log(Abd)) %>%
  filter(Year > (last.year - yrs.do) & Year <= last.year)

# Look at data frame
data.in
ggplot(data.in, aes(Year,logAbd)) + geom_point()



#==============================================================================
# Priors
#==============================================================================

# VERY DIFFUSE PRIORS

# Set priors for intercept, normally distributed: ~N(mean, sig)
intercept_mean <-  median(data.in$logAbd, na.rm=T)
intercept_sig <-  max(data.in$logAbd, na.rm=T)
  # max value in the time-series is an upper bound on plausible sigmas on y-i

# Set priors for slope, normally distributed: ~N(mean, sig)
slope_mean <-  0
slope_sig <-  max(data.in$logAbd, na.rm=T) /
  (max(data.in$Year) - min(data.in$Year))
  # max rise/max run is an upper bound on plausible sigmas on slope

# Set priors for residual variance, inverse gamma distributed:
# ~Inv_Gamma(alpha,beta), where alpha = beta = 'Sig_Gam_Dist'
Sig_Gam_Dist <-  0.00001

#==============================================================================
# Set up data inputs for Stan
#==============================================================================

# test for missing value codes

#data.in$logAbd[c(2)] <- NA

data <- list()
data$Year <- 0:(length(data.in$Year)-1)
data$logAbd <- replace_na(data.in$logAbd, 99) #data.in$logAbd
data$N <- dim(data.in)[1]
data$intercept_mean <- intercept_mean
data$Sig_Gam_Dist <- Sig_Gam_Dist
data$slope_mean <- slope_mean
data$intercept_sig <- intercept_sig
data$slope_sig <- slope_sig

data$logAbd_obs <- as.numeric(!is.na(data.in$logAbd))

                            
#==============================================================================
# Run Stan
#==============================================================================

# stan_fit <- stan(file = 'linear.stan', data = data, iter = 10000,
#                  chains = 6,  control = list(adapt_delta = 0.95))


stan_fit <- stan(file = 'linear-missingValuesv2.stan', data = data, iter = 10000,
                 chains = 6,  control = list(adapt_delta = 0.95))

#==============================================================================
# Pull out parameter estimates
#==============================================================================

# #Plot the distribution for missing values, just to check (only worksfor 1 NA)
# draws_logAbdmis <-
#   as.matrix(stan_fit, pars = "logAbd_Pred")[ , is.na(data.in$logAbd)]
# df<-data.frame(data=draws_logAbdmis)
# ggplot(df, aes(data)) + geom_density()


if (max(bayesplot::rhat(stan_fit)) > 1.02) {
  convergence.flag<- TRUE
  print("Convergence flag")
}


All_Ests <- data.frame(summary(stan_fit)$summary)
All_Ests$Param <- row.names(All_Ests)

logAbd_Fits_Stan <- All_Ests[grepl("logAbd_Fit_out", All_Ests$Param),  ]
logAbd_Preds_Stan <- All_Ests[grepl("logAbd_Pred", All_Ests$Param),  ]


Name <- "Test"
FitsDF <- data.frame(Year = data.in$Year, R = data.in$logAbd, Fit = logAbd_Fits_Stan$X50.,
                     Year = 1:dim(logAbd_Fits_Stan)[1],   Mod = Name,
                     CI_up = logAbd_Fits_Stan$X97.5.,
                     CI_low = logAbd_Fits_Stan$X2.5.,
                     Pred = logAbd_Preds_Stan$X50.,
                     Pred_up = logAbd_Preds_Stan$X97.5.,
                     Pred_low = logAbd_Preds_Stan$X2.5.)


# get  posteriors
fit_values <- rstan::extract(stan_fit)
intercept_Post <- fit_values$intercept
slope_Post <- fit_values$slope


# Return fit and predicted values
out <- list()
out$Fits <- FitsDF
out$Ests <- All_Ests
out$intercept_Post <- intercept_Post
out$slope_Post <- slope_Post

# Calculate % change from the fit at the start and end of the 13 years
mcmc.samples <- data.frame(int=out$intercept_Post, slope=out$slope_Post,
                           Perc_Change = NA,Perc_Change_Raw = NA)
mcmc.samples <- mcmc.samples %>%
  mutate(Fit_Start = int + slope * data$Year[1] ) %>%
  mutate(Fit_End = int + slope * data$Year[yrs.do])

# Identify any years with negative start years (only occurs if abundances are
# scaled prior to log transform)
neg.start.idx <-  mcmc.samples[,"Fit_Start"] < 0


mcmc.samples[,"Perc_Change"][!neg.start.idx] <-
  (exp(mcmc.samples[,"Fit_End"][!neg.start.idx]) - exp(mcmc.samples[,"Fit_Start"][!neg.start.idx])) /
  exp(mcmc.samples[,"Fit_Start"][!neg.start.idx]) * 100

mcmc.samples[,"Perc_Change"][neg.start.idx] <-
  (exp(mcmc.samples[,"Fit_End"][neg.start.idx]) + exp(mcmc.samples[,"Fit_Start"][neg.start.idx])) /
  exp(abs(mcmc.samples[,"Fit_Start"][neg.start.idx])) * 100


pchange <- median(mcmc.samples[,"Perc_Change"])
pchange.df <- data.frame(pchange=mcmc.samples[,"Perc_Change"])
ggplot(pchange.df,aes(pchange)) + geom_density() +
  geom_vline(xintercept=-30, linetype="dashed", colour="yellow") +
  geom_vline(xintercept=-50, linetype="dashed", colour="orange") +
  geom_vline(xintercept=-70, linetype="dashed", colour="red")

probdecl.30 <-
  sum(mcmc.samples[,"Perc_Change"] <= -30) / dim(mcmc.samples)[1] *100
probdecl.50 <-
  sum(mcmc.samples[,"Perc_Change"] <= -50) / dim(mcmc.samples)[1] *100
probdecl.70 <-
  sum(mcmc.samples[,"Perc_Change"] <= -70) / dim(mcmc.samples)[1] *100


# Next,Merge my MetricsCOSEWIC temp branch with main and push
# reload package and make sure it works...!
# then create outlist here like for JAGS, with convergence flag
# then add stan code to R package
# then re-examine priors


#===============================================================================
# End of Stan output
#===============================================================================

#===============================================================================
# Run JAGS and MLE through MetricsCOSEWIC package for the same stock, stk
# JAGS uses same diffuse priors as above
# (previous versions of the code included Rstanarm, which had its own default
# priors)
#===============================================================================

gen <- 4
yrs.do <- (3 * gen) +1


test.df <- SR_Sample %>%
  dplyr::filter(Stock == stk) %>%
  select(Year,Spn)
head(test.df)

fit.out <- comparePercChange(du.label = stk,
                             du.df = test.df,
                             yrs.window = yrs.do ,
                             calc.yr = 2017,
                             samples.out = TRUE,
                             plot.pattern = TRUE,
                             plot.posteriors = TRUE,
                             plot.boxes  = TRUE)


names(fit.out)
fit.out$Summary

#stock 2
# mle -0.00156  10.37837
#jags  -0.00195  10.37111
#rstanarm  -0.01171  10.43364

#stock 3
# -0.01052  11.64572
# -0.01591  11.68371

# stock 18
# mle -0.08297  12.0372
# jags -0.08255  12.03504
# rstanarm -0.05778  11.88134
