#==============================================================================
# Function to run stan to esimate probabiltiy of declines
# Created by: Carrie Holt
# Date last modified: 14 Oct 2022
#==============================================================================



run.stan <- function(du.label,du.df, yrs.window, 
                     perc.change.bm = c(-30,-50,-70), calc.yr, 
                     standardize.data = TRUE, scenario.name = "test", 
                     path.results="", mcmc.plots = TRUE, 
                     prior_sigma_type = "exp", prior_sigma=2.5, H0=FALSE){
                     # samples.out = TRUE, plot.pattern = TRUE, 
                     # plot.fitted = TRUE, plot.posteriors = TRUE, 
                     # plot.boxes  = TRUE){
  
  # rstan_options(auto_write = TRUE)
  # options(mc.cores = parallel::detectCores())
  
  year.scale <- FALSE #standardize x-axis (years) as well. This is not needed
  
  #==============================================================================
  # Create folder for outputs
  #==============================================================================
  out.dir.stock <- here::here(du.label) #here(path.results,scenario.name)
  
  if (file.exists(out.dir.stock) == FALSE){
    dir.create(out.dir.stock)
  }
  
  out.dir <- here::here(du.label, scenario.name) #here(path.results,scenario.name)
  
  if (file.exists(out.dir) == FALSE){
    dir.create(out.dir)
  }
  
  #==============================================================================
  # Standardize log(Abundance) data
  #==============================================================================
 
  if (standardize.data){
    du.df.raw <- du.df# Save raw values
    du.df <- du.df %>% mutate(mean= mean(na.omit(logAbd)), sd= sd(na.omit(logAbd)))
    du.df <- du.df %>% mutate(logAbd = (logAbd-mean)/sd)
  }
  
  #==============================================================================
  # Priors
  #==============================================================================

  
  # ORIGINAL PRIORS- for unstandardized data
  if(!standardize.data){
    # Priors for intercept, normally distributed: ~N(mean, sig)
    intercept_mean <-  median(du.df$logAbd, na.rm=T)
    intercept_sig <-  max(du.df$logAbd, na.rm=T)
    # max value in the time-series is an upper bound on plausible sigmas on y-i
    
    # Priors for slope, normally distributed: ~N(mean, sig)
    slope_mean <-  0
    slope_sig <-  max(du.df$logAbd, na.rm=T) /
      (max(du.df$Year) - min(du.df$Year))
    # max rise/max run is an upper bound on plausible sigmas on slope
    
    
  }
  
  # UPDATED PRIORS- for standardized data
  if(standardize.data){
    # Priors for intercept, normally distributed: ~N(mean, sig)
    intercept_mean <-  0# mean(du.df$logAbd, na.rm=T) #this value is ~0
    intercept_sig <-   prior_sigma# 2.5* sd(du.df$logAbd, na.rm=T) #this value is ~2.5.
    # 2.5 is recommended by Gelman and used in rstanarm as weakly informative
    # https://cran.r-project.org/web/packages/rstanarm/vignettes/priors.html
    # see also https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
    
    # Priors for slope, normally distributed: ~N(mean, sig)
    slope_mean <-  0
    slope_sig <-  prior_sigma
    
    
  }
  
  
  # ORIGINAL Priors for residual variance, inverse gamma distributed:
  # ~Inv_Gamma(alpha,beta), where alpha = beta = 'Sig_Gam_Dist'
  Sig_Gam_Dist <-  0.1#0.00001
  
  # UPDATED priors for residual variance, sigma exponentially distributed:
  # exponential(rate = 1)
  # https://cran.r-project.org/web/packages/rstanarm/vignettes/priors.html
  exp_rate <- 1
  if(year.scale) {exp_rate <- 1/sd(0:(length(du.df$Year)-1))}
  
  #==============================================================================
  # Set up data inputs for Stan
  #==============================================================================
  
  # test for missing value codes
  
  #du.df$logAbd[c(2)] <- NA
  
  data <- list()
  data$Year <- 0:(length(du.df$Year)-1)
  if(year.scale) data$Year <- scale(data$Year)[,1]
  data$logAbd <- replace_na(du.df$logAbd, 99) #du.df$logAbd
  data$N <- dim(du.df)[1]
  data$intercept_mean <- intercept_mean
  data$Sig_Gam_Dist <- Sig_Gam_Dist
  data$exp_rate <- exp_rate
  data$slope_mean <- slope_mean
  data$intercept_sig <- intercept_sig
  data$slope_sig <- slope_sig
  # Which years have observations? (i.e., are not NAs)
  data$logAbd_obs <- as.numeric(!is.na(du.df$logAbd))
  data$slope <- 0# for H0 model
  data$intercept <- 0# for HO model
  
  
  #==============================================================================
  # Run Stan
  #==============================================================================
  
  # stan_fit <- stan(file = 'linear.stan', data = data, iter = 10000,
  #                  chains = 6,  control = list(adapt_delta = 0.95))
  
  if(prior_sigma_type=="exp"){
    if(!H0) stan_fit <- stan(file = 'linear-exp.stan', data = data, 
                             iter = 10000, chains = 4, thin =10,  
                             control = list(adapt_delta = 0.95), refresh = 0)
    if(H0) stan_fit <- stan(file = 'linear-exp-H0.stan', data = data, 
                            iter = 10000, chains = 4, thin =10,  
                            control = list(adapt_delta = 0.95), refresh = 0)
  }
  
  if(prior_sigma_type=="invgamma"){
    stan_fit <- stan(file = 'linear-invgamma.stan', data = data, iter = 10000,
                     chains = 4, thin =10,  control = list(adapt_delta = 0.95), 
                     refresh = 0)
  }
  
  #==============================================================================
  # Pull out parameter estimates
  #==============================================================================
  
  # #Plot the distribution for missing values, just to check (only worksfor 1 NA)
  # draws_logAbdmis <-
  #   as.matrix(stan_fit, pars = "logAbd_Pred")[ , is.na(du.df$logAbd)]
  # df<-data.frame(data=draws_logAbdmis)
  # ggplot(df, aes(data)) + geom_density()
  
  # # Moved this conergence check to output of Rhats in out.list$fit.obj
  # if (max(bayesplot::rhat(stan_fit), na.rm=T) > 1.02) {
  #   convergence.flag<- TRUE
  #   print("Convergence flag")
  # }
  
  slope.converged <- bayesplot::rhat(stan_fit)["slope"] < 1.02
  
  conv.out <- bayesplot::rhat(stan_fit)[c("slope", "intercept", "sigma")]
  
  
  All_Ests <- data.frame(summary(stan_fit)$summary)
  All_Ests$Param <- row.names(All_Ests)
  
  logAbd_Fits_Stan <- All_Ests[grepl("logAbd_Fit_out", All_Ests$Param),  ]
  logAbd_Preds_Stan <- All_Ests[grepl("logAbd_Pred", All_Ests$Param),  ]
  logAbd_PriorPreds_Stan <- All_Ests[grepl("logAbd_PriorPred", All_Ests$Param),  ]  
  
  
  # Compute log marginal likelihood via bridge sampling for H1
  logML <- bridge_sampler(stan_fit, silent = TRUE)
  
  
  # I don't think this is used, though may be helpful for plotting CI and prediction interval (after rescaling predictions)
  FitsDF <- data.frame(Year = du.df$Year, R = du.df$logAbd, Fit = logAbd_Fits_Stan$X50.,
                       Year = 1:dim(logAbd_Fits_Stan)[1],
                       CI_up = logAbd_Fits_Stan$X97.5.,
                       CI_low = logAbd_Fits_Stan$X2.5.,
                       Pred = logAbd_Preds_Stan$X50.,
                       Pred_up = logAbd_Preds_Stan$X97.5.,
                       Pred_low = logAbd_Preds_Stan$X2.5.)
  
  
  if(!H0){
    # get  posteriors
    fit_values <- rstan::extract(stan_fit)
    intercept_Post <- fit_values$intercept 
    slope_Post <- fit_values$slope
    
    if(standardize.data){
      # Unstandardize posteriors
      if(!year.scale){
        sd <- du.df %>% pull(sd) %>% unique()
        mean <- du.df %>% pull(mean) %>% unique()
        intercept_Post <- fit_values$intercept * sd + mean
        slope_Post <- fit_values$slope * sd
      }
      
      if(year.scale){
        sd.y <- du.df %>% pull(sd) %>% unique()
        mean.y <- du.df %>% pull(mean) %>% unique()
        sd.x <- sd(0:(length(du.df$Year)-1))
        mean.x <- mean(0:(length(du.df$Year)-1))
        intercept_Post <- mean.y - fit_values$slope * (sd.y/sd.x) * mean.x
        slope_Post <- fit_values$slope * (sd.y/sd.x)
        
      }
      
    }
    
    # I don't think this list, 'out', is needed, could refer to intercept_Post and slope_Post directly below
    # Return fit and predicted values
    out <- list()
    out$Fits <- FitsDF
    out$Ests <- All_Ests
    out$intercept_Post <- intercept_Post
    out$slope_Post <- slope_Post
    
    # Calculate % change from the fit at the start and end of the 13 years
    mcmc.samples <- data.frame(int=out$intercept_Post, slope=out$slope_Post,
                               Perc_Change = NA,Perc_Change_Raw = NA)
    if(!year.scale) {
      mcmc.samples <- mcmc.samples %>%
        mutate(Fit_Start = int + slope * data$Year[1] ) %>%
        mutate(Fit_End = int + slope * data$Year[yrs.window])
    }
    
    if(year.scale) {
      yrs <-  0:(length(du.df$Year)-1)
      mcmc.samples <- mcmc.samples %>%
        mutate(Fit_Start = int + slope * yrs[1] ) %>%
        mutate(Fit_End = int + slope * yrs[yrs.window])
    }
    
    # Identify any years with negative start years (only occurs if abundances are
    # scaled prior to log transform)
    neg.start.idx <-  mcmc.samples[,"Fit_Start"] < 0
    
    
    mcmc.samples[,"Perc_Change"][!neg.start.idx] <-
      (exp(mcmc.samples[,"Fit_End"][!neg.start.idx]) - exp(mcmc.samples[,"Fit_Start"][!neg.start.idx])) /
      exp(mcmc.samples[,"Fit_Start"][!neg.start.idx]) * 100
    
    mcmc.samples[,"Perc_Change"][neg.start.idx] <-
      (exp(mcmc.samples[,"Fit_End"][neg.start.idx]) + exp(mcmc.samples[,"Fit_Start"][neg.start.idx])) /
      exp(abs(mcmc.samples[,"Fit_Start"][neg.start.idx])) * 100
    
    mcmc.summary <- as.data.frame(summary(stan_fit))
    
    
    pchange <- median(mcmc.samples[,"Perc_Change"])
    pchange.df <- data.frame(pchange=mcmc.samples[,"Perc_Change"])
    probdecl <- data.frame(BM = perc.change.bm,ProbDecl = NA )
    
    for(i in 1:length(perc.change.bm)){
      
      probdecl[i,2] <- sum(mcmc.samples[,"Perc_Change"] <= perc.change.bm[i]) / dim(mcmc.samples)[1] *100
      
    }
    probdecl <- probdecl %>% arrange(BM)
    
    
    # Plot probability distribution of declines relative to thresholds
    p.dist <- ggplot(pchange.df,aes(pchange)) + geom_density() +
      geom_vline(xintercept=-30, linetype="dashed", colour="yellow") +
      geom_vline(xintercept=-50, linetype="dashed", colour="orange") +
      geom_vline(xintercept=-70, linetype="dashed", colour="red")
    ggsave("prob_declines.pdf", p.dist, path=out.dir)
    
    
    if(mcmc.plots){
      posterior <- as.matrix(stan_fit)
      plot_title <- ggtitle("Posterior distributions- slope",
                            "with medians and 80% intervals")
      p1 <- mcmc_areas(posterior,
                       pars = c("slope"),
                       prob = 0.8) + plot_title
      
      plot_title <- ggtitle("Posterior distributions- intercept",
                            "with medians and 80% intervals")
      p2 <- mcmc_areas(posterior,
                       pars = c("intercept"),
                       prob = 0.8) + plot_title
      
      posterior2 <- extract(stan_fit, inc_warmup = FALSE, permuted = FALSE)
      
      color_scheme_set("mix-blue-pink")
      p3 <- mcmc_trace(posterior2,  pars = c("slope", "intercept"), n_warmup = 0,
                       facet_args = list(nrow = 2, labeller = label_parsed)) + 
        facet_text(size = 15)
      
      # Get observed data vector for posterior predictive check
      logAbd_obs <- data$logAbd
      logAbd_obs[which(logAbd_obs==99)]<- NA #Replace 99 with NAs
      #launch_shinystan(stan_fit)
      
      # Get posterior predicted values for ndraws
      pp <- posterior[,grepl("logAbd_Pred", colnames(posterior))]
      ndraws <- 200
      pp <- data.frame(pp) %>% sample_n(ndraws)
      pp <- data.frame(t(pp))
      pp$obs <- logAbd_obs
      
      pp_long <- pivot_longer(pp, cols=1:(ndraws+1))
      pp_long$obs <- 1
      pp_long <- pp_long %>% mutate(obs = replace(obs, name == "obs", 0))
      
      
      # Plot posterior predictive check with observe distribution of log(Abd)
      p4 <- ggplot(pp_long, mapping = aes(x=value)) + 
        geom_density(aes( group=factor(name)), color = "grey", alpha=0.1)  + 
        geom_density(data = subset(pp_long, obs==0), colour="black", size=1.4) + 
        ggtitle("Posterior Predictive Disributions")
      
      # Get prior predictive distribution
      priorp <- posterior[,grepl("logAbd_PriorPred", colnames(posterior))]
      priorp <- data.frame(priorp)
      
      priorp_long <- pivot_longer(priorp, cols=1:(dim(priorp)[2]), names_to="Years", names_prefix="logAbd_PriorPred.")
      # Unstandardize prior predictions
      if (standardize.data){
        mean.y <- mean(du.df.raw$logAbd, na.rm=T) 
        sd.y <- sd(du.df.raw$logAbd, na.rm=T) 
        priorp_long <- priorp_long %>% mutate(rawValue = 
                                                (value/sd.y) + mean.y)
        
      }
      if (!standardize.data){
        du.df.raw<- data.frame(logAbd =  du.df$logAbd)
        priorp_long <- priorp_long %>% mutate(rawValue = value)
      }
      
      # Plot prior predictive check 
      
      priorp_long <- priorp_long %>% mutate(Years=as.factor(as.numeric(Years)))
      
      p5 <- ggplot(priorp_long, mapping = aes(x=rawValue, y=as.factor(Years))) + 
        geom_density_ridges() + 
        ggtitle("Prior Predictive Disributions", "with range of observed data between vertical dashed lines") +
        ylab("Years") + xlab("Log(abundance)") + 
        geom_vline(xintercept=min(du.df.raw$logAbd, na.rm=T), linetype="dotted") +
        geom_vline(xintercept=max(du.df.raw$logAbd, na.rm=T), linetype="dotted") +
        xlim(-50, 100) +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())  
      # xlim(0, 500000) +
      # geom_vline(xintercept=exp(min(du.df.raw$logAbd, na.rm=T)), linetype="dotted") +
      # geom_vline(xintercept=exp(max(du.df.raw$logAbd, na.rm=T)), linetype="dotted")
      
      
      ggsave("slope_posterior.pdf", p1, path=out.dir)
      ggsave("intercept_posterior.pdf", p2, path=out.dir)
      ggsave("trace_plot.pdf", p3, path=out.dir)
      ggsave("posterior_predictive_check.pdf", p4, path=out.dir)
      ggsave("prior_predictive_check.pdf", p5, path=out.dir)
      ggsave("posterior_predictive_check.png", p4, path=out.dir)
      ggsave("prior_predictive_check.png", p5, path=out.dir)
    }# End of mcmc.plots
    
    write.csv(mcmc.samples,file=paste(out.dir,"/mcmc_samples.csv",sep=""))
    write.csv(conv.out,file=paste(out.dir,"/convergenceRhat.csv",sep=""))
    
    
    # Return fit and predicted values
    
    # if(out.type=="short"){ out.list <- list(pchange = pchange,probdecl = probdecl, summary = mcmc.summary,
    #                                         slope.converged = slope.converged, conv.details = conv.out)}
    # 
    out.list <- list(pchange = pchange,probdecl = probdecl, summary = mcmc.summary, 
                     slope.converged = slope.converged, conv.details = conv.out, 
                     samples = mcmc.samples, fit.obj = stan_fit, logML = logML)
    
  }# End of if (!H0)
  
  # For null hypothesis of slope=0, return only marginal likelihood for BF
  if(H0) out.list <- list(logML = logML)
  
  return(out.list)
} #End of runStan




#===============================================================================
# End of Stan output
#===============================================================================
