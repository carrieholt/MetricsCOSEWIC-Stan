#==============================================================================
# Code to run stan code over scenarios with different priors
# Created by: Carrie Holt
# Date last modified: 6 Jan 2023
#==============================================================================


# Package dependencies
library(dplyr)
library(tidyr)
library(rstan)
library(devtools)
library(ggplot2)
library(ggridges)
library(bayesplot)
library(RColorBrewer)
library(gridExtra)
#library(shinystan)
library(bridgesampling)
library(forcats)
install_github("SOLV-Code/MetricsCOSEWIC", dependencies = TRUE,
               build_vignettes = FALSE)
library(MetricsCOSEWIC)
source("runStan.r")

#==============================================================================
# Data
#==============================================================================

# Example data for 22 sockeye stocks are included in MetricsCOSEWIC package
# Quick look at spawners data over the most recent 13 years (3 gens + 1year),
# 2005-2017

SR_Sample %>% filter(Year>(2017-13)) %>%
  ggplot(aes(Year, Spn)) + geom_line(aes(colour=factor(Stock)))


# Read in data
perc.change.bm <- c(-30, -50, -70)
mcmc.plots <- TRUE
out.type <- "long"
standardize.data <- TRUE

du.df.long <- SR_Sample %>% select(Stock,Year,Spn) %>% 
  rename(DU = Stock, Abd = Spn)
gen <- 4

# Decide the window of years over which to calculate trends, and the year at  
# which to calculate it (final or calc.year)
yrs.window <- (3 * gen) +1
calc.year <- 2017

# Choose a stock to look at: "Stock18", "Stock2" and "Stock3" are examples 
# where posteriors were sensitive to choice of prior (i.e., posteriors from 
# previous implementation of RStanarm with informative priors differed from 
# those when we used the diffuse priors described below)
stk <- "Stock3"#"Stock18"#"Stock2"##"Stock18"
du.label <- stk
du.df.long <- du.df.long %>% filter(DU==stk) %>% mutate(logAbd=log(Abd)) 

# Add an NA (for checking)
#du.df.long <- du.df.long %>% mutate(logAbd = replace(logAbd, Year==2008, NA))

# Shorten data set to years used for decline calculations
du.df <- du.df.long %>% filter(Year > (calc.year - yrs.window) & 
                                 Year <= calc.year)


# Look at data frame
du.df
ggplot(du.df, aes(Year,logAbd)) + geom_point()

#=============================================================================
# RUN STAN
#=============================================================================
# Run STAN with standardized data and associated priors
run1 <- run.stan(du.label=stk, du.df=du.df, yrs.window=yrs.window, 
                 standardize.data = TRUE,  scenario.name = "Standardized", 
                 prior_sigma_type = "invgamma")
# Run STAN with unstandardized data and associated priors
run2 <- run.stan(du.label=stk, du.df=du.df, yrs.window=yrs.window, 
                 standardize.data = FALSE,  scenario.name = "Unstandardized",
                 prior_sigma_type = "invgamma")

# Run STAN with standardized data and exponetial(1) prior on residual SD
run3 <- run.stan(du.label=stk, du.df=du.df, yrs.window=yrs.window, 
                 standardize.data = TRUE,  
                 scenario.name = "Standardized-ExpSigma", 
                 prior_sigma_type = "exp")
# # associated null hypothesis with slope = 0, yi = 0 (standardized data)
# run3.H0 <- run.stan(du.label=stk, du.df=du.df, yrs.window=yrs.window, 
#                     standardize.data = TRUE,  
#                     scenario.name = "Standardized-ExpSigma", 
#                     prior_sigma_type = "exp", H0 = TRUE)

# Run STAN with standardized data and various sigma priors on slope and yi
run4 <- run.stan(du.label=stk, du.df=du.df, yrs.window=yrs.window, 
                 standardize.data = TRUE,  
                 scenario.name = "Standardized-PriorSigma1",
                 prior_sigma_type = "exp",  prior_sigma = 1)
# #associated null hypothesis with slope = 0, yi = 0 (standardized data)
# run4.H0 <- run.stan(du.label=stk, du.df=du.df, yrs.window=yrs.window, 
#                  standardize.data = TRUE,  
#                  scenario.name = "Standardized-PriorSigma1",
#                  prior_sigma_type = "exp",  prior_sigma = 1, H0=TRUE)

run5 <- run.stan(du.label=stk, du.df=du.df, yrs.window=yrs.window, 
                 standardize.data = TRUE,  
                 scenario.name = "Standardized-PriorSigma10",
                 prior_sigma_type = "exp",  prior_sigma = 10)
# #associated null hypothesis with slope = 0, yi = 0 (standardized data)
# run5.H0 <- run.stan(du.label=stk, du.df=du.df, yrs.window=yrs.window, 
#                  standardize.data = TRUE,  
#                  scenario.name = "Standardized-PriorSigma10",
#                  prior_sigma_type = "exp",  prior_sigma = 10, H0=TRUE)

# #=============================================================================
# # Bayesian Factor = extent to which the data update the prior odds to the 
# # posterior odds: http://ejwagenmakers.com/2016/LyEtAl2016JMP.pdf
# #=============================================================================
# 
# BayesF.10_3 <- bf(run3$logML, run3.H0$logML) 
# #BayesF.01_3 <- bf(run3.H0$logML, run3$logML)
# #error_measures(run3.H0$logML)
# #error_measures(run3.H0$logML)
# BayesF.10_4 <- bf(run4$logML, run4.H0$logML) 
# BayesF.10_5 <- bf(run5$logML, run5.H0$logML) 

#=============================================================================
# Compile outputs into a single dataframe, run.df
#=============================================================================

run1.df <- data.frame(Value=run1$samples$Perc_Change, PriorSigma=2.5, 
                      std.data=1, VarPrior="InvGamma")
# Remove outliers <0.025 and >0.975 of the distribution for plotting purposes 
# tails distract viewer from the mass of the distribution
run1.df.plot <- run1.df %>% mutate(Percentile = percent_rank(Value)) %>% 
  filter(Percentile >0.025 & Percentile <0.975) %>% select(-Percentile)

run2.df <- data.frame(Value=run2$samples$Perc_Change, PriorSigma=2.5, 
                      std.data=0, VarPrior="InvGamma")
run2.df.plot <- run2.df %>% mutate(Percentile = percent_rank(Value)) %>% 
  filter(Percentile >0.025 & Percentile <0.975) %>% select(-Percentile)

run3.df <- data.frame(Value=run3$samples$Perc_Change, PriorSigma=2.5, 
                      std.data=1, VarPrior="Exp")
run3.df.plot <- run3.df %>% mutate(Percentile = percent_rank(Value)) %>% 
  filter(Percentile >0.025 & Percentile <0.975) %>% select(-Percentile)

run4.df <- data.frame(Value=run4$samples$Perc_Change, PriorSigma=1, 
                      std.data=1, VarPrior="Exp")
run4.df.plot <- run4.df %>% mutate(Percentile = percent_rank(Value)) %>% 
  filter(Percentile >0.025 & Percentile <0.975) %>% select(-Percentile)

run5.df <- data.frame(Value=run5$samples$Perc_Change, PriorSigma=10, 
                      std.data=1, VarPrior="Exp")
run5.df.plot <- run5.df %>% mutate(Percentile = percent_rank(Value)) %>% 
  filter(Percentile >0.025 & Percentile <0.975) %>% select(-Percentile)


#runStd.df <- run1.df.plot %>% bind_rows(run2.df.plot)
runStd.df <- run1.df %>% bind_rows(run2.df)

#runSigType.df <- run1.df.plot %>% bind_rows(run3.df.plot)
runSigType.df <- run1.df %>% bind_rows(run3.df)

# runSigPrior.df <- run3.df.plot %>% bind_rows(run4.df.plot) %>% 
#   bind_rows(run5.df.plot)
runSigPrior.df <- run3.df %>% bind_rows(run4.df) %>% 
  bind_rows(run5.df)

#=============================================================================
# Plot posterior estimates of percent change: standardization of data
#=============================================================================
text.df <- data.frame(std.data=c(0,0,0,1,1,1), Value=rep(c(-70, -50,-30),2), 
                     lab=c(round(run1$probdecl, 1)[,2], 
                             round(run2$probdecl, 1)[,2]))
lab.df <- data.frame(Value=105, std.data=as.factor(0), duName=du.df$DU[1])


p1 <- ggplot(runStd.df, (aes(x=as.factor(std.data), y=Value, 
                             fill = as.factor(std.data)))) + 
  geom_hline(yintercept=-30, linetype="dashed", colour="light grey") + 
  geom_hline(yintercept=-50, linetype="dashed", colour="dark grey") + 
  geom_hline(yintercept=-70, linetype="dashed", colour="black") +  
  geom_violin(alpha=0.3, trim=FALSE) + 
  geom_boxplot(width=0.1) +
  ggtitle("Posterior distributions of percent change", 
          "for original priors on unstandarized data vs. those 
          on standardized data (and probabitlity below 
          thresholds)") +
  theme_minimal() + 
  xlab(element_blank()) +
  scale_x_discrete(labels=c("0"="Unstandardized","1"="Standardized")) +
  ylab("Perent change over three generations") + 
  theme(legend.position="none") + 
  scale_fill_brewer(element_blank(), palette="Dark2") + 
  geom_text(aes(label=lab), data=text.df, colour=grey(0.5), nudge_x=-0.1, 
            nudge_y=4, size=2.8) + 
  ylim(-100,105) + 
  geom_label(data=lab.df, aes(label=duName), nudge_x=-0.4, 
                              label.padding = unit(0.55, "lines"),  
                              fill="white", label.size=0, size=3)

  
p1.cum <- ggplot(runStd.df, (aes(Value, colour=as.factor(std.data)))) + 
  stat_ecdf(geom = "step", pad = TRUE, size = 2) + 
  labs(title=element_blank(), y = "Probability", 
       x="Percent change") +
  scale_colour_brewer(element_blank(), palette="Dark2", labels = c("Unstandardized", "Standardized")) + 
  xlim(-100,20) +
  geom_vline(xintercept=-30, linetype="dashed", colour="light grey") + 
  geom_vline(xintercept=-50, linetype="dashed", colour="dark grey") + 
  geom_vline(xintercept=-70, linetype="dashed", colour="black")  +
  geom_hline(yintercept= run1$probdecl[3,"ProbDecl"]/100, linetype="solid", 
             colour = brewer.pal(3, "Dark2")[1], size = 1) + 
  geom_hline(yintercept= run1$probdecl[2,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[1], size = 1) + 
  geom_hline(yintercept= run1$probdecl[1,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[1], size = 1) +  
  geom_hline(yintercept= run2$probdecl[3,"ProbDecl"]/100, linetype="solid", 
           colour = brewer.pal(3, "Dark2")[2], size = 1) + 
  geom_hline(yintercept= run2$probdecl[2,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[2], size = 1) + 
  geom_hline(yintercept= run2$probdecl[1,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[2], size = 1)  +
  theme(legend.position="bottom")

# Could add transparency to horizontal lines:
#https://stackoverflow.com/questions/55537556/how-to-get-alpha-value-of-a-hexadecimal-color-in-r
#https://www.rdocumentation.org/packages/GISTools/versions/0.7-4/topics/Create%20Transparency

#=============================================================================
# Plot time series with  percent change: standardization of data
#=============================================================================

p2 <- ggplot(du.df.long, aes(Year,logAbd)) + geom_point() + 
  geom_segment(x= calc.year - yrs.window, y = median(run1$samples$Fit_Start),
               xend=calc.year, yend = median(run1$samples$Fit_End), 
               colour=brewer.pal(n = 8, name = "Dark2")[2], size=1.5) + 
  geom_segment(x= calc.year - yrs.window, y = median(run2$samples$Fit_Start),
               xend=calc.year, yend = median(run2$samples$Fit_End), 
               colour=brewer.pal(n = 8, name = "Dark2")[1], size=1.5)+ 
  ggtitle("Time-series of log(abundances)", 
          "with linear declines over most recent three 
          generations for models with the original 
          priors on unstandarized data vs those on 
          standardized data") + 
  ylab("Log(abundances)") + 
  xlim(max(du.df.long$Year) - yrs.window, max(du.df.long$Year))
# note, factors are listed in this order: "Unstandardized", "Standardized"


tabledf <- data.frame(Threshold = c(-70,-50,-30), 
           ProbBelow.Standardized = round(run1$probdecl, 1)[,2], 
           ProbBelow.Unstandardized = round(run2$probdecl, 1)[,2])

print(tabledf)

#=============================================================================
# Plot posterior estimates of percent change: type of sigma prior
#=============================================================================
text.df <- data.frame(VarPrior=c(rep("InvGamma",3),rep("Exp",3)), 
                      Value=rep(c(-70, -50,-30),2), 
                      lab=c(round(run1$probdecl, 1)[,2], 
                            round(run3$probdecl, 1)[,2]))
lab.df <- data.frame(Value=105, VarPrior="Exp", duName=du.df$DU[1])

p3 <- ggplot(runSigType.df, (aes(x=as.factor(VarPrior), y=Value, 
                                 fill=as.factor(VarPrior)))) + 
  # geom_hline(yintercept=-30, linetype="dashed", colour="yellow") + 
  # geom_hline(yintercept=-50, linetype="dashed", colour="orange") + 
  # geom_hline(yintercept=-70, linetype="dashed", colour="red") +  
  geom_hline(yintercept=-30, linetype="dashed", colour="light grey") + 
  geom_hline(yintercept=-50, linetype="dashed", colour="dark grey") + 
  geom_hline(yintercept=-70, linetype="dashed", colour="black") +  
  geom_violin(alpha=0.3, trim=FALSE) + 
  geom_boxplot(width=0.1) +
  ggtitle("Posterior distributions of percent change", 
          "for exponential priors on sigma and inverse 
          Gamma on variance (and probabitlity below 
          thresholds)") +
  theme_minimal() + 
  xlab(element_blank()) + 
  ylab("Perent change over three generations") + 
  scale_fill_brewer(element_blank(), palette="Dark2") +
  theme(legend.position="none") + 
  geom_text(aes(label=lab), data=text.df, colour=grey(0.5), nudge_x=-0.1, 
            nudge_y=4, size=2.8) + 
  ylim(-100,105) + 
  geom_label(data=lab.df, aes(label=duName), nudge_x=-0.4, 
             label.padding = unit(0.55, "lines"),  
             fill="white", label.size=0, size=3)


p3.cum <- ggplot(runSigType.df, (aes(Value, colour=as.factor(VarPrior)))) + 
  stat_ecdf(geom = "step", pad = TRUE, size = 2) + 
  labs(title=element_blank(), y = element_blank(), 
       x="Percent change") +
  scale_colour_brewer(element_blank(), palette="Dark2", labels = c("Exp", "InvGamma")) + 
  xlim(-100,20) +
  geom_vline(xintercept=-30, linetype="dashed", colour="light grey") + 
  geom_vline(xintercept=-50, linetype="dashed", colour="dark grey") + 
  geom_vline(xintercept=-70, linetype="dashed", colour="black")  +
  geom_hline(yintercept= run1$probdecl[3,"ProbDecl"]/100, linetype="solid", 
             colour = brewer.pal(3, "Dark2")[2], size = 1) + 
  geom_hline(yintercept= run1$probdecl[2,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[2], size = 1) + 
  geom_hline(yintercept= run1$probdecl[1,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[2], size = 1) +  
  geom_hline(yintercept= run3$probdecl[3,"ProbDecl"]/100, linetype="solid", 
             colour = brewer.pal(3, "Dark2")[1], size = 1) + 
  geom_hline(yintercept= run3$probdecl[2,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[1], size = 1) + 
  geom_hline(yintercept= run3$probdecl[1,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[1], size = 1)   +
  theme(legend.position="bottom")


#=============================================================================
# Plot time series with  percent change: type of sigma prior
#=============================================================================

p4 <- ggplot(du.df.long, aes(Year,logAbd)) + geom_point() + 
  geom_segment(x= calc.year - yrs.window, y = median(run1$samples$Fit_Start),
               xend=calc.year, yend = median(run1$samples$Fit_End), 
               colour=brewer.pal(n = 8, name = "Dark2")[2], size=1.5) + 
  geom_segment(x= calc.year - yrs.window, y = median(run3$samples$Fit_Start),
               xend=calc.year, yend = median(run3$samples$Fit_End), 
               colour=brewer.pal(n = 8, name = "Dark2")[1], size=1.5) + 
  ggtitle("Time-series of log(abundances)", 
          "with linear declines over most recent three 
          generations for models with exponential 
          priors on residual sigma and inverse Gamma 
          priors on residual variance") + 
  ylab("Log(abundances)") + 
  xlim(max(du.df.long$Year) - yrs.window, max(du.df.long$Year))

# note, factors are listed in this order: "Exp", "InvGamma"

tabledf <- data.frame(Threshold = c(-70,-50,-30), 
                      ProbBelow.InvGamma = round(run1$probdecl, 1)[,2], 
                      ProbBelow.Exponential = round(run3$probdecl, 1)[,2])

print(tabledf)

#=============================================================================
# Plot posterior estimates of percent change: sigma on slope on yi priors
#=============================================================================
text.df <- data.frame(PriorSigma=c(rep(2.5,3),rep(1,3), rep(10,3)), 
                      Value=rep(c(-70, -50,-30),3), 
                      lab=c(round(run3$probdecl, 1)[,2], 
                            round(run4$probdecl, 1)[,2],
                            round(run5$probdecl, 1)[,2]))
lab.df <- data.frame(Value=30, PriorSigma=as.factor(1), duName=du.df$DU[1])

p5 <- ggplot(runSigPrior.df, (aes(x=as.factor(PriorSigma), y=Value, 
                                  fill = as.factor(PriorSigma)))) + 
  # geom_hline(yintercept=-30, linetype="dashed", colour="yellow") + 
  # geom_hline(yintercept=-50, linetype="dashed", colour="orange") + 
  # geom_hline(yintercept=-70, linetype="dashed", colour="red") +  
  geom_hline(yintercept=-30, linetype="dashed", colour="light grey") + 
  geom_hline(yintercept=-50, linetype="dashed", colour="dark grey") + 
  geom_hline(yintercept=-70, linetype="dashed", colour="black") +  
  geom_violin(alpha=0.3, trim=FALSE) + 
  geom_boxplot(width=0.1) +
  ggtitle("Posterior distributions of percent change", 
          "for various sigma values on slope and y-intercept 
          priors (and probabitlity below thresholds)") +
  theme_minimal() + 
  xlab(element_blank()) + 
  ylab("Perent change over three generations") + 
  theme(legend.position="none") +
  scale_fill_brewer(element_blank(), palette="Dark2") +
  scale_x_discrete(labels=c("1"="Sigma=1","2.5"="Sigma=2.5","10"="Sigma=10")) + 
  geom_text(aes(label=lab), data=text.df, colour=grey(0.5), nudge_x=-0.1, 
            nudge_y=4, size=2.8) + 
  ylim(-100,30) + 
  geom_label(data=lab.df, aes(label=duName), nudge_x=-0.4, 
             label.padding = unit(0.55, "lines"),  
             fill="white", label.size=0, size=3)

p5.cum <- ggplot(runSigPrior.df, (aes(Value, colour=as.factor(PriorSigma)))) + 
  stat_ecdf(geom = "step", pad = TRUE, size = 2) + 
  labs(title=element_blank(), y = element_blank(), 
       x="Percent change") +
  scale_colour_brewer(element_blank(), palette="Dark2", labels = 
                        c("1"="Sigma=1","2.5"="Sigma=2.5","10"="Sigma=10")) + 
  xlim(-100,20) +
  geom_vline(xintercept=-30, linetype="dashed", colour="light grey") + 
  geom_vline(xintercept=-50, linetype="dashed", colour="dark grey") + 
  geom_vline(xintercept=-70, linetype="dashed", colour="black")  +
  geom_hline(yintercept= run3$probdecl[3,"ProbDecl"]/100, linetype="solid", 
             colour = brewer.pal(3, "Dark2")[2], size = 1) + 
  geom_hline(yintercept= run3$probdecl[2,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[2], size = 1) + 
  geom_hline(yintercept= run3$probdecl[1,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[2], size = 1) +  
  geom_hline(yintercept= run4$probdecl[3,"ProbDecl"]/100, linetype="solid", 
             colour = brewer.pal(3, "Dark2")[1], size = 1) + 
  geom_hline(yintercept= run4$probdecl[2,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[1], size = 1) + 
  geom_hline(yintercept= run4$probdecl[1,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[1], size = 1) +
  geom_hline(yintercept= run5$probdecl[3,"ProbDecl"]/100, linetype="solid", 
             colour = brewer.pal(3, "Dark2")[3], size = 1) + 
  geom_hline(yintercept= run5$probdecl[2,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[3], size = 1) + 
  geom_hline(yintercept= run5$probdecl[1,"ProbDecl"]/100, linetype="solid", colour = 
               brewer.pal(3, "Dark2")[3], size = 1)  +
  theme(legend.position="bottom")

  
#=============================================================================
# Plot time series with  percent change: sigma on slope on yi priors
#=============================================================================

p6 <- ggplot(du.df.long, aes(Year,logAbd)) + geom_point() + 
  geom_segment(x= calc.year - yrs.window, y = median(run3$samples$Fit_Start),
               xend=calc.year, yend = median(run3$samples$Fit_End), 
               colour=brewer.pal(n = 8, name = "Dark2")[2], size=1.5) + 
  geom_segment(x= calc.year - yrs.window, y = median(run4$samples$Fit_Start),
               xend=calc.year, yend = median(run4$samples$Fit_End), 
               colour=brewer.pal(n = 8, name = "Dark2")[1], size=1.5) + 
  geom_segment(x= calc.year - yrs.window, y = median(run5$samples$Fit_Start),
               xend=calc.year, yend = median(run5$samples$Fit_End), 
               colour=brewer.pal(n = 8, name = "Dark2")[3], size=1.5) + 
  ggtitle("Time-series of log(abundances)", 
          "with linear declines over most recent three 
          generations for models with the various 
          sigmas on normal priors on slope and 
          y-intercept parameters") + 
  ylab("Log(abundances)") + xlim(max(du.df.long$Year) - yrs.window, 
                                 max(du.df.long$Year))


tabledf <- data.frame(Threshold = c(-70,-50,-30), 
                      ProbBelow.Sig1 = round(run4$probdecl, 1)[,2], 
                      ProbBelow.Sig2.5 = round(run3$probdecl, 1)[,2],
                      ProbBelow.Sig10 = round(run5$probdecl, 1)[,2])

print(tabledf)

#=============================================================================
# save/print plots
#=============================================================================

p.all <- grid.arrange(p1,p2,p3,p4,p5,p6, nrow=3)
ggsave(file=paste(stk, "/PriorSens_", du.df$DU[1],".pdf", sep=""), 
       plot=p.all, width=8, height=11)
ggsave(file=paste(stk, "/PriorSens_", du.df$DU[1],".png", sep=""), 
       plot=p.all, width=8, height=11)

p.cum.all <- grid.arrange(p1.cum, p3.cum, p5.cum, nrow=1)
ggsave (file = paste(stk, "/PriorSensCumulative_", du.df$DU[1],".pdf", sep=""), 
       plot = p.cum.all, width=8, height=8)
ggsave (file = paste(stk, "/PriorSensCumulative_", du.df$DU[1],".png", sep=""), 
       plot = p.cum.all, width=11, height=6)

#===============================================================================
# Run JAGS and MLE through MetricsCOSEWIC package for the same stock, stk
# JAGS uses same diffuse priors as above
# (previous versions of the code included Rstanarm, which had its own default
# priors)
#===============================================================================

runJags <- TRUE
if(runJags){
  gen <- 4
  yrs.window <- (3 * gen) +1
  calc.year <- 2017
  stk <- stk
  test.df <- SR_Sample %>%  dplyr::filter(Stock == stk) %>% select(Year,Spn)
  # test.df$Spn[49]<-NA #Add an NA for 2008 (if NA added prior to stan run)
  
  fit.out <- comparePercChange(du.label = stk,
                               du.df = test.df,
                               yrs.window = yrs.window ,
                               calc.yr = calc.year,
                               samples.out = TRUE,
                               plot.pattern = TRUE,
                               plot.posteriors = TRUE,
                               plot.boxes  = TRUE,
                               do.rstanarm = TRUE)
  
  names(fit.out)
  fit.out$Summary

  probdecl.rstanarm <- data.frame(BM = perc.change.bm,ProbDecl = NA )
  for(i in 1:length(perc.change.bm)){
    probdecl.rstanarm[i,2] <- sum(fit.out$samples.rstanarm[,"Perc_Change"] <= 
                                    perc.change.bm[i]) / 
      length(fit.out$samples.rstanarm[,"Perc_Change"]) *100
  }
  probdecl.rstanarm  <- probdecl.rstanarm  %>% arrange(BM)
  
  probdecl.jags <- data.frame(BM = perc.change.bm,ProbDecl = NA )
  for(i in 1:length(perc.change.bm)){
    probdecl.jags[i,2] <- sum(fit.out$samples.jags[,"Perc_Change"] <= 
                                perc.change.bm[i]) / 
      length(fit.out$samples.jags[,"Perc_Change"]) *100
  }
  probdecl.jags  <- probdecl.jags  %>% arrange(BM)
  
  probdecl.rstanarm 
  probdecl.jags 
}

#===============================================================================
# Plot posterior estimates of percent change: stan, JAGS, rstanarm, MLE
# Using same (original) priors for JAGS and stan, and 'null' priors for rstanarm
#===============================================================================

if (stk == "Stock18") maxy <- 105
if (stk != "Stock18") maxy <- 200

run2.df$method <- "stan"
MethodsCompareDF <- run2.df %>% 
  add_row(Value = fit.out$samples.jags[,"Perc_Change"], PriorSigma=0, 
          std.data=0,VarPrior="InvGamma", method="JAGS") %>% 
  add_row(Value = fit.out$samples.rstanarm[,"Perc_Change"], PriorSigma=0, 
          std.data=0,VarPrior="InvGamma", method="rstanarm") %>% 
  add_row(Value = rep(fit.out$Summary["MLE",]["pchange"],500), PriorSigma=0, 
          std.data=0, VarPrior="InvGamma", method="MLE")

# Change order of methods for plotting
MethodsCompareDF <- MethodsCompareDF %>% 
  mutate(method = fct_relevel(as.factor(method), "MLE", "JAGS", "rstanarm", 
                              "stan"))

text.df <- data.frame(Value=rep(c(-70,-50,-30),3), 
                      lab = c(round(fit.out$ProbDecl,1)[,2], 
                              round(fit.out$Probdecl.rstan,1)[,2],
                              round(run2$probdecl, 1)[,2]), 
                      method=c(rep("JAGS",3), rep("rstanarm",3), rep("stan",3)))
lab.df <- data.frame(Value=maxy, method="MLE", duName=du.df$DU[1])

# df <- data.frame(stan=run2.df.plot$Value, method="stan")
# ggplot(df, aes(x=as.factor(method), y=Value))+geom_violin()

p7 <- ggplot(MethodsCompareDF, aes(x=as.factor(method), y=Value, 
                                    fill=as.factor(method))) + 
  geom_hline(yintercept=-30, linetype="dashed", colour="light grey") + 
  geom_hline(yintercept=-50, linetype="dashed", colour="dark grey") + 
  geom_hline(yintercept=-70, linetype="dashed", colour="black") +  
  geom_violin(alpha=0.3, trim=TRUE, scale="width") + 
  geom_boxplot(width=0.1) +
  ggtitle("Posterior distributions of percent change", 
          "for various estimation methods (and probabitlity 
          below thresholds)") +
  theme_minimal() + 
  xlab(element_blank()) + 
  ylab("Perent change over three generations") + 
  theme(legend.position="none") +
  geom_text(aes(label=lab), data=text.df, colour=grey(0.5), nudge_x=-0.2, 
            nudge_y=3, size=2.8) + 
  ylim(-100,maxy) + 
  geom_label(data=lab.df, aes(label=duName), nudge_x=-0.3, 
             label.padding = unit(0.55, "lines"),  
             fill="white", label.size=0, size=4)

ggsave(file=paste(stk,"/MethodSens_", du.df$DU[1],".pdf", sep=""), plot=p7, width=5, height=5)
ggsave(file=paste(stk,"/MethodSens_", du.df$DU[1],".png", sep=""), plot=p7, width=5, height=5)


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

# =============================================================================
# Calculating and plotting BF from Rstanarm fits (exapmle only)
# =============================================================================
BF <- FALSE
if(BF){
  library(rstanarm)
  dataRs.in <- data.frame(logAbd= data$logAbd, Year = data$Year)
  dataRs.in$logAbd[4] <- NA
  
  fit.stan <- stan_glm(logAbd ~ Year, data = dataRs.in, prior=normal(0,2.5), 
                       prior_intercept=normal(0,2.5), prior_aux=exponential(1))
  fit.stan <- stan_glm(logAbd ~ Year, data = dataRs.in, prior=normal(0,1), 
                       prior_intercept=normal(0,1), prior_aux=exponential(1))
  fit.stan$coefficients
  plot(x=dataRs.in$Year, y=dataRs.in$logAbd)
  abline(fit.stan$coefficients)
  library(see)
  bf <- bayesfactor_parameters(fit.stan, null=0) #null model with slope = 0
  bf <- bayesfactor_parameters(fit.stan, null=c(-0.1,Inf)) #null range 
  # change lower limit from -0.1 to slope in standardized log(Abd) = -30%
  bf
  plot(bf)+xlim(c(-1,1))
  # BF for "year" coefficient that is large >3 suggest strong evidence from data 
  # against the null hypothesis that... slope is 0, or slope is >-30%, or ...
  # see my log for refs for this (Oct 23, 2022)
}

# Next
# (1) add cumulative plot like Atlantic salmon
# (2) add stan code to R package



