#==============================================================================
# Code to run stan code for Fraser Sockeye re-assessments (3 DUs)
# Created by: Carrie Holt
# Date last modified: 11 Feb 2023
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
# install_github("Pacific-salmon-assess/MetricsCOSEWIC", dependencies = TRUE,
#                build_vignettes = FALSE)
library(MetricsCOSEWIC)
source("runStan.r")

# Add code to create results directories if they don't yet exist

# Run script for Fraser Sockeye re-assssments
data <- read.csv(here::here("data", "FraserSockeyeReAssessmentData.csv"))
data <- data[,1:4] #Omit notes column
data.long <- pivot_longer(data, cols = c("Harrison..D.S..L", "Kamloops.ES", 
                            "Lillooet.Harrison.L"), 
             names_to = "DU", 
             values_to= "Abd")

# # Not needed if time-series are smoothed
# # Replace zeros with half the min value observed....
# data.long <- data.long %>% group_by(DU) %>% 
#   mutate( Abd = replace(Abd, Abd ==0, NA)) %>% 
#   mutate (min = min(Abd, na.rm=T)/2) %>% ungroup()
# data.long.NA <- data.long %>% filter( is.na(Abd ==TRUE)) %>% mutate (Abd=min/2)
# data.long <- data.long %>% filter(Abd > 0) %>% 
#   add_row(data.long.NA) %>% select(-min)


stk <- "Kamloops.ES"#""Harrison..D.S..L"#"Lillooet.Harrison.L"#"
scenario.name <- "long-time-series" #"short-time-series"
gen <- 4

if (!file.exists(here::here(eval(stk))))  dir.create(here::here(eval(stk)))

if (file.exists(here::here(eval(stk)))){
  if (!file.exists(here::here(eval(stk), "long-time-series"))) {
    dir.create(here::here(eval(stk), "long-times-eries"))
  }
  if (!file.exists(here::here(eval(stk), "short-time-series"))) {
    dir.create(here::here(eval(stk), "short-times-eries"))
  }
}

# Decide the window of years over which to calculate trends, and the year at  
# which to calculate it (final or calc.year)
yrs.window <- (3 * gen) + 1
calc.year <- 2021
du.df.long <- data.long %>% filter(DU==stk) 


# Arithmetic smooth time-series (as in 2017 COSEWIC report on Fr sockeye) 
du.df.long$Abd <- smoothSeries(du.df.long$Abd, gen=gen, filter.sides=1, 
                               log.transform = FALSE, out.exp = FALSE, 
                               na.rm=FALSE)

du.df.long <- du.df.long %>% mutate(logAbd=log(Abd)) 

# Shorten data set to years used to estimate trends
if(scenario.name == "short-time-series"){
  du.df <- du.df.long %>% filter(Year > (calc.year - yrs.window) &
                                   Year <= calc.year)
}
if(scenario.name == "long-time-series"){
  du.df <- du.df.long 
}

#=============================================================================
# RUN STAN
#=============================================================================
# Run STAN with standardized data, exp prior on var, and 2.5 sigma priors on 
# slope and yi (DEFAULT)

stan.out <- run.stan(du.label=stk, du.df=du.df, yrs.window=yrs.window, 
                     standardize.data = TRUE,  
                     scenario.name = scenario.name, 
                     prior_sigma_type = "exp")

out.df <- data.frame(Value=stan.out$samples$Perc_Change, PriorSigma=2.5, 
                     std.data=1, VarPrior="Exp")

#=============================================================================
# Plot cumulative distribution and probbabilty densities
#=============================================================================
  
# Fnd upper x-limit from 90th percentile of probability distribution
y.up <- out.df %>% summarise(quantile = scales::percent(c(0.90)), 
                             Value = quantile(Value, c(0.90)))
y.up <- as.numeric(y.up[[2]])
y.up <- max (0, y.up)
if(y.up >= 10) {
  tick.marks <- seq(10,y.up,20)
  tick.marks <- c(-90, -70, -50, -30, -10, tick.marks)
} 
if(y.up < 10) tick.marks <- c(-90, -70, -50, -30, -10)
if(y.up < 10) axis.text.size <- 18
if(y.up >= 10) axis.text.size <- 12
if (stk == "Harrison..D.S..L") {
  xnudge <- c(0, 0, 15)
  ynudge <- c(0.02, 0.02, 0.01)}
if (stk != "Harrison..D.S..L") {
  xnudge <- 0
  ynudge <- 0.02}

# Text for labeling probabilities at -70%, -50% and -30% declines on plot
text.df <- data.frame(x = c(-95,-95,-95), 
                      Value = round(stan.out$probdecl, 0)[,2]/100+0.02, 
                      lab = round(stan.out$probdecl, 0)[,2]/100,
                      col = brewer.pal(3, "Dark2")[c(2,3,1)])
text.df <- text.df %>% mutate(lab=replace(lab, lab==1.00, ">0.99"))
# Plot cumulative distribution
gg.base1 <- ggplot(out.df, (aes(Value))) +#, colour=as.factor(PriorSigma)))) + 
  stat_ecdf(geom = "step", pad = TRUE, size = 2, colour=grey(0.4)) + 
  labs(title= "Cumulative probability", y = element_blank(), 
       x="Percent change") +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20), 
        plot.title = element_text(size=24)) +
  scale_colour_brewer(element_blank(), palette="Dark2") + 
  coord_cartesian(xlim = c(-100,y.up) )+
  scale_x_continuous(breaks = tick.marks) +
  # geom_text(aes(x= x, y = Value, label=lab), data=text.df,  colour=grey(0.5), nudge_x=0, 
  #           nudge_y=0.02, size=8) + 
  geom_text(aes(x= x, y = Value, label=lab, colour = col), data=text.df, 
            nudge_x = xnudge, 
            nudge_y = ynudge, size = 6) + 
  geom_vline(xintercept=-30, linetype="dashed", 
             colour=brewer.pal(3, "Dark2")[1]) + 
  geom_vline(xintercept=-50, linetype="dashed", 
             colour=brewer.pal(3, "Dark2")[2]) + 
  geom_vline(xintercept=-70, linetype="dashed", 
             colour=brewer.pal(3, "Dark2")[3]) +
  geom_hline(yintercept= stan.out$probdecl[3,"ProbDecl"]/100, linetype="solid", 
             colour = brewer.pal(3, "Dark2")[1], size = 1) + 
  geom_hline(yintercept= stan.out$probdecl[2,"ProbDecl"]/100, linetype="solid", 
             colour = brewer.pal(3, "Dark2")[2], size = 1) + 
  geom_hline(yintercept= stan.out$probdecl[1,"ProbDecl"]/100, linetype="solid", 
             colour = brewer.pal(3, "Dark2")[3], size = 1) +  
  theme(legend.position="none")

# Plot three probability densities for -70, -50, and -30 thresholds
gg.prob.a <- ggplot(out.df, (aes(Value))) + 
  geom_density(colour = grey(0.4), size = 1, fill = brewer.pal(3, "Dark2")[1]) +
  coord_cartesian(xlim = c(-100,y.up) )+
  scale_x_continuous(breaks = tick.marks) +
  scale_y_continuous(breaks = NULL) + 
  theme(axis.text=element_text(size=axis.text.size), axis.title=element_text(size=22), plot.title = element_text(size=24))

p <- ggplot_build(gg.prob.a)$data[[1]]
if(max(p$x, na.rm=T) <= -30) gg.base2a <- gg.prob.a
if(max(p$x, na.rm=T) > -30) gg.base2a <- gg.prob.a + geom_area(data= subset (p, x > -30), aes(x=x, y=y), 
                                   fill= "light grey") +
  labs(title= "Prob. density", x = element_blank(), y = element_blank())  

gg.prob.b <- ggplot(out.df, (aes(Value))) + 
  geom_density(colour = grey(0.4), size = 1, fill = brewer.pal(3, "Dark2")[2]) +
  coord_cartesian(xlim = c(-100,y.up) ) +
  scale_x_continuous(breaks = tick.marks) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.text=element_text(size=axis.text.size), axis.title=element_text(size=22))

p <- ggplot_build(gg.prob.b)$data[[1]]
if(max(p$x, na.rm=T) <= -50) gg.base2b <- gg.prob.b
if(max(p$x, na.rm=T) > -50) gg.base2b <- gg.prob.b + geom_area(data= subset (p, x > -50), aes(x=x, y=y), 
                                   fill= "light grey") +
  labs(title=element_blank(), x = element_blank(), y = element_blank())  

gg.prob.c <- ggplot(out.df, (aes(Value))) + 
  geom_density(colour = grey(0.4), size = 1, fill = brewer.pal(3, "Dark2")[3]) +
  coord_cartesian(xlim = c(-100,y.up) ) +
  scale_x_continuous(breaks = tick.marks) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.text=element_text(size=axis.text.size), axis.title=element_text(size=22))

p <- ggplot_build(gg.prob.c)$data[[1]]
if(max(p$x, na.rm=T) <= -70) gg.base2c <- gg.prob.c
if(max(p$x, na.rm=T) > -70) gg.base2c <- gg.prob.c + geom_area(data= subset (p, x > -70), aes(x=x, y=y), 
                                   fill= "light grey") +
  labs(title=element_blank(), y = element_blank(), 
       x="Percent change")  

g1 <- list(gg.base1, gg.base2a, gg.base2b, gg.base2c)

gg.base <- grid.arrange(grobs = g1, widths = c(2,1), 
                        layout_matrix = rbind(c(1,2), c(1,3), c(1,4)))


# ggsave (file = paste (stk, "/ProbPlots_", du.df$DU[1],".pdf", sep=""), 
#         plot=gg.base, width=8, height=8)
out.dir <- here::here(stk, scenario.name) #here(path.results,scenario.name)

ggsave (file = paste (out.dir, "/ProbPlots_", du.df$DU[1],".png", sep=""), 
        plot=gg.base, width=8, height=8)

