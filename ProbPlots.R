
#==============================================================================
# Code to run stan code over stocks to derive density plots
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
# 
# SR_Sample %>% filter(Year>(2017-13)) %>%
#   ggplot(aes(Year, Spn)) + geom_line(aes(colour=factor(Stock)))


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

# Choose a stock
stk <- "Stock2"#"Stock18"#"Stock2"##"Stock3"
du.label <- stk
du.df.long <- du.df.long %>% filter(DU==stk) %>% mutate(logAbd=log(Abd)) 

# Add an NA (For checking only)
#du.df.long <- du.df.long %>% mutate(logAbd = replace(logAbd, Year==2008, NA))

# Shorten data set to years used to estimate trends
du.df <- du.df.long %>% filter(Year > (calc.year - yrs.window) & 
                                 Year <= calc.year)

#=============================================================================
# RUN STAN
#=============================================================================
# Run STAN with standardized data, exp prior on var, and 2.5 sigma priors on 
# slope and yi (DEFAULT)

stan.out <- run.stan(du.label=stk, du.df=du.df, yrs.window=yrs.window, 
         standardize.data = TRUE,  
         scenario.name = "Standardized-ExpSigma", 
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

# Text for labeling probabilities at -70%, -50% and -30% declines on plot
text.df <- data.frame(x = c(-100,-100,-100), 
                      Value = round(stan.out$probdecl, 0)[,2]/100, 
                      lab = round(stan.out$probdecl, 0)[,2]/100 )

# Plot cumulative distribution
gg.base1 <- ggplot(out.df, (aes(Value))) +#, colour=as.factor(PriorSigma)))) + 
  stat_ecdf(geom = "step", pad = TRUE, size = 2, colour=grey(0.4)) + 
  labs(title= "Cumulative probability", y = element_blank(), 
       x="Percent change") +
  scale_colour_brewer(element_blank(), palette="Dark2") + 
  coord_cartesian(xlim = c(-100,y.up) )+
  scale_x_continuous(breaks = tick.marks) +
  geom_text(aes(x= x, y = Value, label=lab), data=text.df,  colour=grey(0.5), nudge_x=0, 
            nudge_y=0.02, size=4) + 
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
  theme(legend.position="bottom")

# Plot three probability densities for -70, -50, and -30 thresholds
gg.prob.a <- ggplot(out.df, (aes(Value))) ++ 
  geom_density(colour = grey(0.4), size = 1, fill = brewer.pal(3, "Dark2")[1]) +
  coord_cartesian(xlim = c(-100,y.up) )+
  scale_x_continuous(breaks = tick.marks) +
  scale_y_continuous(breaks = NULL) 
p <- ggplot_build(gg.prob.a)$data[[1]]
gg.base2a <- gg.prob.a + geom_area(data= subset (p, x > -30), aes(x=x, y=y), 
                                   fill= "light grey") +
  labs(title= "Probability density", x = element_blank(), y = element_blank())  

gg.prob.b <- ggplot(out.df, (aes(Value))) + 
  geom_density(colour = grey(0.4), size = 1, fill = brewer.pal(3, "Dark2")[2]) +
  coord_cartesian(xlim = c(-100,y.up) ) +
  scale_x_continuous(breaks = tick.marks) +
  scale_y_continuous(breaks = NULL)
p <- ggplot_build(gg.prob.b)$data[[1]]
gg.base2b <- gg.prob.b + geom_area(data= subset (p, x > -50), aes(x=x, y=y), 
                                   fill= "light grey") +
  labs(title=element_blank(), x = element_blank(), y = element_blank())  

gg.prob.c <- ggplot(out.df, (aes(Value))) + 
  geom_density(colour = grey(0.4), size = 1, fill = brewer.pal(3, "Dark2")[3]) +
  coord_cartesian(xlim = c(-100,y.up) ) +
  scale_x_continuous(breaks = tick.marks) +
  scale_y_continuous(breaks = NULL)
p <- ggplot_build(gg.prob.c)$data[[1]]
gg.base2c <- gg.prob.c + geom_area(data= subset (p, x > -70), aes(x=x, y=y), 
                                   fill= "light grey") +
  labs(title=element_blank(), y = element_blank(), 
       x="Percent change")  

g1 <- list(gg.base1, gg.base2a, gg.base2b, gg.base2c)

gg.base <- grid.arrange(grobs = g1, widths = c(2,1), 
                        layout_matrix = rbind(c(1,2), c(1,3), c(1,4)))


ggsave (file = paste (stk, "/ProbPlots_", du.df$DU[1],".pdf", sep=""), 
        plot=gg.base, width=8, height=8)
ggsave (file = paste (stk, "/ProbPlots_", du.df$DU[1],".png", sep=""), 
        plot=gg.base, width=8, height=8)

