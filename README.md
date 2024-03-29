# MetricsCOSEWIC-Stan
Stan code to run probability of decline metrics over 3 generations for COSEWIC, applied to example stocks of Sockeye Salmon

Author: Carrie Holt
Contact: carrie.holt@dfo-mpo.gc.ca


Next steps:
<!--(1) Revise Stan code to estimate missing abundance values in some years

(2) Standardize data so that sigma of normal priors on slope and yi can be specified in standard units. See https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

(3) Set sigmas for priors on slope and yi in those standardized units, and evaluate impact of information in priors on posterior distribution of % change

(4) Try alternative distributions for prior on residual variance, e.g., exponential (possibly uniform, half-cauchy and/or half normal), and evaluate impacts on posterior distribution of % change

(5) Consider setting priors based on the prior predictive distributions (e.g., as described by Wesner and Pomeranz 2021, https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.3739)

(6) Consider robust regression by changing normal distribution in likelihood to a t-distribution in the Stan code
-->
(1) Add stan code to R package (see example here: https://github.com/luke-a-rogers/mmmstan)

(2) Run Stan through the wrapper brms directly from R, and then 

(3) Evaluate impacts of missing values at beginning, throughout, and end of the time-series on posterior distribution of % change  (in simulation)
