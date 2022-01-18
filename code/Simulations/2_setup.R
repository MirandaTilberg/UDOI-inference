#### To run this simulation, first load the functions contained in "1_overlap_fns.R":
# source("../1_overlap_fns.R") # source file manually or with this code

library(ctmm) # (Code built under version 0.6.0)
sessionInfo(package = "ctmm")$otherPkgs$ctmm$Version # print current version


#### 7. Set process parameters ####

reps <- 10              # number of replicates per duration
level <- 0.95           # confidence level for individual HR areas
debias.akde <- TRUE     # whether to debias AKDE area estimate
duration <- 2^5         # number of days for which to simulate data
freq <- 8               # number of observations per day


#### 8. Select model by autocorrelation structure and degree of overlap ####

# Degrees of overlap are coded "LOW" (BC=.15), "MOD" (BC=.5) and "HI" (BC=.85)
OVER_DEG <- "MOD"       # select degree of overlap

# Set mean and variance parameters
sigma <- 1 %#% "km"     # isotropic models used for demonstration
a <- switch(OVER_DEG,   # pre-determined x-coordinate of mean parameters
            "LOW" = 1.947881, 
            "MOD" = 1.17741, 
            "HI" = .570121) * sqrt(sigma)

# Specify pair of OUF models for simulations:
MOD_TYPE <- paste("OUF", OVER_DEG, sep = "-")
MODELS <- list(ctmm(mu = c(-a, 0), sigma = sigma,
                    tau = c(1 %#% "day", .2 %#% "day")),
               ctmm(mu = c(a, 0), sigma = sigma,
                    tau = c(1 %#% "day", .2 %#% "day")))


# # To simulate from OU models, supply only one value for `tau`:
# MOD_TYPE <- paste("OU", OVER_DEG, sep = "-")
# MODELS <- list(ctmm(mu = c(-a, 0), sigma = sigma, tau = c(1 %#% "day")),
#                ctmm(mu = c(a, 0), sigma = sigma, tau = c(1 %#% "day")))


# # To simulate from IID models, do not specify a value for `tau`:
# MOD_TYPE <- paste("IID", OVER_DEG, sep = "-")
# MODELS <- list(ctmm(mu = c(-a, 0), sigma = sigma),
#                ctmm(mu = c(a, 0), sigma = sigma))