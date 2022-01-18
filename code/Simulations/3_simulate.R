#### 9. Set up parallel backend (optional) ####

library(foreach) # (Code built under version 1.5.1) 
library(doParallel) # (Code built under version 1.0.16) 
registerDoParallel(cores = detectCores() - 1)


#### 10. Run simulations ####

# Use current date and time to create pseudo-random seed
timestamp <- gsub(x = gsub(x = Sys.time(), pattern = ":", replacement = "-"),
                  pattern = " ", replacement = "_")
base_seed <- as.numeric(gsub("[^0-9.]", "", timestamp)) %% 1e6
base_seed <- base_seed + 
  as.numeric(paste(rev(strsplit(as.character(base_seed), "")[[1]]), collapse = ""))

# Code may take several minutes to run, depending on the chosen parameters.
# For non-parallel computing, use %do% instead of %dopar%
SIM_RESULTS <- foreach(i = 1:reps, .errorhandling = "pass") %dopar% {
  set.seed(base_seed + i)
  t <- as.numeric(seq(from = ISOdate(2000, 1, 1, 0), 
                      by = paste(60*60*24/freq, "sec"), 
                      length.out = freq * duration))
  
  DATA  <- lapply(MODELS, function(mod) {ctmm::simulate(mod, t = t)})
  GUESS <- lapply(DATA, function(j){ctmm::ctmm.guess(j, interactive = FALSE)})
  FITS  <- lapply(1:2, function(i) ctmm::ctmm.fit(DATA[[i]], GUESS[[i]]))
  UDS   <- tryCatch({ctmm::akde(DATA, FITS, debias = debias.akde)}, 
                    error = function(e) {return(NA)})
  
  if(!is.na(UDS)) {
    OVERLAP <- list(BC = overlap(UDS, method = "BC", debias = FALSE, level = level),
                    BC_DEBIASED = overlap(UDS, method = "BC", debias = TRUE, level = level),
                    APE = overlap(UDS, method = "APE", debias = FALSE, level = level),
                    APE_DEBIASED = overlap(UDS, method = "APE", debias = TRUE, level = level),
                    A12 = overlap(UDS, method = "A12", level = level))
    
    OVERLAP[["UDOI"]] = OVERLAP[["APE"]] * OVERLAP[["A12"]]
    OVERLAP[["UDOI_DEBIASED"]] = OVERLAP[["APE_DEBIASED"]] * OVERLAP[["A12"]]
    
  } else {
    OVERLAP <- NA
  }
  
  list(base_seed = base_seed, rep = i, rep_seed = base_seed + i, 
       duration = duration, freq = freq,  t = t,
       level = level, debias.akde = debias.akde, 
       MOD_TYPE = MOD_TYPE, MODELS = MODELS, DATA = DATA, 
       GUESS = GUESS, FITS = FITS, UDS = UDS, OVERLAP = OVERLAP)
}

stopImplicitCluster()