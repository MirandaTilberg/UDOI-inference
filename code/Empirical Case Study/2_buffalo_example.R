#### To run this example, first load the functions contained in "1_overlap_fns.R":
# source("../1_overlap_fns.R") # source file manually or with this code

library(ctmm) # (Code built under version 0.6.0)
sessionInfo(package = "ctmm")$otherPkgs$ctmm$Version # print current version


# Load data
data("buffalo")


# Initial parameter guess
GUESS_buffalo <- lapply(buffalo, function(j){ctmm.guess(j, interactive = FALSE)})


# Fit models (may take several minutes)
set.seed(32121)
FITS_buffalo <- lapply(1:length(buffalo), function(i){
  ctmm.fit(buffalo[[i]], GUESS_buffalo[[i]])
})

names(FITS_buffalo) <- names(buffalo)


# Fit UDs (may take several minutes)
UDS_buffalo <- akde(buffalo, FITS_buffalo)


# Overlap results with 95% CIs
overlap(UDS_buffalo, method = "UDOI", debias = FALSE, level = .95)
overlap(UDS_buffalo, method = "APE", debias = FALSE, level = .95)
overlap(UDS_buffalo, method = "A12", level = .95)
overlap(UDS_buffalo, method = "BC", debias = TRUE, level = .95)

# Overlap results with 50% CIs
overlap(UDS_buffalo, method = "UDOI", debias = FALSE, level = .50)
overlap(UDS_buffalo, method = "APE", debias = FALSE, level = .50)
overlap(UDS_buffalo, method = "A12", level = .50)
overlap(UDS_buffalo, method = "BC", debias = TRUE, level = .50)

# Note: Debiasing has been shown to improve the BC (Winner, 2018) but does not
#   improve the APE or UDOI. Use 'debias = FALSE' for the APE or UDOI, and 
#   'debias = TRUE' otherwise.