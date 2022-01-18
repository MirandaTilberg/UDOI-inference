library(tidyr) # Code built under version 1.1.2
sessionInfo(package = "tidyr")$otherPkgs$tidyr$Version # print current version

library(purrr) # Code built under version 0.3.4
sessionInfo(package = "purrr")$otherPkgs$purrr$Version # print current version


#### 11. Aggregate overlap results into dataframe ####

OVERLAP_RESULTS <- purrr::map_dfr(SIM_RESULTS, function(x) {
  if(!anyNA(x$OVERLAP)) {
    m <- sapply(x$OVERLAP, function(l) {l[1,2,]})
    df <- as.data.frame(t(m))
    df$duration <- x$duration
    df$rep <- x$rep
    df$level <- x$level
    df <- tibble::rownames_to_column(df, var = "type")
  } else {
    # dummy data frame in case of errors
    df <- data.frame(type = "UDOI_BIASED", low = NA, est = NA, high = NA,
                     duration = x$duration, rep = x$rep)
  }
}) %>% tidyr::pivot_wider(., id_cols = c("duration", "rep"),
                          names_from = type, values_from = ctmm:::NAMES.CI,
                          names_glue = "{type}_{.value}")


#### 12. Aggregate model- and area-based results into dataframe ####

AREA_RESULTS <- purrr::map_dfr(1:length(SIM_RESULTS), function(i){
  SR_F <- SIM_RESULTS[[i]]$FITS
  SR_U <- SIM_RESULTS[[i]]$UDS
  
  df <- purrr::map_dfr(1:2, function(j) {
    x <- summary(SR_F[[j]], units = FALSE)
    data.frame(model = j,
               mean = SR_F[[j]]$mu, 
               sigma.diag = SR_F[[j]]$sigma[1],
               sigma.offdiag = SR_F[[j]]$sigma[2],
               isotropic = SR_F[[j]]$isotropic,
               COV.mu.11 = SR_F[[j]]$COV.mu[1],
               COV.mu.offdiag = SR_F[[j]]$COV.mu[2],
               COV.mu.22 = SR_F[[j]]$COV.mu[4],
               tau.p = SR_F[[j]]$tau[1],
               tau.v = SR_F[[j]]$tau[2],
               method = SR_F[[j]]$method,
               name = x$name,
               DOF.mean = x$DOF[1],
               DOF.area = x$DOF[2],
               DOF.speed = x$DOF[3],
               CI.area.low = x$CI[1,1],
               CI.area.est = x$CI[1,2],
               CI.area.high = x$CI[1,3],
               CI.tau.p.low = x$CI[2,1],
               CI.tau.p.est = x$CI[2,2],
               CI.tau.p.high = x$CI[2,3],
               CI.tau.v.low = x$CI[3,1],
               CI.tau.v.est = x$CI[3,2],
               CI.tau.v.high = x$CI[3,3],
               dr.x = ifelse(all(is.na(SR_U)), NA, SR_U[[j]]$dr["x"]),
               dr.y = ifelse(all(is.na(SR_U)), NA, SR_U[[j]]$dr["y"]),
               dA = ifelse(all(is.na(SR_U)), NA, prod(SR_U[[j]]$dr)))
  })
  df$duration <- SIM_RESULTS[[i]]$duration
  df$rep <- SIM_RESULTS[[i]]$rep
  df$level <- SIM_RESULTS[[i]]$level
  rownames(df) <- NULL
  df
})

# sort columns of AREA_RESULTS for viewing (optional)
first_cols <- c("duration", "rep", "model", "level")
idx <- sapply(first_cols, function(x){which(x == names(AREA_RESULTS))})
AREA_RESULTS <- AREA_RESULTS[,c(idx, (1:ncol(AREA_RESULTS))[-idx])]
