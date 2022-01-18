# Begin by sourcing this file to load `ctmm` and overwrite some `ctmm` functions
#   with local versions.

# Functions are modified from `ctmm` source code, version 0.6.0: 
# (Jan 8th, 2021) https://github.com/ctmm-initiative/ctmm/tree/4c0b7df669269087443c6907b81df6f5d6de5f13/R
# (current) https://github.com/ctmm-initiative/ctmm/tree/master/R 

# Authors' comments from Github source code (link above) are retained below.
# Changes from `ctmm` source code to facilitate UDOI are marked with "--MT"
########################

### 0. Setup ####

if(!require(ctmm)) {install.packages("ctmm")}
library(ctmm) # (Code requires version 0.6.0 or later)
sessionInfo(package = "ctmm")$otherPkgs$ctmm$Version # print current version


#### 1. same.grids() ####

same.grids <- function(UD) {
  n <- length(UD)
  SUB <- ctmm:::grid.intersection(UD)
  for (i in 1:n) {
    UD[[i]]$r$x <- UD[[i]]$r$x[SUB[[i]]$x]
    UD[[i]]$r$y <- UD[[i]]$r$y[SUB[[i]]$y]
    UD[[i]]$PDF <- UD[[i]]$PDF[SUB[[i]]$x, SUB[[i]]$y]
    UD[[i]]$CDF <- UD[[i]]$CDF[SUB[[i]]$x, SUB[[i]]$y] # include CDF --MT
  }
  return(UD)
}


#### 2. APE_D() #### 
# new distance function for APE --MT

# square distance between stationary Gaussian distributions
APE_D <- function(CTMM) {
  CTMM1 <- CTMM[[1]]
  CTMM2 <- CTMM[[2]]
  
  sigma <- (CTMM1$sigma + CTMM2$sigma)/2
  mu <- CTMM1$mu[1,] - CTMM2$mu[1,]
  
  D <- as.numeric(mu %*% ctmm:::PDsolve(sigma) %*% mu)/4 + log(det(4*pi*sigma))/2
  
  return(D)
}


#### 3. CI.UD() ####

CI.UD <- function(object, level.UD = 0.95, level = 0.95, P = FALSE) {
  if (is.null(object$DOF.area) && P) {
    names(level.UD) <- ctmm:::NAMES.CI[2] # point estimate
    return(level.UD)
  }
  
  dV <- prod(object$dr)
  
  # point estimate
  area <- sum(object$CDF <= level.UD) * dV
  names(area) <- ctmm:::NAMES.CI[2] # point estimate
  
  # chi square approximation of uncertainty
  if (!is.null(object$DOF.area)) {
    area <- ctmm:::chisq.ci(area, DOF = 2 * object$DOF.area[1], 
                            alpha = 1 - level)
    names(area) <- ctmm:::NAMES.CI
  }
  
  if (!P) { return(area) }
  
  # probabilities associated with these areas
  P <- round(area/dV)
  
  # fix lower bound for all values of P --MT
  P <- pmax(P, 1)
  # fix upper bound to not overflow
  P[3] <- min(length(object$CDF), P[3])
  if (P[3]==length(object$CDF)) { warning("Outer contour extends beyond raster.") }
  
  P <- sort(object$CDF, method = "quick")[P]
  
  # recorrect point estimate level
  P[2] <- level.UD
  
  names(P) <- ctmm:::NAMES.CI
  return(P)
}


#### 4. overlap.ctmm() #### 

overlap.ctmm <- function(
  object, method = c("BC", "UDOI", "APE", "Bhattacharyya", "Euclidean", "Mahalanobis"),
  debias = NULL, level = 0.95, COV = TRUE, distance = FALSE, ...) 
{
  CTMM1 <- object[[1]]
  CTMM2 <- object[[2]]
  DIM <- length(CTMM1$axes)
  
  # Choose appropriate distance function --MT
  method <- match.arg(method)
  Dfunc <- switch(method, 
                  "Bhattacharyya" = ctmm:::BhattacharyyaD,
                  "BC" = ctmm:::BhattacharyyaD,
                  "Mahanobis" = ctmm:::MahalanobisD,
                  "Euclidean" = ctmm:::EuclideanD,
                  "UDOI" = APE_D,
                  "APE" = APE_D)
  
  STUFF <- ctmm:::gauss.comp(Dfunc, object, COV = COV)
  MLE <- c(STUFF$MLE)
  VAR <- c(STUFF$COV)
  # this quantity is roughly chi-square
  DOF <- 2 * MLE^2/VAR
  
  # approximate debiasing, correct for IID, equal covariance, REML
  ########################
  mu <- CTMM1$mu[1, ] - CTMM2$mu[1, ]
  COV.mu <- CTMM1$COV.mu + CTMM2$COV.mu
  
  if (method == "Euclidean") {
    sigma <- diag(1, DIM)
  } else {
    sigma <- (CTMM1$sigma + CTMM2$sigma)/2 # AM
  }
  
  # trace variances
  s0 <- mean(diag(sigma))
  s1 <- mean(diag(CTMM1$sigma))
  s2 <- mean(diag(CTMM2$sigma))
  
  # approximate average Wishart DOFs
  n1 <- ctmm:::DOF.area(CTMM1)
  n2 <- ctmm:::DOF.area(CTMM2)
  # using mean variance - additive & rotationally invariant
  n0 <- 4 * s0^2/(s1^2/n1 + s2^2/n2)
  # dim cancels out
  
  # hard clamp before soft clamp
  n1 <- ctmm:::clamp(n1, 1, Inf)
  n2 <- ctmm:::clamp(n2, 1, Inf)
  n0 <- ctmm:::clamp(n0, 2, Inf)
  
  # clamp the DOF not to diverge <=DIM+1
  n0 <- ctmm:::soft.clamp(n0, DIM)
  n1 <- ctmm:::soft.clamp(n1, DIM)
  n2 <- ctmm:::soft.clamp(n2, DIM)
  
  # expectation value of log det Wishart
  ElogW <- function(s, n) {
    log(det(s)) + ctmm:::mpsigamma(n/2, dim = DIM) - DIM * log(n/2)
  }
  
  # inverse Wishart expectation value pre-factor
  BIAS <- n0/(n0 - DIM - 1)
  
  if (method == "Euclidean") { BIAS <- 0 } # don't include this term
  
  # mean terms
  BIAS <- sum(diag((BIAS * ctmm:::outer(mu) + COV.mu) %*% ctmm:::PDsolve(sigma)))
  
  if (method == "Bhattacharyya" | method == "BC") {
    BIAS <- BIAS/8
    # AMGM covariance terms
    BIAS <- BIAS + max(ElogW(sigma, n0)/2 - ElogW(CTMM1$sigma, n1)/4 - 
                         ElogW(CTMM2$sigma, n2)/4, 0)
    # this is actually the expectation value?
  } else if (method == "UDOI" | method == "APE") { 
    # compute for APE_D --MT
    BIAS <- BIAS/4
    BIAS <- BIAS + max(ElogW(4*pi*sigma,n0)/2, 0)
  }
  
  # relative bias instead of absolute bias
  BIAS <- BIAS/MLE
  # would subtract off estimate to get absolute bias
  
  # error corrections
  BIAS <- as.numeric(BIAS)
  if(MLE==0) { BIAS <- 1 }
  #####################
  
  if (level) {
    
    # if 'debias = NULL' (new default), choose best option for given metric --MT
    if(is.null(debias)) {
      if(method == "UDOI" | method == "APE") {
        debias <- FALSE
      } else {
        debias <- TRUE
      }
    }
    
    if (debias) { MLE <- MLE/BIAS }
    
    CI <- ctmm:::chisq.ci(MLE, VAR = VAR, alpha = 1 - level)
    if (distance) { return(CI) } # return distance
    
    # transform from (square) distance to overlap measure
    CI <- exp(-rev(CI))
    names(CI) <- ctmm:::NAMES.CI
    
    return(CI)
  }
  else { # return BD or APE_D ingredients --MT
    return(list(MLE = MLE, VAR = VAR, DOF = DOF, BIAS = BIAS))
  }
}


#### 5. overlap.UD() ####

overlap.UD <- function(
  object, method = c("BC", "UDOI", "APE", "A12",
                     "Bhattacharyya", "Euclidean", "Mahalanobis"), 
  debias = NULL, level = 0.95, level.UD = 0.95, ...) 
{
  method <- match.arg(method) # check method argument --MT
  
  CTMM <- list(attr(object[[1]], "CTMM"), attr(object[[2]], "CTMM"))
  type <- c(attr(object[[1]], "type"), attr(object[[2]], "type"))
  type <- type[type != "range"]
  if (length(type)) { stop(type, " overlap is not generally meaningful, biologically.") }
  
  # check resolution and subset to overlapping grid
  object <- same.grids(object)
  # can now be null mass
  
  dr <- object[[1]]$dr
  dA <- prod(dr)
  
  
  # compute overlap and area of overlap (if applicable) --MT
  OVER <- object[[1]]$PDF * object[[2]]$PDF
  if (!is.null(OVER)) {
    if (method == "Bhattacharyya" | method == "BC") {
      
      OVER <- sum(sqrt(OVER)) * dA # BC
      OVER_AREA <- 1 # no scaling necessary
      
    } else if (method == "UDOI" | method == "APE" | method == "A12") {
      
      OVER <- sum(OVER) * dA # APE
      
      # compute area of overlap --MT
      if(length(object[[1]]$PDF)) {
        
        P1 <- CI.UD(object[[1]], level.UD = level.UD, level = level, P = TRUE)
        P2 <- CI.UD(object[[2]], level.UD = level.UD, level = level, P = TRUE)
        
        OVER_AREA <- sapply(1:3, function(l) {
          sum(object[[1]]$CDF <= P1[l] & object[[2]]$CDF <= P2[l]) * dA
        })
      } else {
        OVER_AREA <- rep(0, 3)
      }
    }
  }
  
  if(method == "A12") { return(OVER_AREA) } # return area of overlap --MT
  
  if (!is.null(CTMM)) {
    # calculate Gaussian overlap distance^2 variance, bias, etc.
    CI <- overlap.ctmm(CTMM, method = method, level = FALSE) # use local version --MT
    
    # Bhattacharyya distances or APE_D --MT
    D <- -log(OVER)
    
    
    # if 'debias = NULL' (new default), choose best option for given metric --MT
    if(is.null(debias)) {
      if(method == "UDOI" | method == "APE") {
        debias <- FALSE
      } else {
        debias <- TRUE
      }
    }
    
    # relative debias
    if (debias) {
      D <- D/CI$BIAS 
    }
    
    # calculate new distance^2 with KDE point estimate
    CI <- ctmm:::chisq.ci(D, VAR = CI$VAR, alpha = 1 - level)
    
    # transform from (square) distance to overlap measure
    OVER <- exp(-rev(CI))
  } else {
    OVER <- c(NA, OVER, NA)
  }
  
  names(OVER) <- ctmm:::NAMES.CI
  if(method == "APE") { return(OVER) }
  return(OVER * OVER_AREA)
}


#### 6. overlap() ####

overlap <- function(
  object, method = c("BC", "UDOI", "APE", "A12",
                     "Bhattacharyya", "Euclidean", "Mahalanobis"), 
  level = 0.95, level.UD = 0.95, debias = NULL, ...) 
{
  
  # ctmm:::check.projections(object) # turn off for now to avoid error --MT
  method <- match.arg(method) # check method argument --MT
  
  # if 'debias = NULL' (new default), choose best option for given metric --MT
  if(is.null(debias)) {
    if(method == "UDOI" | method == "APE") {
      debias <- FALSE
    } else {
      debias <- TRUE
    }
  }
  
  # warn against debiasing APE --MT
  if(debias & (method == "UDOI" | method == "APE")) {
    warning("Debiasing is not recommended for the APE or UDOI.")
  }
  
  CLASS <- class(object[[1]])[1]
  if (CLASS == "ctmm") {
    OverlapFun <- overlap.ctmm
  } else if (CLASS == "UD") {
    OverlapFun <- overlap.UD
  } else {
    stop(CLASS, " object class not supported by overlap.")
  }
  
  n <- length(object)
  OVER <- array(0, c(n, n, 3))
  
  # tabulate overlaps
  for (i in 1:n) {
    for (j in i:n) {
      if (j <= n) {
        # make matrices symmetric at assignment --MT
        OVER[i, j, ] <- OVER[j, i, ] <- OverlapFun(object[c(i, j)], method = method, 
                                                   level = level, debias = debias, 
                                                   level.UD = level.UD, ...)
      }
    }
  }
  
  # fix diagonals for BC only --MT
  if(method == "Bhattacharyya" | method == "BC") {
    diag(OVER[, , 1]) <- diag(OVER[, , 2]) <- diag(OVER[, , 3]) <- 1
  }
  
  dimnames(OVER) <- list(names(object), names(object), ctmm:::NAMES.CI)
  return(OVER)
}
