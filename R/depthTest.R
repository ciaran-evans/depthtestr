# sampA = reference sample
# sampB = new sample
# return the depth of each observation in sampB wrt sampA

depthFunHalfspace <- function(sampB, sampA){
  fda.usc::mdepth.TD(sampB, sampA)$dep
}

depthFunHalfspace_self <- function(sampA){
  sapply(1:nrow(sampA), function(i) fda.usc::mdepth.TD(matrix(sampA[i,], nrow=1), sampA[-i,])$dep)
}


depthFunMahalanobis <- function(sampB, sampA){
  fda.usc::mdepth.MhD(sampB, sampA)$dep
}

depthFunMahalanobis_self <- function(sampA){
  sapply(1:nrow(sampA), function(i) fda.usc::mdepth.MhD(matrix(sampA[i,], nrow=1), sampA[-i,])$dep)
}

depthFunPVB <- function(sampB, sampA){
  DepthProc::depthLocal(sampB, sampA)
}

depthFunPVB_self <- function(sampA){
  sapply(1:nrow(sampA), function(i) DepthProc::depthLocal(matrix(sampA[i,], nrow=1), sampA[-i,]))
}

lpDepthObs <- function(sampB_obs, sampA){
  new_dists <- sqrt(rowSums((sweep(sampA, 2, sampB_obs, "-"))^2))
  return(1/(1 + mean(new_dists)))
}

depthFunLp <- function(sampB, sampA){
  apply(sampB, 1, function(x) lpDepthObs(x, sampA))
}

depthFunLp_self <- function(sampA){
  dists <- as.matrix(stats::dist(sampA))
  N1 <- nrow(sampA)
  sapply(1:N1, function(i) 1/(1 + mean(dists[i,-i])))
}

pdDepthObs <- function(sampB_obs, sampA, dists){
  N1 = nrow(sampA)

  new_dists <- sqrt(rowSums((sweep(sampA, 2, sampB_obs, "-"))^2))

  return(1/(N1 * (N1 - 1)) * pdDepthDists(new_dists, dists, N1))
}


depthFunPD <- function(sampB, sampA){
  dists <- as.matrix(stats::dist(sampA))
  apply(sampB, 1, function(x) pdDepthObs(x, sampA, dists))
}

depthFunPD_self <- function(sampA){
  dists <- as.matrix(stats::dist(sampA))
  N1 <- nrow(sampA)
  sapply(1:N1, function(i) {1/((N1-1) * (N1 - 2)) * pdDepthDists(dists[i,-i], dists[-i,-i], N1 - 1)})
}

lcdDepthObs <- function(sampB_obs, sampA, dists){
  N1 = nrow(sampA)

  new_dists <- sqrt(rowSums((sweep(sampA, 2, sampB_obs, "-"))^2))

  return(1/(N1 * (N1 - 1)) * lcdDepthDists(new_dists, dists, N1))
}

depthFunLCD <- function(sampB, sampA){
  dists <- as.matrix(stats::dist(sampA))
  apply(sampB, 1, function(x) lcdDepthObs(x, sampA, dists))
}

depthFunLCD_self <- function(sampA){
  dists <- as.matrix(stats::dist(sampA))
  N1 <- nrow(sampA)
  sapply(1:N1, function(i) {1/((N1-1) * (N1 - 2)) * lcdDepthDists(dists[i,-i], dists[-i,-i], N1 - 1)})
}

#' Test for a difference in distribution
#'
#' @param sample1 A matrix or data frame containing the first sample
#' @param sample2 A matrix or data frame containing the second sample
#' @param method The depth method to use (one of "halfspace", "mahalanobis", "pvb", "lp", "pd", "lcd", or "custom")
#' @param depthFun A depth function which takes two datasets and returns the depths
#'
#' @return The p-value and test statistic for a test of difference in distribution
#' @export
#'
#' @examples
#' library(mvtnorm)
#' samp1 <- rmvnorm(100, mean = rep(0, 10))
#' samp2 <- rmvnorm(100, mean = rep(1, 10))
#' depthTest(samp1, samp2, "mahalanobis")
depthTest <- function(sample1, sample2, method, depthFun = NULL){
  sample1 <- as.matrix(sample1)
  sample2 <- as.matrix(sample2)

  possible_methods <- c("halfspace", "mahalanobis", "pvb",
                        "lp", "pd", "lcd", "custom")
  if(!(method %in% possible_methods)){
    stop("Invalid depth method", paste("", method))
  }
  if(method == "custom"){
    if(is.null(depthFun)){
      stop("Please specify a custom depth function")
    }
  }

  if(method == "halfspace"){
    depthFun <- depthFunHalfspace
    depths_sample2 <- depthFun(sample2, sample1)
    depths_sample1 <- depthFunHalfspace_self(sample1)
  } else if(method == "mahalanobis"){
    depthFun <- depthFunMahalanobis
    depths_sample2 <- depthFun(sample2, sample1)
    depths_sample1 <- depthFunMahalanobis_self(sample1)
  } else if(method == "pvb"){
    depthFun <- depthFunPVB
    depths_sample2 <- depthFun(sample2, sample1)
    depths_sample1 <- depthFunPVB_self(sample1)
  } else if(method == "lp") {
    depthFun <- depthFunLp
    depths_sample2 <- depthFun(sample2, sample1)
    depths_sample1 <- depthFunLp_self(sample1)
  } else if(method == "pd"){
    depthFun <- depthFunPD
    depths_sample2 <- depthFun(sample2, sample1)
    depths_sample1 <- depthFunPD_self(sample1)
  } else if(method == "lcd"){
    depthFun <- depthFunLCD
    depths_sample2 <- depthFun(sample2, sample1)
    depths_sample1 <- depthFunLCD_self(sample1)
  } else {
    depths_sample2 <- depthFun(sample2, sample1)
    depths_sample1 <- depthFun(sample1, sample1)
  }

  # depths_sample2 <- depthFun(sample2, sample1)
  # depths_sample1 <- depthFun(sample1, sample1)

  r_stats <- sapply(depths_sample2, function(x) mean(depths_sample1 <= x))
  q_stat <- mean(r_stats)

  # get the approximate variance of the (unnormalized) test statistic
  N1 = nrow(sample1)
  N2 = nrow(sample2)
  q_var <- 1/((N1 + N2)*12*(N1/(N1 + N2))*(N2/(N1 + N2)))

  teststat <- (q_stat - 0.5)/sqrt(q_var)
  pval <- 2*stats::pnorm(abs(teststat), lower.tail=F)

  output <- list(teststat = teststat,
                 pval = pval)
  return(output)
}
