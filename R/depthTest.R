# sampA = reference sample
# sampB = new sample
# return the depth of each observation in sampB wrt sampA


depthFunHalfspace <- function(sampB, sampA){
  fda.usc::mdepth.TD(sampB, sampA)$dep
}

depthFunMahalanobis <- function(sampB, sampA){
  fda.usc::mdepth.MhD(sampB, sampA)$dep
}

depthFunPVB <- function(sampB, sampA){
  DepthProc::depthLocal(sampB, sampA)
}

lpDepthObs <- function(sampB_obs, sampA){
  new_dists <- sqrt(rowSums((sweep(sampA, 2, sampB_obs, "-"))^2))
  return(1/(1 + mean(new_dists)))
}

depthFunLp <- function(sampB, sampA){
  apply(sampB, 1, function(x) lpDepthObs(x, sampA))
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

lcdDepthObs <- function(sampB_obs, sampA, dists){
  N1 = nrow(sampA)

  new_dists <- sqrt(rowSums((sweep(sampA, 2, sampB_obs, "-"))^2))

  return(1/(N1 * (N1 - 1)) * lcdDepthDists(new_dists, dists, N1))
}

depthFunLCD <- function(sampB, sampA){
  dists <- as.matrix(stats::dist(sampA))
  apply(sampB, 1, function(x) lcdDepthObs(x, sampA, dists))
}


#' Test for a difference in distribution
#'
#' @param sample1 A matrix or data frame containing the first sample
#' @param sample2 A matrix or data frame containing the second sample
#' @param method The depth method to use (one of "halfspace", "mahalanobis", "pvb", "lp", "pd", "lcd", or "custom")
#' @param depthFun A depth function which takes two datasets and returns the depths
#' @param loo_correction Whether to use a leave-one-out correction when computing the depth of a sample to itself
#'
#' @return The p-value and test statistic for a test of difference in distribution
#' @export
#'
#' @examples
#' library(mvtnorm)
#' samp1 <- rmvnorm(100, mean = rep(0, 10))
#' samp2 <- rmvnorm(100, mean = rep(1, 10))
#' depthTest(samp1, samp2, "mahalanobis")
depthTest <- function(sample1, sample2, method,
                      depthFun = NULL, loo_correction = TRUE){
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
  if(method != "custom"){
    message(paste("Using the", method, "depth function"))
  }

  depthFun <- switch(method,
                     "halfspace" = depthFunHalfspace,
                     "mahalanobis" = depthFunMahalanobis,
                     "pvb" = depthFunPVB,
                     "lp" = depthFunLp,
                     "pd" = depthFunPD,
                     "lcd" = depthFunLCD,
                     "custom" = depthFun)

  depths_sample2 <- depthFun(sample2, sample1)

  if(loo_correction){
    depths_sample1 <- sapply(1:nrow(sample1),
                             function(i) depthFun(matrix(sample1[i,], nrow=1),
                                                  sample1[-i,]))
  } else {
    depths_sample1 <- depthFun(sample1, sample1)
  }

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
