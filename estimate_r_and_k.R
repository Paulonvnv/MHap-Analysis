library(Rcpp)
sourceCpp('hmmloglikelihood.cpp')
source('fs_checks.R')

estimate_r_and_k <- function(fs, ds, Ys, epsilon = 0.001, rho = 7.4 * 10 ^ (-7),
                             kinit = 50, rinit = 0.5, warn_fs = TRUE) {
  
  # Convert to a into matrix if not already (loglikelihood_cpp expects a matrix)
  # and perform some checks
  if (!is.matrix(fs)) fs <- as.matrix(fs)
  # If no errors, returns numeric logic for fs zero/non-zero
  # fs_checks_return <- fs_checks(fs, warn = warn_fs, do_return = TRUE)
  # 
  # # Check Ys for NAs
  # if(any(is.na(Ys))) stop("Missing values detected in Ys.\n  Please remove and recompute ds accordingly.")
  # 
  # # Check for alleles that are permissible only if epsilon exceeds zero.
  # # These alleles will break the hmmloglikelihood.cpp code if epsilon is zero.
  # problem <- "Some per-marker allele counts exceed per-marker non-zero allele frequencies."
  # fs_checks_return$fs01 <- 1 * (fs > fs_checks_return$non_zero_fs_lb)
  # Kts <- rowSums(fs_checks_return$fs01)
  # if (any(sapply(1:nrow(Ys), function(i) Ys[i] > (Kts[i]-1)))) {
  #   if (epsilon > 0) warning (paste0(problem, " Data are permissible due to non-zero epsilon."))
  #   else if (epsilon == 0)  stop(paste0(problem, " Data are incompatible with zero-valued epsilon."))
  # }
  
  # Define the function to pass to optim()
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, fs, ds, epsilon, rho)
  
  # Optimise the negative log likelihood
  optimization <- optim(par = c(kinit, rinit), fn = function(x) - ll(x[1], x[2]))
  
  # Extract and name estimates
  rkhats <- c("khat" = optimization$par[1], "rhat" = optimization$par[2])
  
  if (all(rkhats == c(kinit, rinit))) {
    warning("optimization has returned initial parameter values. Data are possibly uniformative.")
  }
  
  # End of function
  return(rkhats)
}
