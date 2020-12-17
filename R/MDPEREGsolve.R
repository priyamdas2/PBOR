#' Minimum Divergence Power Estimation (MDPE) Solver for regression models
#'
#' `MDPEREGsolve` function estimates the minimum divergence power estimator for regression
#' models with assumption that mean of dependent variable being a linear function of
#' independent variables. This function uses a Blackbox optimization technique, which is a
#' variant of Pattern search, known as Prior Based Optimization Routine (PBOR).
#'
#'
#' @param param1 Initial value of beta where distribution mean is given by beta*X.
#' @param param2 Initial value of second parameter : Normal (sigma).
#' @param fixed_param1 First fixed parameter values : Poisson (approx. sum upto).
#'  Default value is 10.
#' @param dist Distribution of sample : "normal", "poisson", "logistic" (with mean logit link)
#' . Default value is "normal".
#' @param Y_array Array of dependent variable (Y) values in sample.
#' @param X_mat Matrix of values of independent variables, each row represents data for
#' a subject.
#' @param alpha Parameter of Power estimator, value must be between 0 and 1.
#' @param center_vec Mean vector of prior normal distribution. It can be set at the
#' neighborhood where the true parameter value is expected.
#' @param sigma_vec The standard deviation vector of the prior normal distribution.
#' @param max_time Maximum time (in seconds) alloted for optimization. Default value is 3600.
#' @param default_ranges 'TRUE' for keeping the default ranges of the beta vector along
#' with sigma paramater (if applicable). In case solution needs to be found within a
#' specific interval, it should be set 'FALSE' and 'l_bound', 'u_bound' should be set
#' accordingly. Default value is TRUE.
#' @param l_bound Lower bound for parameter vector, to activate this bound, user must set
#' 'default_ranges' to be FALSE.
#' @param u_bound Upper bound for parameter vector, to activate this bound, user must set
#' 'default_ranges' to be FALSE.
#' @param rho Step Decay Rate with default value 2
#' @param phi Step size threshold, i.e., lower Bound Of Global Step Size. Default value is \eqn{10^{-6}}
#' @param no_runs Maximum number of runs to be executed, default value is 1000.
#' @param max_iter Max Number Of Iterations in each Run. Default Value is 10000.
#' @param s_init Initial Global Step Size. Default Value is 1.
#' @param tol_fun Termination Tolerance on the function value. Default Value is \eqn{10^{-6}}
#' @param tol_fun_2 Termination Tolerance on the difference of solutions in two consecutive
#' runs. Default Value is \eqn{10^{-6}}.
#'
#' @return A list of length two.
#'\itemize{
#'\item \code{theta:} The solution point where DPE is minimized.
#'\item \code{value:} Value of MDPE.
#'}
#'
#'@references Das, P. and Ghosh, A. (2021). \emph{PBOR paper}. Future Journal, 2, 63-74.
#'@author Priyam Das (maintainer) (\email{priyam_das@@hms.harvard.edu}), Abhik Ghosh.
#' @examples
#' ######### Normal Distribution (mean = X*beta) #################################
#' set.seed(1)
#'
#' p <- 5
#' n <- 100
#' param1 <- rnorm(p)
#' param2 <- 1
#' x0 <- c(param1,param2)
#' lb <- c(array(-Inf,p),0.1)
#' ub <- array(Inf,p+1)
#' beta_true <- rnorm(p)
#' X_mat <- matrix(rnorm(n*p,mean=0,sd=1), n, p)
#' mean_normal <- X_mat%*%beta_true
#' Y_array <- array(NA,n)
#' for(ii in 1:n)
#' {Y_array[ii] <- rnorm(1,mean_normal[ii],1)}
#' alpha <- 0.5
#' center_vec <- array(0,p+1)
#' sigma_vec <- array(1,p+1)
#' alpha <- 0.1
#' PBOR_result <- MDPEREGsolve(param1, param2, fixed_param1 = 10, dist = "normal", Y_array,
#'                             X_mat, alpha = alpha,center_vec = center_vec, sigma_vec = sigma_vec,
#'                             max_time = 60, default_ranges = FALSE, l_bound = lb, u_bound = ub)
#' PBOR_result$value
#' PBOR_result$theta
#'
#' ######### Logistic Distribution (mean = logit(X*beta)) #########################
#'
#' set.seed(1)
#' p <- 10
#' n <- 50
#' param1 <- rnorm(p)
#' x0 <- param1
#' lb <- array(-Inf,p)
#' ub <- array(Inf,p)
#' beta_true <- rnorm(p)
#' X_mat <- matrix(rnorm(n*p,mean=0,sd=1), n, p)
#' pp <- 1/(1+exp(-X_mat%*%beta_true))
#' Y_array <- array(NA,n)
#' for(ii in 1:n)
#' {Y_array[ii] <- rbinom(1,1,pp[ii])}
#' center_vec <- array(0,p)
#' sigma_vec <- array(10,p)
#' alpha <- 0.1
#' PBOR_result <- MDPEREGsolve(param1, dist = "logistic", Y_array = Y_array,
#'                             X_mat = X_mat, alpha = alpha,center_vec = center_vec, sigma_vec = sigma_vec,
#'                             max_time = 60, default_ranges = FALSE, l_bound = lb, u_bound = ub)
#' PBOR_result$value
#' PBOR_result$theta
#'
#' ######### Poisson Distribution (mean = exp(X*beta)) #########################
#'
#' set.seed(1)
#' p <- 5
#' n <- 75
#' param1 <- rnorm(p)
#' x0 <- param1
#' sum_upto <- 1000
#' lb <- array(-Inf,p)
#' ub <- array(Inf,p)
#' beta_true <- rnorm(p)
#' X_mat <- matrix(rnorm(n*p,mean=0,sd=1), n, p)
#' mean_poisson <- exp(X_mat%*%beta_true)
#' Y_array <- array(NA,n)
#' for(ii in 1:n)
#' {Y_array[ii] <- rpois(1,mean_poisson)}
#' center_vec <- array(0,p)
#' sigma_vec <- array(10,p)
#' alpha <- 0.1
#' PBOR_result <- MDPEREGsolve(param1, param2, fixed_param1 = sum_upto, dist = "poisson", Y_array = Y_array,
#'                             X_mat = X_mat, alpha = alpha,center_vec = center_vec, sigma_vec = sigma_vec,
#'                             max_time = 60, default_ranges = FALSE, l_bound = lb, u_bound = ub)
#' PBOR_result$value
#' PBOR_result$theta
#' @export
MDPEREGsolve <- function(param1,
                         param2 = 1,
                         fixed_param1 = 10,
                         dist = "normal",
                         Y_array,
                         X_mat,
                         alpha,
                         center_vec,
                         sigma_vec,
                         max_time = 6000,
                         default_ranges = TRUE,
                         l_bound,
                         u_bound,
                         rho = 2,
                         phi = 10^(-6),
                         no_runs = 1000,
                         max_iter = 500000,
                         s_init = 1,
                         tol_fun = 10^(-6),
                         tol_fun_2 = 10^(-6))
{ M <- length(param1)
if(dist == "normal"){                                ## Normal
  if(default_ranges == TRUE)
  {lb_temp <- array(-Inf,M)
  lb <- c(lb_temp,0)
  ub <- array(Inf,M+1)
  }else{
    lb <- l_bound
    ub <- u_bound
  }
  func_here <- function(z)
  {beta <- z[1:M]
  sigma <- z[M+1]
  mu <- X_mat%*%beta
  part_1 <- 1/(((sqrt(2*pi)*sigma)^alpha)*(sqrt(1+alpha)))
  part_2 <- array(NA,M)
  for(jj in 1:M)
  {part_2[jj] <- dnorm(Y_array[jj],mu[jj],sigma)^alpha
  }
  return(part_1 - ((1+alpha)/alpha) * mean(part_2))
  }
  x0 <- c(param1,param2)
  x_soln <- PBORoptim(x0, func_here, lb = lb, ub = ub, center_vec, sigma_vec, max_time,
                      rho, phi, no_runs, max_iter,s_init, tol_fun, tol_fun_2)
  value_at_x_soln <- func_here(x_soln)
  return(list("theta" = x_soln, "value" = value_at_x_soln))
}else if(dist == "poisson"){                             ## Poisson
  if(default_ranges == TRUE)
  {lb <- array(-Inf, M)
  ub <- array(Inf, M)}
  else
  {lb <- l_bound
  ub <- u_bound}
  sum_upto <- fixed_param1
  func_here <- function(z)
  {beta <- z
  mu <- exp(X_mat%*%beta)
  part_1 <- array(NA,M)
  part_2 <- array(NA,M)
  for(jj in 1:M)
  {part_1 <- sum(dpois(0:sum_upto, mu[jj])^(1+alpha))
  part_2 <-  dpois(Y_array[jj], mu[jj])^alpha
  }
  return(mean(part_1) - ((1+alpha)/alpha) * mean(part_2))
  }
  x0 <- param1
  x_soln <- PBORoptim(x0, func_here, lb = lb, ub = ub, center_vec, sigma_vec, max_time,
                      rho, phi, no_runs, max_iter, s_init, tol_fun, tol_fun_2)
  value_at_x_soln <- func_here(x_soln)
  return(list("theta" = x_soln, "value" = value_at_x_soln))
}else if(dist == "logistic"){                             ## Logistic
  if(default_ranges == TRUE)
  {lb <- array(-Inf, M)
  ub <- array(Inf, M)}
  else
  {lb <- l_bound
  ub <- u_bound}
  func_here <- function(z)
  {beta <- z
  pp <- 1/(1+exp(-X_mat%*%beta))
  part_1 <- array(NA,M)
  part_2 <- array(NA,M)
  for(jj in 1:M)
  {part_1[jj] <- (pp[jj]^(1+alpha)+(1-pp[jj])^(1+alpha))
  part_2[jj] <-  ((pp[jj]^Y_array[jj])*((1-pp[jj])^(1-Y_array[jj])))^alpha
  }
  return(mean(part_1) - ((1+alpha)/alpha) * mean(part_2))
  }
  x0 <- param1
  x_soln <- PBORoptim(x0, func_here, lb = lb, ub = ub, center_vec, sigma_vec, max_time,
                      rho, phi, no_runs, max_iter, s_init, tol_fun, tol_fun_2)
  value_at_x_soln <- func_here(x_soln)
  return(list("theta" = x_soln, "value" = value_at_x_soln))
}
}




