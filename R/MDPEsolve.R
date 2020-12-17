#' Minimum Divergence Power Estimation (MDPE) Solver
#'
#' `MDPEsolve` function estimates the minimum divergence power estimator for list of
#' distributions. This function uses a Blackbox optimization technique, which is a
#' variant of Pattern search, known as Prior Based Optimization Routine (PBOR).
#'
#'
#' @param param1 Initial value of first parameter : Normal (mean), Exponential (mean),
#' Beta (first parameter), Gamma (shape), lognormal (mean of normal), Poisson (mean),
#' Geometric (p), Binomial (p), Negative binomial (p).
#' @param param2 Initial value of second parameter : Normal (sigma), Beta (second parameter),
#' Gamma (scale), lognormal (sigma of normal). Default value is 1.
#' @param fixed_param1 First fixed parameter values : Binomial (n), Poisson (approx. sum upto),
#' Geometric (approx. sum upto), Negative binomial (approx. sum upto). Default value is 10.
#' @param fixed_param2 Second fixed parameter values : Negative binomial (r). Default value
#'  is 3.
#' @param dist Distribution of sample : "normal", "exp", "beta", "gamma", "lognormal",
#' "poisson", "geometric", "binomial", "negbinomial". Default value is "normal".
#' @param X_array Array of values in sample.
#' @param alpha Parameter of Power estimator, value must be between 0 and 1.
#' @param center_vec Mean vector of prior normal distribution. It can be set at the
#' neighborhood where the true parameter value is expected.
#' @param sigma_vec The standard deviation vector of the prior normal distribution.
#' @param max_time Maximum time (in seconds) alloted for optimization. Default value is 3600.
#' @param default_ranges 'TRUE' for keeping the default ranges of the distribution
#' paramaters, e.g., in case of normal distribution, the default ranges of mu and sigma
#' are (-Inf, Inf) and (0, Inf) respectively. In case solution needs to be found within a
#' specific interval, it should be set 'FALSE' and 'l_bound', 'u_bound' should be set
#' accordingly. Default value is TRUE.
#' @param l_bound Lower bound for parameter vector, to activate this bound, user must set
#' 'default_ranges' to be FALSE. Default value is c(0,0).
#' @param u_bound Upper bound for parameter vector, to activate this bound, user must set
#' 'default_ranges' to be FALSE. Default value is c(1,1).
#' @param rho Step Decay Rate with default value 2
#' @param phi Step size threshold, i.e., lower Bound Of Global Step Size. Default value is \eqn{10^{-20}}
#' @param no_runs Maximum number of runs to be executed, default value is 1000.
#' @param max_iter Max Number Of Iterations in each Run. Default Value is 10000.
#' @param s_init Initial Global Step Size. Default Value is 1.
#' @param tol_fun Termination Tolerance on the function value. Default Value is \eqn{10^{-20}}
#' @param tol_fun_2 Termination Tolerance on the difference of solutions in two consecutive
#' runs. Default Value is \eqn{10^{-20}}.
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
#' X_array <- rnorm(10)
#' center_vec <- c(0,0)
#' sigma_vec <- c(10^4, 10^4)
#' alpha <- 0.5
#' MDPEsolve(param1 = 0.1, param2 = 5, dist = "normal", default_ranges = TRUE, X_array = X_array, alpha = 0.5,
#'           center_vec = center_vec, sigma_vec = sigma_vec)
#'
#' X_array <- rexp(5, rate = 1/5)
#' alpha <- 0.7
#' MDPEsolve(0.1, 0, dist = "exp", default_ranges = TRUE, X_array = X_array, alpha = 0.5,
#'           center_vec = center_vec, sigma_vec = sigma_vec)
#'
#' X_array <- rbeta(150, 2, 3)
#' alpha <- 0.5
#' MDPEsolve(2, 2, dist = "beta", default_ranges = FALSE,l_bound = c(1/(1+alpha),1/(1+alpha)),
#'           u_bound = c(10,10),X_array = X_array, alpha = 0.5,center_vec = center_vec,
#'           sigma_vec = sigma_vec)
#'
#' X_array <- rgamma(5, shape = 2, scale = 3) ## k (shape) > alpha/(1+alpha)
#' alpha <- 0.75
#' MDPEsolve(2, 2, dist = "gamma", X_array = X_array, alpha = 0.5,center_vec = center_vec,
#'           sigma_vec = sigma_vec)
#'
#' X_array <- rlnorm(5, meanlog = 5, sdlog = 1) ## MUST BE RUN within a domain
#' MDPEsolve(2, 2.5, dist = "lognormal", l_bound = c(1,2),u_bound = c(10,10), default_ranges = FALSE,
#'           X_array = X_array, alpha = 0.5,center_vec = center_vec, sigma_vec = sigma_vec)
#'
#' X_array <- rpois(50, 3)
#' alpha <- 0.25
#' MDPEsolve(param1 =2, param2 =2, fixed_param1 = 1000, dist = "poisson", X_array = X_array, alpha = 0.5,
#'           center_vec = center_vec, sigma_vec = sigma_vec)
#'
#' X_array <- rgeom(50, 0.6)
#' alpha <- 0.5
#' MDPEsolve(param1 = 0.3, fixed_param1 = 1000, dist = "geometric", X_array = X_array,
#'           alpha = 0.5,center_vec = center_vec, sigma_vec = sigma_vec)
#'
#' X_array <- rbinom(100,8,0.6)
#' alpha <- 0.15
#' MDPEsolve(param1 = 0.3, fixed_param1 = 8, dist = "binomial", X_array = X_array, alpha = 0.5,
#'           center_vec = center_vec, sigma_vec = sigma_vec)
#'
#'
#' X_array <- rnbinom(100,8,0.6) ## k
#' alpha <- 0.65
#' MDPEsolve(param1 = 0.3, fixed_param1 = 1000, fixed_param2 = 30, dist = "negbinomial",
#'           X_array = X_array, alpha = 0.5,center_vec = center_vec, sigma_vec = sigma_vec)
#'
#' @export

MDPEsolve <- function(param1,
                      param2 = 1,
                      fixed_param1 = 10,
                      fixed_param2 = 3,
                      dist = "normal",
                      X_array,
                      alpha,
                      center_vec,
                      sigma_vec,
                      max_time = 3600,
                      default_ranges = TRUE,
                      l_bound = c(0,0),
                      u_bound = c(1,1),
                      rho = 2,
                      phi = 10^(-20),
                      no_runs = 1000,
                      max_iter = 500000,
                      s_init = 1,
                      tol_fun = 10^(-20),
                      tol_fun_2 = 10^(-20))
{  if(dist == "normal"){                                ## Normal
  if(default_ranges == TRUE)
  {lb <- c(-Inf,0)
  ub <- c(Inf,Inf)}
  else
  {lb <- l_bound[1:2]
  ub <- u_bound[1:2]}
  func_here <- function(z)
  {mu <- z[1]
  sigma <- z[2]
  part_1 <- 1/(((sqrt(2*pi)*sigma)^alpha)*(sqrt(1+alpha)))
  part_2 <- ((1+alpha)/alpha) * mean(dnorm(X_array,mu,sigma)^alpha)
  return(part_1 - part_2)
  }
  x0 <- c(param1,param2)
  x_soln <- PBORoptim(x0, func_here, lb = lb, ub = ub, center_vec, sigma_vec, max_time,
                      rho, phi, no_runs, max_iter,s_init = 1, tol_fun, tol_fun_2)
  value_at_x_soln <- func_here(x_soln)
  return(list("theta" = x_soln, "value" = value_at_x_soln))
}else if(dist == "exp"){                             ## Exponential
  if(default_ranges == TRUE)
  {lb <- 0
  ub <- Inf}
  else
  {lb <- l_bound[1]
  ub <- u_bound[1]}
  func_here <- function(z)
  {lambda <- z
  part_1 <- 1/((lambda^alpha)*(1+alpha))
  part_2 <- ((1+alpha)/alpha) * mean(dexp(X_array,rate = 1/lambda)^alpha)
  return(part_1 - part_2)
  }
  x0 <- param1
  # func_here(x0)
  x_soln <- PBORoptim(x0, func_here, lb = lb, ub = ub, center_vec, sigma_vec, max_time,
                      rho, phi, no_runs, max_iter, s_init, tol_fun, tol_fun_2)
  value_at_x_soln <- func_here(x_soln)
  return(list("theta" = x_soln, "value" = value_at_x_soln))
}else if(dist == "beta"){                             ## Beta
  if(default_ranges == TRUE)
  {lb <- c((1+10^(-10))/(alpha+1),(1+10^(-10))/(alpha+1))
  ub <- c(Inf,Inf)}
  else
  {lb <- l_bound[1:2]
  ub <- u_bound[1:2]}
  func_here <- function(z)
  {a <- z[1]
  b <- z[2]
  # integrand <- function(x){return(dbeta(x,a,b)^(1+alpha))}
  # set.seed(1)
  # part_1 <- integrate(integrand, 0, 1)
  part_1 <- beta(a,b)/beta(alpha*(a-1)+a, alpha*(b-1)+b)
  part_2 <- ((1+alpha)/alpha) * mean(dbeta(X_array, a, b)^alpha)
  return(part_1 - part_2)
  }
  x0 <- c(param1,param2)
  # func_here(x0)
  x_soln <- PBORoptim(x0, func_here, lb = lb, ub = ub,center_vec, sigma_vec, max_time,
                      rho, phi, no_runs, max_iter, s_init, tol_fun, tol_fun_2)
  value_at_x_soln <- func_here(x_soln)
  return(list("theta" = x_soln, "value" = value_at_x_soln))
}else if(dist == "gamma"){                           ## Gamma
  if(default_ranges == TRUE)
  {lb <- c(alpha/(1+alpha),0)
  ub <- c(Inf,Inf)}
  else
  {lb <- l_bound[1:2]
  ub <- u_bound[1:2]}
  func_here <- function(z)
  {a <- z[1]     # shape    a must be > alpha/(1+alpha)
  b <- z[2]     # scale
  part_1 <- gamma(a*(1+alpha)-alpha)/((gamma(a)^(1+alpha))*(b^(a*alpha))*((1+alpha)^a))
  part_2 <- ((1+alpha)/alpha) * mean(dgamma(X_array, shape = a, scale = b)^alpha)
  return(part_1 - part_2)
  }
  x0 <- c(param1,param2)
  # func_here(x0)
  x_soln <- PBORoptim(x0, func_here, lb = lb, ub = ub,center_vec, sigma_vec, max_time,
                      rho, phi, no_runs, max_iter, s_init, tol_fun, tol_fun_2)
  value_at_x_soln <- func_here(x_soln)
  return(list("theta" = x_soln, "value" = value_at_x_soln))
}else if(dist == "lognormal"){                           ## Lognormal
  if(default_ranges == TRUE)
  {lb <- c(-Inf,0)
  ub <- c(Inf,Inf)}
  else
  {lb <- l_bound[1:2]
  ub <- u_bound[1:2]}
  func_here <- function(z)
  {mulog <- z[1]
  sigmalog <- z[2]
  integrand <- function(x){return(dlnorm(x, mulog, sigmalog)^(1+alpha))}
  set.seed(1)
  part_1 <- integrate(integrand, 0, Inf)
  part_2 <- ((1+alpha)/alpha) * mean(dlnorm(X_array, mulog, sigmalog)^alpha)
  return(part_1$value - part_2)
  }
  x0 <- c(param1,param2)
  # func_here(x0)
  x_soln <- PBORoptim(x0, func_here, lb = lb, ub = ub, center_vec, sigma_vec, max_time,
                      rho, phi, no_runs, max_iter, s_init, tol_fun, tol_fun_2)
  value_at_x_soln <- func_here(x_soln)
  return(list("theta" = x_soln, "value" = value_at_x_soln))
}else if(dist == "poisson"){                             ## Poisson
  if(default_ranges == TRUE)
  {lb <- 0
  ub <- Inf}
  else
  {lb <- l_bound[1]
  ub <- u_bound[1]}
  sum_upto <- fixed_param1
  func_here <- function(z)
  {part_1 <- sum(dpois(0:sum_upto, z)^(1+alpha))
  part_2 <- ((1+alpha)/alpha) * mean(dpois(X_array, z)^alpha)
  return(part_1 - part_2)
  }
  x0 <- param1
  x_soln <- PBORoptim(x0, func_here, lb = lb, ub = ub, center_vec, sigma_vec, max_time,
                      rho, phi, no_runs, max_iter, s_init, tol_fun, tol_fun_2)
  value_at_x_soln <- func_here(x_soln)
  return(list("theta" = x_soln, "value" = value_at_x_soln))
}else if(dist == "geometric"){                             ## Geometric
  if(default_ranges == TRUE)
  {lb <- 0
  ub <- 1}
  else
  {lb <- l_bound[1]
  ub <- u_bound[1]}
  sum_upto <- fixed_param1
  func_here <- function(z)
  {part_1 <- sum(dgeom(0:sum_upto, z)^(1+alpha))
  part_2 <- ((1+alpha)/alpha) * mean(dgeom(X_array, z)^alpha)
  return(part_1 - part_2)
  }
  x0 <- param1
  x_soln <- PBORoptim(x0, func_here, lb = lb, ub = ub, center_vec, sigma_vec, max_time,
                      rho, phi, no_runs, max_iter, s_init, tol_fun, tol_fun_2)
  value_at_x_soln <- func_here(x_soln)
  return(list("theta" = x_soln, "value" = value_at_x_soln))
}else if(dist == "binomial"){                             ## Binomial
  if(default_ranges == TRUE)
  {lb <- 0
  ub <- 1}
  else
  {lb <- l_bound[1]
  ub <- u_bound[1]}
  num <- fixed_param1
  func_here <- function(z)
  {part_1 <- sum(dbinom(0:num, num, z)^(1+alpha))
  part_2 <- ((1+alpha)/alpha) * mean(dbinom(X_array, num, z)^alpha, max_time)
  return(part_1 - part_2)
  }
  x0 <- param1
  x_soln <- PBORoptim(x0, func_here, lb = lb, ub = ub, center_vec, sigma_vec, max_time,
                      rho, phi, no_runs, max_iter, s_init, tol_fun, tol_fun_2)
  value_at_x_soln <- func_here(x_soln)
  return(list("theta" = x_soln, "value" = value_at_x_soln))
}else if(dist == "negbinomial"){                             ## Negative Binomial
  if(default_ranges == TRUE)
  {lb <- 0
  ub <- 1}
  else
  {lb <- l_bound[1]
  ub <- u_bound[1]}
  num <- fixed_param1
  r <- fixed_param2
  func_here <- function(z)
  {part_1 <- sum(dnbinom(0:num, r, z)^(1+alpha))
  part_2 <- ((1+alpha)/alpha) * mean(dnbinom(X_array, r, z)^alpha)
  return(part_1 - part_2)
  }
  x0 <- param1
  x_soln <- PBORoptim(x0, func_here, lb = lb, ub = ub, center_vec, sigma_vec, max_time,
                      rho, phi, no_runs, max_iter, s_init, tol_fun, tol_fun_2)
  value_at_x_soln <- func_here(x_soln)
  return(list("theta" = x_soln, "value" = value_at_x_soln))
}
}
