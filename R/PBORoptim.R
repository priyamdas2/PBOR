#' Prior Based Optimization Routine
#'
#' `PBOR` function minimizes Black-box function of parameters where each parameter belongs
#' to any one of this kind of intervals : -Inf to Inf, -Inf to F, F to Inf, F to G; where
#' F and G are finite known numbers.
#'
#' @param x0 Initial guess (vector) of solution.
#' @param func Function needs to be minimized.
#' @param lb Lower bound vector of parameter domain.
#' @param ub Upper bound vector of domain.
#' @param center_vec Mean vector of prior normal distribution. It can be set at the
#' neighborhood where the true parameter value is expected. Default value is 0 vector.
#' @param sigma_vec The standard deviation vector of the prior normal distribution.
#' @param max_time Maximum time alloted for optimization.
#' @param rho Step Decay Rate with default value 2
#' @param phi Step size threshold, i.e., lower Bound Of Global Step Size. Default value is \eqn{10^{-20}}
#' @param no_runs Maximum number of runs to be executed, default value is 1000.
#' @param max_iter Max Number Of Iterations in each Run. Default Value is 10000.
#' @param s_init Initial Global Step Size. Default Value is 1.
#' @param tol_fun Termination Tolerance on the function value. Default Value is \eqn{10^{-20}}
#' @param tol_fun_2 Termination Tolerance on the difference of solutions in two consecutive
#' runs. Default Value is \eqn{10^{-20}}.
#'
#' @return The point where func is minimized.
#'
#'@references Das, P. and Ghosh, A. (2021). \emph{PBOR paper}. Future Journal, 2, 63-74.
#'@author Priyam Das (maintainer) (\email{priyam_das@@hms.harvard.edu}), Abhik Ghosh.
#' @examples
#' func <- function(x){return(100*(x[2]-x[1]^2)^2 + (x[1]-1)^2)}
#' PBORoptim(x0 = c(0,0), func, lb = c(-100,-Inf), ub = c(10, Inf), center_vec = c(0,0),
#' sigma_vec = c(10^2, 10^2))
#' @export
PBORoptim <- function(x0,
                      func,
                      lb,
                      ub,
                      center_vec,
                      sigma_vec,
                      max_time = 6000,
                      rho = 2,
                      phi = 10^(-20),
                      no_runs = 10,
                      max_iter = 10000,
                      s_init = 1,
                      tol_fun = 10^(-20),
                      tol_fun_2 = 10^(-20))
{M <- length(ub)
diff <- ub - lb
diff_1 <- x0 - lb
diff_2 <- ub - x0
if(min(diff) <= 0)
{cat("ub must be strictly greater than lb")
   exit()}
if(min(diff_1) < 0)
{cat("x0 must be greater than equal to lb")
   exit()}
if(min(diff_2) < 0)
{cat("ub must be greater than equal to x0")
   exit()}

# Finding the type of interval

type_of_intervals <- interval_types(lb,ub)

# Converting starting point to unit

x0_unit <- Domain_to_unit(x0,type_of_intervals,lb,ub,center_vec,sigma_vec)

run_values <- array(NA,no_runs)

t0 <- Sys.time()

for(run in 1:no_runs)
{s <- s_init
if(run == 1)
{initial_value <- func(x0)
 X <- Domain_to_unit(x0,type_of_intervals,lb,ub,center_vec,sigma_vec)
}
X_domain <- Unit_to_domain(X,type_of_intervals,lb,ub,center_vec,sigma_vec)
run_values[run] <- sqrt(sum(X^2))
for(i in 1:max_iter)
{#print(c(run,i))
 time_elapsed <- Sys.time() - t0
 #cat("time_elapsed is", time_elapsed,"\n")
 if(time_elapsed > max_time)
 {break}
   X_domain <- Unit_to_domain(X,type_of_intervals,lb,ub,center_vec,sigma_vec)
   value_now <- func(X_domain)
   values <- array(NA,2*M)
   s_nows <- array(NA,2*M)
   X_temps <- matrix(NA,2*M,M)
   for(jj in 1:(2*M))
   {X_temp <- X
   s_now <- ((-1)^jj)*s
   j <- ceiling(jj/2)
   x_cord <- X_temp[j]
   temp_value <- x_cord + s_now     #
   if(temp_value > 1)
   {f <- ceiling(log(s_now/(1-x_cord),rho))
   s_now <- s_now/(rho^f)
   temp_value <- x_cord + s_now               #
   }
   if(temp_value < 0)
   {f <- ceiling(log(-s_now/x_cord, rho))
   s_now <- s_now/(rho^f)
   temp_value <- x_cord + s_now               #
   }
   X_temp[j] <- temp_value
   s_nows[jj] <- s_now
   X_temp_domain <- Unit_to_domain(X_temp,type_of_intervals,lb,ub,center_vec,sigma_vec)
   values[jj] <- func(X_temp_domain)
   if(is.null(values[jj]) == TRUE || is.finite(values[jj]) == FALSE)
   {values[jj] <- 999999}
   X_temps[jj,] <- X_temp
   }

   # Finding minimum of all values

   min_of_2M <- min(values)
   min_indexes <- which(values == min_of_2M)
   length_min_indexes <- length(min_indexes)
   if(length_min_indexes > 1)
   {#cord <- rdunif(1, 1, length_min_indexes)
      cord <- 1
   min_index <- min_indexes[cord]
   }
   else
   {min_index <- min_indexes}

   if(values[min_index] < value_now)
   {X <- X_temps[min_index,]
   }
   value_now_2 <- min(min_of_2M,value_now)
   if(abs(value_now_2 - value_now) < tol_fun && s > phi)
   {s <- s/rho
   }

   if(s <= phi)
   {break
   }
}                                           # End of iterations
if(run > 1 && s <= phi)
{new_norm <- sqrt(sum(X^2))
if(abs(run_values[run] - new_norm) < tol_fun_2)
{break
}
}
}
X_domain <- Unit_to_domain(X,type_of_intervals,lb,ub,center_vec,sigma_vec)
return(X_domain)
}
