

#################################################################
exit <- function() {
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}
#################################################################

F_normal <- function(x, center, sigma)
{return(pnorm(x, mean = center, sd = sigma))}

#################################################################

F_normal_inv <- function(F_x, center, sigma)
{return(qnorm(F_x, mean = center, sd = sigma))}

#################################################################

F_trun_normal <- function(x, center, sigma, cut)
{area_at_cut <- pnorm(cut,center,sigma)
area_at_value <- pnorm(x,center,sigma)
if(x >= cut)
{truncated_area <- 1 - area_at_cut
value_to_return <- (area_at_value - area_at_cut)/truncated_area}
if(x < cut)
{truncated_area <- area_at_cut
value_to_return <- (area_at_value)/truncated_area}
return(value_to_return)}

#################################################################

F_trun_normal_inv <- function(F_x, center, sigma, cut, direction)
{if(direction == 1)
{quantile_of_cut <- pnorm(cut,center,sigma)
total_area_right_of_cut <- 1 - pnorm(cut,center,sigma)
actual_area <- F_x*total_area_right_of_cut
area_from_minus_inf <- quantile_of_cut + actual_area
value_to_return <- qnorm(area_from_minus_inf,center, sigma)}
  if(direction == -1)
  {total_area_left_of_cut <- pnorm(cut,center,sigma)
  actual_area <- F_x*total_area_left_of_cut
  area_from_minus_inf <- actual_area
  value_to_return <- qnorm(area_from_minus_inf,center, sigma)}
  return(value_to_return)}


#################################################################

F_trun_normal_both_side <- function(x, center, sigma, cut_1, cut_2)
{area_at_cut_1 <- pnorm(cut_1,center,sigma)
area_at_cut_2 <- pnorm(cut_2,center,sigma)
area_at_value <- pnorm(x,center,sigma)
truncated_area <- area_at_cut_2 - area_at_cut_1
value_to_return <- (area_at_value - area_at_cut_1)/truncated_area
return(value_to_return)}

F_trun_normal_both_side(0, 0, 1, -1.96, 1.96)

##################################################################

F_trun_normal_both_side_inv <- function(F_x, center, sigma, cut_1, cut_2)
{quantile_of_left_cut <- pnorm(cut_1,center,sigma)
total_area_left_cut_to_right_cut <- pnorm(cut_2,center,sigma) - pnorm(cut_1,center,sigma)
actual_area <- F_x*total_area_left_cut_to_right_cut
area_from_minus_inf <- quantile_of_left_cut + actual_area
value_to_return <- qnorm(area_from_minus_inf,center, sigma)
return(value_to_return)}

########### Type of Intervals #####################################

interval_types <- function(lb,ub)
{M <- length(lb)
types <- array(NA, M)
for(i in 1:M)
{if(lb[i] == -Inf && ub[i] == Inf){
  types[i] <- 1
}else if(is.finite(lb[i]) == TRUE && ub[i] == Inf){
  types[i] <- 2
}else if(lb[i] == -Inf && is.finite(ub[i]) == TRUE){
  types[i] <- 3
}else if(is.finite(lb[i]) == TRUE && is.finite(ub[i]) == TRUE){
  types[i] <- 4
}
}
return(types)
}


######### Domain to unit + Unit to domain function ################

Domain_to_unit <- function(x,type_of_intervals,lb,ub,center_vec,sigma_vec)
{M <- length(x)
y <- array(NA,M)
for(i in 1:M)
{if(type_of_intervals[i] == 1){
  y[i] <- F_normal(x[i], center_vec[i], sigma_vec[i])
}else if(type_of_intervals[i] == 2){
  y[i] <- F_trun_normal(x[i], center_vec[i], sigma_vec[i], lb[i])
}else if(type_of_intervals[i] == 3){
  y[i] <- F_trun_normal(x[i], center_vec[i], sigma_vec[i], ub[i])
}else if(type_of_intervals[i] == 4){
  y[i] <- F_trun_normal_both_side(x[i], center_vec[i], sigma_vec[i], lb[i], ub[i])
}
}
return(y)
}



Unit_to_domain <- function(F_x,type_of_intervals,lb,ub,center_vec,sigma_vec)
{M <- length(F_x)
x <- array(NA,M)
for(i in 1:M)
{if(type_of_intervals[i] == 1){
  x[i] <- F_normal_inv(F_x[i], center_vec[i], sigma_vec[i])
}else if(type_of_intervals[i] == 2){
  x[i] <- F_trun_normal_inv(F_x[i], center_vec[i], sigma_vec[i], lb[i], 1)
}else if(type_of_intervals[i] == 2){
  x[i] <- F_trun_normal_inv(F_x[i], center_vec[i], sigma_vec[i], ub[i], -1)
}else if(type_of_intervals[i] == 4){
  x[i] <- F_trun_normal_both_side_inv(F_x[i], center_vec[i], sigma_vec[i], lb[i], ub[i])
}
}
return(x)
}
