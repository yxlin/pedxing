# PACKAGE ---------------------------------------------------
##' Multi-phase decision making 
##' 
##' Version history:
##' \itemize{
##' }
##'
##' @keywords package
##'
##' @name pedxing
##' @docType package
##' @author  Yi-Shin Lin <yishinlin@pm.me>
##' @importFrom Rcpp evalCpp
##' @useDynLib pedxing
NULL

# UM functions----------------
##' @export
ou <- function(xdomain, tdomain, par, condition, add_dw = FALSE) {
  ## time point 0
  nx <- length(xdomain)
  nt <- length(tdomain)
  
  h <- tdomain[2] - tdomain[1]

  s0 <- par[6]
  gamma <- par[7]

  ov <- get_value(par, tdomain, condition, verbose = FALSE)
  sv <- exp(ov)
  
  mu <- matrix(NA, nrow = nt, ncol = nx)
  
  for(i in seq_len(nt)) {
    mu[i, ] <- sv[i] - gamma*xdomain
  }

  idx0 <- .5*(nx+1)
  x0 <- xdomain[idx0]
  dxt0 <- mu[1, idx0] * h  ## 0 to 1
  
  if (add_dw) { 
    dw <- sqrt(h) * s0 * rnorm(nt) 
    dxt0 <- dxt0 + dw[1]
  }

  xtmp <- NULL  ## xtmp is the vector of dXt
  
  
  for(i in 2:nt) {
    
    if (i == 2) {
      xtmp <- c(x0, dxt0)
    } else {
      xtmp <- c(xtmp, dxti) 
    }
    
    xt <- sum(xtmp)  ## accumulated evidence
    ## The difference between A(xt, t) (defined by the evidence space) and 
    ## the moment-to-moment A(xt, t) at the time point i
    tmp <- mu[i, ] - (sv[i] - gamma* xt) 
    
    ## Because the evidence space is also constrained by dx, tmp only closes to 0
    ## - Find the real A(xt, t) on the x domain defined by the user
    idx <- which.min(abs(tmp))
    
    ## calculated based on the defined x domain
    dxti <- ifelse(add_dw, mu[i, idx]*h + dw[i],  mu[i, idx]*h)
  }
  
  return(xtmp)
  
}

##' @export
get_value <- function(par, ts, condition, verbose = TRUE)
{
  # tmp <- pedxing::get_value(tmp_pvec, vehiclets, cond)
  # par <- tmp_pvec
  # ts <- vehiclets
  # condition <- cond
   
  # 1, 2, 3, 9, 10, 11, 12, 13, 14
  # par <- p.vector
  # ts <- t
  # condition
  # verbose = TRUE
  ped_t0 <- par[5]
  
  ped_v <- par[9]
  ped_x0 <- par[8]
  ped_pt <- par[10]
  ped_c <- par[11]
  veh_len <- par[12]
  
  # kdv = 0.5*(kg/walking_speed) + ke
  # ke = 0; originally named Ct; -Ct controls the square of the yaw rate 
  # rotation effectuated for an action. always assume it 0, rendering it unused.
  # Assume the pedestrian walked in a constant speed, so no kda.
  #
  # d^{dot}_{g} is identical to -ped_v when no acceleration is enacted.kg <- p.vector[1]
  kg <- par[1]
  ke <- par[2]
  kc <- par[3]
  kdv <- .5*(kg/ped_v) + ke
  pt <- par[10]

  ## ped_x0, ped_v, ped_pt, ped_t0
  ped_next <- foresee_next(ped_x0, ped_v, ped_pt, ped_t0)
  t2safe <- (ped_c - ped_next) / ped_v
  # par[1]: kg,               par[2]: ke (render unused)
  # par[3]: kc,               par[4]: b0, 
  # par[5]: t0,               par[6]: s0,  
  # par[7]: sigmoid scaling   par[8]: accu_time_scale,
  # par[9]: initial_pedestrian_position  par[10]: walking_speed 
  # par[11]: prediction_time             par[12]: action_duration  
  # par[13]: collision distance          par[14]: vehicle_length 
  veh_v <- condition[1] / condition[2]; # distance and time at t0
  # t2a0 <- get_tta(condition[1], veh_v, ts)
  t2a <- get_tta(condition[1], veh_v, ts, pt)
  # plot(ts, t2a0, ylim = range(t2a0, t2a))
  # points(0, condition[2], cex = 2)
  # lines(ts, t2a)
  
  # t2a <- (condition[1] - veh_speed*(times + par[11])) / veh_speed;
  t2a_copy <- t2a;
  veh_passing_time <- veh_len / veh_v;
  
  # the moment when the vehicle is still far away from the crossing
  # No penalty on the collision utility. kc is unused. The other two
  # parameters speed jointly decide total utility
  t2a_copy[t2a > t2safe] <- Inf  
  if(verbose) cat("Start to suppress at", sum(is.infinite(t2a_copy)), "\n")
  
  # - <0 mean already passed
  # - a small moment when the vehicle's body occupies the crossing zone,
  # t2a is also less than 0, if they are greater than -v_pass_time, it
  # means the car is in the crossing zone. This results in -Inf the 
  # collision utility. Car crash. No go, total # utility = -Inf
  t2a_copy[t2a <= 0 & (t2a >= -veh_passing_time)] <- 0    
  
  # plot(t2a_copy)
  # the moment when the vehicle has moved entirely away from the crossing zone
  # No worries. The vehicle has left the crossing
  t2a_copy[t2a < -veh_passing_time] <- Inf  

  out <- kg*ped_v - kdv*ped_v^2 - kc/t2a_copy
  return(out)
}

##' Foresee a next state
##' 
##' Get a predictive next position for a pedestrian (R function is used 
##' in testing. Another copy of C++ function runs the model fitting.)
##' 
##' Function to get the predictive pedestrian state if beginning to cross. 
##' I assume a pedestrian prepares to walk to the north on the negative
##' side of y axis.
##' 
##' @param x0 an initial (pedestrian) position 
##' @param v pedestrian walking speed  
##' @param pt how far in time a pedestrian can foresee what consequence 
##' its action may result
##' @param t0 how much time is needed to enact an action; assuming
##' instantaneous when t0 = 0; this should map to t0.
##' @import data.table
##' @examples 
##' ## some examples
##' @export
foresee_next <- function(x0, v, pt, t0 = 0)
{
  if (pt < .5*t0) stop("One cannot foresee a duration less than half of its action time")
  x0 + v*(pt - .5*t0)
}

##' @export
get_tta <- function(x0, v, ts, pt = 0)
{
  # get_tta (used for testing only)
  # 
  # x0 = front of the car; initial (vehicle) position
  # v = vehicle speed; NOT pedestrian's speed
  # pt = pedestrian foresees pt seconds into the future 
  # ts = a time vector shared by the two agents
  (x0 - v * (ts + pt)) / v
}

##' @export
get_crossing_onset <- function(ts, evidence, decision_threshold)
{
  above_thresh_samples <- which(evidence >= decision_threshold)
  out <- rep(NA, 2)
  if (length(above_thresh_samples) > 0)
  {
    out[1] <- above_thresh_samples[1]
    out[2] <- ts[out[1]]  
  } 
  names(out) <- c("index", "value")
  out
}
