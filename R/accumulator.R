# LCA functions----------------
##' LCA Process in an R function
##'
##' @param n number of trials
##' @param par par
##' @param I I
##' @param x0 x0
##' @param maxiter maxiter
##' @param nonlinear nonlinear
##' @param random random
##' @param debug debug
##' 
##' @examples 
##' p.vector <- c(kappa = 3, beta = 3, Z = .2, t0 = .3, dtovertau = 0.0001,
##' tau = 10)
##' args(rLCA_R)
##' res <- rLCA_R(n = 1, par = p.vector, I = c(1.2, 1), x0 = c(0, 0))
##' 
##' nacc <- length(I)
##' d1 <- data.table(x = rep(1:nrow(res[[1]]$trace), nacc), 
##'              y1 = as.vector %>% r(res[[1]]$trace), 
##'              y2 = as.vector(res[[1]]$lateral_inhibition),
##'              y3 = as.vector(res[[1]]$leakage),
##'              y4 = as.vector(res[[1]]$dx),
##'              z = factor(rep(1:nacc, each = nrow(res[[1]]$trace))))
##' 
##' 
##' p1 <- ggplot() +
##' geom_line(data = d1, aes(x=x, y=y1, colour =z)) +
##'   geom_hline(yintercept = Z[1]) + xlab("Time step") + ylab("Activation") +
##' theme_bw(base_size = 20) +
##' theme(aspect.ratio=1)
##' p1
##' 
##' @export
rLCA_R <- function(n, par, I, x0, 
                   maxiter=5001, nonlinear = FALSE, random = FALSE, 
                   debug = FALSE) {
  # p.vector <- c(Z = Z, t0 = t0, dtovertau = dtovertau, dt = dt, 
  #               kappa = pmat[j, 2] , beta = pmat[j, 1] )
  
  Z <- par[1]
  T0 <- par[2] 
  dtovertau <- par[3]
  dt <- par[4]
  K <- par[5]
  B <- par[6] 
  sdv <- sqrt(dtovertau)
  
  nacc <- length(I)
  out <- vector("list", n)
  
  for(j in seq_len(n)) {
    output <- numeric(2)
    undone  <- TRUE
    counter <- 0;
    
    dx <- numeric(nacc);   
    leak <- numeric(nacc);   
    inhibition <- numeric(nacc);    
    inhibition_received <- numeric(nacc);   
    f <- numeric(nacc)
    activation <- x0
    f <- x0
    
    ##  initialize accumulator array to record activation value
    ac_mat <- matrix(nrow = maxiter, ncol = nacc)
    dx_mat <- matrix(nrow = maxiter, ncol = nacc)
    in_mat <- matrix(nrow = maxiter, ncol = nacc)
    le_mat <- matrix(nrow = maxiter, ncol = nacc)
    ir_mat <- matrix(nrow = maxiter, ncol = nacc)
    
    ac_mat[counter + 1,] <- activation
    dx_mat[counter + 1, ] <- 0
    
    ##  initialize accumulator array to record activation value
    le_mat[counter + 1,] <- 0
    in_mat[counter + 1,] <- 0
    ir_mat[counter + 1,] <- 0
    
    repeat {
      counter <- counter + 1
      inhibition <- dtovertau * B*f;  
      leak       <- dtovertau * K*activation
      total_inhibition <- sum(inhibition)
      
      for (i in seq_len(nacc)) {
        inhibition_received[i] <- (total_inhibition - inhibition[i])
        dx[i] <- dtovertau * I[i] - leak[i] - inhibition_received[i] 
        activation[i] <- activation[i] + dx[i]
        
        if (random) activation[i] <- activation[i] + rnorm(1, 0, sdv)
        
        if (activation[i] >= Z) { 
          output[1] <- counter * dt - 0.5*dt + T0[i];
          output[2] <- i; 
          undone    <- FALSE; 
          if (debug) message("Activity ", activation[i])
        }
        if (nonlinear && activation[i] < 0) {
          # activation[i] <- 0;
          f[i] <- 0
          
          if (debug) message("Activation value (<0): ", activation[i])
        } else {
          f[i] <- activation[i]
        }
      }
      
      if (counter >= maxiter || !undone) break;
      
      ac_mat[counter + 1, ] <- activation
      dx_mat[counter + 1,]  <- dx
      in_mat[counter + 1, ] <- inhibition
      le_mat[counter + 1, ] <- leak
      ir_mat[counter + 1, ] <- inhibition_received
      
    }
    
    out[[j]] <- list(CRT = output, 
                     trace = ac_mat[1:counter,],
                     leakage = le_mat[1:counter,], 
                     lateral_inhibition = ir_mat[1:counter,],
                     lateral_inhibition_ego = in_mat[1:counter,],
                     dx = dx_mat[1:counter,])    
  }
  
  out
}


##' @export
rBv_R <- function(n, par, I, x0, 
                   maxiter=5001, nonlinear = FALSE, random = FALSE, 
                   debug = FALSE) {
  # n = 1
  # par = p.vector
  # I = c(1.2, .6)
  # x0 = c(0, 0)
  # maxiter = 8001
  # random=F
  # debug = T
  # nonlinear <- T
  
  Z <- par[1]
  T0 <- par[2] 
  dtovertau <- par[3]
  dt <- par[4]
  K <- par[5]
  B <- rep(par[6], 2)

  sdv <- sqrt(dtovertau)
  
  nacc <- length(I)
  out <- vector("list", n)
  
  initial_d <- 40
  initial_v <- 10
  pt <- .5
  tmax <- 8
  # dt <- 10/1000
  t <- seq(0, tmax + .1*dt, dt)
  nt <- length(t)
  t2a <- pedxing::get_tta(initial_d, initial_v, pt, t)
  ped_next <- foresee_next(-4.7, 1.6, .5, .5)
  collision_distance <- 1       
  # self-impose collision tolerance + remain time to reach 
  t2safe <- (collision_distance - ped_next) / 1.6
  veh_len <- 4.96
  veh_v <- 40
  t2a_copy <- t2a;
  veh_passing_time <- veh_len / veh_v;
  t2a_copy[t2a > t2safe] <- Inf  
  t2a_copy[t2a <= 0 & (t2a >= -veh_passing_time)] <- 0    
  t2a_copy[t2a < -veh_passing_time] <- Inf  
  # j <- 1
  for(j in seq_len(n)) {
    output <- numeric(2)
    undone  <- TRUE
    counter <- 0;
    
    dx <- numeric(nacc);   
    leak <- numeric(nacc);   
    inhibition <- numeric(nacc);    
    inhibition_received <- numeric(nacc);   
    f <- numeric(nacc)
    activation <- x0
    f <- x0
    
    ##  initialize accumulator array to record activation value
    ac_mat <- matrix(nrow = maxiter, ncol = nacc)
    dx_mat <- matrix(nrow = maxiter, ncol = nacc)
    in_mat <- matrix(nrow = maxiter, ncol = nacc)
    le_mat <- matrix(nrow = maxiter, ncol = nacc)
    ir_mat <- matrix(nrow = maxiter, ncol = nacc)
    
    ac_mat[counter + 1,] <- activation
    dx_mat[counter + 1, ] <- 0
    
    ##  initialize accumulator array to record activation value
    le_mat[counter + 1,] <- 0
    in_mat[counter + 1,] <- 0
    ir_mat[counter + 1,] <- 0
    # k <- 1
    t2aidx <- is.finite(t2a_copy)
    varyingB <- B[2] + 10 * B[2] / t2a_copy
    
    # ov <- get_value(p.vector, t, condition)
    # sv <- exp(ov)
    
    for(k in seq_len(nt)) {
      counter <- counter + 1
      
      if (t2aidx[k]) {
        B[2] <- varyingB[k]
        if(debug) cat("B", B, "\n")
      }


      inhibition <- dtovertau * B*f;  
      leak       <- dtovertau * K*activation
      total_inhibition <- sum(inhibition)
      
      
      for (i in seq_len(nacc)) {
        inhibition_received[i] <- (total_inhibition - inhibition[i])
        dx[i] <- dtovertau * I[i] - leak[i] - inhibition_received[i] 
        activation[i] <- activation[i] + dx[i]
        
        if (random) activation[i] <- activation[i] + rnorm(1, 0, sdv)
        
        if (activation[i] >= Z) { 
          output[1] <- counter * dt - 0.5*dt + T0[i];
          output[2] <- i; 
          undone    <- FALSE; 
          if (debug) message("Activity ", round(activation[i], 2))
        }
        if (nonlinear && activation[i] < 0) {
          # activation[i] <- 0;
          f[i] <- 0
          
          if (debug) message("Activation value (<0): ", round(activation[i], 2))
        } else {
          f[i] <- activation[i]
        }
      }
      
      ac_mat[counter + 1, ] <- activation
      dx_mat[counter + 1,]  <- dx
      in_mat[counter + 1, ] <- inhibition
      le_mat[counter + 1, ] <- leak
      ir_mat[counter + 1, ] <- inhibition_received
      
      if (counter >= maxiter || !undone) break;
    }
    
    # repeat {
    #   counter <- counter + 1
    #   B[2] <- B[2] + 0.3
    #   cat("B", B, "\n")
    #   inhibition <- dtovertau * B*f;  
    #   leak       <- dtovertau * K*activation
    #   total_inhibition <- sum(inhibition)
    #   
    #   for (i in seq_len(nacc)) {
    #     inhibition_received[i] <- (total_inhibition - inhibition[i])
    #     dx[i] <- dtovertau * I[i] - leak[i] - inhibition_received[i] 
    #     activation[i] <- activation[i] + dx[i]
    #     
    #     if (random) activation[i] <- activation[i] + rnorm(1, 0, sdv)
    #     
    #     if (activation[i] >= Z) { 
    #       output[1] <- counter * dt - 0.5*dt + T0[i];
    #       output[2] <- i; 
    #       undone    <- FALSE; 
    #       if (debug) message("Activity ", activation[i])
    #     }
    #     if (nonlinear && activation[i] < 0) {
    #       # activation[i] <- 0;
    #       f[i] <- 0
    #       
    #       if (debug) message("Activation value (<0): ", activation[i])
    #     } else {
    #       f[i] <- activation[i]
    #     }
    #   }
    # 
    #   ac_mat[counter + 1, ] <- activation
    #   dx_mat[counter + 1,]  <- dx
    #   in_mat[counter + 1, ] <- inhibition
    #   le_mat[counter + 1, ] <- leak
    #   ir_mat[counter + 1, ] <- inhibition_received
    #   
    #   if (counter >= maxiter || !undone) break;
    # 
    # }
    
    out[[j]] <- list(CRT = output, 
                     trace = ac_mat[1:(counter+1),],
                     leakage = le_mat[1:(counter+1),], 
                     lateral_inhibition = ir_mat[1:(counter+1),],
                     lateral_inhibition_ego = in_mat[1:(counter+1),],
                     dx = dx_mat[1:(counter+1),],
                     counter = counter)    
  }
  
  out
}

##' LBA Process in an R function
##' 
##' Reconstruct an LBA in LCA structure
##' 
##' @param n number of trials
##' @param par par
##' @param I I
##' @param x0 x0
##' @param maxiter maxiter
##' @param nonlinear nonlinear
##' @param random random
##' @param debug debug
##' 
##' @examples 
##' ## some examples
##' @export
rInput_R <- function(n, par, I, x0, maxiter=5001, nonlinear = FALSE, 
                   random = FALSE, debug = FALSE) {
  
  Z <- par[1]
  t0 <- par[2]
  dtovertau <- par[3]
  dt <- par[4]
  sdv <- sqrt(dtovertau) 
  
  ## Prepare the process 
  nacc <- length(I)

  counter <- numeric(n);
  out_mat <- matrix(nrow = 2, ncol = n) 
  act_cube <- array(dim = c(nacc, maxiter, n))
  dx_cube  <- array(dim = c(nacc, maxiter, n))
  
  for(j in seq_len(n)) {
    output    <- numeric(2)
    undone    <- TRUE
    counterj  <- 0
    ##  initialize accumulator array to record activation value
    dx <- numeric(nacc);   
    activation <- x0
    act_mat <-  matrix(nrow = nacc, ncol = maxiter)
    dx_mat  <- matrix(nrow = nacc, ncol = maxiter)
    act_mat[,counterj+1] <- activation
    dx_mat[,counterj+1] <- 0
    repeat {
      counterj <- counterj + 1
      for (i in seq_len(nacc)) {
        
        dx[i] <- dtovertau * I[i]
        activation[i] <- activation[i] + dx[i]
        
        if(random) activation[i] <- activation[i] + rnorm(1, 0, sdv)
        
        if (activation[i] >= Z) { 
          output[1] <- counterj * dt - 0.5*dt + t0;
          output[2] <- i; 
          undone    <- FALSE; 
          if (debug) message("Activation value: ", activation[i])
        }
        if (nonlinear && activation[i] < 0) activation[i] <- 0;
      }
      
      act_mat[,counterj+1] <- activation
      dx_mat[,counterj+1]   <- dx
      if (counterj > maxiter || !undone) break;
    }
    
    # idx <- counterj + 1
    act_cube[,,j] <- act_mat
    dx_cube[,,j] <- dx_mat
    out_mat[,j] <- output
    counter[j] <- counterj
    
  }
  
  new("Input", choice_RT = out_mat, activation = act_cube,
      dx = dx_cube, counter = counter)
} 

