library(LaplacesDemon)

distance_fn <- function(p, q, type = 1){
  if (type == 1){
    return( sum(abs(p - q)) )
  } else if (type == 2){
    return( sqrt(sum((p - q)^2)) )
  }
}

log_density_ratio_fn <- function(mu_top, mu_bottom, x, sigma, eps, type = "gauss"){
  if (type == "gauss"){
    return( ( - sum((x - mu_top)^2) + sum((x - mu_bottom)^2) ) >= (2 * sigma^2 * eps) )
  } else if (type == "laplace"){
    return( ( - sum(abs(x - mu_top)) + sum(abs(x - mu_bottom)) ) >= (sigma * eps) )  
  }
}

integral_fn <- function(samp_pts, eps, mu_1, mu_2, sigma, type, int, n){
  
  suppressWarnings(if(is.na(samp_pts)){
    if (int == 1){
      if (type == "gauss"){
        samp_pts <- sapply(mu_1, function(x) stats::rnorm(n, x, sigma))
      } else if (type == "laplace"){
        samp_pts <- sapply(mu_1, function(x) rlaplace(n, location = x, scale = sigma))
      }  
    } else{
      if (type == "gauss"){
        samp_pts <- sapply(mu_2, function(x) stats::rnorm(n, x, sigma))
      } else if (type == "laplace"){
        samp_pts <- sapply(mu_2, function(x) rlaplace(n, location = x, scale = sigma))
      }
    }
  })
  
  count_list <- apply(samp_pts, 1, function(x){
    if ( log_density_ratio_fn(mu_1, mu_2, x, sigma, eps, type) ){
      return(1)
    } else{
      return(0)
    }
  })
  return( sum(count_list)/n )
}


delta_fn <- function(u, eps){
  return( pnorm( - eps / u + u/2 ) - exp(eps) * pnorm( - eps / u - u/2 ) )
}

# u_list <- seq(0.1, 10, 0.1)
# eps_list <- seq(0.1, 10, 0.1)
# input_list <- expand.grid(u_list, eps_list)
# delta_list <- apply(input_list, 1, function(x) delta_fn(x[1], x[2]))
# dat <- data.frame(input_list)
# dat$u <- input_list[,1]
# dat$eps <- input_list[,2]
# dat$delta <- delta_list
# dat <- dat[, -c(1,2)]
# test_plot <- plot_ly(dat, x = ~u, y = ~eps, z = ~delta)
# test_plot <- test_plot %>% add_markers()
# test_plot <- test_plot %>% layout(scene = list(xaxis = list(title = 'mu'),
#                                                yaxis = list(title = 'eps'),
#                                                zaxis = list(title = 'delta')))
# # test_plot <- plot_ly() %>% add_trace(data = dat, x = ~u, y = ~eps, z = ~delta, type="mesh3d" )
# test_plot


intial_interval <- function(f, target, eps){
  k <- 0
  while ( f(2^k, eps) <= target ){
    k <- k + 1
  }
  return( c(0, 2^k) )
}

binary_search <- function(f, target, l, h, precis, eps){
  while ( (h - l) > precis ){
    m <- (l + h) /2
    if (f(m, eps) <= target){
      l <- m
    } else{
      h <- m 
    }
  }
  return(l)
}

gdp_to_dp_2 <- function(u){
  p <- pnorm(- u/2)
  return( log((1 - p) / p) )
}


gdp_u_alt <- function(sensi, sigma, d = 2, n = 3000, eps_max = 3, n_eps = 100, type = "gauss", m){
  
  
  mu_1 <- rep(0, d)
  mu_2 <- c(sensi, rep(0, d-1))

  
  u_list <- c()
  lhs_list <- rep(0, n_eps)
  lhs_list_1 <- rep(0, n_eps)
  lhs_list_2 <- rep(0, n_eps)
  
  eps_list <- seq(eps_max/n_eps, eps_max, eps_max/n_eps)
  for (i in 1:m){
    
    if (type == "gauss"){
      samp_pts_1 <- sapply(mu_1, function(x) stats::rnorm(n, x, sigma))
      samp_pts_2 <- sapply(mu_2, function(x) stats::rnorm(n, x, sigma))
    } else if (type == "laplace"){
      samp_pts_1 <- sapply(mu_1, function(x) rlaplace(n, location = x, scale = sigma))
      samp_pts_2 <- sapply(mu_2, function(x) rlaplace(n, location = x, scale = sigma))
    }
    
   
    lhs_list_tmp <- c()
    lhs_list_1_tmp <- c()
    lhs_list_2_tmp <- c()
    
    for (eps in eps_list){
      lhs_1 <- integral_fn(samp_pts_1, eps, mu_1, mu_2, sigma, type, int = 1, n)
      lhs_2 <- integral_fn(samp_pts_2, eps, mu_1, mu_2, sigma, type, int = 2, n)
      
      suppressWarnings(if (!is.na(samp_pts_1)){
        if (lhs_1 == 0 & lhs_2 == 0){
          break
        } 
      })
      
      lhs <-  max(lhs_1 - exp(eps) * lhs_2, 0)
      lhs_list_tmp <- cbind(lhs_list_tmp, lhs)
      lhs_list_1_tmp <- cbind(lhs_list_1_tmp, lhs_1)
      lhs_list_2_tmp <- cbind(lhs_list_2_tmp, lhs_2)
    }
    
    lhs_list_1_tmp <- c(lhs_list_1_tmp, rep(0, length(eps_list) - length(lhs_list_1_tmp)))
    lhs_list_2_tmp <- c(lhs_list_2_tmp, rep(0, length(eps_list) - length(lhs_list_2_tmp)))
    lhs_list_tmp <- c(lhs_list_tmp, rep(0, length(eps_list) - length(lhs_list_tmp)))
    
    lhs_list_1 <- (lhs_list_1_tmp + lhs_list_1)
    lhs_list_2 <- (lhs_list_2_tmp + lhs_list_2)
    lhs_list <- (lhs_list_tmp + lhs_list)
    
  }
  
  lhs_list_1 <- lhs_list_1 / m
  lhs_list_2 <- lhs_list_2 / m
  lhs_list <- lhs_list / m
  
  i = 1
  for (eps in eps_list){
    lhs <- lhs_list[i]
    if (lhs > 0 & lhs < 1){
      interval <- intial_interval(delta_fn, lhs, eps)
      u <- binary_search(delta_fn, lhs, interval[1], interval[2], 0.001, eps)
      u_list <- cbind(u_list, u)
    } else{
      u_list <- cbind(u_list, 0)
    }
    i <- i + 1
  }
  
  
  res_list <- list()
  res_list$u <- max(u_list)
  res_list$u_list <- c(u_list)
  res_list$eps_list <- eps_list[1:length(u_list)]
  res_list$lhs_list <- c(lhs_list)
  res_list$lhs_list_1 <- c(lhs_list_1)
  res_list$lhs_list_2 <- c(lhs_list_2)
  return(res_list)
}

# sensi <- 1
# sigma <-2
# k <- 2
# eps_max <- k * (sensi / sigma) - (sensi / sigma)^2 / 2
# res_list <- gdp_u_alt(sensi = sensi, sigma = sigma, d = 2, n = 1000, eps_max = eps_max, n_eps = 100, type = "gauss", m = 100); res_list$u # 1.15625
# # # res_list <- gdp_u_alt(sensi = sensi, sigma = sigma, d = 2, n = 100000, eps_max = eps_max, n_eps = 100, type = "gauss", m = 1); res_list$u # 1.15625
# lhs_list_true <- apply(cbind(res_list$eps_list, sensi/sigma), 1, function(x) delta_fn(x[2], x[1]))
# plot(res_list$eps_list, lhs_list_true, xlab = "eps", ylab = "delta", type ="l")
# points(res_list$eps_list, res_list$lhs_list, col = "red")
# plot(res_list$eps_list, lhs_list_true, xlab = "eps", ylab = "delta", ylim = c(0,0.01), type = "l")
# points(res_list$eps_list, res_list$lhs_list, col = "red")
# res_list <- gdp_u(sensi = 0.01, sigma = 0.01, d = 3, n = 5000, eps_max = 10, eps_step = 0.01, type = "gauss"); res_list$u # 1.202148

# lhs_list_true <- apply(cbind(res_list$eps_list, 1), 1, function(x) delta_fn(x[2], x[1]))
# plot(res_list$eps_list, lhs_list_true, xlab = "eps", ylab = "delta", type ="l")
# lhs_list_true <- apply(cbind(res_list$eps_list, 0.5), 1, function(x) delta_fn(x[2], x[1]))
# points(res_list$eps_list, lhs_list_true, xlab = "eps", ylab = "delta", type ="l", col = "red")
# lhs_list_true <- apply(cbind(res_list$eps_list, 2), 1, function(x) delta_fn(x[2], x[1]))
# points(res_list$eps_list, lhs_list_true, xlab = "eps", ylab = "delta", type ="l", col = "blue")

laplace_compare <- function(sensi, sigma, d = 3, n = 3000, eps_max = 10, n_eps = 100, m = 1){
  u_mcmc <- gdp_u_alt(sensi, sigma, d, n, eps_max, n_eps, type = "laplace", m = m)
  eps <- sensi / sigma
  u_true <- -2 * qnorm( 1/(1 + exp(eps)) )
  u_mcmc$u_true <- u_true
  return(u_mcmc)
}

#=== testing laplace_compare() ===
# sensi <- 1
# sigma <-4
# eps_max <- sensi / sigma
# res_list <- laplace_compare(sensi = sensi, sigma = sigma, d = 2, n = 1000, eps_max = eps_max, m = 50); res_list$u; res_list$u_true # 
# plot(res_list$lhs_list, res_list$u_list)

# compare_list <- sapply(1:30, function(x){
#   res <- laplace_compare(sensi = 1, sigma = 1, d=2, n = 5000, eps_max = 10)
#   return(c(res$u, res$u_true))
# }); compare_list
# apply(compare_list, 1, mean) # [1] 1.246527 1.136979
# apply(compare_list, 1, sd) # [1] 0.00000000 0.03805694

# mu_1 <- c(1, 0)
# n <- 500
# sigma <- 1
# samp_pts_1 <- sapply(mu_1, function(x) stats::rnorm(n, x, sigma))
# samp_pts_2 <- sapply(mu_1, function(x) rlaplace(n, location = x, scale = sigma))
# plot(samp_pts_1[,1], samp_pts_1[,2], col = "blue", xlim = c(-5, 5), ylim = c(-5, 5))
# points(samp_pts_2[,1], samp_pts_2[,2], col = "red")

dp_eps <- function(sensi, sigma, d = 2, n = 3000, eps_max = 3, eps_step = 0.01, type = "laplace"){
  
  
  mu_1 <- rep(0, d)
  mu_2 <- c(sensi, rep(0, d-1))
  
  lhs_list <- c()
  eps_res_list <- c()
  
  samp_pts_1 <- sapply(mu_1, function(x) rlaplace(n, location = x, scale = sigma))
  samp_pts_2 <- sapply(mu_2, function(x) rlaplace(n, location = x, scale = sigma))

  eps_list <- seq(eps_step, eps_max, eps_step)
  
  for (i in seq(n, n, 50)){
    
    samp_pts_1_tmp <- samp_pts_1[1:i,]
    samp_pts_2_tmp <- samp_pts_2[1:i,]
    
    eps_list <- seq(eps_step, eps_max, eps_step)
    
    for (eps in eps_list){
      lhs_1 <- integral_fn(samp_pts_1, eps, mu_1, mu_2, sigma, type)
      lhs_2 <- integral_fn(samp_pts_2, eps, mu_1, mu_2, sigma, type)
      if (lhs_1 == 0 & lhs_2 == 0){
        break
      }
      lhs <-  lhs_1 - exp(eps) * lhs_2
      lhs_list <- c(lhs_list, lhs)
    }
    
    eps_res <- eps_list[match(-1, diff(lhs_list > 0)) + 1]
    if(is.na(eps_res)){
      eps_res <- length(lhs_list) + 1
    }
    eps_res_list <- c(eps_res_list, eps_res)
  }
  res_list <- list()
  res_list$eps <- eps_res_list
  res_list$lhs_list <- lhs_list
  return(res_list)
}

# res_list <- dp_eps(1, 1, d = 3, n = 5000, eps_max = 10, eps_step = 0.01)


# ======================== analytic gaussian mechanism ========================

binary_search_2 <- function(f, target, min, l, h, precis, eps){
  if (!min){
    while ( (h - l) > precis ){
      m <- (l + h) /2
      if (f(m, eps) <= target){
        l <- m
      } else{
        h <- m 
      }
    }
    return(l)
  } else{
    while ( (h - l) > precis ){
      m <- (l + h) /2
      if (f(m, eps) <= target){
        h <- m
      } else{
        l <- m 
      }
      # print(c(l, h))
    }
    return(h)
  }
}

intial_interval_2 <- function(f, target, min, eps){
  k <- 0
  if(!min){
    while ( f(2^k, eps) <= target ){
      k <- k + 1
    }
  } else{
    while ( f(2^k, eps) >= target ){
      k <- k + 1
    }
  }
  return( c(0, 2^k) )
}

B_plus <- function(u, eps){
  return( pnorm( sqrt(eps * u) ) - exp(eps) * pnorm(- sqrt(eps * (u + 2)) ) )
}

B_minus <- function(u, eps){
  return( pnorm( - sqrt(eps * u) ) - exp(eps) * pnorm(- sqrt(eps * (u + 2)) ) )
}

analytic_gauss_mechanism <- function(eps, delta, sensi){
  delta_0 <- pnorm(0) - exp(eps) * pnorm(- sqrt(2 * eps))
  if (delta >= delta_0){
    # print("max")
    interval <- intial_interval_2(B_plus, delta, FALSE, eps)
    # print(interval)
    u <- binary_search_2(B_plus, delta, FALSE, interval[1], interval[2], 0.001, eps)
  } else{
    interval <- intial_interval_2(B_minus, delta, TRUE, eps)
    # print(interval)
    u <- binary_search_2(B_minus, delta, TRUE, interval[1], interval[2], 0.001, eps)
    # print(u)
  }
  alpha <- sqrt(1 + u/2) + sqrt(u/2)
  sd_dp <- alpha * sensi / sqrt(2 * eps)
  return(sd_dp)
}

# analytic_gauss_mechanism(eps = 1, delta = 0.01, sensi = 1)

adp_eps_delta <- function(sensi, sigma, d = 2, n = 3000, eps_max = 3, eps_step = 0.01, type = "gauss"){
  
  
  mu_1 <- rep(0, d)
  mu_2 <- c(sensi, rep(0, d-1))
  
  lhs_list <- c()
  eps_res_list <- c()
  
  samp_pts_1 <- sapply(mu_1, function(x) stats::rnorm(n, x, sigma))
  samp_pts_2 <- sapply(mu_2, function(x) stats::rnorm(n, x, sigma))  
  
  eps_list <- seq(eps_step, eps_max, eps_step)
  
  for (i in seq(n, n, 50)){
    
    samp_pts_1_tmp <- samp_pts_1[1:i,]
    samp_pts_2_tmp <- samp_pts_2[1:i,]
    
    for (eps in eps_list){
      lhs_1 <- integral_fn(samp_pts_1, eps, mu_1, mu_2, sigma, type)
      lhs_2 <- integral_fn(samp_pts_2, eps, mu_1, mu_2, sigma, type)
      if (lhs_1 == 0 & lhs_2 == 0){
        break
      }
      lhs <-  lhs_1 - exp(eps) * lhs_2
      lhs_list <- c(lhs_list, lhs)
    }
  }
  res_list <- list()
  res_list$eps_list <- eps_list[1:length(lhs_list)]
  res_list$lhs_list <- lhs_list
  return(res_list)
}

# sigma <- analytic_gauss_mechanism(eps = 0.5, delta = 0.01, sensi = 1)
# res_list <- adp_eps_delta(sensi = 1, sigma = sigma, d = 3, n = 5000, eps_max = 10, eps_step = 0.01)
# cbind(res_list$eps_list[1:length(res_list$lhs_list)], res_list$lhs_list)



