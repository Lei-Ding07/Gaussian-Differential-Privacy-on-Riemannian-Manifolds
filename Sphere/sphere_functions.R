# library(plotly)

norm_fn <- function(v){
  return(sqrt(sum(v^2)))
}

euclid_dist_fn <- function(x, y){
  return(sqrt(sum((x - y)^2)))
}

# from url:https://scicomp.stackexchange.com/questions/27689/numerically-stable-way-of-computing-angles-between-vectors
distance_fn <- function(p, q){
  if (sum(p * q) > 1){
    print(c(p, q))
    warning("sum(p * q) > 1")
    return(0)
  } else{
    dist <- 2 * atan2(norm_fn(p - q), norm_fn(p + q))
    if (dist < 0){
      dist <- dist + pi
    }
    # dist <- acos(sum(p * q))
    return( dist )   
  }
  
}


exp_fn <- function(p, v, r){
  v_norm <- norm_fn(v)
  return( cos( v_norm ) * p + r * sin( v_norm ) / v_norm * v )
}

log_fn <- function(start, end, r){
  if (sum((start - end)^2) == 0){
    n <- length(start)
    return(rep(0, n))
  } else if (sum((start + end)^2) == 0){
    stop("q out of domain")
  } else{
    theta <- distance_fn(start, end)
    return( theta / sin(theta) * ( end - sum(start * end) * start ) )
  }
}


# ======================== sampling from distributions ========================

density_ratio_fn <- function(x, x_new, mu, sigma, type = "gauss"){
  if (type == "gauss"){
    return( exp( ( - distance_fn(x_new, mu)^2 + distance_fn(x, mu)^2 ) / (2 * sigma^2) ) )  
  } else if (type == "laplace"){
    return( exp( ( - distance_fn(x_new, mu) + distance_fn(x, mu) ) / sigma ) )  
  }
}

rSpherical_gauss_single <- function(mu, sigma){
  u <- 1
  r <- 0.5
  v_norm <- 2 * pi
  d <- length(mu)
  while (v_norm > pi | u > r){
    precision <- 1 / sigma^2 + 1 / pi
    sd <- 1 / sqrt(precision)
    v <- rnorm(d, 0, sd)
    v_norm_sq <- sum(v^2) 
    v_norm <- sqrt(v_norm_sq)
    r <- sin(v_norm) / v_norm / exp(- v_norm_sq / (2*pi) ) 
    u <- runif(1)
  }
  v <- v - sum(v * mu) / sum(mu^2) * mu
  return( exp_fn(mu, v, 1) )
}

rSpherical_gauss <- function(n, mu, sigma){
  samp_list <- t(sapply(1:n, function(x){
    rSpherical_gauss_single(mu, sigma)
  }))
}

rS1_gauss <- function(n, mu, sigma){
  samp_list <- truncnorm::rtruncnorm(n, a = -pi, b = pi, mean = 0, sd = sigma)
  samp_list <- sapply(samp_list, function(x){
    c(cos(x), sin(x))
  })
  # print(samp_list)
  rotation_mat <- matrix(c(mu[1], mu[2], -mu[2], mu[1]), 2, 2)
  return(t(rotation_mat %*% samp_list))
}

# test <- rS1_gauss(1000, c(0, -1), 2)
# plot(test[,1], test[,2], xlim = c(-1, 1), ylim = c(-1, 1))

rSpherical_laplace_single <- function(mu, sigma){
  u <- 1
  tau <- 0.5
  iter <- 0
  status = FALSE
  while(!status){
    u <- runif(1)
    y <- rSpherical_gauss_single(mu, sqrt(sigma))
    r <- distance_fn(mu, y)
    tau <- exp( (r^2 - 2*r - pi^2 + 2*pi) / (2 * sigma) ) 
    iter <- iter + 1
    
    if (u <= tau){
      status <- TRUE
    }
    if (iter >= 50){
      break
    }
  }
  if (status){
    return(y)
  } else{
    rSpherical_laplace_single_mh(mu, sigma)
  }
}

rSpherical_laplace_single_mh <- function(mu, sigma){
  n  <- length(mu)
  y0 <- mu + stats::rnorm(n, mean = 0, sd = sigma)
  y0 <- y0 / sqrt(sum(y0^2))
  
  for (i in 1:50){
    y_p <- y0 + stats::rnorm(n, mean = 0, sd = sigma)
    y_p <- y_p / sqrt(sum(y_p^2))
    
    p <- exp( (distance_fn(y0, mu) - distance_fn(y_p, mu)) / sigma )
    if (stats::runif(1) < p){
      y1 <- y_p
    } else {
      y1 <- y0
    }
    y0 <- y1
  }
  return(y0)
}


rSpherical_laplace <- function(n, mu, sigma){
  samp_list <- t(sapply(1:n, function(x){
    rSpherical_laplace_single(mu, sigma)
  }))
}

# === testing for rSpherical_laplace() & rSpherical_gauss() : ===

# sensi <- sensi_fn(10, 1, pi/8) # 0.004048673
# mu_1 <- c(1, 0, 0)
# mu_2 <- c(cos(sensi), sqrt(1 - (cos(sensi))^2), 0)
# test <- rSpherical_gauss(1000, c(1, 0, 0), 10)
# test <- data.frame(test)
# test_plot <- plot_ly(test, x = ~X1, y = ~X2, z = ~X3)
# test_plot <- test_plot %>% add_markers()
# test_plot <- test_plot %>% layout(scene = list(xaxis = list(title = 'X'),
#                                                yaxis = list(title = 'Y'),
#                                                zaxis = list(title = 'Z')))
# test_plot <- test_plot %>% layout(scene = list(xaxis = list(range = c(-1, 1)),
#                                                yaxis = list(range = c(-1, 1)),
#                                                zaxis = list(range = c(-1, 1))))
# test_plot <- test_plot %>% add_trace(x = mu_1[1], y = mu_1[2], z = mu_1[3], name = 'center1', mode = 'markers', color = "red", symbols = c("triangle"))
# test_plot
# test_1 <- rSpherical_laplace(500, mu_1, 5)
# test_2 <- rSpherical_laplace(500, mu_2, 5)
# test <- rbind(test_1, test_2)
# dist <- c(rep(1, 500), rep(0, 500))
# test <- cbind(test, dist)
# test <- data.frame(test)
# colnames(test) <- c("X1", "X2", "X3", "dist")
# # par(mfrow = c(1,1))
# test_plot <- plot_ly(test, x = ~X1, y = ~X2, z = ~X3, color = ~dist, colors = c('#BF382A', '#0C4B8E'))
# test_plot <- test_plot %>% add_markers()
# test_plot <- test_plot %>% layout(scene = list(xaxis = list(title = 'X'),
#                                                yaxis = list(title = 'Y'),
#                                                zaxis = list(title = 'Z')))
# test_plot <- test_plot %>% layout(scene = list(xaxis = list(range = c(-1, 1)),
#                                                yaxis = list(range = c(-1, 1)),
#                                                zaxis = list(range = c(-1, 1))))
# test_plot <- test_plot %>% add_trace(x = mu_1[1], y = mu_1[2], z = mu_1[3], name = 'center1', mode = 'markers', color = "red", symbols = c("triangle"))
# test_plot <- test_plot %>% add_trace(x = mu_2[1], y = mu_2[2], z = mu_2[3], name = 'center2', mode = 'markers', color = "green", symbols = c("triangle"))
# test_plot
# 
# test1 <- t(apply(test, 1, function(x){
#   theta <- acos(x[3])
#   if (x[2] >= 0){
#     phi <- acos(x[1] / sqrt(x[1]^2 + x[2]^2))
#   } else{
#     phi <- - acos(x[1] / sqrt(x[1]^2 + x[2]^2))
#   }
#   return(c(theta, phi))
# }))
# T <- 1000
# # par(mfrow = c(2,1))
# hist(test1[,1])
# hist(test1[,2])


# ======================== monte carlo integration ========================

distance_diff_fn <- function(mu_top, mu_bottom, x, type = "gauss"){
  if (type == "gauss"){
    return( - distance_fn(x, mu_top)^2 + distance_fn(x, mu_bottom)^2 )  
  } else if (type == "laplace"){
    return( - distance_fn(x, mu_top) + distance_fn(x, mu_bottom) )  
  }
}

log_density_ratio_fn <- function(dist_diff, sigma, eps, type = "gauss"){
  if (type == "gauss"){
    return( dist_diff >= (2 * sigma^2 * eps) )  
  } else if (type == "laplace"){
    return( dist_diff >= (sigma * eps) )
  }
}

integral_fn <- function(dist_diff, eps, sigma, type){
  n <- length(dist_diff)
  count_list <- sapply(dist_diff, function(x){
    if ( log_density_ratio_fn(x, sigma, eps, type) ){
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

# === testing for intial_interval_2() and binary_search_2() ===
# interval <- intial_interval_2(delta_fn, 1e-3, 1)
# eps <- binary_search_2(delta_fn, 1e-3, interval[1], interval[2], 0.001, 1)


sensi_fn <- function(n, r = 1, r_0){
  h <- 2 * r_0 / r / tan(2 * r_0 / r)
  return( 2 * r_0 * (2 - h) / (n * h) )
}


gdp_u_alt <- function(sensi, sigma, r = 1, n = 3000, eps_max = 3, n_eps = 100, type = "gauss", m, d = 3){
  
  
  mu_1 <- c(c(r, 0), rep(0, d-2))
  mu_2 <- c(c(cos(sensi) / r, sqrt(r^2 - (cos(sensi) / r)^2)), rep(0, d-2))
  
  u_list <- c()
  lhs_list <- rep(0, n_eps)
  lhs_list_1 <- rep(0, n_eps)
  lhs_list_2 <- rep(0, n_eps)
  
  eps_list <- seq(eps_max/n_eps, eps_max, eps_max/n_eps)
  
  for (i in 1:m){
    
    print(Sys.time())
    print(paste0("============= m: ", i, " ============= "))
    
    if (type == "laplace"){
      samp_pts_1 <- rSpherical_laplace(n, mu = mu_1, sigma = sigma)
      samp_pts_2 <- rSpherical_laplace(n, mu = mu_2, sigma = sigma)
    } else if (type == "gauss"){
      if (d == 2){
        samp_pts_1 <- rS1_gauss(n, mu = mu_1, sigma = sigma)
        samp_pts_2 <- rS1_gauss(n, mu = mu_2, sigma = sigma)  
      } else{
        samp_pts_1 <- rSpherical_gauss(n, mu = mu_1, sigma = sigma)
        samp_pts_2 <- rSpherical_gauss(n, mu = mu_2, sigma = sigma) 
      }
    }
    
    
    dist_diff_1 <- apply(samp_pts_1, 1, function(x){
      distance_diff_fn(mu_1, mu_2, x, type = "gauss")
    })
    
    dist_diff_2 <- apply(samp_pts_2, 1, function(x){
      distance_diff_fn(mu_1, mu_2, x, type = "gauss")
    })
    
    lhs_list_tmp <- c()
    lhs_list_1_tmp <- c()
    lhs_list_2_tmp <- c()
    
    
    for (eps in eps_list){
      lhs_1 <- integral_fn(dist_diff_1, eps, sigma, type)
      lhs_2 <- integral_fn(dist_diff_2, eps, sigma, type)
      
      
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



# sensi <- sensi_fn(10, 1, pi/8)
# sigma <-1 * sensi
# eps_max <- 10
# res_list <- gdp_u_alt(sensi = sensi, sigma = sigma, r = 1, n = 1000, eps_max = eps_max, n_eps = 1000, type = "gauss", m = 100, d = 3); res_list$u # 1.116211
# lhs_list_true <- apply(cbind(res_list$eps_list, sensi/sigma), 1, function(x) delta_fn(x[2], x[1]))
# plot(res_list$eps_list, lhs_list_true, xlab = "eps", ylab = "delta", type ="l")
# points(res_list$eps_list, res_list$lhs_list, col = "red")
# plot(res_list$eps_list, lhs_list_true, xlab = "eps", ylab = "delta", type ="l", ylim = c(0, 0.005))
# points(res_list$eps_list, res_list$lhs_list, col = "red")



dp_eps <- function(sensi, sigma, r = 1, n = 3000, n_burn = 3000, eps_max = 3, eps_step = 0.01, type = "laplace", width = 1, d = 3){
  
  mu_1 <- c(c(r, 0), rep(0, d-2))
  
  mu_2 <- c(c(cos(sensi) / r, sqrt(r^2 - (cos(sensi) / r)^2)), rep(0, d-2))
  
  lhs_list <- c()
  eps_res_list <- c()
  
  samp_pts_1 <- rSpherical_laplace(n, mu = mu_1, sigma = sigma)
  samp_pts_2 <- rSpherical_laplace(n, mu = mu_2, sigma = sigma)
  
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
  return(eps_res_list)
}

# === testing dp_eps() ===
# sensi <- sensi_fn(300, 1, pi/8) # 0.004048673
# lhs_list <- dp_eps(sensi, 0.004, n = 10000, n_burn = 10000, eps_max = 10, eps_step = 0.01, width = 10); lhs_list # 1.02 

# sensi <- sensi_fn(10, 1, pi/8) # 0.1214602
# lhs_list <- dp_eps(sensi, 0.12, n = 10000, n_burn = 10000, eps_max = 10, eps_step = 0.01, width = 10); lhs_list # 1.02
# lhs_list <- dp_eps(sensi, 0.24, n = 10000, n_burn = 10000, eps_max = 10, eps_step = 0.01, width = 10); lhs_list # 0.51


laplace_compare <- function(sensi, sigma, n = 3000, eps_max = 3, eps_step = 0.01){
  u_mcmc <- gdp_u(sensi, sigma, n = n, eps_max = 10, eps_step = eps_step, type = "laplace")
  eps <- sensi / sigma
  u_true <- -2 * qnorm( 1/(1 + exp(eps)) )
  u_mcmc$u_true <- u_true
  return(u_mcmc)
}

#=== testing laplace_compare() ===
# sensi <- sensi_fn(10, 1, pi/8) # 0.12146018
# res_list <- laplace_compare(sensi, 0.12, n = 10000, eps_max = 10, eps_step = 0.005) # 1.2465267 0.7558594
# res_list <- laplace_compare(sensi = 1, sigma = 1, n = 10000, eps_max = 10, eps_step = 0.005) # 1.2465267 0.7558594

# sensi <- sensi_fn(300, 1, pi/8) # 0.004048673
# laplace_compare(sensi, 0.004, n = 10000, eps_max = 10, eps_step = 0.005) # 1.246527 0.781250
# res_list <- laplace_compare(sensi, 0.008, n = 10000, eps_max = 10, eps_step = 0.005) # 0.6314167 0.3867188
# laplace_compare(sensi, 0.002, n = 10000, eps_max = 10, eps_step = 0.005) # 2.383585 1.460938


# compare_list <- sapply(1:30, function(x) laplace_compare(sensi, 0.004, n = 10000, n_burn = 10000, eps_max = 10, width = 10)); compare_list
# apply(compare_list, 1, mean) 
# apply(compare_list, 1, sd)


# ======================== Laplace mechanism vs Gauss mechanism ========================

frechet_mean <- function(data, r = 1, n_iter = 500, step_size = 0.5, lambda = 1e-5){
  
  mean_old <- data[sample(1:nrow(data), 1), ]
  diff <- 1
  iter <- 1
  while (diff > lambda & iter <= n_iter){
    
    v_new <- apply(apply(data, 1, function(x){
      log_fn(start = mean_old, end = x, r) # need to double check
    }), 1, mean)
    
    mean_new <- exp_fn(mean_old, step_size * v_new, r)
    
    diff <- distance_fn(mean_new, mean_old)
    mean_old <- mean_new
    
    iter <- iter + 1
    # print(c(diff, mean_old, iter))
  }
  
  return(mean_old)
}

# === testing for frechet_mean(): ===
# test <- rball(100, 1, pi / 8)
# mean <- frechet_mean(test)
# test_plot <- scatterplot3d::scatterplot3d(test, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1,1))
# test_plot$points3d(mean, col = "red")


rball <- function(n, r, r_0){
  data <- cbind(runif(n, 0, r_0), runif(n, 0, 2 * pi))
  data <- cbind( r * sin(data[, 1]) * cos(data[, 2]), r * sin(data[, 1]) * sin(data[, 2]), r * cos(data[, 1]) )
  return(data)
}

# === testing for rball(): ===
# test <- rball(100, 1, pi / 8)
# scatterplot3d::scatterplot3d(test, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1,1))


gdp_to_dp <- function(u){
  p <- pnorm(- u/2)
  return( log((1 - p) / p) )
}


frechet_mean_compare <- function(n, sigma, r, r_0, u = NA){
  
  data <- rball(n, r, r_0)
  sensi <- sensi_fn(n, r, r_0)
  
  mean_true <- frechet_mean(data)
  
  if (is.na(u)){
    u <- gdp_u_alt(d, sensi, sigma)  
  }
  
  eps <- gdp_to_dp(u)
  mean_gdp <- rSpherical_gauss(1, mu = mean_true, sigma = sigma)
  mean_laplace <- rSpherical_laplace(1, mu = mean_true, sigma = sensi/eps)
  
  sensi_euclid <- 2 * sin(sensi/2)
  mean_gdp_euclid <- mean_true + sapply(1:3, function(x) rnorm(1, 0, sensi_euclid / u))
  
  return(c( distance_fn(mean_gdp, mean_true), distance_fn(mean_laplace, mean_true), euclid_dist_fn(mean_gdp, mean_true), euclid_dist_fn(mean_laplace, mean_true), euclid_dist_fn(mean_gdp_euclid, mean_true) ))
}

# sensi <- sensi_fn(10, 1, pi/8)
# sigma <-1 * sensi
# eps_max <- 10
# res_list <- gdp_u_alt(sensi = sensi, sigma = sigma, r = 1, n = 1000, eps_max = eps_max, n_eps = 1000, type = "gauss", m = 100); res_list$u # [1] 1.118164
# frechet_mean_compare(10, sigma, 1, pi/8, res_list$u)
# compare_list <- t(sapply(1:100, function(x) frechet_mean_compare(10, sigma, 1, pi/8, res_list$u)))
# apply(compare_list, 2, mean) # [1] 0.1530110 0.2852313
# apply(compare_list, 2, sd) # [1] 0.08031061 0.18316671

# ======================== comparsion on the unit circle S1 ========================

lhs_delta_fn <- function(eps, sigma, sensi){
  x1 <- sigma * eps / sensi
  x2 <- sensi / (2 * sigma)
  x3  <- pi / sigma
  eps_exp <- exp(eps)
  
  C <- 1 / (pnorm(x3) - pnorm(-x3))
  
  if (eps <= sensi^2 / (2 * sigma^2)){
    return( C * (pnorm(-x1 + x2) - eps_exp * pnorm(-x1-x2) - pnorm(x1+x2-x3) + eps_exp * pnorm(x1-x2+x3)) - eps_exp )
  } else{
    return( C * (pnorm(-x1 + x2) - eps_exp * pnorm(-x1-x2) - pnorm(x1+x2-x3) + eps_exp * pnorm(x1-x2-x3)))    
  }
}

gdp_u_true <- function(sensi, sigma, r=1, eps_max = 10, n_eps = 1000, type = "gauss", d = 3){
  
  mu_1 <- c(c(r, 0), rep(0, d-2))
  mu_2 <- c(c(cos(sensi) / r, sqrt(r^2 - (cos(sensi) / r)^2)), rep(0, d-2))
  
  eps_max <- min(eps_max, pi * sensi / (2 * sigma^2))
  eps_list <- seq(eps_max/n_eps, eps_max, eps_max/n_eps)
  lhs_list <- c()
  u_list <- c()
  
  for (eps in eps_list){
    lhs <- lhs_delta_fn(eps, sigma, sensi)
    lhs_list <- cbind(lhs_list, lhs)
    if (lhs > 0 & lhs < 1){
      interval <- intial_interval(delta_fn, lhs, eps)
      u <- binary_search(delta_fn, lhs, interval[1], interval[2], 0.001, eps)
      u_list <- cbind(u_list, u)
    } else{
      u_list <- cbind(u_list, 0)
    }
  }
  
  res_list <- list()
  res_list$u <- max(u_list)
  res_list$u_list <- c(u_list)
  res_list$eps_list <- eps_list[1:length(u_list)]
  res_list$lhs_list <- c(lhs_list)
  return(res_list)
}

# sensi <- sensi_fn(10, 1, pi/8)
# 
# sensi <- 1
# sigma <- 5 * sensi
# eps_max <- pi * sensi / (2 * sigma^2)
# res_list_2 <- gdp_u_true(sensi = sensi, sigma = sigma, r = 1, eps_max = eps_max, n_eps = 1000, type = "gauss", d = 2); res_list_2$u # 1.15625
# res_list <- gdp_u_alt(sensi = sensi, sigma = sigma, r = 1, n = 1000, eps_max = eps_max, n_eps = 1000, type = "gauss", m = 100, d = 2); res_list$u # 1.15625
# lhs_list_true <- res_list_2$lhs_list
# # lhs_list_true <- apply(cbind(res_list$eps_list, sensi/sigma), 1, function(x) delta_fn(x[2], x[1]))
# plot(res_list_2$eps_list, lhs_list_true, xlab = "eps", ylab = "delta", type ="l")
# points(res_list$eps_list, res_list$lhs_list, col = "red")
# # plot(res_list$eps_list, res_list$lhs_list, xlab = "eps", ylab = "delta", col ="red")
# #points(res_list$eps_list, lhs_list_true, col = "black", type ="l")
# plot(res_list$eps_list, lhs_list_true, xlab = "eps", ylab = "delta", type ="l", ylim = c(0, 0.005))
# points(res_list$eps_list, res_list$lhs_list, col = "red")

 
