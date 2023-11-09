rm(list =ls())
path <- "/Users/yangdi/Library/Mobile Documents/com~apple~CloudDocs/My Documents/Graduate School Materials/Graduate Projects/GDP on Manifold/sphere/"
# path <- "C:\\Users\\jiang\\iCloudDrive\\My Documents\\Graduate School Materials\\Graduate Projects\\GDP on Manifold\\sphere\\"
source(paste0(path, "sphere_functions.R"))


n_trial <- 20

sensi <- 1
n_eps <- 1000
n <- 1000
m <- 100

for (z in seq(0.25, 4, 0.25)){
  sigma <- z * sensi
  eps_max <- min(pi * sensi / (2 * sigma^2), 10)
  # parallelization 
  no_cores <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(no_cores)
  parallel::clusterExport(cl, c("path", "sensi", "sigma", "eps_max", "n_eps", "n", "m", "n_trial", "z"))
  parallel::clusterEvalQ(cl, source(paste0(path, "sphere_functions.R")))
  
  
  print(Sys.time())
  # gdp_u_alt(sensi = sensi, sigma = sigma, r = 1, n = n, eps_max = eps_max, n_eps = n_eps, type = "gauss", m = 2, d = 2)$u
  # print(Sys.time())
  
  res_list <- parallel::parSapply(cl, 1:n_trial, function(i){
    # dist_list <- sapply(1:100, function(i){  
    gdp_u_alt(sensi = sensi, sigma = sigma, r = 1, n = n, eps_max = eps_max, n_eps = n_eps, type = "gauss", m = m, d = 2)$u
  })
  
  k <- 1
  
  u_true <- gdp_u_true(sensi = sensi, sigma = sigma, r = 1, eps_max = eps_max, n_eps = n_eps, type = "gauss", d = 2)$u
  res_list <- cbind(rep(u_true, n_trial), res_list)
  write.csv(res_list, file = paste0(path, "s1_result_", k, "_", sensi, "_", z, ".csv"))
  parallel::stopCluster(cl)
  print(Sys.time())
}