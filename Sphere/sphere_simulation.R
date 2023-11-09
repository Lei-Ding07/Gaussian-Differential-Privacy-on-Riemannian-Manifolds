rm(list =ls())
path <- "/Users/yangdi/Library/Mobile Documents/com~apple~CloudDocs/My Documents/Graduate School Materials/Graduate Projects/GDP on Manifold/sphere/"
# path <- "C:\\Users\\jiang\\iCloudDrive\\My Documents\\Graduate School Materials\\Graduate Projects\\GDP on Manifold\\sphere\\"
source(paste0(path, "sphere_functions.R"))


print(Sys.time())
n <- 10
d <- 2
r <- pi / 8
sensi <- sensi_fn(10, 1, r)
z_list <- seq(0.25, 4, 0.25)

no_cores <- parallel::detectCores() - 2
cl <- parallel::makeCluster(no_cores)

parallel::clusterExport(cl, c("path", "n", "d", "r", "sensi"))
parallel::clusterEvalQ(cl, source(paste0(path, "sphere_functions.R")))

print(Sys.time())
u_list <- t(parallel::parSapply(cl, z_list, function(z){
  
  sigma <- z * sensi
  eps_max <- max(5 * sensi / sigma + sensi^2 / (2 * sigma^2), 10)
  res_list <- gdp_u_alt(sensi = sensi, sigma = sigma, r = 1, n = 1000, eps_max = eps_max, n_eps = 1000, type = "gauss", m = 100)
  u <- res_list$u
}))
parallel::stopCluster(cl)
print(Sys.time())

n_trial <- 1000
i <- 1
for (z in z_list){
  
  sigma <- z * sensi
  eps_max <- max(5 * sensi / sigma + sensi^2 / (2 * sigma^2), 10)
  
  # parallelization 
  no_cores <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(no_cores)
  
  parallel::clusterExport(cl, c("path", "n", "n_trial",  "d", "sigma", "r", "u_list", "i"))
  parallel::clusterEvalQ(cl, source(paste0(path, "sphere_functions.R")))
  
  print(Sys.time())
  dist_list <- t(parallel::parSapply(cl, 1:n_trial, function(x){
    # dist_list <- t(sapply(1:100, function(i){  
    frechet_mean_compare(n, sigma, 1, r, u_list[i])
  }))
  dist_list <- cbind(rep(u_list[i], n_trial), dist_list)
  write.csv(dist_list, file = paste0(path, "sphere_compare_result_", z, ".csv"))
  parallel::stopCluster(cl)
  print(Sys.time())
  i = i + 1
}

# res <- c()
# for (i in z_list){
#   res_list <- read.csv(paste0(path, "sphere_compare_result_", i, ".csv"))
#   res_list <- res_list[, -1]
#   res <- rbind(res, res_list)
#   print(apply(res, 2, mean))
#   print(apply(res, 2, sd))
# }


