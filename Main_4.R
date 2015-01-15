####**********************************************************************
####  Written and Developed by: 
####**********************************************************************
####
####    Julien Duvanel, Copyright 2014
####    email: duvanel@stanford.edu
####
####**********************************************************************
####**********************************************************************

####**********************************************************************
####
####  Main file
####
####**********************************************************************

suppressPackageStartupMessages(library(devtools)) 

### Clean workspace, set directory and load functions
rm(list = ls())

setwd("/Users/julien/Dropbox/Ecole/EPFL/5eme annee/MA4 - PDM/PackageCH/")
#setwd("/home/duvanel/git/PackageCH/")

# Main R file
with_debug(load_all(pkg = "PackageCH"))

### Settings
options(digits.secs=10)
Project <- SetupProject()


####################
### Comparison DCOR
####################
param.list <- list()

n_Estimation <- 5
n_MAX <- 500

##########
## Compare 
##########
param.list[[1]] <- list(n = rep(500, times = n_Estimation),
                        s = seq(from = 50, to = n_MAX, by = n_MAX/n_Estimation),
                        t = rep(100, times = n_Estimation),
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),
                        variable = c("s"))

param.list[[2]] <- list(n = seq(from = 50, to = n_MAX, by = n_MAX/n_Estimation),
                        s = rep(500, times = n_Estimation),
                        t = rep(100, times = n_Estimation),
                        snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),
                        variable = c("n"))

# param.list[[1]] <- list(n = rep(500, times = n_Estimation),
#                         s = rep(500, times = n_Estimation),
#                         t = seq(from = 10, to = 100, by = 100 / n_Estimation),
#                         m = rep(1, times = n_Estimation),
#                         snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),
#                         variable = c("t"))
# 
# param.list[[2]] <- list(n = rep(500, times = n_Estimation),
#                         s = rep(500, times = n_Estimation),
#                         t = seq(from = 10, to = 100, by = 100 / n_Estimation),
#                         m = rep(2, times = n_Estimation),
#                         snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),
#                         variable = c("t"))
# 
# param.list[[3]] <- list(n = rep(500, times = n_Estimation),
#                         s = rep(500, times = n_Estimation),
#                         t = seq(from = 10, to = 100, by = 100 / n_Estimation),
#                         m = rep(3, times = n_Estimation),
#                         snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),
#                         variable = c("t"))
# 
# param.list[[4]] <- list(n = rep(500, times = n_Estimation),
#                         s = rep(500, times = n_Estimation),
#                         t = seq(from = 10, to = 100, by = 100 / n_Estimation),
#                         m = rep(4, times = n_Estimation),
#                         snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),
#                         variable = c("t"))


# param.list[[1]] <- list(n = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
#                         s = rep(500, times = n_Estimation),
#                         t = rep(10, times = n_Estimation),
#                         m = rep(1, times = n_Estimation),
#                         snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),
#                         variable = c("n"))                   

# param.list[[2]] <- list(n = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
#                         s = rep(500, times = n_Estimation),
#                         t = rep(10, times = n_Estimation),
#                         m = rep(2, times = n_Estimation),
#                         snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),
#                         variable = c("n"))                   
# 
# param.list[[3]] <- list(n = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
#                         s = rep(500, times = n_Estimation),
#                         t = rep(10, times = n_Estimation),
#                         m = rep(3, times = n_Estimation),
#                         snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),
#                         variable = c("n"))  
# 
# param.list[[4]] <- list(n = seq(from = 10, to = n_MAX, by = n_MAX/n_Estimation),
#                         s = rep(500, times = n_Estimation),
#                         t = rep(10, times = n_Estimation),
#                         m = rep(4, times = n_Estimation),
#                         snps_value = matrix(c(0,1,1,1,1,1,1,2), nrow = n_Estimation, ncol = 8, byrow = T),
#                         variable = c("n"))  



res <- list()
# Do estimation and generate parameters
for(i in 1:length(param.list)) {
    
    cat(paste0("= Param list: ", i, "/", length(param.list), "\n"))
    
    res[[i]] <- compare_dcor_multi(n = param.list[[i]]$n, 
                                   s = param.list[[i]]$s,
                                   t = param.list[[i]]$t, 
                                   snps_value = param.list[[i]]$snps_value,
                                   variable = param.list[[i]]$variable)
    
    cat(paste0("= End of ", i, "/", length(param.list), "\n"))
    
}

