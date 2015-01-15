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

n_Estimation <- 2
n_MAX <- 50

##########
## Compare variation of importance of delta dom and epi in a specific setup
##########
                

param.list[[1]] <- list(n = rep(n_MAX, times = n_Estimation),
                      s = floor(seq(from = 20, to = n_MAX, length.out = n_Estimation)),
                      u = rep(20, times = n_Estimation),
                      b = c(0),
                      noise.sd = rep(0, times = n_Estimation),
                      delta_add = rep(1, times = n_Estimation),
                      delta_dom = rep(0, times = n_Estimation),
                      delta_epi = rep(0, times = n_Estimation),
                      snps_value = matrix(c(0,1,2), nrow = n_Estimation, ncol = 3, byrow = T),
                      variable = c("s"))


#param.list[[2]] <- list(n = rep(n_MAX, times = n_Estimation),
#                      s = rep(n_MAX, times = n_Estimation),
#                      u = rep(20, times = n_Estimation),
#                      b = c(0),
#                      noise.sd = floor(seq(from = 0, to = 25, length.out = n_Estimation)),
#                      delta_add = rep(1, times = n_Estimation),
#                      delta_dom = rep(0, times = n_Estimation),
#                      delta_epi = rep(0, times = n_Estimation),
#                      snps_value = matrix(c(0,1,2), nrow = n_Estimation, ncol = 3, byrow = T),
#                      variable = c("noise.sd"))  

# 
# param.list[[1]] <- list(n = rep(n_MAX, times = n_Estimation),
#                        s = rep(n_MAX, times = n_Estimation),
#                        u = rep(20, times = n_Estimation),
#                        b = c(0),
#                        noise.sd = rep(0, times = n_Estimation),
#                        delta_add = rep(1, times = n_Estimation),
#                        delta_dom = seq(from = 0.0, to = 0.3, length.out = n_Estimation),
#                        delta_epi = seq(from = 0, to = 0, length.out = n_Estimation),   
#                        snps_value = matrix(c(0,1,2), nrow = n_Estimation, ncol = 3, byrow = T),
#                        variable = c("delta_dom"))  
# 
# # number of individuals
#param.list[[4]] <- list(n = floor(seq(from = 100, to = n_MAX, length.out = n_Estimation)),
#                       s = rep(n_MAX, times = n_Estimation),
#                       u = rep(20, times = n_Estimation),
#                       b = c(0),
#                       noise.sd = rep(0, times = n_Estimation),
#                       delta_add = rep(1, times = n_Estimation),
#                       delta_dom = rep(0, times = n_Estimation),
#                       delta_epi = rep(0, times = n_Estimation),
#                       snps_value = matrix(c(0,1,2), nrow = n_Estimation, ncol = 3, byrow = T),
#                       variable = c("n"))   


res <- list()
# Do estimation and generate parameters
for(i in 1:length(param.list)) {
    
    cat(paste0("= Param list: ", i, "/", length(param.list), "\n"))
    
    res[[i]] <- compare_dcor(n = param.list[[i]]$n, 
                             s = param.list[[i]]$s,
                             u = param.list[[i]]$u, 
                             b = param.list[[i]]$b, 
                             noise.sd = param.list[[i]]$noise.sd,
                             delta_add = param.list[[i]]$delta_add,
                             delta_dom = param.list[[i]]$delta_dom,
                             delta_epi = param.list[[i]]$delta_epi,
                             snps_value = param.list[[i]]$snps_value,
                             variable = param.list[[i]]$variable)
    
    cat(paste0("= End of ", i, "/", length(param.list), "\n"))
                       
}

