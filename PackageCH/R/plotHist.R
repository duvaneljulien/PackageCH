# #setwd("/home/duvanel/PackageCH/")
# method <- c("GCTA", "TheoKin", "IBD", "IBS")
# #load("~/PackageCH/results/02102014_164924_results_IBD.RData")
# 
# for(i in 1:length(method)) {
#     
#     pdf(file = paste0("results/Hist", method[i], ".pdf"), 
#                       width = 17, height = 7)
#     
#         par(mfrow=c(1,2))
#         hist(as.numeric(as.vector(full.res[[i]][,2])), 
#              breaks = 20, 
#              xlab = "Value of estimate",
#              main = "Histogram of heritability estimates")
#     
#         hist(as.numeric(as.vector(full.res[[i]][,3])), 
#          breaks = 20, 
#          xlab = "Value of estimate",
#          main = "Histogram of heritability estimates dcov")
#     
#     dev.off()
#     
# }