library(doMC)
library(foreach)

build_matrix_G_manually <- function(path.snps = "data/SNPs_14102014_145630.raw") {
    
    N <- 20
    
    # Get results (stored into SNPS_datetime.stamp.raw)
    Phenotypes <- read.table(path.snps, header = TRUE)
    # Remove rows with NA's 
    Phenotypes <- Phenotypes[complete.cases(Phenotypes), ]
    # We must have "famid" and "id" as the two first columns
    colnames(Phenotypes) <- c("famid", "id")
    
    # Value of SNPs are stored in the last N columns
    snps <- Phenotypes[, (ncol(Phenotypes)-N+1):ncol(Phenotypes)]
    
    # Way (!) faster than doing two for loop
    G <- dist(snps)
} 

HeritabilityEstimation <- function(P, G) {
    
    P <- as.matrix(P)
    
    ### We keep one phenotype and normalize its
    V <- matrix(P, 
                nrow = nrow(P),
                ncol = ncol(P))
    
    V <- t(t(V) - colMeans(V))
    V.sd <- colSds(V)
    
    # Return
    V <- t(t(V) * (1/V.sd))
    
    # Get data for the specific phenotype   
    Y.distance.carre <- matrix(rep(V^2,nrow(V)),
                               ncol = nrow(V)) + 
        t(matrix(rep(V^2,nrow(V)),
                 ncol = nrow(V))) - 
        2 * V %*%t(V)
    
    Y.distance <- sqrt(Y.distance.carre)
    
    dcov(Y.distance, G) / (sqrt(dcov(G, G)) * sqrt(dcov(Y.distance, Y.distance)))
    
}
