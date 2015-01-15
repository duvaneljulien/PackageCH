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
####  Tools
####
####**********************************************************************
library(matrixStats)
library(energy)
#' Generate a correlation matrix with parameters in argument.
#' Notice that corr can be a vector.
#'
#' @title Get the row ID of P according to the full_id
#' @param full_id the id that we want to match
#' @param P the matrix we want to find the row (must contains P$full_id)
#' @return row number of P$full_id == full_id
#' @author Julien Duvanel
#' @export
get_P_row_from_full_id <- function (full_id, P) {
  which(P$full_id == full_id)
}

#' Obtain a list of availale phenotypes
#'
#' @title List of available phenotypes
#' @param P the matrix that contain phenotypes
#' @return row id and phenotype name
#' @author Julien Duvanel
#' @export
get_available_phenotypes <- function(P) {
  cbind(4:ncol(P), colnames(P)[4:ncol(P)]) 
}

#' Filter P to keep the correct phenotypes (and remove NA)
#'
#' @title Filter P
#' @param P the matrix that contain phenotypes
#' @param phenotype ID
#' @return P filtered and without NA
#' @author Julien Duvanel
#' @export
filter_P_by_id <- function(P, id) {
  
  # Keep the correct column
  P <- as.data.frame(P[, c(which(colnames(P) == "full_id"), id)])
  
  # Return without NA
  as.data.frame(P[complete.cases(P),])
  
}

#' Filter P to keep the correct phenotypes (and remove NA)
#'
#' @title Filter P
#' @param P the matrix that contain phenotypes
#' @param phi.matrix the matrix containing the values for phi_pq
#' @return P and phi filtered and without NA
#' @author Julien Duvanel
#' @export
filter_data <- function(P, phi.matrix) {
  
  P.unique <- unique(P$full_id)
  
  # We expect that the number of people didn't change
  # i.e. that this number is unique. If not, this means we lost information
  expect_that(nrow(P), equals(length(P.unique)))
  
  phi.unique <- unique(colnames(phi.matrix))
  
  # Keep data where we have rows for P and K
  index <- intersect(P.unique, phi.unique)
  
  # Filter K.matrix and P
  phi.matrix.filtered <- phi.matrix[as.vector(index), 
                                    as.vector(index)]
  P.filtered <- P[as.vector(index),]
  
  expect_that(dim(phi.matrix.filtered), equals(c(nrow(P.filtered), nrow(P.filtered))))
  
  list(P.filtered = P.filtered[, c(which(colnames(P.filtered) != "full_id"))], 
       phi.matrix.filtered = phi.matrix.filtered)
  
}

#' Build the matrix V
#'
#' @title Build matrix V
#' @param P the matrix that contain one phenotype
#' @return matrix V
#' @author Julien Duvanel
#' @export
build_matrix_V <- function(P) {
#   
#   P <- as.matrix(P)
#   
#   ### We keep one phenotype and normalize its
#   V <- matrix(P, 
#               nrow = nrow(P),
#               ncol = ncol(P))
#   
#   V <- t(t(V) - colMeans(V))
#   V.sd <- colSds(V)
#   
#   # Return
#   t(t(V) * (1/V.sd))
  as.matrix(P)
}

#' Build the matrix K
#'
#' @title Build matrix K
#' @param K the matrix that contain theo kinship
#' @param P the matrix that contain phenotypes
#' @return matrix K with all relationship
#' @author Julien Duvanel
#' @export
build_matrix_K <- function(K) {
  
  K.unique <- unique(union(K$full_ID1, K$full_ID2))
  
  K.matrix <- matrix(data = 0, nrow = length(K.unique),
                     ncol = length(K.unique))
  
  diag(K.matrix) <- 1
  
  colnames(K.matrix) <- K.unique
  rownames(K.matrix) <- K.unique
  
  print(" = Building K, row by row")
  
  for(i in 1:nrow(K)) {
    
    if((i %% floor(nrow(K) / 10)) == 0) {
      print(paste0(" == ", i, " / ", nrow(K)))  
    }
    
    # Because this matrix is symmetric, we have to set 2 values
    K.matrix[K$full_ID1[i], 
             K$full_ID2[i]] <- K$value[i]
    
    K.matrix[K$full_ID2[i], 
             K$full_ID1[i]] <- K$value[i]
  }
  
  print(" = K built")
  
  list(K.matrix = K.matrix)
  
}

#' Wrapper to easily use one of the viewport option
#' 
#' @title Wrapper for viewport
#' @param x number of rows
#' @param y number of columns
#' @author Julien Duvanel
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)


#' R script to read the GRM binary file
#' 
#' @title R script to read the GRM binary file
#' @param prefix name of the files
#' @param allN .
#' @param size .
#' @author plink website
ReadGRMBin <- function(prefix, AllN = F, size = 4) {
  
  BinFileName <- paste(prefix,".grm.bin",sep="")
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  id <- read.table(IDFileName)
  
  n <- dim(id)[1]
  BinFile <- file(BinFileName, "rb")
  grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile <- file(NFileName, "rb")
  
  if(AllN == T){
    N <- readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  } else {
    N <- readBin(NFile, n=1, what=numeric(0), size=size)
  }
  
  i=sapply(1:n, function(i) sum(1:i))
  
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
  
}


#' Build G matrix 
#' 
#' @title G matrix construction
#' @param snps a matrix containing n rows and p snps
#' @return return a matrix of size n \times n
#' @author Julien Duvanel
build_matrix_G <- function(snps) {
        
    G <- matrix(0, 
                ncol = nrow(snps),
                nrow = nrow(snps))
    
    for(i in 1:nrow(snps)) {
        for(j in 1:nrow(snps)) {
            if(j %% 10 == 0) {
                cat("i = ", i, "/", nrow(snps), " and j = ", j, "/", nrow(snps), "\n")
            }
            G[i,j] <- sqrt(sum((snps[i,] - snps[j,]) ^ 2))
        }
    }
    
    G
}


#' Get frequency of the minor allele for the SNP snp
#' 
#' @title Get frequency of minor allele
#' @param G genomic matrix
#' @param snp a snp (index)
#' @return frequency of the minor allele
#' @author Julien Duvanel
get_frequency <- function(G, snp) {
    # Since we are interested by the minor allele
    # it is very easy to obtain it because it's just the sum of the values 
    # for this snp across all individuals divided by 2 times the sample size
    sum(G[,snp]) / (2 * nrow(G))
}

#' Build W a normalizd version of the genomic matrix
#' 
#' @title Build a normalized version of the genomic matrix
#' @param G a genomic matrix
#' @return W a normalized version of the genomic matrix
#' @author Julien Duvanel
compute_W <- function(G) {
    
    # Create W a matrix that has the same size as G
    W <- matrix(0, ncol = ncol(G), nrow = nrow(G))
    
    # For each SNPs
    for(i in 1:ncol(G)) {
        # Get the minor allele frequency of this SNP
        p_i <- get_frequency(G, i)
        
        # For each individuals
        for(j in 1:nrow(G)) {
            # Normalize it using minor allele frequency
            W[j,i] <- (G[j,i] - 2 * p_i) / sqrt(2*p_i*(1 - p_i))
        }    
    }
    
    # Return the matrix
    W
}

