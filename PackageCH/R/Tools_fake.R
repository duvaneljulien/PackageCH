library(energy)
library(matrixStats)
library(reshape2)
library(plyr)

#' Build alpha s.t. X = G * alpha
#'
#' @title Build alpha
#' @param u number of non_zero coeff
#' @param s number of SNPS
#' @param b value of runif(b[1],b[2])
#' @return alpha
#' @author Julien Duvanel
#' @export
build_alpha <- function(u, s, b) {
    
    # build alpha
    if(length(b) == 2) {
        alpha <- sample(c(runif(ifelse(u == 1, 2, u), b[1], b[2]), 
                          rep(x = 0, times = s - u)), size = s)
    } else {
        alpha <- sample(c(rnorm(ifelse(u == 1, 2, u), 0, 1), 
                          rep(x = 0, times = s - u)), size = s)
    }
    
}

#' Build alpha multi s.t. X = G * alpha 
#'
#' @title Build alpha
#' @param s number of SNPS
#' @param t number of traits
#' @param mode 1 = random, 2 = disjoint, 3 = half random / half noise, 4 = half disjoint / half noise
#' @return alpha
#' @author Julien Duvanel
#' @export
build_alpha_multi <- function(s, t, mode = 1) {
    
    if(t %% 2 != 0) stop("t has to be an even number !")
    
    if(mode == 1) {
        # Completly random mode
        alpha <- sapply(1:t, function(x) runif(n = s, min = -1, max = 1))
    } else if(mode == 2) {
        # We can have floor ( s / t ) number of SNPs for each traits (so that they are disjoint)
        alpha <- sapply(1:t, function(x) {
            c(rep(0, times = (x-1)*floor(s/t)),
              runif(n = floor(s/t), min = -1, max = 1),
              rep(0, times = s - (x*floor(s/t))))
        })
    } else if(mode == 3) {
        # half random / half noise (i.e. 0 and noise is added later)
        alpha <- cbind(sapply(1:(t/2), function(x) runif(n = s, min = -1, max = 1)),
                       sapply(1:(t/2), function(x) rep(0, times = s)))
    } else if(mode == 4) {
        # half disjoint / half noise
        t_half <- t/2
        alpha <- cbind(sapply(1:t_half, function(x) {
            c(rep(0, times = (x-1)*floor(s/t_half)),
              runif(n = floor(s/t_half), min = -1, max = 1),
              rep(0, times = s - (x*floor(s/t_half))))}),
                       sapply(1:t_half, function(x) rep(0, times = s)))
    } else {
        # Only noise
        alpha <- sapply(1:t, function(x) rep(0, times = s))
    }
    
    alpha
    
}


#' Build alpha s.t. X = G * alpha but alpha(g_i)
#'
#' @title Build alpha
#' @param u number of non_zero coeff
#' @param s number of SNPS
#' @param b value of runif(b[1],b[2])
#' @return alpha with dominant effect
#' @author Julien Duvanel
#' @export
build_alpha_dominant <- function(u, s, b) {
    
    # We create an alpha for the dominant effect
    alpha_d <- matrix(0, nrow = 3, ncol = s)
    alpha_d[1, ] <- sample(c(runif(ifelse(u == 1, 2, u), -1, 1), 
                             rep(x = 0, times = s - u)), size = s)
    alpha_d[2, ] <- sample(c(runif(ifelse(u == 1, 2, u), 3, 4), 
                             rep(x = 0, times = s - u)), size = s)
    alpha_d[3, ] <- sample(c(runif(ifelse(u == 1, 2, u), 15, 20), 
                             rep(x = 0, times = s - u)), size = s)
    
    alpha_d
    
}

#' Build alpha s.t. X = G * alpha but alpha(g_i,g_j)
#'
#' @title Build alpha
#' @param u number of non_zero coeff
#' @param s number of SNPS
#' @param b value of runif(b[1],b[2])
#' @return alpha with dominant effect
#' @author Julien Duvanel
#' @export 
build_alpha_epistatic <- function(u, s, b) {
    
    # We create an alpha for the dominant effect
    alpha_d <- matrix(0, nrow = 4, ncol = s)
    
    # If product of g_i and g_j = 0
    alpha_d[1, ] <- sample(c(runif(ifelse(u == 1, 2, u), -1, 1), 
                             rep(x = 0, times = s - u)), size = s)
    # If product of g_i and g_j = 1
    alpha_d[2, ] <- sample(c(runif(ifelse(u == 1, 2, u), 3, 4), 
                             rep(x = 0, times = s - u)), size = s)
    # If product of g_i and g_j = 2
    alpha_d[3, ] <- sample(c(runif(ifelse(u == 1, 2, u), 13, 15), 
                             rep(x = 0, times = s - u)), size = s)
    # If product of g_i and g_j = 4
    alpha_d[4, ] <- sample(c(runif(ifelse(u == 1, 2, u), 15, 20), 
                             rep(x = 0, times = s - u)), size = s)    
    alpha_d
}

#' Build X using alpha and M with multi traits and mode
#'
#' @title Build X
#' @param M genome matrix
#' @param alpha matrix of traits
#' @param m the mode
#' @return X 
#' @author Julien Duvanel
#' @export 
build_X_multi <- function(M, alpha) {
    
    # Build fake trait X
    X <- product_snps_alpha(M, alpha)
    
    # We are more interested to have X + noise 
    # noise is 10% of sd(X)
    noise <- sapply(1:ncol(alpha), function(x) rnorm(nrow(M), 
                                                     mean = 0, 
                                                     sd = 0.1 * ifelse(sd(X) == 0, 1, sd(X))))  
    
    # Normalize X 
    X_norm_plus_noise <- t(t(X + noise) - colMeans(X + noise))
    X.sd <- colSds(X + noise)
    
    # what we return
    t(t(X_norm_plus_noise) * (1/X.sd))
    
}

#' mat * alpha but alpha(mat)
#'
#' @title Product between mat and alpha for additive effect
#' @param M a matrix with s columns
#' @param alpha a vector 
#' @return product of mat * alpha
#' @author Julien Duvanel
#' @export
product_snps_alpha <- function(M, alpha) { M %*% alpha }

#' mat * alpha but alpha(mat)
#'
#' @title Product between mat and alpha for dominant effect
#' @param M a matrix with s columns
#' @param alpha a matrix with 3 rows and s columsn
#' @return product of mat * alpha
#' @author Julien Duvanel
#' @export
product_snps_alpha_dominant <- function(M, alpha) {
    # Usually, as.matrix(mat) %*% alpha does the job
    # but here, the multiplication differs if mat[i,j] has 
    # a specific value.
    
    # For every row, we do the "special product"
    res <- apply(M, 1, function(x) {
        ret <- 0
        for(i in 1:length(x)) {
            ret <- ret + x[i]*alpha[x[i]+1, i]
        }
        ret
    })
    
    as.matrix(res)
}

#' mat * alpha but alpha(mat)
#'
#' @title Product between mat and alpha for epistatic effect
#' @param M a matrix with s columns
#' @param alpha a matrix with 3 rows and s columsn
#' @return product of mat * alpha
#' @author Julien Duvanel
#' @export
product_snps_alpha_epistatic <- function(M, alpha) {
    # Usually, as.matrix(mat) %*% alpha does the job
    # but here, the multiplication differs if mat[i,j] has 
    # a specific value.
    
    f <- function(x) {
        if(x == 0) 1
        else if(x == 1) 2
        else if(x == 2) 3
        else if(x == 4) 4
        else 1
    }
    # For every row, we do the "special product"
    res <- apply(M, 1, function(x) {
        ret <- 0
        for(i in 1:(length(x)-1)) {
            if(length(x) != 1) {
                ret <- ret + x[i]*x[i+1]*alpha[f(x[i]*x[i+1]), i]
            }
        }
        ret
    })
    
    as.matrix(res)
}

#' Build a random Snps matrix
#'
#' @title Build a random Snps matrix
#' @param n sample size
#' @param s number of SNPS
#' @param snps_value value that SnPs can take
#' @return a random SnPs matrix
#' @author Julien Duvanel
#' @export
build_SNPs_matrix <- function(n, s, snps_value = c(0,1,2)) {
    
    # build a completly random matrix
    M <- matrix(sample(snps_value, size = n * s, replace = T), 
                nrow = n)
    
}
