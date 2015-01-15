library(energy)
library(matrixStats)
library(reshape2)
library(plyr)

dcor.perso <- function (x, y, index = 1) 
{
    #     if (!(class(x) == "dist")) 
    #         x <- dist(x)
    #     if (!(class(y) == "dist")) 
    #         y <- dist(y)
    #     x <- as.matrix(x)
    #     y <- as.matrix(y)
    n <- nrow(x)
    m <- nrow(y)
    if (n != m) 
        stop("Sample sizes must agree")
    if (!(all(is.finite(c(x, y))))) 
        stop("Data contains missing or infinite values")
    if (index < 0 || index > 2) {
        warning("index must be in [0,2), using default index=1")
        index = 1
    }
    stat <- 0
    dims <- c(n, ncol(x), ncol(y))
    Akl <- function(x) {
        d <- as.matrix(x)^index
        rowM <- rowMeans(d)
        colM <- colMeans(d)
        M <- mean(d)
        
        a <- sweep(d, 1, rowM)
        b <- sweep(a, 2, colM)
        
        return(b + M)
    }
    A <- Akl(x)
    B <- Akl(y)
    dCov <- (mean(A * B))
    dVarX <- (mean(A * A))
    dVarY <- (mean(B * B))
    V <- sqrt(dVarX * dVarY)
    if (V > 0) 
        dCor <- dCov/V
    else dCor <- 0
    return(list(dCov = dCov, dCor = dCor, dVarX = dVarX, dVarY = dVarY))
}

#' Compare DCOR with different kind of settings
#'
#' @title Compare DCOR and LM
#' @param n sample size
#' @param s number of SNPS
#' @param u number of SnPs that really explain the phenotype
#' @param b X = G * runif(, -b, b)
#' @param noise noise
#' @param delta_add importance of additive effect
#' @param delta_dom importance of dominant effect
#' @param delta_epi importance of epistatic effect
#' @param snps_values the value randomly chosen for the G matrix
#' @param variable the variable that vary
#' @return export a pdf file
#' @author Julien Duvanel
#' @export
compare_dcor <- function(n, 
                         s, 
                         u, 
                         b,
                         noise.sd,
                         delta_add,
                         delta_dom,
                         delta_epi,
                         snps_value,
                         variable = "") {
    
    # We have to check that dimensions agree
    if (length(n) != length(s) | 
        length(n) != length(u) |
        length(n) != nrow(snps_value)) stop("Problem, length of n has to be the same as s, delta_add/dom/epi and snps_value.")    
    
    # We loop through length(n) (but they all have same length at this point)
    res <- c()
    for(i in 1:length(n)) {
        
        # Build two fake genome matrices
        # the second one is only used to compare with random results
        M <- build_SNPs_matrix(n[i], s[i], snps_value[i,])
        M_tilde <- build_SNPs_matrix(n[i], s[i])
        
        # Prepare all alpha's (additive/dominant/epistatic effect)
        alpha_add <- build_alpha(u = u[i], 
                                 s = s[i], 
                                 b = b)
        
        alpha_dom <- build_alpha_dominant(u = u[i], 
                                          s = s[i], 
                                          b = b)
        
        alpha_epi <- build_alpha_epistatic(u = u[i], 
                                           s = s[i], 
                                           b = b)
        # Build fake trait X
        X <- delta_add[i] * product_snps_alpha(M, alpha_add) + 
             delta_dom[i] * product_snps_alpha_dominant(M, alpha_dom) +
             delta_epi[i] * product_snps_alpha_epistatic(M, alpha_epi)
        
        # We are more interested to have X + noise 
        # noise is 10% of sd(X)
        noise <- rnorm(n[i], mean = 0, sd = noise.sd[i])   
        X.noise <- X + noise
        
        # Compute in advance distance matrices for
        # linear regression estimate of heritability
        # Actually, we use GRM as stated in one of the papers
        M_W <- compute_W(G = M)
        A_M <- M_W %*% t(M_W) / ncol(M_W)
        A_M.tri <- A_M[lower.tri(A_M)]
        
        M_tilde_W <- compute_W(G = M_tilde)
        A_M_tilde <- M_tilde_W %*% t(M_tilde_W) / ncol(M_tilde_W)
        A_M_tilde.tri <- A_M_tilde[lower.tri(A_M_tilde)]
        
        Z <- X.noise %*% t(X.noise)
        Z.tri <- Z[lower.tri(Z)]
        
        dist_M <- dist(M, p = 1)
        dist_M.tri <- as.matrix(dist_M)[lower.tri(as.matrix(dist_M))]
        dist_X <- dist(X.noise, p = 1)
        dist_X.tri <- as.matrix(dist_X)[lower.tri(as.matrix(dist_X))]
        
        dcor_XX_GRM <- dcor.perso(Z, A_M)
        dcor_distX_GRM <- dcor.perso(as.matrix(1-dist_X), A_M)
        
        cat(dcor_XX_GRM$dVarX, "\n")
        cat(dcor_distX_GRM$dVarX, "\n")
        cat(var(X.noise), "\n")
        
        # Do estimates
        res <- rbind(res, 
                     c(get(variable[1])[i], 
                       mean(abs(dist(X.noise, p = 1))) / (sqrt(n[i]-1)*dcov(X.noise,X.noise)),

                       # Correlation estimate
                       dcor_XX_GRM$dCor,
                       dcor_distX_GRM$dCor,
                       cor(as.vector(Z.tri), as.vector(A_M.tri)),

                       # Heritability estimate
                       (dcor_XX_GRM$dCor * dcor_XX_GRM$dVarX^0.5) / (dcor_XX_GRM$dVarY^0.5 * var(X.noise)),                       
                       dcor_distX_GRM$dCov / (dcor_distX_GRM$dVarY^0.5 * var(X.noise)),                       
                       lm( Z.tri ~ A_M.tri )$coefficients[2] / var(X.noise),
                                              
                       var(delta_add[i] * product_snps_alpha(M, alpha_add)) / var(X.noise),
                       var(X) / var(X.noise)
                     ))
        cat("-> i = ", i, "/", length(n), "\n")
    }
    
    # Gather data into a dataframe
    res <- data.frame(var = res[,1],
                      lim = res[,2],
                      
                      dcor_XX_GRM = res[,3],
                      dcor_distX_GRM = res[,4],
                      cor_XX_GRM = res[,5],
                      
                      dcor_H2 = res[,6],
                      dcor_distX_H2 = res[,7],
                      lm_H2 = res[,8],
                      
                      h2 = res[,9],
                      H2 = res[,10])
    
    # melt data to be able plot group into ggplots
    data.melt <- melt(res, measure.vars = c(
                                            "dcor_XX_GRM",
                                            "dcor_distX_GRM",
                                            "cor_XX_GRM",
                                            
                                            "dcor_H2",
                                            "dcor_distX_H2",
                                            "lm_H2",
                                            
                                            "h2",
                                            "H2"))
    
    data.melt$var <- jitter(data.melt$var, factor = 0.1)
    data.melt$value <- abs(data.melt$value)

    label_graph <- c(
                     expression(dcor(XX^T, GRM)),
                     expression(dcor(distX, GRM)),
                     expression(cor(XX^T, GRM)),
                     
                     expression(frac(dcor(XX^T, GRM), var(X)) %*% sqrt(frac(dVar(XX^T),dVar(GRM)))),
                     expression(frac(dcor(distX, GRM), var(X)) %*% sqrt(frac(dVar(distX),dVar(GRM)))),
                     expression(frac(lm(XX^T, GRM), var(X))),
                     
                     expression(paste("true ", h^2)),
                     expression(paste("true ", H^2)))
    

    # Export a pdf file
    p <- ggplot(data = data.melt,
                aes_string(x = "var" , y = "value")) +
            geom_point(aes_string(color = "variable", shape = "variable"), 
                       size = 3) +
            geom_line(aes_string(color = "variable"),
                      alpha = 0.6, size = 0.8) +
            geom_line(data = res, aes(x = var, y = lim)) +
            xlab(paste0("Value of ", paste(variable, collapse=", "))) +
            ylab("Heritability estimate") + 
            scale_shape_manual(values = 0:length(unique(data.melt$variable)),
                               labels=label_graph) +
#             scale_y_continuous(limits = c(-0.1, 1)) +
            scale_colour_discrete(labels=label_graph) + 
            theme(legend.text.align = 0) +  
            GetCustomGgplotTheme()
            
        
    ggsave(plot = p,
           filename = paste0("results/plots/dcor_vs_lm_", 
                             "n", min(n), "-", max(n), "_",
                             "s", min(s), "-", max(s), "_",
                             expression(delta_add), min(delta_add), "-", max(delta_add), "_",
                             "delta_dom", min(delta_dom), "-", max(delta_dom), "_",
                             "delta_epi", min(delta_epi), "-", max(delta_epi), "_",
                             "u", min(u), "-", max(u), "_",
                             "noise.sd", min(noise.sd), "-", max(noise.sd), "_",                             
                             format(Sys.time(), "%d%m%Y_%H%M%S"),
                             ".pdf"),
           width = 11, height = 7)
    
    save(list = c("res",
                  "n",
                  "s",
                  "delta_add",
                  "delta_dom",
                  "delta_epi",
                  "u",
                  "noise.sd"), 
         file = paste0("results/data/dcor_vs_lm_", 
                        as.character(substitute(get_alpha)), "_",
                        format(Sys.time(), "%d%m%Y_%H%M%S"),
                        ".RData"))
    
    # return
    list(res = res, X = X, M = M)    
    
}

#' Compare DCOR with different kind of settings and multi-traits
#'
#' @title Compare DCOR and LM
#' @param n sample size
#' @param s number of SNPS
#' @param t the number of traits
#' @param m mode (see build_alpha_multi)
#' @param snps_values the value randomly chosen for the G matrix
#' @param variable the variable that vary
#' @return export a pdf file
#' @author Julien Duvanel
#' @export
compare_dcor_multi <- function(n, 
                               s, 
                               t,
                               snps_value,
                               variable = "") {
    
    # We have to check that dimensions agree
    if (length(n) != length(s) | 
            length(n) != length(t) |
            length(n) != nrow(snps_value)) stop("Problem, length of n has to be the same as s, delta_add/dom/epi and snps_value.")    
    
    # We loop through length(n) (but they all have same length at this point)
    res <- data.frame()
    res_tilde <- data.frame()
    for(i in 1:length(n)) {
        
        # Build two fake genome matrices
        # the second one is only used to compare with random results
        M <- build_SNPs_matrix(n[i], s[i], snps_value[i,])
        M_tilde <- build_SNPs_matrix(n[i], s[i])
        
        X <- list()
        
        # Prepare all multi alpha
        for(j in 1:5) {
            X[[j]] <- build_X_multi(M = M,
                                    build_alpha_multi(s = s[i],
                                                      t = t[i],
                                                      mode = j))
        }

        # Compute in advance distance matrices for
        # linear regression estimate of heritability
        # dist_M <- 1 - as.vector(as.matrix(dist(M)))
        # dist_M_tilde <- 1 - as.vector(as.matrix(dist(M_tilde)))
                
        # Do estimates
        res_est <- vector(mode = "list", length = 5)
        res_est_tilde <- vector(mode = "list", length = 5)
        for(l in 1:5) {
            for(j in 1:ncol(X[[l]])) {
                if(j %% 2 == 0) {
                    res_est[[l]] <- c(res_est[[l]], dcor(X[[l]][,1:j], M))
                    res_est_tilde[[l]] <- c(res_est_tilde[[l]], dcor(X[[l]][,1:j], M_tilde))                
                }
            }        
        }

        for(l in 1:5) {
            res <- rbind.fill(res, 
                              data.frame(var = get(variable[1])[i] + (l-1)*0.02*max(s), 
                                         est = t(res_est[[l]])))
            res_tilde <- rbind.fill(res_tilde, 
                                    data.frame(var = get(variable[1])[i] + (l-1)*0.02*max(s), 
                                               est = t(res_est_tilde[[l]])))
        }

        
        cat("-> i = ", i, "/", length(n), "\n")
    }
    
    # Gather data into a dataframe
    res <- data.frame(var = res[,1],
                      est = res[,-1])
    res_tilde <- data.frame(var = res_tilde[,1],
                            est = res_tilde[,-1])
    
    # melt data to be able plot group into ggplots
    data.melt <- melt(res, measure.vars = names(res)[-1])
    data.melt.tilde <- melt(res_tilde, measure.vars = names(res_tilde)[-1])
    
    
    # Export a pdf file
    p <- ggplot(data = data.melt,
                aes_string(x = "var" , y = "value")) +
        geom_point(aes_string(colour = "variable"), 
                   size = 3, 
                   position = position_jitter(w = 0.005 * sd(data.melt$value), 
                                              h = 0)) +
        xlab(paste0("Value of ", paste(variable, collapse=", "))) +
        ylab("Heritability estimate") + 
        scale_y_continuous(limits = c(0, 1)) +     
        GetCustomGgplotTheme() +
        theme(legend.position="none")
    
    p_tilde <- ggplot(data = data.melt.tilde,
                aes_string(x = "var" , y = "value")) +
        geom_point(aes_string(color = "variable"), size = 3, position = position_jitter(w = 0.005 * sd(data.melt.tilde$value), h = 0)) +
        xlab(paste0("Value of ", paste(variable, collapse=", "))) +
        ylab("Heritability estimate w/ rand. gen. matrix") + 
        scale_y_continuous(limits = c(0, 1)) +          
        GetCustomGgplotTheme() +
        theme(legend.position="none")
    
    # Create a grid with 2 rows, length(p) columns
    # Plot results, 3 columns (= 3 methods)
    pdf(file = paste0("results/plots/dcor_multi_", 
                      "n", min(n), "-", max(n), "_",
                      "s", min(s), "-", max(s), "_",
                      "t", min(t), "-", max(t), "_",
                      format(Sys.time(), "%d%m%Y_%H%M%S"),
                      ".pdf"), 
        width = 17, 
        height = 7)
    
        grid.newpage()
            pushViewport(viewport(layout = grid.layout(1, 2)))   
    
            print(p, vp = vplayout(1,1))
            print(p_tilde, vp = vplayout(1,2))        
    
        upViewport(0)
    
    dev.off()
    
    save(list = c("data.melt", "data.melt.tilde"), file = paste0("results/data/dcor_multi_", 
                                     "n", min(n), "-", max(n), "_",
                                     "s", min(s), "-", max(s), "_",
                                     "t", min(t), "-", max(t), "_",
                                     format(Sys.time(), "%d%m%Y_%H%M%S"),
                                     ".RData"))
    
    # return
    list(res = res, X = X, M = M)    
    
}

