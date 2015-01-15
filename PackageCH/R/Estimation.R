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
####  Data simulation
####
####**********************************************************************

#' Do an estimation using the method described 
#' in the function get_estimate (which is an argument here)
#'
#' @title Do estimation
#' @param get_estimate This argument has to be a function and must return estimate's in $h
#' @param ... arguments if needed
#' @return estimate h
#' @author Julien Duvanel
#' @export
DoEstimation <- function(get_estimate, V, phi) {
  
    # Get the beta's
    tryCatch({

        res <- get_estimate(V, phi)
        
        # We need to store the auc of all simulation
        # to be able to find the quantile (0.1 and 0.9 for example)
        # For each simulation, get the ROC Curves' points (as well as the MSPE's points)
        list(h = res$h)
      
      }, error = function(e) print(e))
  
}

#' Do the simulations for each element of param.list
#' 
#' @title doing estimations in batch
#' 
#' @param param.list List of parameters (containing n,p,g, etc.)
#' @param rfunction from which distribution do we want to simulate data ?
#' @param ... arguments used by the GenerateDatasets() function.
#' @return Nothing but export pdf in results/time_stamp[...].pdf
#' @author Julien Duvanel
#' @export
DoBatchEstimation <- function(param.list = NULL, ...) {
  
  # We need to have parameters, if not we stop.
  if(is.null(param.list) || length(param.list) < 1) {
    stop("A list of parameters is needed.")
  }
  
  # These informations are used to track the processus
  cat("Date time : ", date(), "\n\n")
  ptm <- proc.time()
  
  print(param.list)
  
  # For each parameters in param.list
  for(k in 1:length(param.list)) {
    
    cat("========================================================\n")
    cat("== Running parameters ", k, "/", length(param.list), "\n")
    cat("========================================================\n")
        
    full.res <- list()
    
    for(j in 1:length(param.list[[k]]$method.to.test)) {
      
      cat("== starting with method ", param.list[[k]]$method.to.test[j], "\n")
      
      # Get the parameters
      name <- param.list[[k]]$method.to.test[j]
      phenotypes.id <- param.list[[k]]$phenotypes.id
      
      # Timing (we wanna know how long it takes)
      datetime.stamp <- format(Sys.time(), "%d%m%Y_%H%M%S")
      ptm_precise <- proc.time()

      # The DoEstimation function needs a get_estimate function as an argument.
      # This function must return $h which contains the matrix of the heritability estimate 
      
      # One of the functions needs this...
      pdf(file = "/dev/null")
      plot.new()
      #########################
      ## Processing phenotypes given in param.list[[k]]$phenotypes.id
      #########################
      cat("\n")
      
      data.res <- c()
      
      for(l in 1:length(phenotypes.id)) {
        
        cat("   Processing ", phenotypes.id[l])
        
          method_data <- match.fun(paste0("load_data_", param.list[[k]]$method.to.test[j]))
        
          # Filter data before estimating
          P <- filter_P_by_id(P = P.raw, phenotypes.id[l])
          
          # This action is rquired so that their dimension match
          data.filtered <- filter_data(P, method_data()$phi)
          
          V <- build_matrix_V(data.filtered$P.filtered)
          colnames(V) <- colnames(P.raw)[phenotypes.id[l]]
        
          # Since we have to filter data first, we expect V to have dimension 1
          expect_that(ncol(V), equals(1))
        
          res <- c()
          for(m in 1:length(param.list[[k]]$model.to.test)) {

              est <- DoEstimation(get_estimate=match.fun(paste0("get_estimate_", 
                                                                param.list[[k]]$method.to.test[j],
                                                                param.list[[k]]$model.to.test[m])),
                                  V = V, 
                                  phi = data.filtered$phi.matrix.filtered)$h
              
              res <- cbind(res, est)
         
          }
        
        data.res <- rbind(data.res, cbind(colnames(P)[2], res))
        colnames(data.res) <- c("Phenotypes", param.list[[k]]$model.to.test)
        
        cat("   -> processed !\n")
                
      }
      
      full.res[[param.list[[k]]$method.to.test[j]]] <- data.res
  
      ###############################
      # Ending j loop (method.to.test)
      ###############################
      #         final_precise <- proc.time() - ptm_precise
      #         cat("Time for ", k, "/", length(param.list), " : ", final_precise, "\n\n")
      cat("\n== ending with method ", param.list[[k]]$method.to.test[j], "\n")
      
    }
    
    ###############################
    # Save datasets
    ###############################
    # We have to save the datasets in a file
    cat("== Saving results for method ", k, "/", length(param.list))
    
    # data
    save(full.res, file=paste0("results/data/", 
                               datetime.stamp, 
                               "_results_", name, ".RData"))
    
    cat("   > saved ! \n")
    
    # For each param.list entry, we wanna plot the graph
    dev.off()
    
    # Ending k loop
    cat("========================================================\n")
    cat("== End of Running with parameters ", k, "/", length(param.list), "\n")
    cat("========================================================\n\n")
    
  }
  
  final.time <- proc.time() - ptm
  cat("Total time : ", final.time, "\n")
  cat("Date time : ", date(), "\n")
  
}