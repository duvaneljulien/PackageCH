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
####  Diagnostic plots
####
####**********************************************************************

#' Plot the heatmap of X's correlation matrix
#'
#' @param X original matrix
#' @param correlation choosen groups
#' @return Plot the heatmap
#' @author Julien Duvanel
#' @export
PlotHeatmapCorrelationMatrix <- function(X, 
                                         correlation = FALSE) {
  
  if(correlation == TRUE) {
    # Get the correlation matrix of X
    z <- cor(X)
  } else {
    z <- X
  }

  # Get colors between white and black
  col <- colorRampPalette(c("white", "black"), space="rgb")
  
  breaks <- seq(min(X), max(X), length.out=10)
  
  image(z, col=col(length(breaks)-1), breaks=breaks, xaxt="n", yaxt="n", ylab="", xlab="", useRaster=T)

}
