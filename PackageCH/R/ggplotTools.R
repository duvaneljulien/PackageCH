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
####  ggplot2 Tools
####
####**********************************************************************

#' Get legend from a ggplot's object
#' 
#' @title Get ggplot's legend
#' @param a.gplot a ggplot object
#' @return legend of a ggplot's object
#' @author Julien Duvanel
GetLegendFromGgplot2 <- function(a.gplot) {
  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  
  return(legend)
  
}

#' Get a custom ggplot theme
#' 
#' @title Custom ggplot theme
#' @return ggplot's theme
#' @author Julien Duvanel
GetCustomGgplotTheme <- function() {
  #theme_bw() + 
  theme(title = element_text(size = rel(1.5)),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16, face = "bold"),
        plot.margin = grid::unit(c(1,1,0,0), "cm"),
        legend.justification = c(1,0), 
        #legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_rect(fill = "transparent"),
        legend.key.width = grid::unit(0.5, "line"),
        legend.key.height = grid::unit(1.8,"line"),
        legend.key.size = grid::unit(3, "line"))
}

#' This function gives the final color (including transparency)
#'
#' @title Alpha color
#' @param colour the colour we wanna get transparent
#' @param alpha the degree of transparency
#' @return a new colour
#' @author Julien Duvanel
alpha <- function (colour, alpha = NA) {
  
  col <- col2rgb(colour, TRUE)/255
  
  if (length(colour) != length(alpha)) {
    if (length(colour) > 1 && length(alpha) > 1) {
      stop("Only one of colour and alpha can be vectorised")
    }
    if (length(colour) > 1) {
      alpha <- rep(alpha, length.out = length(colour))
    }
    else if (length(alpha) > 1) {
      col <- col[, rep(1, length(alpha)), drop = FALSE]
    }
  }
  
  alpha[is.na(alpha)] <- col[4, ][is.na(alpha)]
  new_col <- rgb(col[1, ], col[2, ], col[3, ], alpha)
  new_col[is.na(colour)] <- NA
  new_col
  
}