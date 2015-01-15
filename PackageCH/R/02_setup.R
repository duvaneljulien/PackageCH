#' setup the project
#' 
#' create required directories, and return constants
#' 
#' @return list
#' @export
SetupProject <- function() {
  
  options(digits.secs=10)
  GetProjectConstants()

}

# This is used to avoid having "NOTE" when R CMD check'ing the project
utils::globalVariables(c("Project"))