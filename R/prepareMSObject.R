#' @title Prepare MS Object
#' @description Prepares a mass spectrometry experiment (as stored in 
#' the Bioconductor package ``pRolocdata``) for modelling, extracting the 
#' numerical data, the classes and the indicator matrix for which labels are 
#' observed.
#' @param MS_object A mass spectrometry experiment such as ``tan2009r1`` from 
#' ``pRolocdata``.
#' @return A list of ``X``, the fracitonation data from a LOPIT experiment, 
#'  ``fixed``, the matrix indicating which labels are observed,  
#'  ``initial_labels``, the matrix of the initial labels to be input into 
#'  ``runMCMCChains`` (note that ``"unknown"`` organelles are represented 
#'  arbitrarily with a ``1`` as these will be sampled again within the wrapper 
#'  of ``callMDI``) and ``class_key`` which maps the numeric representation of 
#'  organelles back to the original naming.
#'  @importFrom Biobase exprs fData
#' @export
prepareMSObject <- function(MS_object) {
  
  # Extract the LOPIT data and the organelles
  X <- Biobase::exprs(MS_object)
  organelles <- Biobase::fData(MS_object)[, "markers"]
  
  # Create a data frame of the classes present and their associated number;\
  # this can be used to map the numeric representation of the classes back to
  # an organelle
  organelles_present <- pRoloc::getMarkerClasses(MS_object)
  class_key <- data.frame(Organelle = organelles_present, 
    Key = seq(1, length(organelles_present))
  )
  
  # Number of components modelled
  K <- length(organelles_present)
  
  # Number of views modelled
  V <- 1
  
  # Number of samples modelled
  N <- nrow(X)
    
  # Prepare initial labels
  initial_labels <- fixed <- matrix(0, nrow = N, ncol = V)
    
  # Fix training points, allow test points to move component
  fix_vec <- (organelles != "unknown") * 1
  fixed[, 1] <- fix_vec
  
  # Assign initial labels
  initial_labels[, 1] <- class_key$Key[match(organelles, class_key$Organelle)]
  
  # Any unknown labels are given an arbitrary value which will be reset in the 
  # model call function.
  initial_labels[is.na(initial_labels)] <- 1
  
  data_modelled <- list(
    X
  )
  
  # Return the prepared objects
  list(X = X, 
    fixed = fixed,
    initial_labels = initial_labels,
    class_key = class_key
  )
}
