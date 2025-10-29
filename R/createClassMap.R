#' @title Create Class Label Map
#' @description Create bidirectional mapping between class names and numeric 
#' representations, handling non-contiguous class presence across views
#' @param class_names Character vector of unique class names
#' @param start_index Integer. Starting index for numeric labels (default: 1)
#' @returns List with:
#' * `name_to_num`: Named vector mapping class names to numeric labels
#' * `num_to_name`: Named vector mapping numeric labels to class names  
#' * `n_classes`: Total number of classes
#' @export
#' @examples
#' class_map <- createClassMap(c("A", "B", "C", "D", "E"))
#' class_map$name_to_num["C"] # Returns 3
createClassMap <- function(class_names, start_index = 1) {
  
  if (!is.character(class_names)) {
    stop("class_names must be a character vector")
  }
  
  unique_classes <- unique(class_names)
  n_classes <- length(unique_classes)
  
  if (n_classes == 0) {
    stop("No classes provided")
  }
  
  numeric_labels <- seq(start_index, start_index + n_classes - 1)
  
  name_to_num <- setNames(numeric_labels, unique_classes)
  num_to_name <- setNames(unique_classes, numeric_labels)
  
  list(
    name_to_num = name_to_num,
    num_to_name = num_to_name,
    n_classes = n_classes
  )
}

#' @title Convert Class Labels to Numeric
#' @description Convert character class labels to numeric using a class map
#' @param labels Character vector of class labels
#' @param class_map Output from `createClassMap()`
#' @param allow_missing Logical. If TRUE, NA values are preserved; if FALSE, 
#' they trigger an error
#' @returns Integer vector of numeric labels
#' @export
#' @examples  
#' class_map <- createClassMap(c("A", "B", "C"))
#' labelsToNumeric(c("A", "C", "B"), class_map)
labelsToNumeric <- function(labels, class_map, allow_missing = FALSE) {
  
  if (!is.character(labels)) {
    stop("labels must be a character vector")
  }
  
  numeric_labels <- class_map$name_to_num[labels]
  
  if (!allow_missing && any(is.na(numeric_labels))) {
    unknown <- unique(labels[is.na(numeric_labels)])
    stop(paste(
      "Unknown class labels found:",
      paste(unknown, collapse = ", ")
    ))
  }
  
  as.integer(numeric_labels)
}

#' @title Convert Numeric Labels to Class Names
#' @description Convert numeric labels back to character class names
#' @param labels Integer vector of numeric labels
#' @param class_map Output from `createClassMap()`
#' @returns Character vector of class names
#' @export
numericToLabels <- function(labels, class_map) {
  
  if (!is.numeric(labels)) {
    stop("labels must be numeric")
  }
  
  char_labels <- class_map$num_to_name[as.character(labels)]
  
  if (any(is.na(char_labels))) {
    unknown <- unique(labels[is.na(char_labels)])
    warning(paste(
      "Unknown numeric labels found:",
      paste(unknown, collapse = ", ")
    ))
  }
  
  char_labels
}

#' @title Prepare Multi-View Labels for MDI
#' @description Convert class labels to contiguous numeric format for MDI,
#' handling non-contiguous class representations across views. This function
#' solves the problem where different views may have different subsets of  
#' classes (e.g., view 1 has classes 1-5, view 2 has only classes 1,3,5).
#' 
#' @param labels_list List of N x V matrices or list of N-vectors, one per view.
#' Each contains class labels (character or numeric)
#' @param class_map Optional. Output from `createClassMap()`. If NULL, created 
#' from unique labels across all views
#' @param fixed_list Optional list of N x V matrices indicating fixed labels. 
#' If provided, validates that observed classes match the global class set
#' @returns List with:
#' * `labels_matrix`: N x V matrix of contiguous integer labels (1:K)
#' * `class_map`: The class mapping used
#' * `original_K`: Vector of K values per view before remapping
#' * `global_K`: Total number of unique classes across all views
#' * `view_class_sets`: List showing which global classes appear in each view
#' @export
#' @examples
#' # View 1 has classes A-E, View 2 has only A, C, E
#' labels_v1 <- sample(c("A", "B", "C", "D", "E"), 100, replace = TRUE)
#' labels_v2 <- sample(c("A", "C", "E"), 100, replace = TRUE)
#' 
#' result <- prepareMDILabels(list(labels_v1, labels_v2))
#' # result$labels_matrix contains contiguous 1:5 for both views
#' # result$view_class_sets shows which classes appear in each view
prepareMDILabels <- function(labels_list,
                             class_map = NULL,
                             fixed_list = NULL) {
  
  # Input validation
  if (!is.list(labels_list)) {
    stop("labels_list must be a list")
  }
  
  V <- length(labels_list)
  N <- length(labels_list[[1]])
  
  # Validate dimensions
  for (v in seq_len(V)) {
    if (length(labels_list[[v]]) != N) {
      stop(paste("View", v, "has different number of samples"))
    }
  }
  
  # Create class map if not provided
  if (is.null(class_map)) {
    all_labels <- unlist(labels_list)
    unique_labels <- unique(all_labels[!is.na(all_labels)])
    class_map <- createClassMap(as.character(unique_labels))
  }
  
  global_K <- class_map$n_classes
  
  # Convert each view to numeric
  labels_matrix <- matrix(NA_integer_, nrow = N, ncol = V)
  view_class_sets <- vector("list", V)
  original_K <- integer(V)
  
  for (v in seq_len(V)) {
    view_labels <- as.character(labels_list[[v]])
    
    # Track which classes appear in this view
    view_classes <- unique(view_labels[!is.na(view_labels)])
    view_class_sets[[v]] <- sort(class_map$name_to_num[view_classes])
    original_K[v] <- length(view_classes)
    
    # Convert to numeric using global class map
    labels_matrix[, v] <- labelsToNumeric(
      view_labels,
      class_map,
      allow_missing = TRUE
    )
    
    # Validate with fixed labels if provided
    if (!is.null(fixed_list)) {
      fixed_v <- fixed_list[[v]]
      observed_idx <- which(fixed_v == 1)
      
      if (length(observed_idx) > 0) {
        observed_classes <- unique(labels_matrix[observed_idx, v])
        
        # Check all observed classes are valid
        invalid <- observed_classes[!(observed_classes %in% seq_len(global_K))]
        if (length(invalid) > 0) {
          stop(paste(
            "View", v, "has observed classes not in global class set:",
            paste(invalid, collapse = ", ")
          ))
        }
      }
    }
  }
  
  list(
    labels_matrix = labels_matrix,
    class_map = class_map,
    original_K = original_K,
    global_K = global_K,
    view_class_sets = view_class_sets
  )
}