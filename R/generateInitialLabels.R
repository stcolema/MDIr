#' @title Generate initial labels
#' @description For V views, generate initial labels allowing for both
#' semi-supervised and unsupervised views.
#' @param labels $N x V$ matrix of initial labels. The actual values passed only
#' matter for semi-supervised views (i.e. the views for which some labels are
#' observed).
#' @param fixed $N x V$ matrix indicating which labels are observed and hence
#' fixed. If the vth column has no 1's then this view is unsupervised.
#' @param alpha The concentration parameter (vector).
#' @param K The number of components modelled in each  view.
#' @param labels_as_intended Logical indicating if ``labels`` is the intended
#' output.
#' @return An N vector of labels.
#' @export
generateInitialLabels <- function(labels, fixed, K, alpha,
                                  labels_as_intended = FALSE) {
  V <- ncol(fixed)
  N <- nrow(fixed)

  view_indices <- seq(1, V)

  is_semisupervised <- apply(fixed, 2, function(x) any(x == 1))

  not_a_matrix <- !is.matrix(labels)

  if (not_a_matrix) {
    stop("``initial_labels`` must be a matrix of integers.")
  }

  wrong_number_of_views <- ncol(labels) != V
  if (wrong_number_of_views) {
    stop("Number of columns in labels must equal number of datasets/views in ``X``")
  }

  wrong_number_of_samples <- nrow(labels) != N
  if (wrong_number_of_samples) {
    stop("Number of rows in labels must equal number in each dataset/view in ``X``")
  }

  for (v in view_indices) {
    labels_v <- labels[, v]
    fixed_v <- fixed[, v]
    K_v <- K[v]
    alpha_v <- alpha[v]

    if (labels_as_intended) {
      checkLabels(labels_v, K_v)
    } else {
      if (is_semisupervised[v]) {

        # Assign the unobserved labels to an arbitrary class, in this case 1.
        labels_v[which(fixed_v == 0)] <- 1

        # Only need to check semisueprvised labels closely as all of the
        # unsupervised labels are generated
        checkLabels(labels_v, K_v)
        labels_v <- generateInitialSemiSupervisedLabels(labels_v, fixed_v)

        # Check the labels are all contiguous, i.e. all values in the range 0 to K-1
        # are represented
        non_contiguous_labels <- max(labels_v) != (length(unique(labels_v)))

        if (non_contiguous_labels) {
          stop("initial labels are not all contiguous integers.")
        }
      } else {
        labels_v <- generateInitialUnsupervisedLabels(N, alpha_v, K_v)
      }
    }

    # Check that the initial labels starts at 0, if not remedy this.
    no_zero_label <- !any(labels_v == 0)

    if (no_zero_label) {
      labels_v <- as.numeric(as.factor(labels_v)) - 1
    }

    # Update the labels matrix with the vector for the current view
    labels[, v] <- labels_v
  }
  labels
}
