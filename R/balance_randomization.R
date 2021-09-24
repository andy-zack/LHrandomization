#' Randomize a dataset, checking for balance with Lawley-Hotelling statistic
#'
#' This function is a wrapper around the randomizr functions. It runs a
#' randomization several times, selecting the one that is most balanced based
#' on the specified balance variables. The best balance is measured using the
#' Lawley-Hotelling statistic, which was added with the help of Brian Osserman.
#'
#' The L-W statistic is a measure of balance which is independent of affine
#' transformations of covariates (in particular, it does not depend on units).
#' However, it cannot distinguish more important covariates from less important
#' ones, so it is advisable not to include unnecessary covariates among the
#' balance variables.
#'
#' The L-W statistic has some strange behavior when some covariates are highly
#' correlated, and in the extreme of 100% correlation it will fail due to being
#' unable to invert a matrix.
#'
#' @param data a dataframe or tibble to be randomized
#' @param block_vars a vector of quoted column names for blocking variables
#' @param cluster_vars an optional vector of quoted column names to cluster on
#' @param balance_vars a vector of quoted column names for checking balance
#' @param n_treatment the number of experimental groups to randomize into
#' @param attempts the number of randomizations to conduct. the most balanced of these randomizations will be returned
#' @param condition_names a vector of quoted names of the experimental groups
#' @param treat_probs a vector of proportion to be placed into each group. If only two experimental groups, you may specifcy the proportion for one group and the remaining will be in the other group.
#' @param print_status boolean, will print status if set to TRUE
#' @param ...

#' @importFrom magrittr "%>%"
#' @importFrom dplyr select any_of mutate count
#' @importFrom tidyr unite
#' @importFrom randomizr block_ra block_and_cluster_ra
#'
#' @return a vector of treatment assignments
#' @export
#'
#' @examples
#' iris$assignments <- balance_randomization(data = iris,
#' block_vars = "Species",
#' balance_vars = c("Sepal.Length", "Sepal.Width"))
#'
#' #'\dontrun{library(nycflights13)
#' randomized_data <- flights[complete.cases(flights), ] %>%
#' filter(carrier %in% c("AA", "DL")) %>%
#' head(10000) %>%
#' balance_randomization(n_treatment = 3,
#'                         block_vars = c("carrier"),
#'                         cluster_vars = c("carrier", "tailnum"),
#'                         balance_vars = c("air_time", "distance"),
#'                         treat_probs = c(.2, .2, .6))}

balance_randomization <- function(data,
                                  block_vars,
                                  cluster_vars = NULL,
                                  balance_vars = NULL,
                                  n_treatment=2,
                                  attempts = 10,
                                  condition_names = NULL,
                                  treat_probs = .5,
                                  print_status = FALSE,
                                  ...) {


  # Stop Conditions ---------------------------------------------------------
  # If user declares condition names, it's length should equal n_treatment
  if(! length(condition_names) %in% c(0, n_treatment) ) {
    stop("If you specify condition_names, its length must equal n_treatment. n_treatment defaults to 2")
  }

  # N treatment should equal length of treat_probs unless n_treatment = 2 then you can specificy 1 treat prob
  if(n_treatment != 2 & length(treat_probs) != n_treatment) {
    stop("treat_probs must be a vector with length equal to n_treatment. The only exception is that if n_treatment = 2, then treat_probs can be one number, and the second prob is 1 - treat_probs.")
  }

  # make sure treat_probs add to one
  if(length(treat_probs) > 1 & sum(treat_probs) != 1) {
    stop("The probabilities in treat_prob must add up to 1.")
  }

  # Prepare for Randomization -----------------------------------------------
  # Create condition names if user did not specify
  if (n_treatment == 2 & is.null(condition_names)) {
    condition_names <- c("control", "treatment")
  } else if (is.null(condition_names)) {
    condition_names <- c("control", paste0("T", 2:n_treatment))
  }

  # Deal with situation where user specifices 1 treat_prob and n_treatment=2
  if(n_treatment == 2 & length(treat_probs) == 1) {
    treat_probs <- c(1 - treat_probs, treat_probs)
  }

  # make skinnier DF with just data that is necessary
  randomization_df <- data %>%
    select(any_of(balance_vars),
           any_of(cluster_vars),
           any_of(block_vars)) %>%
    unite("blocking_var",
          any_of(block_vars),
          remove = FALSE)

  if (!is.null(cluster_vars)){
    randomization_df <- unite(randomization_df,
                              "cluster_var",
                              any_of(cluster_vars), remove = FALSE)
  }

  # create formula to check balance
  model_formula <- paste("assignment ~ ", paste(balance_vars, collapse = " + ")) %>%
    as.formula()

  #initialize variables to store lowest Lawley-Hotelling and randomization
  lowest_LH <- -1
  best_randomization_df <- NULL


  # Run Randomizations, Checking for Balance --------------------------------
  # Do this part n times and pick version with lowest Lawley-Hotelling
  for (i in 1:attempts){

    if (n_treatment == 2) {

      # run randomization
      if (is.null(cluster_vars)) {
        this_randomization_df <- randomization_df %>%
          mutate(assignment = block_ra(randomization_df$blocking_var,
                                       prob =  treat_probs[2],
                                       conditions = condition_names))
      } else { # if we're doing cluster randomization
        this_randomization_df <- randomization_df %>%
          mutate(assignment = block_and_cluster_ra(randomization_df$blocking_var,
                                                   prob = treat_probs[2],
                                                   conditions = condition_names,
                                                   clusters = randomization_df$cluster_var))
      }

    } else { # if there are more than 2 treatment groups
      if (is.null(cluster_vars)) {
        this_randomization_df <- randomization_df %>%
          mutate(assignment = block_ra(randomization_df$blocking_var,
                                       num_arms = n_treatment,
                                       prob_each = treat_probs,
                                       conditions = condition_names))
      } else {
        this_randomization_df <- randomization_df %>%
          mutate(assignment = block_and_cluster_ra(randomization_df$blocking_var,
                                                   num_arms = n_treatment,
                                                   prob_each = treat_probs,
                                                   conditions = condition_names,
                                                   clusters = randomization_df$cluster_var))
      }

    }


    # Compute Lawley-Hotelling statistic, a multigroup-generalization of
    # Mahalanobis distance

    # wmat is the matrix of means for each treatment group and covariate
    wmatbig <- aggregate(as.formula(paste("cbind(", paste(balance_vars, collapse = ", "), ") ~ assignment", sep = "")), data = this_randomization_df, mean)
    #wmat <- select(wmatbig, any_of(balance_vars))
    wmat <- wmatbig[,balance_vars]
    group_names <- wmatbig$assignment
    rownames(wmat) <- group_names

    # v is the vector of means (of the balance variables) for each treatment group
    v <- this_randomization_df %>%
      select(any_of(balance_vars)) %>%
      colMeans(na.rm = TRUE)

    # nums is the number of observations in each treatment group
    nums = count(this_randomization_df, assignment)

    # B is the "between groups sum of squares" matrix
    B <- matrix(0L, nrow = length(balance_vars), ncol = length(balance_vars))
    for (j in 1:n_treatment) {
      w <- as.numeric(wmat[j,])
      n <- as.numeric(nums[j,2])
      B <- B + n * (w - v) %*% t(w - v)
    }

    # E is the "error sum of squares" (or "within groups sum of squares") matrix
    E <- matrix(0L, nrow = length(balance_vars), ncol = length(balance_vars))
    for (j in 1:nrow(this_randomization_df)) {
      u <- as.numeric(this_randomization_df[j,balance_vars])
      w <- as.numeric(wmat[this_randomization_df$assignment[j],])
      E <- E + (u - w) %*% t(u - w)
    }

    # The Lawley-Hotelling statistic is obtained as the trace of the matrix E^{-1} B
    this_LH <- sum(diag(solve(E) %*% B))

    # show means and Lawley-Hotelling statistic for each randomization
    if(print_status) {
      #print(aggregate(as.formula(paste("cbind(", paste(balance_vars, collapse = ", "), ") ~ assignment", sep = "")), data = this_randomization_df, mean))
      print(paste0("Randomization ", i, " of ", attempts, ". L-H stat: ", round(this_LH, 6)))
    }

    # if this randomization is better (as measured by lowest Lawley-Hotelling), overwrite the last version
    if ((this_LH < lowest_LH) || (lowest_LH == -1)){
      lowest_LH <- this_LH
      best_randomization_df <- this_randomization_df
    }
  }

  return(best_randomization_df$assignment)

}
