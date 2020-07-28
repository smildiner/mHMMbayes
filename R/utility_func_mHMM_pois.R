#' @keywords internal
is.mHMM_cont <- function(x) {
  inherits(x, "mHMM_pois")
}

#' Obtain the AIC for a fitted multilevel HMM with a Poisson conditional
#' distribution
#'
#' \code{get_aic_pois} obtains the Akaike information criterion (AIC) for a
#' fitted multilevel hidden Markov model with a Poisson state-dependent
#' conditional distribution, averaged across subjects.
#'
#' @param x An object of class \code{mHMM_pois}, generated
#'   by the function \code{\link{mHMM_pois}} or \code{\link{mHMM_pois_rm}},
#'   respectively.
#' @param burn_in An integer which specifies the number of iterations to discard
#'   when obtaining the model parameter summary statistics. When left
#'   unspecified (\code{burn_in = NULL}), the burn in period specified when
#'   creating the \code{mHMM_pois} or \code{mHMM_pois_rm} object will be used.
#'
#' @return \code{get_aic_pois} returns the AIC for the fitted model averaged
#' across subjects, based on the number of parameters used in the model and the
#' loglikelihood. When a repeated measures design has been used, the AIC is
#' calculated based on the loglikelihood summed over repeated measures of a
#' same subject, and the AIC averaged across subjects.
#'
#' @export

get_aic_pois <- function(x, burn_in = NULL){

  # Retrieve model parameters
  n_subj  <- x[["input"]][["n_subj"]]
  m       <- x[["input"]][["m"]]
  n_dep   <- x[["input"]][["n_dep"]]
  J       <- x[["input"]][["J"]]

  if(is.null(burn_in)) {
    burn_in <- x[["input"]][["burn_in"]]
  }

  # Extract MAP LL per subject
  LL      <- numeric(n_subj)
  for(i in 1:n_subj){
    LL[i] <- median(x$PD_subj[[i]][((burn_in + 1): J), ((n_dep * m) + m*m + 1)])
  }

  # Calculate and return AIC
  AIC<-2*(sum(n_dep*m)+(m-1)*m) - (2*LL)

  # Return average AIC over subjects
  return(mean(AIC))

}


#' Burn-in the MC chain of each state forward probabilities
#'
#' @keywords internal

burn_fw <- function(x, burn_in = NULL) {

  # Number of burn_in samples
  if(is.null(burn_in)) {
    burn_in <- x[["input"]][["burn_in"]]
  }
  J <- x$input$J
  burn_me <- x$forward_prob

  # Burn
  # Lowest level repeated measures or subjects?
  if(mode(burn_me[[1]][[1]]) == "list"){
    x <- lapply(burn_me,
                function(burn_day) lapply(burn_day,
                                          function(burn_trial) lapply(burn_trial,
                                                                      function(burn_fw, burn_in) burn_fw[,(burn_in+1):ncol(burn_fw)],
                                                                      burn_in)))
  } else {
    x <- lapply(burn_me,
                function(burn_day) lapply(burn_day,
                                          function(burn_fw, burn_in) burn_fw[,(burn_in+1):ncol(burn_fw)],
                                          burn_in))
  }

  # Return
  return(x)
}


#' Gather forward probabilities from a nested list into a matrix
#'
#' @keywords internal

gather_fw <- function(x, target){

  if(mode(x[[1]][[1]][[1]]) == "list") {
    do.call(rbind,
            lapply(seq_along(x),
                   function(n_subj) do.call(rbind,
                                            lapply(seq_along(x[[n_subj]]),
                                                   function(n_rm) {
                                                     out <- do.call(cbind, lapply(x[[n_subj]][[n_rm]], function(state) state[[target]]))
                                                     out <- cbind("subj" = n_subj, "rm" = n_rm, out)
                                                   }
                                            )
                   )
            )
    )
  } else {
    do.call(rbind,
            lapply(seq_along(x),
                   function(n_subj) {
                     out <- do.call(cbind, lapply(x[[n_subj]], function(state) state[[target]]))
                     out <- cbind("subj" = n_subj, out)
                   }
            )
    )
  }

}



#' Extract the maxima a posteriori (MAPs) for the state-specific forward
#' probabilities from a fitted multilevel HMM with a Poisson conditional
#' distribution
#'
#' \code{get_map_fw} extracts the \code{mean}, \code{median} and \code{sd} of
#' each state-specific forward probability for the full sequence length of all
#' the lowest level units in a fitted multilevel hidden Markov model with a
#' Poisson state-dependent conditional distribution. When
#' \code{\link{mHMM_pois}} has been used to train a model the lowest level
#' corresponds to subjects, whereas for \code{\link{mHMM_pois_rm}},
#' repeated measurements.
#'
#' @param x An object of class \code{mHMM_pois}, generated
#'   by the function \code{\link{mHMM_pois}} or \code{\link{mHMM_pois_rm}},
#'   respectively.
#' @param burn_in An integer which specifies the number of iterations to discard
#'   when obtaining the model parameter summary statistics. When left
#'   unspecified (\code{burn_in = NULL}), the burn in period specified when
#'   creating the \code{mHMM_pois} or \code{mHMM_pois_rm} object will be used.
#' @param target A character specifying the statistic to be extracted. It
#'   accepts \code{"mean"}, \code{"median"}, or \code{"se"} as inputs. If
#'   no entry is sepcified, it defaults to \code{"median"}.
#'
#' @return \code{get_map_fw} depending on the \code{target} chosen, the
#' \code{mean}, \code{median} or \code{sd} of each state-specific forward
#' probability for the full sequence length of all the lowest level units in
#' a fitted multilevel hidden Markov model with a Poisson state-dependent
#' conditional distribution is returned.
#'
#' @export

get_map_fw <- function(x, burn_in = NULL, target = "median") {

  # Extract parameters
  if(is.null(burn_in)) {
    burn_in <- x[["input"]][["burn_in"]]
  }
  J       <- x$input$J
  burn_me <- x$forward_prob

  # Initialize elements
  map_out <- vector("list", length(burn_me))

  # Burn
  burned  <- burn_fw(x, burn_in = burn_in)

  # Get MAPs on burned data
  for(n_subj in seq_along(burned)) {

    # Lowest level repeated measures or subjects?
    if(mode(burned[[1]][[1]]) == "list") {
      for(n_rm in seq_along(burned[[n_subj]])) {
          map_out[[n_subj]][[n_rm]] <- lapply(burned[[n_subj]][[n_rm]], function(x) {
            list(
              "mean" = unname(apply(x, 1, mean)),
              "median" = unname(apply(x, 1, median)),
              "se" = unname(apply(x, 1, sd))
            )}
          )
      }
    } else {
      map_out[[n_subj]] <- lapply(burned[[n_subj]], function(x) {
        list(
          "mean" = unname(apply(x, 1, mean)),
          "median" = unname(apply(x, 1, median)),
          "se" = unname(apply(x, 1, sd))
        )}
      )
    }
  }

  # Unlist into a data frame
  map_out_df <- gather_fw(map_out, target = target)

  # Return MAPs
  return(map_out_df)
}
