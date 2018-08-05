#' Computing the best sc model for each electrode
#'
#' This function computes the best sc model for each electrode, following three
#' steps:
#' 1. start with the initial model at the best uniform lag;
#' 2. for each electrode in turn, evaluate the risks of the 3*13 models
#' that correspond to using the 3*13 SC equations for that electrode, and replace
#' 	the current equation with the best equation.
# 3. repeat until the model no longer changes
#' @param ztrain neural data in the training set
#' @param ztest neural data in the testing set
#' @param model starting model, one per electrode
#' @param n_units number of units in the dataset
#' @param n_transf number of sc transforamtion used, by default 3
#' @param n_lags number of lags used, by default 13
#' @param eqs_df data frame for equation reference
#' @param max_iter max iteration of the algorithm if it doesn't converge first
#' @param mle whether mle is true or false
#' @param save_name file name where saving the objects
#' @param random_order whether to iterate randomly over the electrodes
#' @param int_dec vector of integer on when the trial starts
#' @param time_dec vector of time of the intervals
#' @param vel velocity in the training set
#' @param veltest velocity in the testing set
#' @param is_sc whether the best model should run on spike counts, default equal TRUE
#' @param wf_mom wf moment to consider, only if is_sc is FALSE
#' @param wf_key wf key to consider, only if is_sc is FALSE
#' @param transf_selected transformation selected 1 linear, 2 ace, 3 sqrt, NULL all
#' @return the final model
get_best_model <- function(ztrain, ztest, model, n_units, n_transf = 3,
                            n_lags = 13, eqs_df, max_iter = 10, mle = FALSE,
                            save_name, random_order = TRUE, int_dec, time_dec,
                            vel, veltest, is_sc = TRUE, wf_mom = NULL,
                            wf_feat = NULL, transf_selected = NULL){


  # Set the initial tests
  if(is.null(save_name)) stop("file name has to be passed")

  if(eqs_df[model,]$key %>% unique %>% length != n_units)
    stop("The starting model needs to be one equation per electrode")

  if(is_sc) wf_mom <- wf_feat <- NULL

  if((!is_sc & is.null(wf_mom)) | (!is_sc & is.null(wf_mom)))
    stop("When sc is FALSE, wf_mom and wk_key need to be specified")


  # Initializing useful quantities
  forward_mse <- new_model <- final_model <-
    lapply(vector('list', max_iter), function(x) vector('list', n_units))
  all_risk <- vector('list', max_iter)

  iter <- 1
  final_model[[iter]][[1]] <- model

  # Compute current error with current model
  current_model <- get_current_model_info(z = ztrain[,model], ztest = ztest[,model],
                                        encoding=TRUE, decoding=TRUE, mle = mle,
                                        decoding_train = FALSE)
  old_risk <- current_model$mse %>% mean
  Qm1_tot <- NULL
  if(mle) Qm1_tot <- solve(current_model$Q)

  # Generate random order of the electordes
  electrodes <- 1:n_units
  if(random_order) electrodes <- sample(1:n_units, n_units, replace=FALSE)
  condition_valid <- TRUE

  while(condition_valid & (iter <= max_iter)){

    cat("***** iteration -", iter, "\n")

    for(ii in 1:n_units){
      t0 <- proc.time()
      cat("--- unit -", electrodes[ii], "\n")
      cat('old equation: ', model[electrodes[ii]], '\n')
      obj_by_electrode <- get_best_model_by_electrode(ii, electrodes, all_risk,
                              iter, old_risk, Qm1_tot, eqs_df, mle, model,
                              new_model, int_dec, time_dec, forward_mse,
                              ztrain, ztest, final_model, is_sc, wf_mom, wf_feat,
                              transf_selected)

      cat('time = ', (proc.time() - t0)[[3]], "\n")


      # Assign object of a  list to single variable to update them in next iteration
      all_risk  <- obj_by_electrode$all_risk
      old_risk  <- obj_by_electrode$old_risk
      Qm1_tot   <- obj_by_electrode$Qm1_tot
      model     <- obj_by_electrode$model
      new_model <- obj_by_electrode$new_model
      forward_mse <- obj_by_electrode$forward_mse
      final_model <- obj_by_electrode$final_model

      cat('new equation: ', model[electrodes[ii]], '\n')
      # Run this check only if not spike counts
      if(!is_sc){

        if((eqs_df[model,]$wf %>% unique != wf_mom) |
           (eqs_df[model,]$wf_key %>% unique != wf_feat))
          stop("Model is proposing wrong wf key or feature")

        cat('wf: ', unique(eqs_df[model, ]$wf), ', wf-feat: ',
            unique(eqs_df[model, ]$wf_key), '\n')
      } else{
        if(eqs_df[model,]$wf %>% unique != 0)
          stop("This should be a spike count model")
      }


      # Save the iteration
      save(model, new_model, forward_mse, final_model, all_risk, electrodes,
           file = save_name)

    }

    # Update quantities for while
    if(length(unique(all_risk[[iter]])) == 1) condition_valid <- FALSE
    iter <- iter + 1
    final_model[[iter]][[1]] <- model
  }

  return(model)
}

#' Computing the best sc model for a single electrodee
#'
#' This function computes the best sc model for a single electrode, trying to
#' substute the equation for a specificic electrode with all the 13*3 equations
#' that belong to the same electrode
#' @param ii iterator for the electrode vector
#' @param electrodes vector of the electrodes randomly ordered
#' @param all_risk vector with all risks computed so far
#' @param iter iteration
#' @param old_risk risk at the previous stage, used to compare new equations
#' @param Qm1_tot inverse of full model, if mle
#' @param eqs_df data frame with equation info
#' @param mle whether mle or kf
#' @param model current model
#' @param new_model list with history of all models tried
#' @param int_dec vector of integer on when the trial starts
#' @param time_dec vector of time of the intervals
#' @param forward_mse is history with all past forward mse values
#' @param ztrain neural data in the training set
#' @param ztest neural data in the testing set
#' @param final_model list with history of all successfull models iter by iter
#' @param is_sc whether the best model should run on spike counts, default = TRUE
#' @param wf_mom wf moment to consider, only if is_sc is FALSE
#' @param wf_key wf key to consider, only if is_sc is FALSE
#' @param transf_selected transformation selected 1 linear, 2 ace, 3 sqrt, NULL all
#' @return a list of vectors update and needed for next iteration
get_best_model_by_electrode <- function(ii, electrodes, all_risk, iter,
                                           old_risk, Qm1_tot, eqs_df, mle, model,
                                           new_model, int_dec, time_dec,
                                           forward_mse, ztrain, ztest, final_model,
                                           is_sc, wf_mom, wf_feat, transf_selected){

  all_risk[[iter]] <- c(all_risk[[iter]], old_risk)
  unit <- electrodes[ii]

  Qm1 <- NULL
  if(mle) Qm1 <- get_inverse_matrix_m1(Qm1_tot, unit)

  # Extract equations for the same electrode
  if(is_sc){
    candidate_rows <- (eqs_df %>%
                        subset(key == unit) %>%
                        subset(wf == 0))$eq
  } else {
    candidate_rows <- (eqs_df %>%
                        subset(key == unit) %>%
                        subset(wf == wf_mom) %>%
                        subset(wf_key == wf_feat))$eq
  }
  candidate_eqs <-  candidate_rows %>% setdiff(model[unit])

  # If only a specific transformation wants to be selected
  if(!is.null(transf_selected) & is_sc)
    candidate_eqs <- subset(eqs_df[candidate_rows,],
                            transf == transf_selected)$eq  %>% setdiff(model[unit])

  if(!is.null(transf_selected) & !is_sc)
    candidate_eqs <- subset(eqs_df[candidate_rows,],
                            wf_transf == transf_selected)$eq  %>% setdiff(model[unit])

  # Check the model is selecting the right marginal equations
  if(eqs_df[candidate_eqs, ]$wf %>% unique !=  eqs_df[model,]$wf %>% unique |
     eqs_df[candidate_eqs, ]$wf_key %>% unique != eqs_df[model,]$wf_key %>% unique)
      stop("Either type of waveform or waveform feature doesn't match ")



  # Delete current equation and compute length of current model and model to add
  new_model[[iter]][[ii]] <- c(model[-unit], candidate_eqs)
  current_model <- new_model[[iter]][[ii]]
  l <- length(model) - 1
  L <- length(current_model)
  ztrain_selected <- ztrain[,current_model]
  ztest_selected  <- ztest[,current_model]

  # Compute add-one-in score
  encoding_obj <- get_current_model_info(z = ztrain_selected, ztest = ztest_selected,
                                      encoding = TRUE, decoding = FALSE, mle = mle,
                                      decoding_train = FALSE)
  forward_mse[[iter]][[ii]] <- get_score(ztrain = ztrain_selected,
                                        ztest = ztest_selected,
                                        A = encoding_obj$A, W = encoding_obj$W,
                                        H = encoding_obj$H, Q = encoding_obj$Q,
                                        Qm1 = Qm1, method = "forward",
                                        mle = mle, full_size = L, current_size = l,
                                        decoding_train = FALSE)


  # Select the equation to replace and re-order by electrode
  which_lower_risk <- (forward_mse[[iter]][[ii]]< old_risk) %>% which
  wwhich <-	which_lower_risk[forward_mse[[iter]][[ii]][which_lower_risk] %>%
                               which.min]

  # Only if there is something that beats the risk update it
  if(length(wwhich) > 0) {
    updated_obj <- update_model(model, eqs_df, unit, wwhich, candidate_eqs, mle,
                                ztrain, ztest, Qm1)
    model   <- updated_obj$model
    Qm1_tot <- updated_obj$Qm1_tot
  }

  # Assign model
  final_model[[iter]][[ii]] <- model

  # Update risk
  cat('old risk: ', old_risk, "\n")

  old_risk <- ifelse(length(wwhich) > 0, forward_mse[[iter]][[ii]][wwhich], old_risk)
  cat("new risk: ", old_risk, "\n")

  return(list(all_risk = all_risk, old_risk = old_risk, Qm1_tot = Qm1_tot,
              model = model, new_model = new_model, forward_mse = forward_mse,
              final_model = final_model))
}

#' Updates the model if risk is lowered
#'
#' This function udates the model if the risk is lowered by one of the 13*3 equations
#' that belong to the same electrode
#' @param model current model
#' @param eqs_df data frame with equation info
#' @param unit unit selected
#' @param wwhich which are the candidate to replace that given equation
#' @param candidate_eqs list of all candidate eqs
#' @param mle whether mle or kf
#' @param ztrain neural data in the training set
#' @param ztest neural data in the testing set
#' @param Qm1 inverse of the covariance matrix
#' @return a list of vectors update and needed for next iteration
update_model <- function(model, eqs_df, unit, wwhich, candidate_eqs, mle, ztrain,
                         ztest, Qm1){


  model <- c(model[-unit], candidate_eqs[wwhich])
  Qm1_tot <- NULL

  # If mle I need the new Qm1
  if(mle){
    Q <- get_current_model_info(z = ztrain[,model], ztest = ztest[,model],
                                encoding = TRUE, decoding = FALSE, mle = mle)$Q
    # Add the eq in the cov matrix
    Qm1_tot <- get_inverse_matrix_p1(Q, Qm1)

    # Replace the order of the column correctly
    Qm1_tot <- Qm1_tot[eqs_df[model, ]$key %>% order,
                       eqs_df[model, ]$key %>% order]
  }

  # Update the final model
  model <- model[eqs_df[model, ]$key %>% order]

  return(list(Qm1_tot = Qm1_tot, model = model))
}




