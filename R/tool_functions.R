#' This functions checks whether a code produces error/warning/message
#'
#' @param code code to test
#' @return string with either error/warning/message
show_condition <- function(code) {
  tryCatch(code,
           error = function(c) "error",
           warning = function(c) "warning",
           message = function(c) "message"
  )
}


#' Leave-one-out/Add-one-in score wrapper function
#'
#' This functions wrapper function that computes leave-one-out/add-one-in score
#' in either MLE or KF settings calling the right functions depending on the settings
#' @param ztrain neural data in the training set
#' @param ztest neural data in the testing set
#' @param A the state matrix
#' @param W the state covariance matrix
#' @param H the encoding matrix
#' @param Q the encoding covariance matrix
#' @param Qm1 inverse of the encoding covariance matrix
#' @param int_train vector of integer on when the trial starts in the training set
#' @param time_train vector of time of the intervals in the training set
#' @param int_test vector of integer on when the trial starts in the testing set
#' @param time_test vector of time of the intervals in the testing set
#' @param mle boolean TRUE if mle, FALSE if KF
#' @param full_size size of the full model
#' @param current_size size of the current model
#' @param vel velocity in the training set
#' @param method forward for add one in and backward for leave one out
#' @return the average mse score for each trial left-out or add-in
get_score <- function(ztrain, ztest = NULL,
                      A = NULL, W = NULL, H = NULL, Q = NULL, Qm1 = NULL,
                      decoding_train = FALSE,
                      int_train = int.Train, time_train = Time.Train,
                      int_test = int.Test, time_test = Time.Test,
                      vel = velocity.Train, veltest = velocity.Test,
                      full_size = NULL, current_size = NULL,
                      mle, method="forward"){

  # by defualt decoding testing set
  int_dec  <- int_test
  time_dec <- time_test
  zdec     <- ztest
  vdec     <- veltest

  # decide on which set decoding
  if(decoding_train)	{
    int_dec  <- int_train
    time_dec <- time_train
    zdec     <- ztrain
    vdec     <- vel
  }

  # rescale the vars
  spike_center <- apply(ztrain, 2, mean)
  v_center	   <- apply(velocity.Train, 2, mean)
  zcenter      <- scale(zdec, center = spike_center, scale = FALSE)
  velcenter    <- scale(vdec, center = v_center, scale = FALSE)

  # wrap the correspondent function based on mle or kf opton
  score_fun <- ifelse(mle, match.fun(get_mle_score), match.fun(get_kf_score))
  mse <- score_fun(zcenter = zcenter, A = A, W = W, H = H, Q = Q, Qm1 = Qm1,
                 current_size = current_size, full_size = full_size,
                 velcenter = velcenter, int = int_dec, Time = time_dec, method = method)

  return(mse)
}

#' Computes the mle score
#'
#' This function computes the mle score, based on the method chosen
#' @param zcenter centered neural data
#' @param A the state matrix
#' @param W the state covariance matrix
#' @param H the encoding matrix
#' @param Q the encoding covariance matrix
#' @param Qm1 inverse of the encoding covariance matrix
#' @param int vector of integer on when the trial starts
#' @param Time vector of time of the intervals
#' @param full_size size of the full model
#' @param current_size size of the current model
#' @param velcener velocity in the training set
#' @param method forward for add one in and backward for leave one out
#' @return the average mle score for each trial left-out or add-in
get_mle_score <- function(zcenter, A, W, H, Q, Qm1, current_size, full_size,
                          velcenter, int, Time, method){

  # pick a function depending on the method
  FUN_mle <- ifelse(method == "forward", match.fun(get_mle_aoo_score),
                    match.fun(get_mle_loo_score))
  mle_score <- FUN_mle(zcenter = zcenter, H = H, Q = Q, Qm1 = Qm1,
                       current_size = current_size, full_size = full_size,
                       velcenter = velcenter, int = int, Time = Time)
  return(mle_score)
}


#' Computes the kf score
#'
#' This function computes the kf score, based on the method chosen
#' @param zcenter centered neural data
#' @param A the state matrix
#' @param W the state covariance matrix
#' @param H the encoding matrix
#' @param Q the encoding covariance matrix
#' @param Qm1 inverse of the encoding covariance matrix
#' @param int vector of integer on when the trial starts
#' @param Time vector of time of the intervals
#' @param full_size size of the full model
#' @param current_size size of the current model
#' @param velcener velocity in the training set
#' @param method forward for add one in and backward for leave one out
#' @return the average kf score for each trial left-out or add-in
get_kf_score <- function(zcenter, A, W, H, Q, Qm1, current_size, full_size,
                         velcenter, int, Time, method){

  # pick a function depending on the method
  mse_size <- ifelse(method == "forward", full_size - current_size, current_size)
  FUN_kf   <- ifelse(method == "forward", match.fun(kalman_decode_aoo),
                     match.fun(kalman_decode_loo))
  kf_score <- matrix(NA, ncol = mse_size, nrow = length(int))

  for(j in 1:length(int))

    kf_score[j,] <- get_kf_score_bytrial(zcenter, A, W, H, Q, Qm1, current_size,
                                         full_size, velcenter, int, Time, method,
                                         FUN_kf, j)

  kf_score <- apply(kf_score, 2, mean)
  return(kf_score)
}


#' Computes the kf score by trial
#'
#' This function computes the kf score by trial, based on the method chosen
#' @param zcenter centered neural data
#' @param A the state matrix
#' @param W the state covariance matrix
#' @param H the encoding matrix
#' @param Q the encoding covariance matrix
#' @param Qm1 inverse of the encoding covariance matrix
#' @param int vector of integer on when the trial starts
#' @param Time vector of time of the intervals
#' @param full_size size of the full model
#' @param current_size size of the current model
#' @param velcener velocity in the training set
#' @param method forward for add one in and backward for leave one out
#' @return the kf score by trial for each equations left-out or add-in
get_kf_score_bytrial <- function(zcenter, A, W, H, Q, Qm1, current_size, full_size,
                                 velcenter, int, Time, method, FUN_kf, trial_j){

  # for each trial computes the prediction
  ix   <- c(int[trial_j]:ifelse(trial_j == length(int), length(Time),
                                int[trial_j + 1] - 1))
  pred <- FUN_kf(A = A, W = W, H = H, Q = Q, length_current = current_size,
                 z = zcenter[ix,], initial_state = velcenter[ix[1],])

  # for each average mse for leave-out/add-in
  kf_score_bytrial <- vector(mode = 'numeric', length = length(pred))
  for(ii in 1:length(pred))
    kf_score_bytrial[ii] = get_mse_trial(pred[[ii]], velcenter[ix,])

  return(kf_score_bytrial)
}




#' Computes the mle leave-one-out mean error for each equation left out
#'
#' This function computes mle leave-one-out mean error for each equation left out
#' @param zcenter centered neural data
#' @param H the encoding matrix
#' @param Q the encoding covariance matrix
#' @param Qm1 inverse of the encoding covariance matrix
#' @param int vector of integer on when the trial starts
#' @param Time vector of time of the intervals
#' @param mle boolean TRUE if mle, FALSE if KF
#' @param full_size size of the full model
#' @param current_size size of the current model
#' @param velcenter centered  velocity in the training set,
#' if zcenter is from the training set
#' @return a vector of mean error for each equation left out
get_mle_loo_score <- function(zcenter, H, Q, Qm1, current_size, full_size,
                              velcenter, int, Time){

  # compute leave-one out mean error for each equation left out
  mse_mi <- vector(mode = 'numeric', length = nrow(H))

  for(ii in 1:nrow(H)){

    # leave-one-out parameters
    Hi    <- H[-ii,]
    Qi_m1 <- get_inverse_matrix_m1(Qm1, ind = ii)

    # leave-one-out prediction
    Mole <- get_Mole(H = Hi, Qm1 = Qi_m1)
    vhat <- vhat_fun(Mole, int, Time, zcenter[,-ii], velcenter, mle = TRUE)
    err  <- get_tot_mse(int, Time, velcenter, vhat)
    mse_mi[ii] <- err %>% mean
  }

  return(mse_mi)
}


#' Computes the mle add-one-in mean error for each equation included
#'
#' This function computes mle add-one-in mean error for each equation included
#' @param zcenter centered neural data
#' @param H the encoding matrix
#' @param Q the encoding covariance matrix
#' @param Qm1 inverse of the encoding covariance matrix
#' @param int vector of integer on when the trial starts
#' @param Time vector of time of the intervals
#' @param mle boolean TRUE if mle, FALSE if KF
#' @param full_size size of the full model
#' @param current_size size of the current model
#' @param velcenter centered  velocity in the training set,
#' if zcenter is from the training set
#' @return a vector of mean error for each equation included
get_mle_aoo_score <- function(zcenter, H, Q, Qm1, current_size, full_size,
                              velcenter, int, Time){

  # size of the big model
  n_incl = current_size
  n_left = full_size - n_incl

  # compute leave-one out
  mse_pi <- vector(mode = 'numeric', length = n_left)

  for(jj in 1:n_left){

    # add-one-in parameters
    index_sel <- c(1:n_incl, n_incl + jj)
    Hj <- H[index_sel,]
    Qj <- Q[index_sel, index_sel]
    Qjm1_p1 <- get_inverse_matrix_p1(Qj, Qm1)

    # add-one-in prediction
    Mole <- get_Mole(H = Hj, Qm1 = Qjm1_p1)
    vhat <- vhat_fun(Mole, int, Time, zcenter[,index_sel], velcenter, mle = TRUE)
    err  <- get_tot_mse(int, Time, velcenter, vhat)
    mse_pi[jj] <- err %>% mean
  }

  return(mse_pi)
}

#' Prune the model
#'
#' This function prunes the current model
#' @param ztrain training neural data
#' @param ztest testing neural data
#' @param int vector of integer on when the trial starts
#' @param Time vector of time of the intervals
#' @param vel is the reference velocity to use, need to be compatible with int
#' and Time
#' @param mle whether mle or kf
#' @param condition whether to prune only on some equation
#' @return a list with pruned model and current risk
prune_model <- function(ztrain, ztest, model, int, Time, vel, mle,
                        condition = 0){

  if(ncol(ztrain) <= length(model)) stop("Training set should no be evaluated
                                         in model. Please pass the whole training
                                         and testing set.")

  # compute the encoding model and current risk
  backward <- TRUE
  encoding_obj <- get_current_model_info(z = ztrain[,model], ztest = ztest[,model],
                                          encoding = TRUE, decoding = TRUE,
                                          decoding_train = FALSE, mle = mle)

  A <- encoding_obj$A
  W <- encoding_obj$W
  H <- encoding_obj$H
  Q <- encoding_obj$Q
  old_risk <- mean(encoding_obj$mse)

  while(backward){

    cat(' # equations =', length(model), '\n')

    # compute backward list
    backward_obj <- do_backward(new_reference = old_risk,
                                ztrain = ztrain[,model], ztest = ztest[,model],
                                A = A, W = W, H = H, Q = Q, model = model,
                                mle = mle, int = int, Time = Time, vel = vel,
                                condition = condition)

    # update quantities based on backward ouput
    old_risk <- backward_obj$reference
    backward <- backward_obj$backward
    model    <- backward_obj$model
    H <- backward_obj$H
    Q <- backward_obj$Q
  }

  return(list(backward_model = model, risk = old_risk))

}


#' Backward mse
#'
#' This function computes leave one out mse for a given model
#' @param new_reference error of the current model
#' @param ztrain training neural data
#' @param ztest testing neural data
#' @param A the state matrix
#' @param W the state covariance matrix
#' @param H the encoding matrix
#' @param Q the encoding covariance matrix
#' @param model given model
#' @param mle whether mle or kf
#' @param int vector of integer on when the trial starts
#' @param Time vector of time of the intervals
#' @param vel is the reference velocity to use, need to be compatible with int
#' and Time
#' @param condition whether to prune only on some equation
#' @return a list with pruned model and current risk
do_backward <- function(new_reference, ztrain, ztest, A, W, H, Q, model,
                        mle, int, Time, vel, condition){

  # initializing
  forward  = FALSE
  backward = TRUE

  # invert covariance matrix only in the mle case
  if(mle) Qm1 = solve(Q)

  # get leave_one_out score and updating the reference
  backward_mse <- get_score(ztrain = ztrain, ztest = ztest,
                       A = A, W = W, H = H, Q = Q, Qm1 = Qm1,
                       mle = mle, full_size = NULL, current_size = length(model),
                       method="backward", decoding_train = FALSE)


  # update references
  old_reference <- new_reference
  new_reference <- min(backward_mse)
  cat('old error', old_reference, "\n")
  cat('new error', new_reference, "\n")


  # select the equation to delete
  new_indeces <- !(backward_mse == min(backward_mse))
  index <- which(backward_mse == min(backward_mse))

  # update the model if we are deleting equations and the eq to delete
  # meets the condition
  if((new_reference <= old_reference) & (model[index] > condition)){
    cat("deleting equation", model[index], "\n")

    model  <- model[new_indeces]
    ztrain <- ztrain[, new_indeces]
    ztest  <- ztest[, new_indeces]
    H <- H[new_indeces, ]
    Q <- Q[new_indeces, new_indeces]

    # If mle update inverse of covariance matrix as well
    if(mle) Qm1 = get_inverse_matrix_m1(Qm1, index)

  }
  else{
    backward <- FALSE
    forward <- TRUE
  }

  return(list(ztrain = ztrain, ztest = ztest, model = model, H = H, Q = Q,
              backward = backward, forward = forward,
              reference = new_reference, model = model))
}

