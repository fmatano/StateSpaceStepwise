#' Combine and prune models
#'
#' This function combines and prunes the marginal models
#' @param load_name file to all with all data
#' @param sc_fname name of the file where to save spike count optimal model
#' @param wf_fname name of the file where to save waveform optimal model
#' @param save_name_partial name of the name where to store partial results
#' @param save_name_final name of the name where to store final results
#' @param mle whether is mle or kf
#' @param general_path what's the path for sc_fname and wf_fname
combine_and_prune <- function(load_name, sc_fname, wf_fname, save_name_partial,
                              save_name_final, mle = NULL, general_path,
                              pick_only_transformation = NULL){

  if(is.null(mle)) stop("mle is NULL, it needs to be TRUE/FALSE")

  # Get the data and set the path where to get the data from
  if(is.null(load_name)) cat("You are not loading the files, the algorithm
                             is using whatever data is loaded in your
                             working directory \n ")
  if(!is.null(load_name)) source(load_name)

  # Subselect only data needed here, if for instance you're running only ace
  if(!is.null(pick_only_transformation)){

    eqs_sorted_eq <- eqs_sorted_eq %>%
      subset((transf == pick_only_transformation) |
               wf_transf == pick_only_transformation)

    Ztrain_tot <- Ztrain_tot[,eqs_sorted_eq$eq]
    Ztest_tot  <- Ztest_tot[,eqs_sorted_eq$eq]
    eqs_sorted_eq$eq <- 1:length(eqs_sorted_eq$eq)
    rownames(eqs_sorted_eq) <- 1:length(eqs_sorted_eq$eq)
  }


  # Load spike count optimal model
  load(paste(general_path, sc_fname, ".RData", sep = ""))

  cat("-- Pruning spike count model --", "\n")

  # prune sc model
  sc_model <- prune_model(ztrain = Ztrain_tot, ztest = Ztest_tot,
                          model = model, mle = mle,
                          condition = 0, int = int.Test, Time = Time.Test,
                          vel = velocity.Test)$backward_model
  testthat::expect_equal(eqs_sorted_eq[model,]$wf %>% unique, 0)

  # Load waveform model
  augmented_model <- wf_model_after <- wf_model_before <- list()
  for(wf_mom in 1:2){
    for(wf_feat in 1:4){

      # Load wf model
      f_name <- paste(general_path, wf_fname, "_wfmom", wf_mom,"-feat",
                      wf_feat, ".RData", sep = "")
      load(f_name)

      # Test the files contain the right information
      testthat::expect_equal(eqs_sorted_eq[model,]$wf %>% unique, wf_mom)
      testthat::expect_equal(eqs_sorted_eq[model,]$wf_key %>% unique, wf_feat)

      cat("-- Pruning wf model: moment", wf_mom, ", feature", wf_feat, "--", "\n")

      # Step1: prune wf model
      wf_model_before  <- append(wf_model_before, list(wf_model = model))
      current_wf_model <- prune_model(ztrain = Ztrain_tot, ztest = Ztest_tot,
                                      model = model, mle = mle,
                                      condition = n_units*n_transf*n_lags,
                                      int = int.Test, Time = Time.Test,
                                      vel = velocity.Test)$backward_model
      wf_model_after <- append(wf_model_after, list(wf_model = current_wf_model))


      # Step2: combine the wf and sc model (pruned)
      combined_model <- c(sc_model, current_wf_model)

      cat("--Combining and pruning spike count and wf model: moment", wf_mom,
          ", feature", wf_feat, "--", "\n")


      # Step3: prune the combined model on everything
      tmp <- prune_model(ztrain = Ztrain_tot, ztest = Ztest_tot, mle = mle,
                         model = combined_model, int=int.Test, Time=Time.Test,
                         vel = velocity.Test,
                         condition = 0)$backward_model
                         # condition = n_units*n_transf*n_lags)$backward_model
      augmented_model <- append(augmented_model, list(wf_model = tmp))
      print(augmented_model)

      # Save the model
      save(augmented_model, sc_model, wf_model_after, wf_model_before,
           file = save_name_partial)

    }
  }


  # If stage 1 is already finished, do only this following part
  # Combine all the models and prune
  # source(load_file)
  load(save_name_partial)

  final_scwf_model <- augmented_model %>% unlist %>% unique
  # tmp <- sample(final_scwf_model, 500)
  final_scwf_model_pruned <- prune_model(Ztrain_tot, Ztest_tot, final_scwf_model,
  																		  int = int.Test, Time = Time.Test,
  																		  vel = velocity.Test, mle = mle,
                            						condition = 0)$backward_model

  save(augmented_model, sc_model, wf_model_after, wf_model_before,
       final_scwf_model, final_scwf_model_pruned,
  					 file = save_name_final, tmp)
}
