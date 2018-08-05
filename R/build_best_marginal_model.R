#' Build marginal model
#'
#' This function builds the marginal model for spike count model and all
#' waveoforms
#' @param load_file file to all with all data
#' @param f_name name of the file where to save spike count optimal model
#' @param mle whether is mle or kf
#' @param wf_mom waveform moment, if 0 computes for spike counts
#' @param wf_feat waveform feature, if 0 computes for spike counts
#' @param max_iter number of max iteration before stopping
#' @param transformation transformation selected 1 linear, 2 ace, 3 sqrt, NULL all
#' @param pick_only_transformation in case you only want to select a subset of
#' the data, indicate which transformation in particular
#' @param general_path path for saving the file
build_marginal_model <- function(load_name, f_name, mle = NULL,
                                 wf_mom = 0, wf_feat = 0, max_iter = 50,
                                 transformation = NULL,
                                 pick_only_transformation = NULL, general_path,
                                 eqs_df = NULL, starting_model = NULL){

  # Initial condition check
  is_sc <- FALSE
  if(is.null(mle)) stop("mle is NULL, it needs to be TRUE/FALSE")
  if((wf_mom == 0) |(wf_feat == 0)) is_sc = TRUE


  # Get the data and set the path where to get the data from
  if(is.null(load_name)) cat("You are not loading the files, the algorithm
                             is using whatever data is loaded in your
                             working directory \n ")
  if(!is.null(load_name)) source(load_name)

  if(is.null(eqs_df)) eqs_df <- eqs_sorted_eq

  # Best sc model, starting from optimal uniform lag
  # general_path <- "inst/data/lindsay_unsorted/robustness/"

  starting_transf <- 1
  if(!is.null(transformation)) starting_transf <- transformation

  # If pick_only_transformation is not null, select only that part of dataset
  if(!is.null(pick_only_transformation)){

    eqs_df <- eqs_df %>%
      subset((transf == pick_only_transformation) |
               wf_transf == pick_only_transformation)

    Ztrain_tot <- Ztrain_tot[,eqs_df$eq]
    Ztest_tot  <- Ztest_tot[,eqs_sdf$eq]
    starting_transf  <- pick_only_transformation
    eqs_df$eq <- 1:length(eqs_sdf$eq)
    rownames(eqs_df) <- 1:length(eqs_df$eq)
  }

  # Define Starting model
  if(is_sc){

    cat("******* BUILD BEST SC MODEL ******* \n")
    if(is.null(starting_model))
      starting_model <- (eqs_df %>%
                           subset(lag == opt_lag) %>%
                           subset(transf == starting_transf) %>%
                           subset(wf == 0))$eq %>% sort

    save_name <- paste(general_path, f_name, ".RData", sep = "")
  } else{

    cat("*** COMPUTE BEST MODEL WF MOMENT", wf_mom, "- FEATURE", wf_feat, "\n")
    if(is.null(starting_model))
      starting_model <- (eqs_df %>%
                           subset(lag == opt_lag) %>%
                           subset(wf_key == wf_feat) %>%
                           subset(wf == wf_mom) %>%
                           subset(wf_transf == starting_transf))$eq %>% sort
    save_name <- paste(general_path, f_name, "_wfmom", wf_mom,"-feat",
                       wf_feat, ".RData", sep = "")
  }


  if(length(starting_model) != n_units) stop("Starting model should have same
                                             length as number of units!")

  get_best_model(ztrain = Ztrain_tot,  ztest = Ztest_tot, eqs_df = eqs_df,
                 model = starting_model, n_units = n_units, mle = mle,
                 save_name = save_name, int_dec = int.Test, time_dec = Time.Test,
                 vel = velocity.Train, veltest = velocity.Test, is_sc = is_sc,
                 wf_mom = wf_mom, wf_feat = wf_feat, max_iter = max_iter,
                 transf_selected = starting_transf)
}

