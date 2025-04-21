#' Make spatial partitions (using k-means clustering)
#'
#' Divides a spatial sample of polygons into irregularly shaped spatial partitions
#'     using k-means clustring. It uses the spatial coordinates defined by the polygon centroids.
#' @param x A SpatialPolygonsDataframe
#' @param n_folds The number of folds (i.e. partitions)
#' @param seed_cv The random seed defined to reproduce results
#' @return A vector containing numeric labels defining each fold


spcvFoldsPoly <- function(x, n_folds = 5, seed_cv = 11082017){
  # x is a polygon spatial object
  x$objectid <- 1:nrow(x)
  centroid_poly <- rgeos::gCentroid(x, byid=TRUE, id = x$objectid)
  set.seed(seed_cv)
  poly_crds <- sp::coordinates(centroid_poly)
  sp_resamp <- sperrorest::partition_kmeans(poly_crds, nfold = n_folds)

  cv_labels <- rep(NA, nrow(x))
  # Spatial cv labels
  for(i in 1:n_folds){
    test_label <- sp_resamp[[1]][[i]]$test
    cv_labels[test_label] <- i
  }

  return(cv_labels)
}




#' Perform (repeated) spatial cross-validation of random walk optimization
#'
#' @param x Random walk grid search results (list) from multiple runouts
#' @param slide_plys A SpatialPolygonsDataframe
#' @param n_folds The number of folds (i.e. partitions)
#' @param repetitions Number of repetitions
#' @param from_save (Logical) if TRUE, will load save files from current working directory
#' @return A vector containing numeric labels defining each fold

rwSPCV <- function(x, slide_plys, n_folds, repetitions, from_save = FALSE){

  if(from_save){

    x <- list()

    files <- list.files(pattern = "result_rw_roc")

    for(i in 1:length(files)){
      res_nm <- paste("result_rw_roc_", i, ".Rd", sep="")
      res_obj_nm <- load(res_nm)
      roc_result <- get(res_obj_nm)
      x[[i]] <- roc_result
    }

  }


  rwslp_vec <- as.numeric(dimnames(x[[1]])[[1]])
  rwexp_vec <- as.numeric(dimnames(x[[1]])[[2]])
  rwper_vec <- as.numeric(dimnames(x[[1]])[[3]])

  rep_spcv_rw <- list()

  for(n in 1:repetitions){

    rep_seed <- 11082017 + n

    fold_labels <- spcvFoldsPoly(slide_plys, n_folds = n_folds, seed_cv = rep_seed)

    spcv_results <- list()

    # Folds ###############
    for(k in 1:n_folds){

      roc_median <- list()
      roc_iqr <- list()
      roc_mean <- list()
      roc_sd <- list()

      polyid_train.vec <- which(fold_labels != k)
      polyid_test.vec <- which(fold_labels == k)

      for(PER in 1:length(rwper_vec)){

        per_list <- list()

        for(i in 1:length(polyid_train.vec)){
          obj.id <- polyid_train.vec[i]
          res <- x[[obj.id]]
          per_list[[i]] <- as.matrix(res[,,PER])

        }

        Y <- do.call(cbind, per_list)
        Y <- array(Y, dim=c(dim(per_list[[1]]), length(per_list)))

        Y.median <- apply(Y, c(1, 2), median, na.rm = TRUE)
        Y.iqr<- apply(Y, c(1, 2), IQR, na.rm = TRUE)
        Y.mean <- apply(Y, c(1, 2), mean, na.rm = TRUE)
        Y.sd<- apply(Y, c(1, 2), sd, na.rm = TRUE)

        roc_median[[PER]] <- Y.median
        roc_iqr[[PER]] <- Y.iqr
        roc_mean[[PER]] <- Y.mean
        roc_sd[[PER]] <- Y.sd

      }

      rocMedian <- do.call(cbind, roc_median)
      rocMedian <- array(rocMedian, dim=c(dim(roc_median[[1]]), length(roc_median)))

      rocIQR <- do.call(cbind, roc_iqr)
      rocIQR <- array(rocIQR, dim=c(dim(roc_iqr[[1]]), length(roc_iqr)))

      rocMean <- do.call(cbind, roc_mean)
      rocMean <- array(rocMean, dim=c(dim(roc_mean[[1]]), length(roc_mean)))

      rocSD <- do.call(cbind, roc_sd)
      rocSD <- array(rocSD, dim=c(dim(roc_sd[[1]]), length(roc_sd)))


      # Find which combination of RW parameters had the highest auroc
      ROCwh <- which(rocMedian==max(rocMedian), arr.ind=T)
      ROCwh

      # In ties, take first record
      if(length(ROCwh > 3)) {ROCwh <- ROCwh[1,]}

      #train_auc <- rocMedian[ROCwh]

      train_auc_median <- rocMedian[ROCwh[1], ROCwh[2], ROCwh[3]]
      train_auc_iqr <- rocIQR[ROCwh[1], ROCwh[2], ROCwh[3]]
      train_auc_mean <- rocMean[ROCwh[1], ROCwh[2], ROCwh[3]]
      train_auc_sd <- rocSD[ROCwh[1], ROCwh[2], ROCwh[3]]

      slp_opt <- rwslp_vec[ROCwh[1]]
      exp_opt <- rwexp_vec[ROCwh[2]]
      per_opt <- rwper_vec[ROCwh[3]]

      # to get test results
      roc_median <- list()
      roc_iqr <- list()
      roc_mean <- list()
      roc_sd <- list()

      for(PER in 1:length(rwper_vec)){

        per_list <- list()

        for(i in 1:length(polyid_test.vec)){
          obj.id <- polyid_test.vec[i]
          res <- x[[obj.id]]
          per_list[[i]] <- as.matrix(res[,,PER])

        }

        Y <- do.call(cbind, per_list)
        Y <- array(Y, dim=c(dim(per_list[[1]]), length(per_list)))

        Y.median <- apply(Y, c(1, 2), median, na.rm = TRUE)
        Y.iqr<- apply(Y, c(1, 2), IQR, na.rm = TRUE)
        Y.mean <- apply(Y, c(1, 2), mean, na.rm = TRUE)
        Y.sd<- apply(Y, c(1, 2), sd, na.rm = TRUE)

        roc_median[[PER]] <- Y.median
        roc_iqr[[PER]] <- Y.iqr
        roc_mean[[PER]] <- Y.mean
        roc_sd[[PER]] <- Y.sd

      }

      test_rocMedian <- do.call(cbind, roc_median)
      test_rocMedian <- array(test_rocMedian, dim=c(dim(roc_median[[1]]), length(roc_median)))

      test_rocIQR <- do.call(cbind, roc_iqr)
      test_rocIQR <- array(test_rocIQR, dim=c(dim(roc_iqr[[1]]), length(roc_iqr)))

      test_rocMean <- do.call(cbind, roc_mean)
      test_rocMean <- array(test_rocMean, dim=c(dim(roc_mean[[1]]), length(roc_mean)))

      test_rocSD <- do.call(cbind, roc_sd)
      test_rocSD <- array(test_rocSD, dim=c(dim(roc_sd[[1]]), length(roc_sd)))

      #test_roc <- do.call(cbind, roc_median)
      #test_roc<- array(test_roc, dim=c(dim(roc_median[[1]]), length(roc_median)))

      test_auc_median <- test_rocMedian[ROCwh[1], ROCwh[2], ROCwh[3]]
      test_auc_iqr <- test_rocIQR[ROCwh[1], ROCwh[2], ROCwh[3]]
      test_auc_mean <- test_rocMean[ROCwh[1], ROCwh[2], ROCwh[3]]
      test_auc_sd <- test_rocSD[ROCwh[1], ROCwh[2], ROCwh[3]]

      # result for test dataset
      result_list <- list(
        n_train = length(polyid_train.vec),
        n_test = length(polyid_test.vec),
        rw_slp_opt = slp_opt,
        rw_exp_opt = exp_opt,
        rw_per_opt = per_opt,
        test_auroc_median = test_auc_median,
        train_auroc_median = train_auc_median
      )

      spcv_results[[k]] <- result_list

    }


    spcv_rw_folds <- data.frame(matrix(unlist(spcv_results), nrow = length(spcv_results), byrow = T))
    names(spcv_rw_folds) <- names(spcv_results[[1]])

    rep_spcv_rw[[n]] <- spcv_rw_folds

  }

  rep_spcv_rw$settings <- list(rwslp_vec = rwslp_vec, rwexp_vec = rwexp_vec, rwper_vec = rwper_vec)


  return(rep_spcv_rw)

}





#' Pool random walk spatial cross-validation results
#'
#' Creates a summary of the frequency and performance (AUROC) of optimal parameter
#'     sets across spatial cross-validation repetitions.
#' @param x The (list) result of spcvRndWalk()
#' @param plot_freq (Logical) if TRUE, will a create bubble plot of optimal parameter
#'      set frequencies across grid search space
#' @return A data frame summarizing the frequency and performance of each
#'     optimal parameter sets


rwPoolSPCV<- function(x, plot_freq = FALSE){

  pool_rw <- do.call(rbind, x[1:(length(x) - 1)])

  rwslp_vec <- x$settings$rwslp_vec
  rwexp_vec <- x$settings$rwexp_vec
  rwper_vec <- x$settings$rwper_vec

  # summarize by frequency and median score for optimal parameter sets
  opt_sets <- data.frame(slp = pool_rw$rw_slp_opt, per = pool_rw$rw_per_opt, exp = pool_rw$rw_exp_opt)
  freq_sets <- table(opt_sets)

  slp_nm <- as.numeric(attributes(freq_sets)$dimnames$slp)
  per_nm <- as.numeric(attributes(freq_sets)$dimnames$per)
  exp_nm <- as.numeric(attributes(freq_sets)$dimnames$exp)

  # Get array index of pairs
  set_ind <- which(freq_sets !=0, arr.ind = TRUE)
  sets <- data.frame(slp = slp_nm[set_ind[,1]],
                     per = per_nm[set_ind[,2]],
                     exp = exp_nm[set_ind[,3]],
                     freq = freq_sets[set_ind])

  sets$rel_freq <- sets$freq/nrow(pool_rw) * 100

  #Find median AUROC values
  set_id <- paste(sets$slp, sets$per, sets$exp)
  pool_rw$set_id <- paste(pool_rw$rw_slp_opt, pool_rw$rw_per_opt, pool_rw$rw_exp_opt)

  for(i in 1:length(set_id)){

    auroc_i <- pool_rw$test_auroc_median[which(pool_rw$set_id == set_id[i])]
    sets$median_auroc[i] <- median(auroc_i, na.rm = TRUE)
    sets$iqr_auroc[i] <- IQR(auroc_i, na.rm = TRUE)

  }

  if(plot_freq){
    gg <- ggplot2::ggplot(sets, ggplot2::aes(x=sets$per, y=sets$exp)) +
      # Can improve by making different colors for different slope thresholds...

      ggplot2::geom_point(alpha=0.5, ggplot2::aes(col = as.factor(sets$slp), size = sets$rel_freq)) +

      ggplot2::labs(size = "Relative\nfrequency (%)", col = "Slope threshold") +

      ggplot2::scale_x_continuous(expression(paste("Persistence factor")),
                         limits = c(min(rwper_vec), max = max(rwper_vec))) +
      ggplot2:: scale_y_continuous(expression(paste("Exponent of divergence")),
                         limits = c(min(rwexp_vec), max = max(rwexp_vec)+.1)) +
      ggplot2::theme_light() +
      ggplot2::theme(text = ggplot2::element_text(family = "Arial", size = 8), axis.title = ggplot2::element_text(size = 9),
            axis.text = ggplot2::element_text(size = 8))

    print(gg)
  }

  return(sets)

}




