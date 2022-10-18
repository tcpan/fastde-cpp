
#' @rdname FastFoldChange
#' @param object Any object that's not a seurat
#' @param cells.clusters cell labels, integer cluster ids for each cell.
#' @param features Subset a particular features.
#' @param slot Slot to pull data from; note that if \code{test.use} is "negbinom", "poisson", or "DESeq2",
#' \code{slot} will be set to "counts"
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @param base base for log function. does not apply when slot is "scaled.data"
#' @param fc.name Name of the fold change, average difference, or custom function column
#' in the output data.frame
#' @param return.dataframe if TRUE, return a dataframe, else return a 2D matrix.
#' @param verbose Print report verbosely
#' @return If return.dataframe == FALSE, returns a list of matrices, first is fold change values, 
#' next 2 are percent above thresholds in the 2 classes.
#' Else return a dataframe with the matrices reshaped as columns, values grouped by gene.
#'
#' If NULL, use all features
#' 
#' @importFrom methods is
#'
#' @concept differential_expression
#' @name FastFoldChange
#' @export
#' @method FastFoldChange default
FastFoldChange.default <- function(
  object,
  cells.clusters = NULL,
  features = NULL,
  slot = "data",
  pseudocount.use = TRUE,
  base = 2,
  fc.name = NULL,
  return.dataframe = TRUE,
  verbose = FALSE,
  ...
) {
  if (verbose) { print("TCP FASTDE: FastFoldChange.default") }
  if (verbose) { tictoc::tic("FastFoldChange.default setup") }

  log.use <- ifelse(
    test = slot == "scale.data",
    yes = FALSE,
    no = TRUE
  )
  expm1.use <- ifelse(
    test = slot == "data",
    yes = TRUE,
    no = FALSE
  )

  if ( is.null(x = cells.clusters) ) {
    clusters <- Seurat::Idents(object = object)
  } else {
    clusters <- cells.clusters
  }
  if ( is.null(clusters) ) {
    message("ERROR:  null cluster")
  } else if (verbose) {
      message("using cluster names")
      print(levels(clusters))
  }
  # clusters <- cells.clusters %||% Seurat::Idents(object = object)

  # Omit the decimal value of e from the column name if base == exp(1)
  base.text <- ifelse(
    test = base == exp(1),
    yes = "",
    no = base
  )
  if ( is.null(x = fc.name) ) {
    fc.name = ifelse(
      test = slot == "scale.data",
      yes = "avg_diff",
      no = paste0("avg_log", base.text, "FC")
    )
  }
  # fc.name <- fc.name %||% ifelse(
  #   test = slot == "scale.data",
  #   yes = "avg_diff",
  #   no = paste0("avg_log", base.text, "FC")
  # )
  
  if (! is.null(features)) {
    data <- object[features, , drop=FALSE]
  } else {
    data <- object
  }
  if (verbose) { tictoc::toc() }
  
  tictoc::tic("FastFoldChange.default FastPerformFC")
  
  PerformFCFunc <- if (is(data, 'dgCMatrix') | is(data, 'dgCMatrix64') )  {
    FastPerformSparseFC
  } else {
    FastPerformFC
  }

  fc.results <- PerformFCFunc(data, clusters,
    features_as_rows = TRUE,
    calc_percents = TRUE, fc_name = fc.name, 
    use_expm1 = expm1.use, min_threshold = 0.0, 
    use_log = log.use, log_base = base, use_pseudocount = pseudocount.use, 
    as_dataframe = return.dataframe,
    threads = get_num_threads())
  tictoc::toc()

  if (verbose) { print("TCP FASTDE: FastFoldChange.default DONE") }
  return(fc.results)
}

#' @rdname FastFoldChange
#' @param object A Seurat assay object
#' @param cells.clusters cell labels, integer cluster ids for each cell.
#' @param features Subset a particular features.
#' @param slot Slot to pull data from; note that if \code{test.use} is "negbinom", "poisson", or "DESeq2",
#' \code{slot} will be set to "counts"
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @param base base for log function. does not apply when slot is "scaled.data"
#' @param fc.name Name of the fold change, average difference, or custom function column
#' in the output data.frame
#' @param return.dataframe if TRUE, return a dataframe, else return a 2D matrix.
#' @param verbose Print report verbosely
#' @return If return.dataframe == FALSE, returns a list of matrices, first is fold change values, 
#' next 2 are percent above thresholds in the 2 classes.
#' Else return a dataframe with the matrices reshaped as columns, values grouped by gene.
#'  
#' @concept differential_expression
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat Idents
#' @name FastFoldChange
#' @export
#' @method FastFoldChange Assay
FastFoldChange.Assay <- function(
  object,
  cells.clusters = NULL,
  features = NULL,
  slot = "data",
  pseudocount.use = 1,
  base = 2,
  fc.name = NULL,
  return.dataframe = TRUE,
  verbose = FALSE,
  ...
) {
  if (verbose) { print("TCP FASTDE: FastFoldChange.Assay") }
  if (verbose) { tictoc::tic("FastFoldChange.Assay setup") }
  data <- Seurat::GetAssayData(object = object, slot = slot)

  # using the passed in seurat object to get cluster ids.
  if ( is.null(x = cells.clusters) ) {
    clusters <- Seurat::Idents(object = object)
  } else {
    clusters <- cells.clusters
  }
  # clusters <- cells.clusters %||% Seurat::Idents(object = object)
  
  if (verbose) { tictoc::toc() }
  if (verbose) { tictoc::tic("FastFoldChange.Assay dispatch FastFoldChange") }
  fc.results <- FastFoldChange(
    object = data,
    cells.clusters = clusters,
    features = features,
    slot = slot,
    pseudocount.use = pseudocount.use,
    base = base,
    fc.name = fc.name,
    return.dataframe = return.dataframe,
    verbose = verbose,
    ...
  )
  if (verbose) { tictoc::toc() }
  if (verbose) { print("TCP FASTDE: FastFoldChange.Assay DONE") }
  return(fc.results)
}

#' @rdname FastFoldChange
#' @param object A Seurat dim reduced object
#' @param cells.clusters cell labels, integer cluster ids for each cell.
#' @param features Subset a particular features.
#' @param fc.name Name of the fold change, average difference, or custom function column
#' in the output data.frame
#' @param return.dataframe if TRUE, return a dataframe, else return a 2D matrix.
#' @param verbose Print report verbosely
#' @return If return.dataframe == FALSE, returns a list of 1 matrix, fold change values.
#' Else return a dataframe with the matrix reshaped as a column, values grouped by gene.
#' 
#' @concept differential_expression
#' @importFrom Seurat Embeddings
#' @importFrom Seurat Idents
#' @importFrom methods is
#'
#' @name FastFoldChange
#' @export
#' @method FastFoldChange DimReduc
FastFoldChange.DimReduc <- function(
  object,
  cells.clusters = NULL,
  features = NULL,
  fc.name = NULL,
  return.dataframe = TRUE,
  verbose = FALSE,
  ...
) {
  if (verbose) { print("TCP FASTDE: FastFoldChange.DimReduc") }
  if (verbose) { tictoc::tic("FastFoldChange.DimReduc init") }

  if ( is.null(fc.name) ) {
    fc.name <- "avg_diff"
  }
  # fc.name <- fc.name %||% "avg_diff"
  data <- Seurat::Embeddings(object = object)  # transpose as the reduction appears to be transposed.
  if (! is.null(x = features) ) {
    data <- data[, features, drop=FALSE]
  }

  # does the size match?
  if ( is.null(x = cells.clusters) ) {
    clusters <- Seurat::Idents(object = object)
  } else {
    clusters <- cells.clusters
  }
  # clusters <- cells.clusters %||% Seurat::Idents(object = object)

  if (verbose) { tictoc::toc() }
  tictoc::tic("FastFoldChange.DimReduc FastPerformFC")
  # Calculate avg difference.  This is just rowMeans.
  
  PerformFCFunc <- if (is(data, 'dgCMatrix') | is(data, 'dgCMatrix64') )  {
    FastPerformSparseFC
  } else {
    FastPerformFC
  }

  fc.results <- PerformFCFunc(data, clusters,
    features_as_rows = FALSE,
    calc_percents = FALSE, fc_name = fc.name, 
    use_expm1 = FALSE, min_threshold = 0.0, 
    use_log = FALSE, log_base = 2, use_pseudocount = FALSE, 
    as_dataframe = return.dataframe,
    threads = get_num_threads())
  tictoc::toc()
  if (verbose) { print("TCP FASTDE: FastFoldChange.DimReduc DONE") }
  return(fc.results)
}

#' @rdname FastFoldChange
#' 
#' @param object A Seurat object
#' @param cells.clusters cell labels, integer cluster ids for each cell.
#' @param features Subset a particular features.
#' @param group.by Regroup cells into a different identity class prior to
#' calculating fold change (see example in \code{\link{FastFindMarkers}})
#' @param subset.ident Subset a particular identity class prior to regrouping.
#' Only relevant if group.by is set (see example in \code{\link{FastFindMarkers}})
#' @param assay Assay to use in fold change calculation.  mutually exclusive with reduction
#' @param reduction Reduction to use - will calculate average difference on cell embeddings
#' @param slot Slot to pull data from; note that if \code{test.use} is "negbinom", "poisson", or "DESeq2",
#' \code{slot} will be set to "counts"
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @param base The base with respect to which logarithms are computed.
#' @param fc.name Name of the fold change, average difference, or custom function column
#' in the output data.frame
#' @param return.dataframe if TRUE, return a dataframe, else return a 2D matrix.
#' @param verbose Print report verbosely
#' @return If return.dataframe == FALSE, returns a list of matrices, first is fold change values, 
#' next 2 are percent above thresholds in the 2 classes (not for dim reduc).
#' Else return a dataframe with the matrices reshaped as columns, values grouped by gene.
#' 
#' 
#' @concept differential_expression
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat Idents
#' @name FastFoldChange
#' @export
#' @method FastFoldChange Seurat
FastFoldChange.Seurat <- function(
  object,
  cells.clusters = NULL,
  features = NULL,
  group.by = NULL,
  subset.ident = NULL,
  assay = NULL,
  reduction = NULL,
  slot = 'data',
  pseudocount.use = 1,
  base = 2,
  fc.name = NULL,
  return.dataframe = TRUE,
  verbose = FALSE,
  ...
) {
  if (verbose) { print("TCP FASTDE: FastFoldChange.Seurat") }
  if (verbose) { tictoc::tic("FastFoldChange.Seurat init") }
  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Seurat::Idents(object = object) <- group.by
  }
  if (!is.null(x = assay) && !is.null(x = reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  # select which data to use
  if (is.null(x = reduction)) {
    if ( is.null(x = assay) ) {
      assay = Seurat::DefaultAssay(object = object)
    }
    # assay <- assay %||% Seurat::DefaultAssay(object = object)
    data.use <- object[[assay]]
    cellnames.use <-  colnames(x = data.use)
  } else {
    data.use <- object[[reduction]]
    cellnames.use <- rownames(data.use)
  }

  if ( is.null(x = cells.clusters) ) {
    clusters <- Seurat::Idents(object = object)
  } else {
    clusters <- cells.clusters
  }
  # clusters <- cells.clusters %||% Seurat::Idents(object = object)

  if (verbose) { tictoc::toc() }
  if (verbose) { tictoc::tic("FastFoldChange.Seurat dispatch FastFoldChange") }

  fc.results <- FastFoldChange(
    object = data.use,
    cells.clusters = clusters,
    features = features,
    slot = slot,
    pseudocount.use = pseudocount.use,
    base = base,
    fc.name = fc.name,
    return.dataframe = return.dataframe,
    verbose = verbose,
    ...
  )
  if (verbose) { tictoc::toc() }
  if (verbose) { print("TCP FASTDE: FastFoldChange.Seurat DONE") }
  return(fc.results)
}

#' Fold Change
#'
#' Calculate log fold change and percentage of cells expressing each feature
#' for different identity classes.
#'
#' If the slot is \code{scale.data} or a reduction is specified, average difference
#' is returned instead of log fold change and the column is named "avg_diff".
#' Otherwise, log2 fold change is returned with column named "avg_log2_FC".
#'
# @examples
# data("pbmc_small")
# FastFoldChange(pbmc_small, ident.1 = 1)
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @rdname FastFoldChange
#' @export FastFoldChange
#' @return Returns a data.frame, or a list of 1 or 3 matrices.
#' @seealso \code{FastFindMarkers}
FastFoldChange <- function(object, ...) {
  UseMethod(generic = 'FastFoldChange', object = object)
}

FastPerformSparseFC <- function(data, clusters, 
  features_as_rows = FALSE,
  calc_percents = FALSE, 
  fc_name = NULL, 
  use_expm1 = FALSE, 
  min_threshold = 0.0, 
  use_log = FALSE, 
  log_base = 2, 
  use_pseudocount = FALSE,
  as_dataframe = TRUE,
  threads = 1) {

  # Calculate fold change incrementally - if data is too big.
  # print("dims data.use")
  # print(dim(data))
  print("PERFORM SPARSE MAT FC")

  if (! is.factor(clusters)) {
    clusters <- as.factor(clusters)
  }

  # no need to break up the data - can't deal with larger datasets than have non-zero entries.

  # # need to put features into columns.
  # if (features_as_rows == TRUE) {
  #   # slice and transpose
  #   # if (is(data, 'dgCMatrix64') )  {
  #   #   dd <- fastde::sp_transpose(data)
  #   # } else {
  #   #   dd <- fastde::sp_transpose(data)
  #   # }
  #   dd <- fastde::sp_transpose(data, threads=threads)
  # } else {
  #   # slice the data
  #   dd <- data
  # }

  # print('IN FC')
  # str(dd)

  # PerformFCFunc <- ComputeFoldChangeSparse


  # output has features in columns, and clusters in rows
  fc.results <- ComputeFoldChangeSparse(data, as.integer(clusters),
    features_as_rows = features_as_rows,
    calc_percents = calc_percents, fc_name = fc_name, 
    use_expm1 = use_expm1, min_threshold = min_threshold, 
    use_log = use_log, log_base = log_base, use_pseudocount = use_pseudocount, 
    as_dataframe = as_dataframe,
    threads = threads)
  if (as_dataframe == TRUE) {
    fc.results$cluster <- factor(as.numeric(fc.results$cluster), labels=levels(clusters))
  } else {
    # return data the same way we got it
    if (features_as_rows == TRUE) {
      fc.results[[fc_name]] <- t(fc.results[[fc_name]])   # put features back in rows
      if (calc_percents == TRUE) {
        fc.results$pct.1 <- t(fc.results$pct.1)   # now features are in rows
        fc.results$pct.2 <- t(fc.results$pct.2)   # now features are in rows
      }
    }
  }
  return(fc.results)
}


# incrementally perform FoldChange.  DENSE ONLY
FastPerformFC <- function(data, clusters, 
  features_as_rows = FALSE,
  calc_percents = FALSE, 
  fc_name = NULL, 
  use_expm1 = FALSE, 
  min_threshold = 0.0, 
  use_log = FALSE, 
  log_base = 2, 
  use_pseudocount = FALSE,
  as_dataframe = TRUE, 
  threads = 1) {

  print("PERFORM DENSE MAT FC")

  # Calculate fold change incrementally - if data is too big.
  # print("dims data.use")
  # print(dim(data))

  if (! is.factor(clusters)) {
    clusters <- as.factor(clusters)
  }

  # ======== break up the data if needed.
  if (features_as_rows == TRUE) {
    # features in rows
    nfeatures = nrow(data)
    # samples in columns
    nsamples = ncol(data)
  } else {
    # features in columns already
    nfeatures = ncol(data)
    # samples in rows.
    nsamples = nrow(data)
  }
  # print("nfeatures, nsamples")
  # print(nfeatures)
  # print(nsamples)

  # get the number of features to process at a time.
  max_elem <- 1024*1024*1024
  block_size <- pmin(max_elem %/% nsamples, nfeatures)
  nblocks <- (nfeatures + block_size - 1) %/% block_size
  # print("block size, nblocks")
  # print(block_size)
  # print(nblocks)

  # need to put features into columns.
  if (features_as_rows == TRUE) {
    # slice and transpose
    dd <- t(data[1:block_size, ])
  } else {
    # slice the data
    dd <- data[, 1:block_size]
  }
  # output has features in columns, and clusters in rows
  fc.results <- ComputeFoldChange(dd, as.integer(clusters),
    calc_percents = calc_percents, fc_name = fc_name, 
    use_expm1 = use_expm1, min_threshold = min_threshold, 
    use_log = use_log, log_base = log_base, use_pseudocount = use_pseudocount, 
    as_dataframe = as_dataframe,
    threads = threads)
  if (as_dataframe == FALSE) {
    # return data the same way we got it
    if (features_as_rows == TRUE) {
      fc.results[[fc_name]] <- t(fc.results[[fc_name]])   # put features back in rows
      if (calc_percents == TRUE) {
        fc.results$pct.1 <- t(fc.results$pct.1)   # now features are in rows
        fc.results$pct.2 <- t(fc.results$pct.2)   # now features are in rows
      }
    }
  } # if dataframe, already in the right format.
  # print("dims fc.results 1")
  # print(dim(fc.results))

  if (nblocks > 1) {
    for (i in 1:(nblocks - 1)) {
      # compute bounds 
      start <- i * block_size + 1
      end <- pmin(nfeatures, (i + 1) * block_size )
      # slice the data
      if (features_as_rows == TRUE) {
        # slice and transpose
        dd <- t(data[start:end, ])
      } else {
        # slice the data
        dd <- data[, start:end]
      }

      fcr <- ComputeFoldChange(dd, as.integer(clusters),
        calc_percents = calc_percents, fc_name = fc_name, 
        use_expm1 = use_expm1, min_threshold = min_threshold, 
        use_log = use_log, log_base = log_base, use_pseudocount = use_pseudocount, 
        as_dataframe = as_dataframe,
        threads = threads)

      if (as_dataframe == TRUE)  {
        fc.results <- rbind(fc.results, fcr)
      } else {
        # return data the same way we got it
        if (features_as_rows == TRUE) {
          fc.results[[fc_name]] <- rbind(fc.results[[fc_name]], t(fcr[[fc_name]]))
          if (calc_percents == TRUE) {
            fc.results$pct.1 <- rbind(fc.results, t(fcr$pct.1))   # now features are in rows
            fc.results$pct.2 <- rbind(fc.results, t(fcr$pct.2))   # now features are in rows
          }
        } else {
          fc.results[[fc_name]] <- cbind(fc.results[[fc_name]], fcr[[fc_name]])
          if (calc_percents == TRUE) {
            fc.results$pct.1 <- cbind(fc.results, fcr$pct.1)   # features are in cols
            fc.results$pct.2 <- cbind(fc.results, fcr$pct.2)   # features are in cols
          }
        }
      }
      # print("dims fc.results i")
      # print(dim(fc.results))

    }
  }
  # print("dims fc.results final")
  # print(dim(fc.results))
  if (as_dataframe == TRUE) {
    fc.results$cluster <- factor(as.numeric(fc.results$cluster), labels=levels(clusters))
  }
  return(fc.results)

}
