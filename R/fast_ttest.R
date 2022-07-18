library(tictoc)

#' Differential expression using t-test, Sparse matrix.
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a Wilcoxon Rank Sum test. Makes use of limma::rankSumTestWithCorrelation for a
#' more efficient implementation of the wilcoxon test. Thanks to Yunshun Chen and
#' Gordon Smyth for suggesting the limma implementation.
#'
#' @param data.use Data matrix to test.  rows = features, columns = samples
#' @param cells.clusters cell labels, integer cluster ids for each cell.  array of size same as number of samples
#' @param features.as.rows controls transpose
#' @param verbose Print report verbosely
#' @param return.dataframe if TRUE, return a dataframe, else return a 2D matrix.
#' @param ... Extra parameters passed to wilcox.test
#'
#' @return If return.dataframe == FALSE, returns a p-value matrix of putative differentially expressed
#' features, with genes in rows and clusters in columns.  else return a dataframe with "p-val" column, 
#' with results for the clusters grouped by gene.
#'
#' @export
#'
# @examples
# data("pbmc_small")
# pbmc_small
# FastSparseDiffTTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#             cells.2 = WhichCells(object = pbmc_small, idents = 2))
#
FastSparseDiffTTest <- function(
  data.use,
  cells.clusters,
  features.as.rows,
  verbose = FALSE,
  return.dataframe = TRUE,
  ...
) {
  # input has each row being a gene, and each column a cell (sample).  -
  #    based on reading of the naive wilcox.test imple in Seurat.
  # fastde input is expected to : each row is a gene, each column is a sample
  message("USING FastSparseTTestDE")
  if (verbose) {
    print(rownames(data.use[1:20, ]))
    print(cells.clusters[1:20]) 
  }
  if (! is.factor(cells.clusters)) {
    cells.clusters <- as.factor(cells.clusters)
  }
  labels <- levels(cells.clusters)
  if (verbose) {
    print(as.numeric(cells.clusters[1:20]))
    print(cells.clusters[1:20])
  }

  # print("dims data.use")
  # print(dim(data.use))

  # NOTE: label's order may be different - cluster ids are stored in unordered_map so order may be different.
  # IMPORTANT PART IS THE SAME ID is kep.

  # two sided : 2
  # print(head(data.use))
  tic("FasTTestDETest fast_t")


  # need to put features into columns.
  if (features.as.rows == TRUE) {
    # slice and transpose
    if (is(data.use, 'dgCMatrix64') )  {
      dd <- fastde::sp_transpose(data.use)
    } else {
      dd <- fastde::sp_transpose(data.use)
    }
  } else {
    # slice the data
    dd <- data.use
  }

  PerformDEFunc <- if (is(dd, 'dgCMatrix64') )  {
    sparse64_ttest_fast
  } else {
    sparse_ttest_fast
  }


  # default type is 2 (2 sided)
  p_val <- PerformDEFunc(dd, as.integer(cells.clusters),
          as_dataframe = return.dataframe, threads = get_num_threads(),
          alternative = as.integer(2), var_equal = FALSE)
  if (return.dataframe == FALSE) {
    # return data the same way we got it
    if (features.as.rows == TRUE) {
      p_val <- t(p_val)   # put features back in rows
    }
  } # if dataframe, already in the right format.
  # print("dims p_val 1")
  # print(dim(p_val))

  toc()
  # print("dims pval final")
  # print(dim(p_val))

  if ( return.dataframe == TRUE ) {
    if (verbose) {
      print("head of p_val orig")
      print(p_val[1:20, ])
    }
    p_val$cluster <- factor(as.numeric(p_val$cluster), labels = labels)
    # gene names are same in order, but replicated nlabels times.
    if (features.as.rows == TRUE) {
      genesf <- as.factor(rownames(data.use))
    } else {
      genesf <- as.factor(colnames(data.use))
    }
    gene_labels <- levels(genesf)  # 
    genes <- rep(as.numeric(genesf), each = length(labels))
    p_val$gene <- factor(genes, labels = gene_labels)      
         
    if (verbose) {
      print("head of p_val")
      print(p_val[1:20, ])
    }
  } else {
    if (features.as.rows == TRUE ) {
      colnames(p_val) <- labels
      rownames(p_val) <- rownames(data.use)
    } else {
      rownames(p_val) <- labels
      colnames(p_val) <- colnames(data.use)
    }
  }
  
  return(p_val)
}



#' Differential expression using Wilcoxon Rank Sum
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a Wilcoxon Rank Sum test. Makes use of limma::rankSumTestWithCorrelation for a
#' more efficient implementation of the wilcoxon test. Thanks to Yunshun Chen and
#' Gordon Smyth for suggesting the limma implementation.
#'
#' @param data.use Data matrix to test.  rows = features, columns = samples
#' @param cells.clusters cell labels, integer cluster ids for each cell.  array of size same as number of samples
#' @param features.as.rows controls transpose
#' @param verbose Print report verbosely
#' @param return.dataframe if TRUE, return a dataframe, else return a 2D matrix.
#' @param ... Extra parameters passed to wilcox.test
#'
#' @return If return.dataframe == FALSE, returns a p-value matrix of putative differentially expressed
#' features, with genes in rows and clusters in columns.  else return a dataframe with "p-val" column, 
#' with results for the clusters grouped by gene.
#'
#' @export
#'
# @examples
# data("pbmc_small")
# pbmc_small
# FastDiffTTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#             cells.2 = WhichCells(object = pbmc_small, idents = 2))
#
FastDiffTTest <- function(
  data.use,
  cells.clusters,
  features.as.rows,
  verbose = FALSE,
  return.dataframe = TRUE,
  ...
) {
  # input has each row being a gene, and each column a cell (sample).  -
  #    based on reading of the naive wilcox.test imple in Seurat.
  # fastde input is expected to : each row is a gene, each column is a sample
  message("USING FastTTestDE")
  if (verbose) {
    print(rownames(data.use[1:20, ]))
    print(cells.clusters[1:20]) 
  }
  if (! is.factor(cells.clusters)) {
    cells.clusters <- as.factor(cells.clusters)
  }
  labels <- levels(cells.clusters)
  if (verbose) {
    print(as.numeric(cells.clusters[1:20]))
    print(cells.clusters[1:20])
  }

  # print("dims data.use")
  # print(dim(data.use))

  # NOTE: label's order may be different - cluster ids are stored in unordered_map so order may be different.
  # IMPORTANT PART IS THE SAME ID is kep.

  # two sided : 2
  # print(head(data.use))
  tic("FastTTestDETest fast_t")
  if (features.as.rows == TRUE) {
    # features in rows
    nfeatures = nrow(data.use)
    # samples in columns
    nsamples = ncol(data.use)
  } else {
    # features in columns already
    nfeatures = ncol(data.use)
    # samples in rows.
    nsamples = nrow(data.use)
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
  if (features.as.rows == TRUE) {
    # slice and transpose
    dd <- t(data.use[1:block_size, ])
  } else {
    # slice the data
    dd <- data.use[, 1:block_size]
  }
  # default type is 2 (2 sided)
  p_val <- ttest_fast(dd, as.integer(cells.clusters),
          as_dataframe = return.dataframe, threads = get_num_threads(),
          alternative =  as.integer(2), var_equal = FALSE)

  if (return.dataframe == FALSE) {
    # return data the same way we got it
    if (features.as.rows == TRUE) {
      p_val <- t(p_val)   # put features back in rows
    }
  } # if dataframe, already in the right format.
  # print("dims p_val 1")
  # print(dim(p_val))

  if (nblocks > 1) {
    for (i in 1:(nblocks - 1)) {
      # compute bounds 
      start <- i * block_size + 1
      end <- pmin(nfeatures, (i + 1) * block_size )
      # slice the data
      if (features.as.rows == TRUE) {
        # slice and transpose
        dd <- t(data.use[start:end, ])
      } else {
        # slice the data
        dd <- data.use[, start:end]
      }

      # default type is 2 (2 sided)
      pv <- ttest_fast(dd, as.integer(cells.clusters),
          as_dataframe = return.dataframe, threads = get_num_threads(),
          alternative =  as.integer(2), var_equal = FALSE)
      if (return.dataframe == TRUE)  {
        p_val <- rbind(p_val, pv)
      } else {
        # return data the same way we got it
        if (features.as.rows == TRUE) {
          p_val <- rbind(p_val, t(pv))
        } else {
          p_val <- cbind(p_val, pv)
        }
      }
      # print("dims pval")
      # print(dim(p_val))
    }
  }
  toc()
  # print("dims pval final")
  # print(dim(p_val))

  if ( return.dataframe == TRUE ) {
    if (verbose) {
      print("head of p_val orig")
      print(p_val[1:20, ])
    }
    p_val$cluster <- factor(as.numeric(p_val$cluster), labels = labels)
    # gene names are same in order, but replicated nlabels times.
    if (features.as.rows == TRUE) {
      genesf <- as.factor(rownames(data.use))
    } else {
      genesf <- as.factor(colnames(data.use))
    }
    gene_labels <- levels(genesf)  # 
    genes <- rep(as.numeric(genesf), each = length(labels))
    p_val$gene <- factor(genes, labels = gene_labels)      
         
    if (verbose) {
      print("head of p_val")
      print(p_val[1:20, ])
    }
  } else {
    if (features.as.rows == TRUE ) {
      colnames(p_val) <- labels
      rownames(p_val) <- rownames(data.use)
    } else {
      rownames(p_val) <- labels
      colnames(p_val) <- colnames(data.use)
    }
  }
  
  return(p_val)
}

