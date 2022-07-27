

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
#' @importFrom BioQC wmwTest
#' @importFrom methods is
#'
#' @export
#'
# @examples
# data("pbmc_small")
# pbmc_small
# WilcoxDETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#             cells.2 = WhichCells(object = pbmc_small, idents = 2))
#
BioQCDETest <- function(
  data.use,
  cells.clusters,
  features.as.rows,
  verbose = FALSE,
  return.dataframe = TRUE,
  ...
) {
  bioqc.check <- fastDEPackageCheck("BioQC", error = FALSE)
  message("USING BioQC")
  if (verbose) {
    print(rownames(data.use[1:20, ]))
    print(cells.clusters[1:20])
  }
  # input has each row being a gene, and each column a cell (sample).  -
  #    based on reading of the naive wilcox.test imple in Seurat.
  if (! is.factor(cells.clusters)) {
    cells.clusters <- as.factor(cells.clusters)
  }
  LL <- levels(cells.clusters)
  if (verbose) {
    print(as.numeric(cells.clusters[1:20]))
    print(LL)
  }
  if ( bioqc.check[1] ) {

    tictoc::tic("BioQCDETest BioQC::wmwTest")
    labels = list()
    for (i in 1:length(LL)) {
      labels[[i]] <- as.integer(cells.clusters) %in% i
    }

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
    # get the number of features to process at a time.
    max_elem <- 1024*1024*1024
    block_size <- pmin(max_elem %/% nsamples, nfeatures)
    nblocks <- (nfeatures + block_size - 1) %/% block_size

    # need to put features into columns.
    tictoc::tic("BioQC sparse to dense")
    if (is(data.use, 'dgCMatrix')) {
      if (features.as.rows == TRUE) {
        # slice and transpose
        dd <- fastde::sp_to_dense_transposed(data.use[1:block_size, ])
      } else {
        # slice the data
        dd <- fastde::sp_to_dense(data.use[, 1:block_size])
      }
    } else if (is(data.use, 'dgCMatrix64') ) {

      if (features.as.rows == TRUE) {
        # slice and transpose
        dd <- fastde::sp_to_dense_transposed(data.use[1:block_size, ])
      } else {
        # slice the data
        dd <- fastde::sp_to_dense(data.use[, 1:block_size])
      }
    } else {
      if (features.as.rows == TRUE) {
        # slice and transpose
        dd <- t(data.use[1:block_size, ])
      } else {
        # slice the data
        dd <- data.use[, 1:block_size]
      }
    }
    tictoc::toc()

    pv <- BioQC::wmwTest(dd, labels, valType = "p.two.sided")
    # return data the same way we got it
    if ((return.dataframe == FALSE) && (features.as.rows == TRUE)) {
      pv <- t(pv)   # put features back in rows
    }
    
    if (nblocks > 1) {
      for (i in 1:(nblocks - 1)) {
        # compute bounds 
        start <- i * block_size + 1
        end <- pmin(nfeatures, (i + 1) * block_size )
        # slice the data
        tictoc::tic("BioQC sparse to dense")
        if (is(data.use, 'dgCMatrix')) {
          if (features.as.rows == TRUE) {
            # slice and transpose
            dd <- fastde::sp_to_dense_transposed(data.use[start:end, ])
          } else {
            # slice the data
            dd <- fastde::sp_to_dense(data.use[, start:end])
          }
        } else if (is(data.use, 'dgCMatrix64') )  {

          if (features.as.rows == TRUE) {
            # slice and transpose
            dd <- fastde::sp_to_dense_transposed(data.use[start:end, ])
          } else {
            # slice the data
            dd <- fastde::sp_to_dense(data.use[, start:end])
          }
        } else { # dense matrix
          if (features.as.rows == TRUE) {
            # slice and transpose
            dd <- t(data.use[start:end, ])
          } else {
            # slice the data
            dd <- data.use[, start:end]
          }
        }
        tictoc::toc()
        pvi <- BioQC::wmwTest(dd, labels, valType = "p.two.sided")

        # return data the same way we got it
        if ((return.dataframe == FALSE) && (features.as.rows == TRUE)) {
          pv <- rbind(pv, t(pvi))
        } else {
          pv <- cbind(pv, pvi)
        }
      }
    }
    # p_val HERE has clusters in rows and features/genes in columns
    tictoc::toc()

    if (return.dataframe  == TRUE ) {
      # NOTE that bioqc is iterating the cluster labels in numeric order, so no need for factor here.
      clusters <- rep(1:length(LL), times=nfeatures)  # clusters repeated nfeatures number of times
      if (verbose) {
        print(clusters[1:20])
      }
      # use rep to make the genes vector, which need numeric.
      # so need to convert rownames to factor
      # but now also need to get levels, as the mapping may be different.

      if (verbose) {
        print(rownames(data.use)[1:20])
      }

      if (features.as.rows == TRUE) {
        genesf <- as.factor(rownames(data.use))
      } else {
        genesf <- as.factor(colnames(data.use))
      }
      gene_labels <- levels(genesf)  # 
      if (verbose) {
        print(as.numeric(genesf[1:20]))
        print(genesf[1:20])
      }
      # genenames <- levels(genes)
      genes <- rep(as.numeric(genesf), each = length(LL))
      # genes <- as.factor(genes)
      # levels(genes) <- c(genenames)
      if (verbose) {
        print(genes[1:20])
      }
      # print(length(genes))

      p_val <- data.frame(cbind(clusters, genes, as.vector(pv)))
      colnames(p_val) <- c("cluster", "gene", "p_val")
      if (verbose) {
        print(p_val[1:20, , drop=FALSE])
      }

      p_val$cluster <- factor(as.numeric(clusters), labels = LL)
      p_val$gene <- factor(genes, labels = gene_labels)      
      if (verbose) {
        print("head of p_val")
        print(p_val[1:20, , drop=FALSE])
      }
    } else {
      # just a matrix returned.
      if (features.as.rows == TRUE) {
        colnames(p_val) <- LL
        rownames(p_val) <- rownames(data.use)
      } else {
        rownames(p_val) <- LL
        colnames(p_val) <- colnames(data.use)
      }

    }
  } else {
    message("Please install the BioQC package")
  }
  return(p_val)
}
