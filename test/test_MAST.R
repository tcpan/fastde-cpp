#' @include generics.R
#'
NULL

library(tictoc)

# Differential expression using MAST
#
# Identifies differentially expressed genes between two groups of cells using
# a hurdle model tailored to scRNA-seq data. Utilizes the MAST package to run
# the DE testing.
#
# @references Andrew McDavid, Greg Finak and Masanao Yajima (2017). MAST: Model-based
# Analysis of Single Cell Transcriptomics. R package version 1.2.1.
# https://github.com/RGLab/MAST/
#
# @param data.use Data to test
# @param labels cell labels as a list
# @param latent.vars Confounding variables to adjust for in DE test. Default is
# "nUMI", which adjusts for cellular depth (i.e. cellular detection rate). For
# non-UMI based data, set to nGene instead.
# @param verbose print output
# @param \dots Additional parameters to zero-inflated regression (zlm) function
# in MAST
# @details
# To use this method, please install MAST, using instructions at https://github.com/RGLab/MAST/
#
# @return Returns a p-value ranked matrix of putative differentially expressed
# genes.
#
#' @importFrom stats relevel
#
# @export
#
# @examples
# \dontrun{
#   data("pbmc_small")
#   pbmc_small
#   MASTDETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, idents = 1),
#               cells.2 = WhichCells(object = pbmc_small, idents = 2))
# }
#
FastMASTDETest <- function(
  data.use,
  labels,
  latent.vars = NULL,
  verbose = TRUE,
  ...
) {
  # Check for MAST
  if (!PackageCheck('MAST', error = FALSE)) {
    stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
  }
  group.info <- data.frame(row.names = c(cells.1, cells.2))   # cells in groups 1 and 2.  would this an empty data frame?
  latent.vars <- latent.vars %||% group.info                # %||% is a conditional assignment if null.
  # latent.vars should be columns (nGene, nUMI), where as group.info is rows.
  group.info[cells.1, "group"] <- "Group1"  # this is basically the labels for the 2 grops
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])  # then change that column to factors.
  latent.vars.names <- c("condition", colnames(x = latent.vars))  # not sure about this...
  latent.vars <- cbind(latent.vars, group.info)  # add the group info to latent vars.
  latent.vars$wellKey <- rownames(x = latent.vars)   # these should be either passed in, or cell names, 
  fdat <- data.frame(rownames(x = data.use))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  # exprsArray:     col == cells, row == genes. 
  # cData is column are cell ids.  wellKey
  # fData is rows: primerid
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(x = data.use),
    check_sanity = FALSE,
    cData = latent.vars,
    fData = fdat
  )
  cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  cond <- relevel(x = cond, ref = "Group1")
  SummarizedExperiment::colData(sca)$condition <- cond
  fmla <- as.formula(
    object = paste0(" ~ ", paste(latent.vars.names, collapse = "+"))
  )
  zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
  summaryCond <- summary(object = zlmCond, doLRT = 'conditionGroup2')
  summaryDt <- summaryCond$datatable
  # fcHurdle <- merge(
  #   summaryDt[contrast=='conditionGroup2' & component=='H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
  #   summaryDt[contrast=='conditionGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid'
  # ) #logFC coefficients
  # fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  p_val <- summaryDt[summaryDt[, "component"] == "H", 4]
  genes.return <- summaryDt[summaryDt[, "component"] == "H", 1]
  # p_val <- subset(summaryDt, component == "H")[, 4]
  # genes.return <- subset(summaryDt, component == "H")[, 1]
  to.return <- data.frame(p_val, row.names = genes.return)
  return(to.return)
}
