library(Seurat)
library(rhdf5)
library(fastde)   # loading library can take a while (3 s) if not preloaded here.
library(tictoc)
library(Matrix)

comparemat <- function(name, A, B) {
    diff <- A - B
    maxdiff <- max(diff)
    mindiff <- min(diff)
    mediandiff <- median(diff)
    meandiff <- mean(diff * diff)
    stdevdiff <- sd(diff * diff)
    cat(sprintf("%s : diff range [%f, %f], median %f, mean %f, sd %f\n", 
        name, mindiff, maxdiff, mediandiff, meandiff, stdevdiff))
}

# tic("read")

# # read input
# samplenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis0")
# genenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis1")

# input <- rsparsematrix(length(samplenames), length(genenames), 0.1)

# labels <- as.vector(h5read("/home/tpan/build/wave/labels.h5", "array/block0_values"))

# colnames(input) <- genenames
# rownames(input) <- samplenames
# toc()

cat(sprintf("NOTE: rounding of negative number are different in R and C++:\n\tR: round(%f) = %f.  C would return %f\n\tMATRIX COMPARISON RESULTS ARE EXPECTED TO BE DIFFERENT BUT WITH ZERO MEAN.\n\tFASTDE FoldChange results are also not rounded.\n", -0.5, round(-0.5), round(-0.6)))


tic("generate")
# max is 2B.
ncols = 2000  # features
nrows = 2000  # samples
nclusters = 30

clusters = 1:nclusters
labels <- sample(clusters, nrows, replace = TRUE)

input <- rsparsematrix(nrows, ncols, 0.05)

samplenames <- paste(1:nrows);
genenames <- paste(1:ncols);
colnames(input) <- genenames
rownames(input) <- samplenames
toc()


tic("get unique labels")
# count number of labels
#num_labels = nlevels(labels)

cat(sprintf("input size:  r %d X c %d\n", nrow(input), ncol(input)))
cat(sprintf("Labels rows: %d \n", length(labels)))

L <- unique(sort(labels))

cat(sprintf("Labels unique: %d \n", length(L)))

toc()



# time and run wilcox.test
tic("Seurat builtin")

x <- t(as.matrix(input))
colnames(x) <- samplenames
rownames(x) <- genenames

seuratfc <- matrix(c(0), ncol = ncol(input), nrow = length(L) )
seuratperc1 <- matrix(c(0), ncol = ncol(input), nrow = length(L) )
seuratperc2 <- matrix(c(0), ncol = ncol(input), nrow = length(L) )
cat(sprintf("dim seuratfc %d %d\n", dim(seuratfc)[1], dim(seuratfc)[2]))

rownames(seuratfc) <- L
rownames(seuratperc1) <- L
rownames(seuratperc2) <- L

for ( c in L ) {
    cells.1 <- which(labels %in% c)
    cells.2 <- which(! (labels %in% c))

    cat(sprintf("cells.1 %d. cells.2 %d\n", length(cells.1), length(cells.2)))

    v <- Seurat::FoldChange(x, cells.1, cells.2, mean.fxn=rowMeans, fc.name="fc")
    # cat(colnames(v))
    # cat(sprintf(" size %d x %d, cluster %d\n", dim(v)[1], dim(v)[2], c))
    # cat(sprintf(" output %d %f, %f\n", length(as.vector(v$fc)), as.vector(v$fc)[11], v[c,11]))

    # cat(sprintf("R wilcox %f\n", v))
    pos = which(L == c)
    seuratfc[pos, ] <- as.vector(v$fc)
    seuratperc1[pos, ] <- as.vector(v$pct.1)
    seuratperc2[pos, ] <- as.vector(v$pct.2)

    # seuratfc
}
toc()

# seuratfc[, 1]
# seuratperc1[, 1]
# seuratperc2[, 1]


tic("fastde df")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdefc_df <- fastde::ComputeFoldChange(as.matrix(input), labels, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = TRUE, threads = as.integer(4))
toc()

tic("fastde")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdefc <- fastde::ComputeFoldChange(as.matrix(input), labels, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = FALSE, threads = as.integer(4))
toc()
cat(sprintf("output %d X %d\n", nrow(fastdefc$fc), ncol(fastdefc$fc)))

tic("Ordering by cluster num")
x <- as.integer(row.names(fastdefc$fc))
ord <- order(x)
# x
# ord

fastdefcsorted = fastdefc$fc[ord, ]
fastdefc_pct1_sorted = fastdefc$pct.1[ord, ]
fastdefc_pct2_sorted = fastdefc$pct.2[ord, ]
#fastdefc
toc()

# fastdefcsorted[, 1]
# fastdefc_pct1_sorted[, 1]
# fastdefc_pct2_sorted[, 1]


# print(fastdefc_df)


## compare by calculating the residuals.

comparemat("R vs fastde fc", seuratfc, fastdefcsorted)
comparemat("R vs fastde pct1", seuratperc1, fastdefc_pct1_sorted)
comparemat("R vs fastde pct2", seuratperc2, fastdefc_pct2_sorted)

diff = which(abs(seuratperc1 - fastdefc_pct1_sorted) >= 0.0005)
cat(sprintf("different 1: seurat %f, fastde %f \n", seuratperc1[diff[1]], fastdefc_pct1_sorted[diff[1]]))
cat(sprintf("different 2: seurat %f, fastde %f \n", seuratperc1[diff[2]], fastdefc_pct1_sorted[diff[2]]))
cat(sprintf("different 3: seurat %f, fastde %f \n", seuratperc1[diff[3]], fastdefc_pct1_sorted[diff[3]]))


tic("sparse fastde df")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdesfc_df <- fastde::ComputeSparseFoldChange(input, labels, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = TRUE, threads = as.integer(4))
toc()


tic("sparse fastde")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdesfc <- fastde::ComputeSparseFoldChange(input, labels, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = FALSE, threads = as.integer(4))
toc()
cat(sprintf("output %d X %d\n", nrow(fastdefc$fc), ncol(fastdefc$fc)))



tic("Ordering by cluster num")
x <- as.integer(row.names(fastdesfc$fc))
ord <- order(x)
# x
# ord

fastdesfcsorted = fastdesfc$fc[ord, ]
fastdesfc_pct1_sorted = fastdesfc$pct.1[ord, ]
fastdesfc_pct2_sorted = fastdesfc$pct.2[ord, ]
#fastdefc
toc()

# fastdefcsorted[, 1]
# fastdefc_pct1_sorted[, 1]
# fastdefc_pct2_sorted[, 1]



# print(fastdefc_df)



## compare by calculating the residuals.

comparemat("R vs sparse fastde fc", seuratfc, fastdesfcsorted)
comparemat("R vs sparse fastde pct1", seuratperc1, fastdesfc_pct1_sorted)
comparemat("R vs sparse fastde pct2", seuratperc2, fastdesfc_pct2_sorted)

