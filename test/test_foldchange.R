# library(Seurat)

library(rhdf5)

# library(fastde)

library(tictoc)

# read input
input <- h5read("/home/tpan/build/wave/input.h5", "array/block0_values")
samplenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis0")
genenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis1")

labels <- as.vector(h5read("/home/tpan/build/wave/labels.h5", "array/block0_values"))

colnames(input) <- genenames
rownames(input) <- samplenames

# count number of labels
#num_labels = nlevels(labels)

cat(sprintf("input size:  r %d X c %d\n", nrow(input), ncol(input)))
cat(sprintf("Labels rows: %d \n", length(labels)))

L <- unique(sort(labels))

cat(sprintf("Labels unique: %d \n", length(L)))

# typeof(labels)
# for ( i in labels) {
#     cat(sprintf("%d, ", i))
# }
# cat(sprintf("\n"))
# labels

tic("fastde")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdefc <- fastde::ComputeFoldChange(input, labels, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = FALSE, threads = as.integer(4))
toc()

x <- as.integer(row.names(fastdefc$fc))
ord <- order(x)
x
ord

fastdefcsorted = fastdefc$fc[ord, ]
fastdefc_pct1_sorted = fastdefc$pct.1[ord, ]
fastdefc_pct2_sorted = fastdefc$pct.2[ord, ]
#fastdefc


fastdefcsorted
# fastdefc$pct.1[, 1]
# fastdefc$pct.2[, 1]


tic("fastde df")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdefc_df <- fastde::ComputeFoldChange(input, labels, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = TRUE, threads = as.integer(4))
toc()

print(fastdefc_df)


# time and run wilcox.test
tic("Seurat builtin")

x <- t(input)

colnames(x) <- samplenames
rownames(x) <- genenames

seuratfc <- matrix(c(0), ncol = ncol(input), nrow = length(L) )
seuratperc1 <- matrix(c(0), ncol = ncol(input), nrow = length(L) )
seuratperc2 <- matrix(c(0), ncol = ncol(input), nrow = length(L) )
cat(sprintf("dim seuratfc %d %d\n", dim(seuratfc)[1], dim(seuratfc)[2]))


for ( c in L ) {
    cells.1 <- which(labels %in% c)
    cells.2 <- which(! (labels %in% c))

    cat(sprintf("cells.1 %d. cells.2 %d\n", length(cells.1), length(cells.2)))

    v <- Seurat::FoldChange(x, cells.1, cells.2, mean.fxn=rowMeans, fc.name="fc")
    # cat(colnames(v))
    # cat(sprintf(" size %d x %d, cluster %d\n", dim(v)[1], dim(v)[2], c))
    # cat(sprintf(" output %d %f, %f\n", length(as.vector(v$fc)), as.vector(v$fc)[11], v[c,11]))

    # cat(sprintf("R wilcox %f\n", v))
    seuratfc[c, ] <- as.vector(v$fc)
    seuratperc1[c, ] <- as.vector(v$pct.1)
    seuratperc2[c, ] <- as.vector(v$pct.2)

    # seuratfc

}
toc()

seuratfc
# seuratperc1[, 1]
# seuratperc2[, 1]



## compare by calculating the residuals.
res = seuratfc - fastdefcsorted
residual = sqrt(mean(res * res))

cat(sprintf("R naive vs fastde residual fc = %f\n", residual))


res = seuratperc1 - fastdefc_pct1_sorted
residual = sqrt(mean(res * res))

cat(sprintf("R naive vs fastde residual pct.1 = %f\n", residual))


res = seuratperc2 - fastdefc_pct2_sorted
residual = sqrt(mean(res * res))

cat(sprintf("R naive vs fastde residual pct.2 = %f\n", residual))
