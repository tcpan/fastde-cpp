unlink(".RData")
library("devtools")

library(rhdf5)
library(fastde)
library(tictoc)



tic("READ pbmc3k")

#f <- "~/data/SingleCell/pbmc3k/filtered_gene_bc_matrices/hg19/"
f <- "~/scgc/data/pbmc3k/filtered_gene_bc_matrices/hg19/"
sobject <- Seurat::Read10X(f)
message("PBMC3K:  x size ", length(sobject@x), " i size ", length(sobject@i), " p size", length(sobject@p))
head(sobject@x)
head(sobject@i)
head(sobject@p)
head(sobject@Dimnames[[1]][1:10])
head(sobject@Dimnames[[2]][1:10])
message("PBMC3K dim")
sobject@Dim
class(sobject@Dim)
nrows <- sobject@Dim[1]
toc()

message("i = 0")
print(which(sobject@i == 0))
message("end")

str(sobject)
nclusters = 30

# generate fake class labels.
tic("GEN clusters")
clusters = 1:nclusters
labels <- sample(clusters, nrows, replace = TRUE)
toc()

tic("sparse fastde FC")
# time and run BioQC
fastdefc <- fastde::ComputeFoldChangeSparse(sobject, labels, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = FALSE, threads = as.integer(1))
toc()
str(fastdefc)



tic("wilcox pbmc3k")
# run wilcox
wilcox_de <- fastde::FastFindAllMarkers64(sobject, idents.clusters = labels, test.use = 'fastwmw')
toc()

# run ttest
tic("ttest pbmc3k")
ttest_de <- fastde::FastFindAllMarkers64(sobject, idents.clusters = labels, test.use = 'fast_t')
toc()

message("FINISHED 1")

# so2 <- new('dgCMatrix64', x = sobject@x, i = sobject@i, p = sobject@p, dimension = sobject@Dim, 
#     Dim = sobject@Dim, Dimnames = sobject@Dimnames)



tic("convert to 64bit")
so64 <- as.dgCMatrix64(sobject)

nrows <- so64@dimension[1]
toc()



tic("sparse fastde 64 FC")
# time and run BioQC
fastdefc2 <- fastde::ComputeFoldChangeSparse64(so64, labels, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = FALSE, threads = as.integer(1))
toc()
str(fastdefc2)

tic("wilcox pbmc3k 64")
# run wilcox
wilcox_de <- fastde::FastFindAllMarkers64(so64, idents.clusters = labels, test.use = 'fastwmw')
toc()

# run ttest
tic("ttest pbmc3k 64")
ttest_de <- fastde::FastFindAllMarkers64(so64, idents.clusters = labels, test.use = 'fast_t')
toc()

message("FINISHED 2")


print("======write to h5 and read in as dgCMatrix64")

tic("WRITE pbm3k h5")
#f <- "~/data/SingleCell/pbmc3k.h5"
f <- "~/scgc/data/pbmc3k.h5"
fastde::Write10X_h5(sobject, f)
toc()

tic("READ pbmc3k h5")
sobject4 <- Seurat::Read10X_h5(f)
str(sobject4)
toc()

tic("READ pbmc3k h5")
sobject3 <- fastde::Read10X_h5_big(f)
str(sobject3)
toc()

message("FINISHED 3")

print("==========READING Big H5")

tic("READ brain 1.3M")
# f2 <- "~/data/SingleCell/1M_neurons_filtered_gene_bc_matrices_h5.h5"
f2 <- "~/scgc/data/1M_neurons_filtered_gene_bc_matrices_h5.h5"

sobject2 <- fastde::Read10X_h5_big(f2)
#sobject2
toc()

nrows <- sobject2@dimension[1]
str(sobject2)

nclusters = 30

tic("GEN clusters")
# generate fake class labels.
clusters = 1:nclusters
labels2 <- sample(clusters, nrows, replace = TRUE)
toc()


tic("sparse fastde 64 FC")
# time and run BioQC
fastdefc3 <- fastde::ComputeFoldChangeSparse64(sobject2, labels2, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = FALSE, threads = as.integer(1))
toc()
str(fastdefc3)

# run wilcox
tic("wilcox 1.3M")
wilcox_de2 <- fastde::FastFindAllMarkers64(sobject2, idents.clusters = labels2, test.use = 'fastwmw')
toc()
str(wilcox_de2)

# run ttest
tic("ttest 1.3M")
wilcox_de2 <- fastde::FastFindAllMarkers64(sobject2, idents.clusters = labels2, test.use = 'fast_t')
toc()
str(wilcox_de2)

message("FINISHED 4")
