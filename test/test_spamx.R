unlink(".RData")
library("devtools")

library(rhdf5)
library(spam)
library(spam64)
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
head(sobject@Dimnames)
message("PBMC3K dim")
sobject@Dim
nrows <- sobject@Dim[1]
toc()

nclusters = 30

# generate fake class labels.
tic("GEN clusters")
clusters = 1:nclusters
labels <- sample(clusters, nrows, replace = TRUE)
toc()

tic("wilcox pbmc3k")
# run wilcox
wilcox_de <- fastde::FastFindAllMarkers64(sobject, idents.clusters = labels, test.use = 'fastwmw')
toc()

# run ttest
tic("ttest pbmc3k")
ttest_de <- fastde::FastFindAllMarkers64(sobject, idents.clusters = labels, test.use = 'fast_t')
toc()



tic("READ brain 1.3M")
# f2 <- "~/data/SingleCell/1M_neurons_filtered_gene_bc_matrices_h5.h5"
f2 <- "~/scgc/data/1M_neurons_filtered_gene_bc_matrices_h5.h5"

sobject2 <- fastde::Read10X_h5_big(f2)
message("Brain 1M:  x size ", length(sobject2@entries), " i size ", length(sobject2@colindices), " p size", length(sobject2@rowpointers))
head(sobject2@entries)
head(sobject2@colindices)
head(sobject2@rowpointers)
head(sobject2@Dimnames)
message("Brain 1M dim")
sobject2@dimension
#sobject2
toc()

nrows <- sobject@dimension[1]

nclusters = 30

tic("GEN clusters")
# generate fake class labels.
clusters = 1:nclusters
labels <- sample(clusters, nrows, replace = TRUE)
toc()

# run wilcox
tic("wilcox 1.3M")
wilcox_de <- fastde::FastFindAllMarkers64(sobject, idents.clusters = labels, test.use = 'fastwmw')
toc()

# run ttest
tic("ttest 1.3M")
wilcox_de <- fastde::FastFindAllMarkers64(sobject, idents.clusters = labels, test.use = 'fast_t')
toc()

# res <- PCA_big(sobject2)
# res


# fid <- H5Fopen(f)
## open up the dataset to add attributes to, as a class
# did <- H5Dopen(fid, "aNEONSite/temperature")
