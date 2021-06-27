library(BioQC)
library(rhdf5)
library(fastde)
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

# read input
# tic("read")
# labels <- as.vector(h5read("/home/tpan/build/wave/labels.h5", "array/block0_values"))
# genenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis1")
# samplenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis0")
# wilcox <- h5read("/home/tpan/build/wave/test-wilcox.h5", "array/block0_values")

# input <- rsparsematrix(length(samplenames), length(genenames), 0.1)

# colnames(input) <- genenames
# rownames(input) <- samplenames
# toc()

tic("generate")
# max is 2B.
ncols = 10  # features  - test with a small number. this is totally parallel
nrows = 200000  # samples
nclusters = 30

clusters = 1:nclusters
labels <- sample(clusters, nrows, replace = TRUE)

input <- rsparsematrix(nrows, ncols, 0.05)

samplenames <- as.character(1:nrows);
genenames <- as.character(1:ncols);
colnames(input) <- genenames
rownames(input) <- samplenames
toc()


tic("get unique labels")
# count number of labels
#num_labels = nlevels(labels)

# cat(sprintf("test size:  r %d X c %d.\n", nrow(wilcox), ncol(wilcox)))
cat(sprintf("input size:  r %d X c %d\n", nrow(input), ncol(input)))
cat(sprintf("Labels: %d \n", length(labels)))

L <- unique(sort(labels))

cat(sprintf("Labels unique: %d \n", length(L)))

toc()

# typeof(labels)
# for ( i in labels) {
#     cat(sprintf("%d, ", i))
# }
# cat(sprintf("\n"))


# time and run wilcox.test
tic("R builtin")
Rwilcox <- matrix(, ncol = ncols, nrow = length(L) )
for ( gene in 1:ncol(input) ) {
    x <- matrix(as.vector(input[, gene]), ncol=1)
    # cat(sprintf("gene %d \n", gene))
    # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    i <- 1
    for ( c in L ) {
        lab <- labels %in% c

        v <- wilcox.test(x ~ lab, alternative="two.sided", correct=TRUE)$p.value
        # cat(sprintf("R wilcox %f\n", v))
        Rwilcox[i, gene] <- v
        i <- i + 1
    }
}
toc()

tic("Limma")
# time and run LIMMA
Limmawilcox <- matrix(, ncol = ncols, nrow = length(L) )
for ( gene in 1:ncol(input) ) {
    x <- matrix(as.vector(input[, gene]), ncol=1)
    # cat(sprintf("gene %d \n", gene))
    # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    
    i <- 1
    for ( c in L ) {
        ind <- which(labels %in% c)
        
        ## two sided
        v <- min(2 * min(limma::rankSumTestWithCorrelation(index = ind, statistics = x)), 1)  # two sides.
        
        ## less
        # v <- limma::rankSumTestWithCorrelation(index = ind, statistics = x)[1]  # less, left tail
        ## greater
        # v <- limma::rankSumTestWithCorrelation(index = ind, statistics = x)[2]    # greater, right tail
        # cat(sprintf("limma %f\n", v))
        Limmawilcox[i, gene] <- v
        i <- i + 1
    }
}
toc()
# warnings()





tic("fastde")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdewilcox <- fastde::wmwfast(as.matrix(input), labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
toc()


tic("fastde_df")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdewilcox_df <- fastde::wmwfast(as.matrix(input), labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(4))
toc()

# print(fastdewilcox_df)
# fastdewilcox_df$p_val
# fastdewilcox_df$cluster
# fastdewilcox_df$genes



# tic("BioQC2")
# # time and run BioQC
# BioQCwilcox2 <- matrix(, ncol = ncols, nrow = length(L) )
# inds = list()
# for ( c in L ) {
#     inds[[c]] <- labels %in% c
# }
# # cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
# BioQCwilcox2 <- BioQC::wmwTest(as.matrix(input), inds, valType = "p.two.sided")
# toc()



# tic("BioQC3")
# # time and run BioQC
# BioQCwilcox3 <- matrix(, ncol = ncols, nrow = length(L) )
# inds = list()
# for ( c in L ) {
#     inds[[c]] <- labels %in% c
#     cat(sprintf("label %d:  %d\n", c, length(which(inds[[c]]))))
# }
# cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
# for ( gene in 1:ncol(input) ) {
#     x <- matrix(as.vector(input[, gene]), ncol=1)
    
#     BioQCwilcox3[, gene] <- BioQC::wmwTest(x, inds, valType = "p.two.sided")
# }
# toc()





# tic("BioQC")
# # time and run BioQC
# BioQCwilcox <- matrix(, ncol = ncols, nrow = length(L) )
# for ( gene in 1:ncol(input) ) {
#     x <- matrix(as.vector(input[, gene]), ncol=1)
#     # cat(sprintf("gene %d \n", gene))
#     # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    
#     i <- 1
#     for ( c in L ) {
#         ind <- which(labels %in% c)
        
#         v <- BioQC::wmwTest(x, ind, valType = "p.two.sided")
#         # cat(sprintf("bioQC %f\n", v))
#         BioQCwilcox[i, gene] <- v
#         i <- i + 1
#     }
# }
# toc()

# BioQCwilcox



# fastdewilcox[, 1]
# BioQCwilcox2[, 1]
# Limmawilcox[, 1]
# Rwilcox[, 1]

# fastdewilcox[1, ]
# BioQCwilcox2[1, ]
# Limmawilcox[1, ]
# Rwilcox[1, ]


# # time and run wilcox.test
# tic("R builtin 2sided")
# Rwilcox2 <- matrix(, ncol = ncols, nrow = length(L) )
# for ( gene in 1:ncol(input) ) {
#     x <- matrix(as.vector(input[, gene]), ncol=1)
#     # cat(sprintf("gene %d \n", gene))
#     # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
#     i <- 1
#     for ( c in L ) {
#         lab <- labels %in% c

#         v <- wilcox.test(x ~ lab, alternative="two.sided")$p.value
#         # cat(sprintf("R wilcox %f\n", v))
#         Rwilcox2[i, gene] <- v
#         i <- i + 1
#     }
# }
# toc()

# Rwilcox2[, 1]


## compare by calculating the residuals.
comparemat("R vs fastde", Rwilcox, fastdewilcox)
# comparemat("bioqc vs fastde", BioQCwilcox2, fastdewilcox)
comparemat("limma vs fastde", Limmawilcox, fastdewilcox)



#==================


tic("fastde")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdewilcox <- fastde::sparsewmwfast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
toc()


tic("fastde_df")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdewilcox_df <- fastde::sparsewmwfast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(4))
toc()

# print(fastdewilcox_df)
# fastdewilcox_df$p_val
# fastdewilcox_df$cluster
# fastdewilcox_df$genes

# fastdewilcox[, 1]
# BioQCwilcox2[, 1]
# Limmawilcox[, 1]
# Rwilcox[, 1]

# fastdewilcox[1, ]
# BioQCwilcox2[1, ]
# Limmawilcox[1, ]
# Rwilcox[1, ]


## compare by calculating the residuals.
comparemat("R vs sparse fastde", Rwilcox, fastdewilcox)
# comparemat("bioqc vs sparse fastde", BioQCwilcox2, fastdewilcox)
comparemat("limma vs sparse fastde", Limmawilcox, fastdewilcox)




