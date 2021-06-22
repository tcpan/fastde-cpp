library(BioQC)

library(rhdf5)

library(fastde)

library(tictoc)

# read input
input <- h5read("/home/tpan/build/wave/input.h5", "array/block0_values")
labels <- as.vector(h5read("/home/tpan/build/wave/labels.h5", "array/block0_values"))
genenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis1")
samplenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis0")
wilcox <- h5read("/home/tpan/build/wave/test-wilcox.h5", "array/block0_values")


colnames(input) <- genenames
rownames(input) <- samplenames

# count number of labels
#num_labels = nlevels(labels)

cat(sprintf("test size:  r %d X c %d.\n", nrow(wilcox), ncol(wilcox)))
cat(sprintf("input size:  r %d X c %d\n", nrow(input), ncol(input)))
cat(sprintf("Labels rows: %d \n", length(labels)))

L <- unique(sort(labels))

cat(sprintf("Labels unique: %d \n", length(L)))

# typeof(labels)
# for ( i in labels) {
#     cat(sprintf("%d, ", i))
# }
# cat(sprintf("\n"))

tic("fastde")
# time and run BioQC
fastdewilcox <- matrix(, ncol = ncol(wilcox), nrow = nrow(wilcox) )
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdewilcox <- fastde::wmwfast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
toc()

fastdewilcox[, 1]


tic("fastde_df")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdewilcox_df <- fastde::wmwfast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(4))
toc()

print(fastdewilcox_df)
# fastdewilcox_df$p_val
# fastdewilcox_df$cluster
# fastdewilcox_df$genes



tic("BioQC2")
# time and run BioQC
BioQCwilcox2 <- matrix(, ncol = ncol(wilcox), nrow = nrow(wilcox) )
inds = list()
for ( c in L ) {
    inds[[c]] <- labels %in% c
}
# cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
BioQCwilcox2 <- BioQC::wmwTest(input, inds, valType = "p.two.sided")
toc()

BioQCwilcox2[, 1]


# tic("BioQC3")
# # time and run BioQC
# BioQCwilcox3 <- matrix(, ncol = ncol(wilcox), nrow = nrow(wilcox) )
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

# BioQCwilcox3



tic("Limma")
# time and run LIMMA
Limmawilcox <- matrix(, ncol = ncol(wilcox), nrow = nrow(wilcox) )
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

Limmawilcox[, 1]
# warnings()



# tic("BioQC")
# # time and run BioQC
# BioQCwilcox <- matrix(, ncol = ncol(wilcox), nrow = nrow(wilcox) )
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


# time and run wilcox.test
tic("R builtin")
Rwilcox <- matrix(, ncol = ncol(wilcox), nrow = nrow(wilcox) )
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

Rwilcox[, 1]


# # time and run wilcox.test
# tic("R builtin 2sided")
# Rwilcox2 <- matrix(, ncol = ncol(wilcox), nrow = nrow(wilcox) )
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
res = Rwilcox - fastdewilcox
residual = sqrt(mean(res * res))

cat(sprintf("R naive vs fastde residual tail = %f\n", residual))

res = BioQCwilcox2[, 0:5] - fastdewilcox[, 0:5]
residual = sqrt(mean(res * res))

cat(sprintf("BioQC vs fastde residual = %f\n", residual))

# res = BioQCwilcox2[, 0:5] - Rwilcox[, 0:5]
# residual = sqrt(mean(res * res))

# cat(sprintf("BioQC vs R naive residual = %f\n", residual))


# res = Limmawilcox - Rwilcox
# residual = sqrt(mean(res * res))

# cat(sprintf("Limma vs R naive residual = %f\n", residual))

res = Limmawilcox - fastdewilcox
residual = sqrt(mean(res * res))

cat(sprintf("Limma vs fastde residual = %f\n", residual))

