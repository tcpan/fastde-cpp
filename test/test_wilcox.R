library(BioQC)
library(rhdf5)
library(fastde)
library(tictoc)

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
input <- h5read("/home/tpan/build/wave/input.h5", "array/block0_values")
labels_all <- as.vector(h5read("/home/tpan/build/wave/labels.h5", "array/block0_values"))
genenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis1")
samplenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis0")
# wilcox <- h5read("/home/tpan/build/wave/test-wilcox.h5", "array/block0_values")


colnames(input) <- genenames
rownames(input) <- samplenames


# count number of labels
#num_labels = nlevels(labels)

# cat(sprintf("test size:  r %d X c %d.\n", nrow(wilcox), ncol(wilcox)))
cat(sprintf("input size:  r %d X c %d\n", nrow(input), ncol(input)))

labels <- labels_all[1:nrow(input)]
cat(sprintf("Labels: %d \n", length(labels)))

L <- unique(sort(labels))

cat(sprintf("Labels unique: %d \n", length(L)))

# typeof(labels)
# for ( i in labels) {
#     cat(sprintf("%d, ", i))
# }
# cat(sprintf("\n"))

# time and run BioQC
cat(sprintf("warm up\n"));
fastdewilcox <- fastde::wmw_fast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))

# # print(fastdewilcox_df)
# # fastdewilcox_df$p_val
# # fastdewilcox_df$cluster
# # fastdewilcox_df$genes
# comparemat("c++ vs fastde tstat2", wilcox, fastdewilcox)


tic("fastde")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdewilcox <- fastde::wmw_fast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
toc()


tic("fastde_df")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdewilcox_df <- fastde::wmw_fast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(4))
toc()

# print(fastdewilcox2_df)



# tic("BioQC2")
# # time and run BioQC
# BioQCwilcox2 <- matrix(, ncol = ncol(input), nrow = length(L) )
# inds = list()
# for ( c in L ) {
#     inds[[c]] <- labels %in% c
# }
# # cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
# BioQCwilcox2 <- BioQC::wmwTest(input, inds, valType = "p.two.sided")
# toc()



# tic("BioQC3")
# # time and run BioQC
# BioQCwilcox3 <- matrix(, ncol = ncol(input), nrow = length(L) )
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




tic("Limma")
# time and run LIMMA
Limmawilcox <- matrix(, ncol = ncol(input), nrow = length(L) )
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



# tic("BioQC")
# # time and run BioQC
# BioQCwilcox <- matrix(, ncol = ncol(input), nrow = length(L) )
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
Rwilcox <- matrix(, ncol = ncol(input), nrow = length(L) )
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


# fastdewilcox[, 1]
# BioQCwilcox2[, 1]
# Limmawilcox[, 1]
# Rwilcox[, 1]

# fastdewilcox[1, ]
# BioQCwilcox2[1, ]
# Limmawilcox[1, ]
# Rwilcox[1, ]

# comparemat("c++ vs R wilcox", wilcox, Rwilcox)
# comparemat("c++ vs fastde wilcox", wilcox, fastdewilcox)
comparemat("c++ vs fastde wilcox2", Rwilcox, fastdewilcox)
# comparemat("c++ vs BioQC wilcox", Rwilcox, BioQCwilcox2)
comparemat("c++ vs Limma wilcox", Rwilcox, Limmawilcox)


# # time and run wilcox.test
# tic("R builtin 2sided")
# Rwilcox2 <- matrix(, ncol = ncol(input), nrow = length(L) )
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


# ## compare by calculating the residuals.
# res = Rwilcox - fastdewilcox
# residual = sqrt(mean(res * res))

# cat(sprintf("R naive vs fastde residual tail = %f\n", residual))

# res = BioQCwilcox2[, 0:5] - fastdewilcox[, 0:5]
# residual = sqrt(mean(res * res))

# cat(sprintf("BioQC vs fastde residual = %f\n", residual))

# res = BioQCwilcox2[, 0:5] - Rwilcox[, 0:5]
# residual = sqrt(mean(res * res))

# cat(sprintf("BioQC vs R naive residual = %f\n", residual))


# res = Limmawilcox - Rwilcox
# residual = sqrt(mean(res * res))

# cat(sprintf("Limma vs R naive residual = %f\n", residual))

# res = Limmawilcox - fastdewilcox
# residual = sqrt(mean(res * res))

# cat(sprintf("Limma vs fastde residual = %f\n", residual))

