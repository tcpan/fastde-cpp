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
    cat(sprintf("%s : diff range [%.17g, %.17g], median %.17g, mean %.17g, sd %.17g\n", 
        name, mindiff, maxdiff, mediandiff, meandiff, stdevdiff))

        # mxpos = which(diff == maxdiff, arr.ind = TRUE)
    # mnpos = which(diff == mindiff, arr.ind = TRUE)

    # if ( abs(maxdiff) > .Machine$double.eps)
    #    cat(sprintf("%s : max diff at pos %d:  A %.17g - B %.17g = DIFF %.17g.\n", 
    #         name, mxpos, A[mxpos], B[mxpos], diff[mxpos]))
    # if  ( abs(mindiff) > .Machine$double.eps)
    #     cat(sprintf("%s : min diff at pos %d:  A %.17g - B %.17g = DIFF %.17g.\n", 
    #         name, mnpos, A[mnpos], B[mnpos], diff[mnpos]))
    
}


# read input
input <- h5read("/home/tpan/data/gnw2000/gnw2000.h5", "array/block0_values")
labels_all <- as.vector(h5read("/home/tpan/data/gnw2000/gnw2000_truth.h5", "array/block0_values"))
genenames <- h5read("/home/tpan/data/gnw2000/gnw2000.h5", "array/axis1")
samplenames <- h5read("/home/tpan/data/gnw2000/gnw2000.h5", "array/axis0")
# wilcox <- h5read("/home/tpan/build/wave/test-wilcox.h5", "array/block0_values")


colnames(input) <- genenames
rownames(input) <- samplenames

input <- input[, 1:100]


# count number of labels
#num_labels = nlevels(labels)

# cat(sprintf("test size:  r %d X c %d.\n", nrow(wilcox), ncol(wilcox)))
cat(sprintf("input size:  r %d X c %d\n", nrow(input), ncol(input)))

labels <- as.integer(1 - labels_all[1:nrow(input)])  # inputs are 0 an 1.   wilcox treat class 0 as other in a formula.   fastde treat 0 as first class.   We need to flip this...
cat(sprintf("Labels: %d \n", length(labels)))

L <- unique(sort(labels))

cat(sprintf("Labels unique: %d \n", length(L)))

# typeof(labels)
# for ( i in labels) {
#     cat(sprintf("%d, ", i))
# }
# cat(sprintf("\n"))

# time and run BioQC

# # print(fastdewilcox_df)
# # fastdewilcox_df$p_val
# # fastdewilcox_df$cluster
# # fastdewilcox_df$genes
# comparemat("c++ vs fastde tstat2", wilcox, fastdewilcox)


tic("fastde")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdewilcox <- fastde::wmw_fast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
toc()

tic("fastde stat")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdewilcox_stat <- fastde::wmw_fast(input, labels, rtype=as.integer(3), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
toc()

tic("fastde stat")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdewilcox_z <- fastde::wmw_fast(input, labels, rtype=as.integer(4), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
toc()

tic("fastde_df")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdewilcox_df <- fastde::wmw_fast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(1))
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




# tic("Limma")
# # time and run LIMMA
# Limmawilcox <- matrix(, ncol = ncol(input), nrow = length(L) )
# for ( gene in 1:ncol(input) ) {
#     x <- matrix(as.vector(input[, gene]), ncol=1)
#     # cat(sprintf("gene %d \n", gene))
#     # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    
#     i <- 1
#     for ( c in L ) {
#         ind <- which(labels %in% c)
        
#         ## two sided
#         v <- min(2 * min(limma::rankSumTestWithCorrelation(index = ind, statistics = x)), 1)  # two sides.
        
#         ## less
#         # v <- limma::rankSumTestWithCorrelation(index = ind, statistics = x)[1]  # less, left tail
#         ## greater
#         # v <- limma::rankSumTestWithCorrelation(index = ind, statistics = x)[2]    # greater, right tail
#         # cat(sprintf("limma %f\n", v))
#         Limmawilcox[i, gene] <- v
#         i <- i + 1
#     }
# }
# toc()
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
Rwilcox_stat <- matrix(, ncol = ncol(input), nrow = length(L) )
Rwilcox_z <- matrix(, ncol = ncol(input), nrow = length(L) )
for ( gene in 1:ncol(input) ) {
    x <- matrix(as.vector(input[, gene]), ncol=1)
    # cat(sprintf("gene %d \n", gene))
    # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    i <- 1
    for ( c in L ) {
        lab <- labels %in% c

        xx <- x[which(labels == c)]
        yy <- x[which(labels != c)]

        # v <- wilcox.test(x ~ lab, alternative="two.sided", correct=TRUE)
        v <- wilcox.test(x = xx, y= yy, alternative="two.sided", correct=TRUE)
        # cat(sprintf("R wilcox %f\n", v))
        Rwilcox[i, gene] <- v$p.value
        # cat(sprintf("R wilcox %f\n", v))
        Rwilcox_stat[i, gene] <- v$statistic

        nx = length(xx)
        ny = length(yy)
        z <- v$statistic - nx * ny * 0.5
        r <- rank(c(xx, yy))
        NTIES <- table(r)
        NTIES2 <- NTIES[which(NTIES > 1)]
        SIGMA <- sqrt((nx * ny / 12) *
                          ((nx + ny + 1)
                           - sum(NTIES^3 - NTIES)
                           / ((nx + ny) * (nx + ny - 1))))
            CORRECTION <- 0 # sign(z) * 0.5
	    z <- (z - CORRECTION) / SIGMA
        Rwilcox_z[i, gene] <- z

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

# head(Rwilcox[, 11:20])
# head(fastdewilcox[ , 11:20])

# comparemat("c++ vs R wilcox", wilcox, Rwilcox)
# comparemat("c++ vs fastde wilcox", wilcox, fastdewilcox)
comparemat("R vs fastde wilcox2", Rwilcox, fastdewilcox)
# comparemat("c++ vs BioQC wilcox", Rwilcox, BioQCwilcox2)
# comparemat("R vs Limma wilcox", Rwilcox, Limmawilcox)


# head(Rwilcox_stat[, 11:20])
# head(fastdewilcox_stat[ , 11:20])

# comparemat("c++ vs R wilcox", wilcox, Rwilcox)
# comparemat("c++ vs fastde wilcox", wilcox, fastdewilcox)
comparemat("R vs fastde wilcox stat", Rwilcox_stat, fastdewilcox_stat)

# head(Rwilcox_z[, 11:20])
# head(fastdewilcox_z[ , 11:20])

# comparemat("c++ vs R wilcox", wilcox, Rwilcox)
# comparemat("c++ vs fastde wilcox", wilcox, fastdewilcox)
comparemat("R vs fastde wilcox z", Rwilcox_z, fastdewilcox_z)

# pos <- which((Rwilcox_z[1, ] - fastdewilcox_z[1, ]) != 0)
# Rwilcox_z[1, pos[1:10]]
# fastdewilcox_z[1, pos[1:10]]
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

