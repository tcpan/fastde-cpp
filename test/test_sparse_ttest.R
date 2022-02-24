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
ncols = 2000  # features  - test with a small number. this is totally parallel
nrows = 2000  # samples
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


tic("R builtin, t val")
Rmean <- matrix(, ncol = ncols, nrow = length(L) )
dense_input = as.matrix(input);
print(typeof(input))
print(typeof(dense_input))
print(dim(input))
print(dim(dense_input))
for ( gene in 1:ncol(input) ) {
    dat <- matrix(as.vector(dense_input[, gene]), ncol = 1)
    # cat(sprintf("gene %d \n", gene))
    # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    i <- 1
    for ( c in L ) {
        lab <- labels %in% c

        # v <- t.test(x = dat[lab], y = dat[!lab])$p.value
        v <- t.test(x = dat[lab], y = dat[!lab])
        # print(v$estimate[['mean of x']])
        # cat(sprintf("R ttest %f\n", v))
        Rmean[i, gene] <- v$statistic
        i <- i + 1
    }
}
toc()

print(typeof(Rmean))
print(dim(Rmean))


# print(Rttest)
# fastdettest_df$p_val
# fastdettest_df$cluster
# fastdettest_df$genes

# Rmean[, 1]

# Rmean[1, ]

tic("fastde t stat")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdetval <- fastde::ttest_fast(as.matrix(input), labels, as_dataframe = FALSE, threads = as.integer(4), alternative =  as.integer(3), var_equal = FALSE)
toc()

comparemat("R vs sparse tstat", Rmean, fastdetval)


tic("sparse fastde t stat")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdesparsetval <- fastde::sparse_ttest_fast(input, labels, as_dataframe = FALSE, threads = as.integer(4), alternative = as.integer(3), var_equal = FALSE)
toc()

comparemat("R vs sparse tstat", Rmean, fastdesparsetval)






# time and run ttest.test
tic("R builtin, pval")
Rttest <- matrix(, ncol = ncols, nrow = length(L) )
dense_input = as.matrix(input);
print(typeof(input))
print(typeof(dense_input))
print(dim(input))
print(dim(dense_input))
for ( gene in 1:ncol(input) ) {
    dat <- matrix(as.vector(dense_input[, gene]), ncol = 1)
    # cat(sprintf("gene %d \n", gene))
    # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    i <- 1
    for ( c in L ) {
        lab <- labels %in% c

        v <- t.test(x = dat[lab], y = dat[!lab])$p.value
        # print(v$statistic)
        # cat(sprintf("R ttest %f\n", v))
        Rttest[i, gene] <- v
        i <- i + 1
    }
}
toc()

print(typeof(Rttest))
print(dim(Rttest))


# print(Rttest)
# fastdettest_df$p_val
# fastdettest_df$cluster
# fastdettest_df$genes

# Rttest[, 1]

# Rttest[1, ]

tic("R builtin sparse, pval")
Rttest2 <- matrix(, ncol = ncols, nrow = length(L) )
for ( gene in 1:ncol(input) ) {
    dat <- as.vector(input[, gene, drop = FALSE])
    # cat(sprintf("gene %d \n", gene))
    # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    i <- 1
    for ( c in L ) {
        lab <- labels %in% c

        v <- t.test(x = dat[lab], y = dat[!lab])$p.value
        # cat(sprintf("R ttest %f\n", v))
        Rttest2[i, gene] <- v
        i <- i + 1
    }
}
toc()

# Rttest2[, 1]

# Rttest2[1, ]
comparemat("R vs R sparse", Rttest, Rttest2)




tic("fastde ttest")
# time and run ttest two sided (2)
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdettest <- fastde::ttest_fast(as.matrix(input), labels, as_dataframe = FALSE, threads = as.integer(4), alternative =  as.integer(2), var_equal = FALSE)
toc()


tic("fastde_df ttest")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdettest_df <- fastde::ttest_fast(as.matrix(input), labels,  as_dataframe = TRUE, threads = as.integer(4), alternative =  as.integer(2), var_equal = FALSE)
toc()



# print(fastdettest_df)
# fastdettest_df$p_val
# fastdettest_df$cluster
# fastdettest_df$genes



# fastdettest[, 1]
# # Rttest[, 1]

# fastdettest[1, ]
# # Rttest[1, ]




## compare by calculating the residuals.
comparemat("R vs fastde", Rttest, fastdettest)



tic("sparse fastde")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdesparsettest <- fastde::sparse_ttest_fast(input, labels, as_dataframe = FALSE, threads = as.integer(4), alternative =  as.integer(2), var_equal = FALSE)
toc()


tic("sparse fastde_df")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdesparsettest_df <- fastde::sparse_ttest_fast(input, labels,  as_dataframe = TRUE, threads = as.integer(4), alternative =  as.integer(2), var_equal = FALSE)
toc()


# print(fastdesparsettest_df)
# fastdettest_df$p_val
# fastdettest_df$cluster
# fastdettest_df$genes



# fastdesparsettest[, 1]
# # Rttest[, 1]

# fastdesparsettest[1, ]
# # Rttest[1, ]




## compare by calculating the residuals.
comparemat("R vs fastde", Rttest, fastdesparsettest)

