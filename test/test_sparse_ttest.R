library(rhdf5)
library(fastde)
library(tictoc)
library(Matrix)

print_mat <- function(mat, count) {
    print(head(mat, n = count)[, 1:count])
    print("...")
    print(tail(mat, n = count)[, 1:count])
}

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
# tic("read")
# labels <- as.vector(h5read("/home/tpan/build/wave/labels.h5", "array/block0_values"))
# genenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis1")
# samplenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis0")
# wilcox <- h5read("/home/tpan/build/wave/test-wilcox.h5", "array/block0_values")

# input <- rsparsematrix(length(samplenames), length(genenames), 0.1)

# colnames(input) <- genenames
# rownames(input) <- samplenames
# toc()

nthreads <- as.integer(1)

tic("generate")
# max is 2B.
ncols = 1000  # features  - test with a small number. this is totally parallel
nrows = 10000  # samples
nclusters = 15

clusters = 1:nclusters
labels <- as.integer(sample(clusters, nrows, replace = TRUE))

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
print_mat(Rmean, 5)


# print(Rttest)
# fastdettest_df$p_val
# fastdettest_df$cluster
# fastdettest_df$genes

# Rmean[, 1]

# Rmean[1, ]

tic("fastde t stat")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdetval <- fastde::ttest_fast(as.matrix(input), labels, as_dataframe = FALSE, threads = nthreads, alternative =  as.integer(3), var_equal = FALSE)
toc()

comparemat("R vs dense tstat", Rmean, fastdetval)
print_mat(fastdetval, 5)

tic("sparse fastde t stat")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdesparsetval <- fastde::sparse_ttest_fast(input, labels, features_as_rows = FALSE, as_dataframe = FALSE, threads = nthreads, alternative = as.integer(3), var_equal = FALSE)
toc()

comparemat("R vs sparse tstat", Rmean, fastdesparsetval)
print_mat(fastdesparsetval, 5)




tic("R builtin, dof")
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
        Rmean[i, gene] <- v$parameter
        i <- i + 1
    }
}
toc()

print(typeof(Rmean))
print(dim(Rmean))
print_mat(Rmean, 5)

# print(Rttest)
# fastdettest_df$p_val
# fastdettest_df$cluster
# fastdettest_df$genes

# Rmean[, 1]

# Rmean[1, ]

tic("fastde t dof")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdetval <- fastde::ttest_fast(as.matrix(input), labels, as_dataframe = FALSE, threads = nthreads, alternative =  as.integer(4), var_equal = FALSE)
toc()

comparemat("R vs dense dof", Rmean, fastdetval)
print_mat(fastdetval, 5)


tic("sparse fastde t dof")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdesparsetval <- fastde::sparse_ttest_fast(input, labels, features_as_rows = FALSE,  as_dataframe = FALSE, threads = nthreads, alternative = as.integer(4), var_equal = FALSE)
toc()

comparemat("R vs sparse dof", Rmean, fastdesparsetval)
print_mat(fastdesparsetval, 5)






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
print_mat(Rttest, 5)



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
print_mat(Rttest2, 5)




tic("fastde ttest")
# time and run ttest two sided (2)
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdettest <- fastde::ttest_fast(as.matrix(input), labels, as_dataframe = FALSE, threads = nthreads, alternative =  as.integer(2), var_equal = FALSE)
toc()


tic("fastde_df ttest")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdettest_df <- fastde::ttest_fast(as.matrix(input), labels,  as_dataframe = TRUE, threads = nthreads, alternative =  as.integer(2), var_equal = FALSE)
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
comparemat("R vs dense ttest", Rttest, fastdettest)
print_mat(fastdettest, 5)



tic("sparse fastde")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdesparsettest <- fastde::sparse_ttest_fast(input, labels,features_as_rows = FALSE,  as_dataframe = FALSE, threads = nthreads, alternative =  as.integer(2), var_equal = FALSE)
toc()


tic("sparse fastde_df")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdesparsettest_df <- fastde::sparse_ttest_fast(input, labels, features_as_rows = FALSE,  as_dataframe = TRUE, threads = nthreads, alternative =  as.integer(2), var_equal = FALSE)
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
comparemat("R vs sparse ttest", Rttest, fastdesparsettest)
print_mat(fastdesparsettest, 5)



print("Convert to sparse 64bit")

input64 <- as.dgCMatrix64(input)

tic("sparse64 fastde")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input64), ncol(input64)))
fastdesparsettest <- fastde::sparse_ttest_fast(input64, labels, features_as_rows = FALSE, as_dataframe = FALSE, threads = nthreads, alternative =  as.integer(2), var_equal = FALSE)
toc()


tic("sparse64 fastde_df")
# time and run BioQC
cat(sprintf("input %d X %d\n", nrow(input64), ncol(input64)))
fastdesparsettest_df <- fastde::sparse_ttest_fast(input64, labels, features_as_rows = FALSE,  as_dataframe = TRUE, threads = nthreads, alternative =  as.integer(2), var_equal = FALSE)
toc()


## compare by calculating the residuals.
comparemat("R vs sparse64 ttest", Rttest, fastdesparsettest)
print_mat(fastdesparsettest, 5)
