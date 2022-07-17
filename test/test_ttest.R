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
#wilcox <- h5read("/home/tpan/build/wave/test-wilcox.h5", "array/block0_values")

nthreads <- as.integer(1)

colnames(input) <- genenames
rownames(input) <- samplenames

# count number of labels
#num_labels = nlevels(labels)

#cat(sprintf("test size:  r %d X c %d.\n", nrow(wilcox), ncol(wilcox)))
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

# two sides = 2.

cat(sprintf("warm up\n"));
fastdettest2 <- fastde::ttest_fast(input, labels, as_dataframe = FALSE, threads = nthreads, alternative = as.integer(2), var_equal = FALSE)

tic("fastde dof")
# time and run dense test not as dataframe.
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdettest2 <- fastde::ttest_fast(input, labels, as_dataframe = FALSE, threads = nthreads, alternative = as.integer(4), var_equal = FALSE)
toc()

tic("fastde dof df")
# time and run dense ttest
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdettest_df2 <- fastde::ttest_fast(input, labels, as_dataframe = TRUE, threads = nthreads, alternative = as.integer(4), var_equal = FALSE)
toc()


# time and run ttest.test
tic("R builtin")
Rttest <- matrix(, ncol = ncol(fastdettest2), nrow = nrow(fastdettest2) )
for ( gene in 1:ncol(input) ) {
    dat <- as.vector(input[, gene])
    # if ( (gene %% 100) == 0 ) cat(sprintf("R building t.test.  gene %d \n", gene))
    # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    i <- 1
    for ( c in L ) {
        lab <- labels %in% c

        v <- t.test(x = dat[lab], y = dat[!lab])$parameter
        # cat(sprintf("R wilcox %f\n", v))
        Rttest[i, gene] <- v
        i <- i + 1
    }
}
toc()


# Rttest[, 1]

# Rttest[1, ]

comparemat("R vs fastde df", Rttest, fastdettest2)


## compare by calculating the residuals.

res2 = Rttest - fastdettest2
residual2 = sqrt(mean(res2 * res2))

cat(sprintf("R naive vs fastde2 residual tail = %f\n", residual2))





tic("fastde_stat")
# time and run dense ttest
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdettest2 <- fastde::ttest_fast(input, labels, as_dataframe = FALSE, threads = nthreads, alternative = as.integer(3), var_equal = FALSE)
toc()


tic("fastde_stat")
# time and run dense ttest
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdettest_df2 <- fastde::ttest_fast(input, labels, as_dataframe = TRUE, threads = nthreads, alternative = as.integer(3), var_equal = FALSE)
toc()


# time and run ttest.test
tic("R builtin")
Rttest <- matrix(, ncol = ncol(fastdettest2), nrow = nrow(fastdettest2) )
for ( gene in 1:ncol(input) ) {
    dat <- as.vector(input[, gene])
    # if ( (gene %% 100) == 0 ) cat(sprintf("R building t.test.  gene %d \n", gene))
    # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    i <- 1
    for ( c in L ) {
        lab <- labels %in% c

        v <- t.test(x = dat[lab], y = dat[!lab])$statistic
        # cat(sprintf("R wilcox %f\n", v))
        Rttest[i, gene] <- v
        i <- i + 1
    }
}
toc()


# Rttest[, 1]

# Rttest[1, ]

comparemat("R vs fastde tstat", Rttest, fastdettest2)


## compare by calculating the residuals.

res2 = Rttest - fastdettest2
residual2 = sqrt(mean(res2 * res2))

cat(sprintf("R naive vs fastde2 residual tail = %f\n", residual2))



tic("fastde_pval")
# time and run dense ttest
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdettest2 <- fastde::ttest_fast(input, labels, as_dataframe = FALSE, threads = nthreads, alternative = as.integer(2), var_equal = FALSE)
toc()

tic("fastde_pval")
# time and run dense ttest
cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
fastdettest_df2 <- fastde::ttest_fast(input, labels, as_dataframe = TRUE, threads = nthreads, alternative = as.integer(2), var_equal = FALSE)
toc()


# time and run ttest.test
tic("R builtin")
Rttest <- matrix(, ncol = ncol(fastdettest2), nrow = nrow(fastdettest2) )
for ( gene in 1:ncol(input) ) {
    dat <- as.vector(input[, gene])
    # if ( (gene %% 100) == 0 ) cat(sprintf("R building t.test.  gene %d \n", gene))
    # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    i <- 1
    for ( c in L ) {
        lab <- labels %in% c

        v <- t.test(x = dat[lab], y = dat[!lab])$p.value
        # cat(sprintf("R wilcox %f\n", v))
        Rttest[i, gene] <- v
        i <- i + 1
    }
}
toc()


# Rttest[, 1]

# Rttest[1, ]

comparemat("R vs fastde pval", Rttest, fastdettest2)


## compare by calculating the residuals.

res2 = Rttest - fastdettest2
residual2 = sqrt(mean(res2 * res2))

cat(sprintf("R naive vs fastde2 residual tail = %f\n", residual2))
