
library(fastde)
library(tictoc)
#library(jointprof)
# make a random matrix

library(Matrix)

# rows = 15052
# cols = 38000
rows = 15052
cols = 2800

M <- rsparsematrix(rows, cols, 0.04)
rownames(M) <- paste0("r", 1:rows)
colnames(M) <- paste0("c", 1:cols)
# M <- rsparsematrix(100, 150, 0.05)
# rownames(M) <- paste0("r", 1:100)
# colnames(M) <- paste0("c", 1:150)
# str(M)

tic("rcpp overhead")
M2 <- rttest_dgCMatrix(M)
toc()
rm(M2)

tic("fastde SPARSE MAT TRANSPOSE sexp")
tM3 <- fastde::rc_sp_transpose(M)
toc()

tic("fastde SPARSE MAT TRANSPOSE")
tM <- fastde::sp_transpose(M)
toc()
# str(tM)

tic("R sparse mat transpose")
tM2 <- t(M)
toc()
# str(tM2)
print(all.equal(tM2, tM))
print(all.equal(tM2, tM3))
rm(tM2)
rm(tM3)

ttM <- fastde::sp_transpose(tM)
# str(ttM)
print(all.equal(M, ttM))

rm(tM)
rm(ttM)


tic("R convert sparse to dense")
dM <- as.matrix(M)
toc()



tic("fastde convert sparse to dense")
cdM <- fastde::sp_to_dense(M)
toc()
print(all.equal(dM, cdM))
rm(cdM)

out_file <- "sp_to_dense.prof" #tempfile("fastde_sp_to_dense", fileext = ".out")
print(out_file)
fastde::start_profiler(out_file)
for (i in 1:20) cdM2 <- fastde::sp_to_dense(M)
fastde::stop_profiler()

out_file <- "rc_sp_to_dense.prof" #tempfile("fastde_sp_to_dense", fileext = ".out")
print(out_file)
fastde::start_profiler(out_file)
for (i in 1:20) cdM2 <- fastde::rc_sp_to_dense(M)
fastde::stop_profiler()
rm(cdM2)

# tic("fastde convert sparse to dense 2")
# cdM4 <- fastde::sp_to_dense2(M)
# toc()
# print(all.equal(dM, cdM4))
# rm(cdM4)

tic("fastde convert sparse to dense sexp")
cdM2 <- fastde::rc_sp_to_dense(M)
toc()
print(all.equal(dM, cdM2))
rm(cdM2)

# str(dM)
# str(cdM)
rm(dM)

tic("R convert sparse to dense then transpose")
tdM <- t(as.matrix(M))
toc()


tic("fastde convert sparse to dense then R transpose")
ctdM <- t(fastde::sp_to_dense(M))
toc()

print(all.equal(tdM, ctdM))
rm(ctdM)

tic("R transpose than convert to dense")
tdM <- as.matrix(t(M))
toc()


tic("fastde transpose sparse then convert to dense")
cdtM <- fastde::sp_to_dense(fastde::sp_transpose(M))
toc()
# str(tdM)
# str(cdtM)
print(all.equal(tdM, cdtM))
rm(cdtM)

tic("fastde to transposed dense")
cdtM2 <- fastde::sp_to_dense_transposed(M)
toc()

print(all.equal(tdM, cdtM2))
rm(cdtM2)


tic("fastde to transposed dense sexp")
cdtM3 <- fastde::rc_sp_to_dense_transposed(M)
toc()

# str(tdM)
# str(cdtM3)
print(all.equal(tdM, cdtM3))
rm(cdtM3)


rm(tdM)

rm(M)
