
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

print("TESTING dgCMatrix functions")



tic("R sparse mat transpose")
tM2 <- t(M)
toc()
rm(tM2)

tic("fastde SPARSE MAT TRANSPOSE sexp")
tM <- fastde::sp_transpose(M, threads = 1)
toc()

tic("R sparse mat transpose")
tM2 <- t(M)
toc()

tic("fastde SPARSE MAT TRANSPOSE sexp")
tM <- fastde::sp_transpose(M, threads = 1)
toc()

tic("R sparse mat transpose")
tM2 <- t(M)
toc()

tic("fastde SPARSE MAT TRANSPOSE sexp")
tM <- fastde::sp_transpose(M, threads = 8)
toc()
# str(tM2)
print("R and c++ transpose the same?")
print(all.equal(tM2, tM))
rm(tM2)

ttM <- fastde::sp_transpose(tM, threads = 8)

# str(ttM)
print("c++ transpose invertible?")
print(all.equal(M, ttM))
rm(tM)
rm(ttM)


tic("R convert sparse to dense")
dM <- as.matrix(M)
toc()

tic("fastde convert sparse to dense")
cdM <- fastde::sp_to_dense(M, threads = 8)
toc()
print("c++ to dense the same as R?")
print(all.equal(dM, cdM))
rm(dM)
rm(cdM)


tic("R convert sparse to dense transposed")
dM <- as.matrix(t(M))
toc()

tic("fastde convert sparse to dense transposed")
cdM <- fastde::sp_to_dense_transposed(M, threads = 8)
toc()
print("c++ to dense transposed the same as R?")
print(all.equal(dM, cdM))
rm(dM)
rm(cdM)


# out_file <- "sp_to_dense.prof" #tempfile("fastde_sp_to_dense", fileext = ".out")
# print(out_file)
# fastde::start_profiler(out_file)
# for (i in 1:20) cdM2 <- fastde::sp_to_dense(M)
# fastde::stop_profiler()

# out_file <- "rc_sp_to_dense.prof" #tempfile("fastde_sp_to_dense", fileext = ".out")
# print(out_file)
# fastde::start_profiler(out_file)
# for (i in 1:20) cdM2 <- fastde::rc_sp_to_dense(M)
# fastde::stop_profiler()
# rm(cdM2)

# tic("fastde convert sparse to dense 2")
# cdM4 <- fastde::sp_to_dense2(M)
# toc()
# print(all.equal(dM, cdM4))
# rm(cdM4)



tic("convert to 64")
M64 <- as.dgCMatrix64(M)
toc()


print("TESTING dgCMatrix64 functions")

tic("R convert sparse to dense")
dM <- as.matrix(M)
toc()

tic("fastde convert sparse to dense 64")
cdM <- fastde::sp_to_dense(M64, threads = 8)
toc()
print("c++ to dense 64 the same as R?")
print(all.equal(dM, cdM))
rm(dM)
rm(cdM)

tic("R sparse mat transpose")
tM642 <- as.dgCMatrix64(t(M))
toc()

tic("fastde SPARSE 64 MAT TRANSPOSE sexp")
tM64 <- fastde::sp_transpose(M64, threads = 8)
toc()

# str(tM2)
print("R and c++ transpose 64 the same?")
print(all.equal(fastde::sp_to_dense(tM642), fastde::sp_to_dense(tM64)))
rm(tM642)

ttM64 <- fastde::sp_transpose(tM64, threads = 8)

# str(ttM)
print("c++ transpose 64 invertible?")
print(all.equal(fastde::sp_to_dense(M64), fastde::sp_to_dense(ttM64)))
rm(tM64)
rm(ttM64)




tic("R convert sparse to dense transposed")
dM <- as.matrix(t(M))
toc()

tic("fastde convert sparse to dense transposed 64")
cdM64 <- fastde::sp_to_dense_transposed(M64, threads = 8)
toc()
print("c++ to dense transposed 64 the same as R?")
print(all.equal(dM, cdM64))
rm(dM)
rm(cdM64)

# test cbind and rbind that can produce output that's 64bit.

Ms  <- list(M, M, M)

rMR <- rbind(M, M, M)
cMR <- cbind(M, M, M)



tic("R rbind")
rMR <- rbind(M, M, M)
toc()

tic("fastde rbind")
rM <- fastde::sp_rbind(Ms, 1)
toc()

identical(rMR@x, rM@x)
identical(rMR@i, rM@i)
identical(rMR@p, rM@p)



tic("fastde rbind 4t")
rM4 <- fastde::sp_rbind(Ms, 4)
toc()


identical(rMR@x, rM4@x)
identical(rMR@i, rM4@i)
identical(rMR@p, rM4@p)



tic("R cbind")
cMR <- cbind(M, M, M)
toc()

tic("fastde cbind")
cM <- fastde::sp_cbind(Ms, 1)
toc()

identical(cMR@x, cM@x)
identical(cMR@i, cM@i)
identical(cMR@p, cM@p)

tic("fastde cbind 4t")
cM4 <- fastde::sp_cbind(Ms, 4)
toc()

identical(cMR@x, cM4@x)
identical(cMR@i, cM4@i)
identical(cMR@p, cM4@p)
