
library(fastde)
library(tictoc)
# make a random matrix

library(Matrix)

# rows = 15052
# cols = 38000
rows = 15052
cols = 28000

M <- rsparsematrix(rows, cols, 0.04)
rownames(M) <- paste0("r", 1:rows)
colnames(M) <- paste0("c", 1:cols)
# M <- rsparsematrix(100, 150, 0.05)
# rownames(M) <- paste0("r", 1:100)
# colnames(M) <- paste0("c", 1:150)
str(M)

tic("SPARSE MAT TRANSPOSE")
tM <- fastde::sp_transpose(M)
toc()
str(tM)

ttM <- fastde::sp_transpose(tM)
str(ttM)

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
rm(dM)
rm(cdM)

tic("R convert sparse to dense then transpose")
tdM <- t(as.matrix(M))
toc()

tic("fastde convert sparse to dense then transpose")
ctdM <- t(fastde::sp_to_dense(M))
toc()

print(all.equal(tdM, ctdM))
rm(ctdM)

tic("fastde transpose sparse then convert to dense")
cdtM <- fastde::sp_to_dense(fastde::sp_transpose(M))
toc()

print(all.equal(tdM, cdtM))
rm(tdM)
rm(cdtM)

rm(M)
