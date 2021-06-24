
library(fastde)
library(tictoc)
# make a random matrix

library(Matrix)
M <- Matrix(c(0, 0,  0, 2,
              6, 0, -1, 5,
              0, 4,  3, 0,
              0, 0,  5, 0),
            byrow = TRUE, nrow = 4, sparse = TRUE)
rownames(M) <- paste0("r", 1:4)
colnames(M) <- paste0("c", 1:4)
M
str(M)

tic("SPARSE MAT TRANSPOSE")
tM <- fastde::sp_transpose(M)
toc()
tM
str(tM)

ttM <- fastde::sp_transpose(tM)
ttM
str(ttM)

print(all.equal(M, ttM))

