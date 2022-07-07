library(Matrix)
library(spam)
library(spam64)
library(fastde)

print("small sparse mat")
# small test.
x <- rsparsematrix(1000, 2000, 0.05)
test_x <- rttest_dgCMatrix(x)  # testd okay.
all.equal(x, test_x)

print("convert to spam32")
dx <- as.matrix(x)
sx <- as.spam(dx)
class(sx@dimension)
rm(dx)
test_x <- rttest_spam32(sx)
class(test_x@dimension)
all.equal(sx@dimension, test_x@dimension)
all.equal(sx@entries, test_x@entries)
all.equal(sx@colindices, test_x@colindices)
all.equal(sx@rowpointers, test_x@rowpointers)
rm(sx)
rm(x)
rm(test_x)

print("test large dgCMatrix")
# test create a dgCMatrix
x <- rsparsematrix(2^30, 1, 0.05)
test_x <- rttest_dgCMatrix(x)  # testd okay.
all.equal(x, test_x)
rm(x)
rm(test_x)
# dx <- as.matrix(x)
# sx <- as.spam(dx)

print("test large spam32")
# test create a spam object
y <- spam_random(2^30, 1, 0.05)
test_y <- rttest_spam32(y)
all.equal(y@dimension, test_y@dimension)
all.equal(y@entries, test_y@entries)
all.equal(y@colindices, test_y@colindices)
all.equal(y@rowpointers, test_y@rowpointers)
rm(y)
rm(test_y)


print("test large spam64")
# test create a spam64 object
z <- spam_random(2^31, 1, 0.05)
test_z <- rttest_spam64(z)
all.equal(z@dimension, test_z@dimension)
all.equal(z@entries, test_z@entries)
all.equal(z@colindices, test_z@colindices)
all.equal(z@rowpointers, test_z@rowpointers)
rm(z)
rm(test_z)

print("test large spam64 2")
# test create a spam64 object
w <- spam_random(2^30, 2, 0.05)
test_w <- rttest_spam64(w)
all.equal(w@dimension, test_w@dimension)
all.equal(w@entries, test_w@entries)
all.equal(w@colindices, test_w@colindices)
all.equal(w@rowpointers, test_w@rowpointers)
rm(w)
rm(test_w)




