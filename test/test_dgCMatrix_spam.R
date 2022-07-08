library(Matrix)
library(dotCall64)
library(spam)
library(spam64)
library(fastde)

print("small sparse mat")
# small test.
x <- rsparsematrix(1000, 2000, 0.05)
test_x <- rttest_dgCMatrix(x)  # testd okay.
all.equal(x, test_x)
rm(test_x)

dx <- as.matrix(x)

print("test small spam from R")
sx <- as.spam(dx)
class(sx)
class(sx@dimension)
test_x <- rttest_spam32(sx)
class(test_x)
class(test_x@dimension)
all.equal(sx, test_x)
rm(sx)
rm(test_x)

print("test small spamx from R")
sx <- as.spamx(dx)
class(sx)
class(sx@dimension)
test_x <- rttest_spamx32(sx)
class(test_x)
class(test_x@dimension)
all.equal(sx, test_x)
rm(sx)
rm(test_x)

rm(dx)
rm(x)


print("test small spam32")
options("spam.force64" = FALSE)

# test create a spam object
y <- spam(1, 1000, 2000)
class(y)
class(y@dimension)
test_y <- rttest_spam32(y)
class(test_y)
class(test_y@dimension)
all.equal(y, test_y)
rm(test_y)

z <- as.spamx(y)
class(z)
class(z@dimension)
test_z <- rttest_spamx32(z)
class(test_z)
class(test_z@dimension)
all.equal(z, test_z)
rm(z)
rm(test_z)
rm(y)

y <- spamx(1, 1000, 2000)
class(y)
class(y@dimension)
test_y <- rttest_spamx32(y)
class(test_y)
class(test_y@dimension)
all.equal(y, test_y)
rm(y)
rm(test_y)


print("test small spam64")
options("spam.force64" = TRUE)
y <- spam(1, 1000, 2000)
class(y)
class(y@dimension)
test_y <- rttest_spam64(y)
class(test_y)
class(test_y@dimension)
all.equal(y, test_y)
rm(test_y)

z <- as.spamx(y)
class(z)
class(z@dimension)
test_z <- rttest_spamx64(z)
class(test_z)
class(test_z@dimension)
all.equal(z, test_z)
rm(z)
rm(test_z)
rm(y)

y <- spamx(1, 1000, 2000)
class(y)
class(y@dimension)
test_y <- rttest_spamx64(y)
class(test_y)
class(test_y@dimension)
all.equal(y, test_y)
rm(y)
rm(test_y)


options("spam.force64" = TRUE)

print("test large spamx32")
# test create a spam object
y <- spam(1, 2^30, 1)
test_y <- rttest_spam32(y)
all.equal(y, test_y)
rm(test_y)

z <- as.spamx(y)
test_z <- rttest_spamx32(z)
all.equal(z, test_z)
rm(z)
rm(test_z)
rm(y)


print("test large spamx64")
# test create a spamx64 object
y <- spam(1, 2^31, 1)
test_y <- rttest_spam64(y)
all.equal(y, test_y)
rm(test_y)

z <- as.spamx(y)
test_z <- rttest_spamx64(z)
all.equal(z, test_z)
rm(z)
rm(test_z)
rm(y)



print("test large spamx64 2")
# test create a spamx64 object
y <- spamx(1, 2^30, 2)
test_y <- rttest_spam64(y)
all.equal(y, test_y)
rm(test_y)

z <- as.spamx(y)
test_z <- rttest_spamx64(z)
all.equal(z, test_z)
rm(z)
rm(test_z)
rm(y)



