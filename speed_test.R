library(xgboost)

set.seed(1)

n <- 10^3
train <- data.frame(x1=rnorm(n), x2=rnorm(n))
train$y <- 1*train$x1^2 + 2*train$x2^2 + runif(n)
train <- as.matrix(train)

n <- 5*10^5
test <- data.frame(x1=rnorm(n), x2=rnorm(n))
test <- as.matrix(test)


system.time(xgb.train(
  params=list(eta=0.05, max_depth=4, nthread=1),
  data=xgb.DMatrix(train[,1:2], label=train[,3]),
  nrounds=100000, verbose=FALSE))


system.time(xgb.train(
  params=list(eta=0.05, max_depth=4, nthread=15),
  data=xgb.DMatrix(train[,1:2], label=train[,3]),
  nrounds=100000, verbose=FALSE))


require(xgboost)
x <- matrix(rnorm(100 * 10000), 10000, 100)
y <- x %*% rnorm(100) + rnorm(1000)
system.time({
  bst <- xgboost(data = x, label = y, nthread = 1, nround = 100, verbose = F)
})

system.time({
  bst <- xgboost(data = x, label = y, nthread = 10, nround = 100, verbose = F)
})

