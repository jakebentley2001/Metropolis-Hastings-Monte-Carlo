
#defining the gamma function up to a constant
gammafunc <- function(x, alpha, beta){
  if (any(x < 0)) return (0)
  stopifnot(alpha > 0, beta > 0)
  return((x^(alpha-1)*exp(-beta*x)))
}
#Runs the algorithm
MHS <- function(alpha, beta, n=30000){
  x <- numeric(n)
  #sets the initial state
  x[1] <- rexp(1,rate = 1)
  u <- runif(n)
  #this defines the rest of the states
  for (i in 2:n) {
    xt <- x[i-1]
    Y <- rexp(1,rate = xt)
    #defines the acceptance rate
    num <- gammafunc(Y,alpha,beta) * dexp(xt, rate = Y)
    den <- gammafunc(xt,alpha,beta) * dexp(Y, rate = xt)
    #checks if it will be accepted
    if (u[i] <= num/den) {
      #if accepted changed to new accepted state
      x[i] <- Y
    } else {
      #if rejected the state is the same as the previous on and we try again
      x[i] <- xt
    }
  }
  return(x)
}
#Here we display a generated histogram that should follow a gamma pdf with alpha = 3 and beta = 2 
n <- 100000
trials <- MHS(3,2,n)
#Gets rid of the beginning data because it does not approach a gamma distribution at the very beginning of the chain
trials <- trials[round(n/10):n]
hist(trials, xlab = 'X', main = "Method1",
     probability = TRUE,freq=F,col = "red", breaks = 100)
x <- seq(0,1,length=100)
curve(dgamma(x,3,2),col = "blue",lwd=2,add=TRUE)





set.seed(100)
chifunc <- function(x,k) {
  if (any(x < 0)) return (0)
  stopifnot(k > 0, k > 0)
  return(x^((k/2)-1)*exp(-x/2))
}
MHS2 <- function(k, n=30000) {
  x <- numeric(n)
  x[1] <- rgamma(1,shape = k, scale = k)
  u <- runif(n)
  for (i in 2:n) {
    xt <- x[i-1]
    Y <- rgamma(1,shape = xt, scale = xt/2)
    num <- chifunc(Y,k) * dgamma(xt, shape = Y, scale = Y/2)
    den <- chifunc(xt,k) * dgamma(Y, shape = xt, scale = xt/2)
    if (u[i] <= num/den) {
      x[i] <- Y
    } else {
      x[i] <- xt
    }
  }
  return(x)
}
n <- 300000
trials <- MHS2(2,n)
trials <- trials[round(n/4):n]
hist(trials, xlab = 'X', main = "Method2",
     probability = TRUE,freq=F,col = "red", breaks = 100)
x <- seq(0,1,length=100)
curve(dchisq(x,2),col = "blue",lwd=2,add=TRUE)

