# Load in the required libraries for data manipulation
# and multivariate normal distribution
require(MASS)
require(plyr)
require(reshape2)
require(ggplot2)

# Set a seed for repeatable plots
set.seed(12345)

# Calculates the covariance matrix sigma using a
# simplified version of the squared exponential function.
#
# Although the nested loops are ugly, I've checked and it's about
# 30% faster than a solution using expand.grid() and apply()
#
# Parameters:
# X1, X2 = vectors
# l = the scale length parameter
# Returns:
# a covariance matrix
calcSigma <- function(X1,X2,l=1) {
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
    }
  }
  return(Sigma)
}

# 1. Plot some sample functions from the Gaussian process
# as shown in Figure 2.2(a)

# Define the points at which we want to define the functions
x.star <- seq(-5,5,len=50)

# Calculate the covariance matrix
sigma <- calcSigma(x.star,x.star)

# Generate a number of functions from the process
n.samples <- 3
values <- matrix(rep(0,length(x.star)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  values[,i] <- mvrnorm(1, rep(0, length(x.star)), sigma)
}
values <- cbind(x=x.star,as.data.frame(values))
values <- melt(values,id="x")

# Plot the result
fig2a <- ggplot(values,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable)) +
  theme_bw() +
  scale_y_continuous(lim=c(-2.5,2.5), name="output, f(x)") +
  xlab("input, x")
plot(fig2a)
# 2. Now let's assume that we have some known data points;
# this is the case of Figure 2.2(b). In the book, the notation 'f'
# is used for f$y below. I've done this to make the ggplot code
# easier later on.
f <- data.frame(x=c(-4,-3,-1,0,2),
                y=c(-2,0,1,2,-1))

# Calculate the covariance matrices
# using the same x.star values as above
x <- f$x
k.xx <- calcSigma(x,x)
k.xxs <- calcSigma(x,x.star)
k.xsx <- calcSigma(x.star,x)
k.xsxs <- calcSigma(x.star,x.star)

# These matrix calculations correspond to equation (2.19)
# in the book.
f.star.bar <- k.xsx%*%solve(k.xx)%*%f$y
cov.f.star <- k.xsxs - k.xsx%*%solve(k.xx)%*%k.xxs

# This time we'll plot more samples. We could of course
# simply plot a +/- 2 standard deviation confidence interval
# as in the book but I want to show the samples explicitly here.
n.samples <- 50
values <- matrix(rep(0,length(x.star)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  values[,i] <- mvrnorm(1, f.star.bar, cov.f.star)
}
values <- cbind(x=x.star,as.data.frame(values))
values <- melt(values,id="x")

# Plot the results including the mean function
# and constraining data points
fig2b <- ggplot(values,aes(x=x,y=value)) +
  geom_line(aes(group=variable), colour="grey80") +
  geom_line(data=NULL,aes(x=x.star,y=f.star.bar),colour="red", size=1) +
  geom_point(data=f,aes(x=x,y=y)) +
  theme_bw() +
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
plot(fig2b)
# 3. Now assume that each of the observed data points have some
# normally-distributed noise.

# The standard deviation of the noise
sigma.n <- 0.1

# Recalculate the mean and covariance functions
f.bar.star <- k.xsx%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%f$y
cov.f.star <- k.xsxs - k.xsx%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%k.xxs

# Recalulate the sample functions
values <- matrix(rep(0,length(x.star)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  values[,i] <- mvrnorm(1, f.bar.star, cov.f.star)
}
values <- cbind(x=x.star,as.data.frame(values))
values <- melt(values,id="x")

# Plot the result, including error bars on the observed points
gg <- ggplot(values, aes(x=x,y=value)) +
  geom_line(aes(group=variable), colour="grey80") +
  geom_line(data=NULL,aes(x=x.star,y=f.bar.star),colour="red", size=1) +
  geom_errorbar(data=f,aes(x=x,y=NULL,ymin=y-2*sigma.n, ymax=y+2*sigma.n), width=0.2) +
  geom_point(data=f,aes(x=x,y=y)) +
  theme_bw() +
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")

plot(gg)


####################################################
####Start of Mark's new code
####################################################



####################################################
####Set up generalized K function
####################################################

K <- function(X1,X2,l=1,sigma.f=1) {
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      k1<- sigma.f*exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
      Sigma[i,j] <- k1
    }
  }
  return(Sigma)
}

####################################################
####Sample some data from this function to check fit
####################################################
# Define the seed points at which we want to define the functions
x <- seq(-5,5,len=50)

# Calculate the covariance matrix
l=1
sigma.f=10##this is the variance of the function
sigma.n=0.1##this is the variance of the noise
sigma <- K(x,x,l,sigma.f)

# Generate an observed data vector from the process
truth <- data.frame(x=x, y=mvrnorm(1, rep(0, length(x)), sigma))


# Plot the result, including error bars on the observed points
gg <- ggplot(truth, aes(x=x,y=y)) +
  geom_line(colour="red", size=1) +
  geom_ribbon(aes(ymin=y-1.96*sqrt(sigma.n), ymax=y+1.96*sqrt(sigma.n)),alpha=0.2) +
  geom_point() +
  theme_bw() +
  scale_y_continuous( name="output, f(x)") +
  xlab("input, x")

plot(gg)


#################################
####Now generate some sample data
#################################
x.star <- seq(-5,5,len=100)
# Calculate the covariance matrices
# using the same x.star values as above
k.xx <- K(truth$x,truth$x,l,sigma.f)
k.xx <-k.xx+sigma.n*diag(1,ncol(k.xx))
k.xxs <- K(truth$x,x.star,l,sigma.f)
k.xsx <- K(x.star,truth$x,l,sigma.f)
k.xsxs <- K(x.star,x.star,l,sigma.f)
k.xsxs <-k.xsxs+sigma.n*diag(1,ncol(k.xsxs))

# These matrix calculations correspond to equation (2.19)
# in the book.

f.bar.star <- k.xsx%*%solve(k.xx)%*%truth$y
cov.f.star <- k.xsxs - k.xsx%*%solve(k.xx)%*%k.xxs

# This time we'll plot more samples. We could of course
# simply plot a +/- 2 standard deviation confidence interval
# as in the book but I want to show the samples explicitly here.
n.samples <- 5
values <- matrix(rep(0,length(x.star)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  values[,i] <- mvrnorm(1, f.bar.star, cov.f.star)
}
values <- cbind(x=x.star,as.data.frame(values))
values <- melt(values,id="x")
colnames(values)[3]="y"

# Plot the results including the mean function
# and constraining data points
fig2b <- ggplot(values,aes(x=x,y=y)) +
  geom_line(aes(group=variable), colour="grey80") +
  geom_line(data=NULL,aes(x=x.star,y=f.bar.star),colour="red", size=1) +
  geom_ribbon(data=truth,aes(ymin=y-1.96*sqrt(sigma.n), ymax=y+1.96*sqrt(sigma.n)),alpha=0.2) +
  geom_point(data=truth,aes(x=x,y=y)) +
  theme_bw() +
  scale_y_continuous(name="output, f(x)") +
  xlab("input, x")
plot(fig2b)

####################################################
####Use mlegp package to estimate GP parameters
####################################################
library(mlegp)
x<-values$x
y<-values$y

fit = mlegp(x,y, nugget=0.1 ,nugget.known = 1)


summary(fit)
plot(fit)

fit_preds<-predict(fit,se.fit = TRUE)
fit_preds$x<-fit$X
fit_preds$se.fit<-sqrt(fit_preds$se.fit^2+fit$nugget)
fit_dat<-data.frame(fit_preds)

fig<- ggplot(values,aes(x=x,y=y)) +
  geom_line(aes(group=variable), colour="grey80") +
  geom_line(data=fit_dat,aes(x=x,y=fit),colour="red", size=1) +
  geom_ribbon(data=fit_dat,aes(x=x,y=NULL,ymin=fit-1.96*se.fit, ymax=fit+1.96*se.fit), alpha=0.2) +
  geom_point(data=truth,aes(x=x,y=y)) +
  theme_bw() +
  scale_y_continuous( name="output, f(x)") +
  xlab("input, x")
plot(fig)





fit = mlegp(x,y,nugget.known = 1)

summary(fit)
plot(fit)

fit_preds<-predict(fit,se.fit = TRUE)
fit_preds$x<-fit$X
fit_preds$se.fit<-sqrt(fit_preds$se.fit^2+fit$nugget)
fit_dat<-data.frame(fit_preds)

fig<- ggplot(values,aes(x=x,y=y)) +
  geom_line(aes(group=variable), colour="grey80") +
  geom_line(data=fit_dat,aes(x=x,y=fit),colour="red", size=1) +
  geom_ribbon(data=fit_dat,aes(x=x,y=NULL,ymin=fit-1.96*se.fit, ymax=fit+1.96*se.fit), alpha=0.2) +
  geom_point(data=truth,aes(x=x,y=y)) +
  theme_bw() +
  scale_y_continuous( name="output, f(x)") +
  xlab("input, x")
plot(fig)


fit = mlegp(x,y,nugget.known = 0)

summary(fit)
plot(fit)

fit_preds<-predict(fit,se.fit = TRUE)
fit_preds$x<-fit$X
fit_preds$se.fit<-sqrt(fit_preds$se.fit^2+fit$nugget)
fit_dat<-data.frame(fit_preds)

fig<- ggplot(values,aes(x=x,y=y)) +
  geom_line(aes(group=variable), colour="grey80") +
  geom_line(data=fit_dat,aes(x=x,y=fit),colour="red", size=1) +
  geom_ribbon(data=fit_dat,aes(x=x,y=NULL,ymin=fit-1.96*se.fit, ymax=fit+1.96*se.fit), alpha=0.2) +
  geom_point(data=truth,aes(x=x,y=y)) +
  theme_bw() +
  scale_y_continuous( name="output, f(x)") +
  xlab("input, x")
plot(fig)


####################################################
####Set up the marginal likelihood to estimate the parameters
####################################################
marginal_lik<-function(par,data){
  l=par[1]
  sigma.f=par[2]
  sigma.n=par[3]
  k.xx <- K(data$x,data$x,l,sigma.f)
  k.xx<-k.xx+sigma.n^2*diag(1,ncol(k.xx))
  lik<-(-0.5)*t(data$y)%*%solve(k.xx)%*%data$y-(0.5)*log(det(k.xx))-(length(data$y)/2)*log(2*pi)
  return(-lik)###return the negative for minimization
}

result<-optim(par=c(1,1,1),marginal_lik,data=values)

result
