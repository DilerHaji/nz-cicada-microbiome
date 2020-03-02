library(R2jags)
library(lme4)

# 1. 
load("../ProblemSet10.rdata")
attach(ps.data)

# 2. 
par(mfrow = c(1,3))
plot(elev, count[,1])
plot(forest, count[,1])
plot(wind[,1], count[,1])

# 4. 
sink("jags_ps.R")
cat("
    model{

    # random effects 
    for (k in 1:n.obs) {
      for (p in 1:n.visit) {
        wind[k,p] ~ dnorm(0, tau.bw)
      }
    }
    for (w in 1:n.visit) {
    b[w] ~ dnorm(0, tau.b)
    }

    # process model - likelihood
    for (i in 1:n.obs) {
      for (j in 1:n.visit) {
        log(mu[i,j]) <- beta0 + beta1*elev[i] + beta2*forest[i] + beta3*elev[i]*forest[i] + wind[i,j] + b[j]
        count[i,j] ~ dpois(mu[i,j])
      }
    }

    # priors
    tau.b ~ dgamma(0.001, 0.001)
    sigma.b <-  pow(tau.b, -0.5)
    tau.bw ~ dgamma(0.001, 0.001)
    sigma.bw <- pow(tau.bw, -0.5)
    beta0 ~ dnorm(0, 0.001)
    beta1 ~ dnorm(0, 0.001)
    beta2 ~ dnorm(0, 0.001)
    beta3 ~ dnorm(0, 0.001)
    
  }
  ") ; sink()

# 5. 
jags.data <- list(
  elev = elev,
  forest = forest,
  wind = wind,
  count = count,
  n.obs = length(elev),
  n.visit = dim(count)[2]
)

jags.model <- jags(data = jags.data,
                   parameters.to.save = c("beta0", "beta1", "beta2", "beta3", "sigma.b", "sigma.bw"), 
                   model.file = "jags_ps.R",
                   n.chains = 10, 
                   n.iter = 20000,
                   n.burnin = 1000)

traceplot(jags.model)


# 6.

# Elevation
par(mfrow = c(1,2))
mcount <- apply(count, 1, mean)
plot(elev, mcount, xlab = "elevation", ylab = "count", col = "grey")
mod.post <- array(dim = c(length(jags.model$BUGSoutput$sims.list$beta1), length(mcount)))
for (i in 1:length(mcount)) {
  mod.post[,i] <- exp(jags.model$BUGSoutput$sims.list$beta0 + jags.model$BUGSoutput$sims.list$beta1*elev[i])
}
mod.mean <- apply(mod.post, 2, mean)
points(elev, mod.mean, pch = 16, cex = 0.5)
mod.random <- array(dim = c(length(jags.model$BUGSoutput$sims.list$sigma.b), length(mcount)))
for (i in 1:length(mcount)) {
  mod.random[,i] <- jags.model$BUGSoutput$sims.list$beta0 + jags.model$BUGSoutput$sims.list$beta1*elev[i] + jags.model$BUGSoutput$sims.list$sigma.bw + jags.model$BUGSoutput$sims.list$sigma.b
}
ranup <- apply(mod.random, 2, mean)
ranup <- exp(ranup)
mod.random <- array(dim = c(length(jags.model$BUGSoutput$sims.list$sigma.b), length(mcount)))
for (i in 1:length(mcount)) {
  mod.random[,i] <- jags.model$BUGSoutput$sims.list$beta0 + jags.model$BUGSoutput$sims.list$beta1*elev[i] - jags.model$BUGSoutput$sims.list$sigma.bw - jags.model$BUGSoutput$sims.list$sigma.b
}
randown <- apply(mod.random, 2, mean)
randown <- exp(randown)
arrows(elev, ranup, elev, randown, length=0.05, angle=90, code=3, lwd = 0.1)

# Forest Cover 
mcount <- apply(count, 1, mean)
plot(forest, mcount, xlab = "forest cover", ylab = "count", col = "grey")
mod.post <- array(dim = c(length(jags.model$BUGSoutput$sims.list$beta1), length(mcount)))
for (i in 1:length(mcount)) {
  mod.post[,i] <- exp(jags.model$BUGSoutput$sims.list$beta0 + jags.model$BUGSoutput$sims.list$beta2*forest[i])
}
mod.mean <- apply(mod.post, 2, mean)
points(forest, mod.mean, pch = 16, cex = 0.5)
mod.random <- array(dim = c(length(jags.model$BUGSoutput$sims.list$sigma.b), length(mcount)))
for (i in 1:length(mcount)) {
  mod.random[,i] <- jags.model$BUGSoutput$sims.list$beta0 + jags.model$BUGSoutput$sims.list$beta2*forest[i] + jags.model$BUGSoutput$sims.list$sigma.bw + jags.model$BUGSoutput$sims.list$sigma.b
}
ranup <- apply(mod.random, 2, mean)
ranup <- exp(ranup)
mod.random <- array(dim = c(length(jags.model$BUGSoutput$sims.list$sigma.b), length(mcount)))
for (i in 1:length(mcount)) {
  mod.random[,i] <- jags.model$BUGSoutput$sims.list$beta0 + jags.model$BUGSoutput$sims.list$beta2*forest[i] - jags.model$BUGSoutput$sims.list$sigma.bw - jags.model$BUGSoutput$sims.list$sigma.b
}
randown <- apply(mod.random, 2, mean)
randown <- exp(randown)
arrows(forest, ranup, forest, randown, length=0.05, angle=90, code=3, lwd = 0.1)


# 7.
count <- as.vector(count)
wind <- as.vector(wind)
forest <- rep(forest, 3)
elev <- rep(elev, 3)
block <- c(rep(1, 267), rep(2, 267), rep(3, 267))

glmer(count ~ elev*forest + (1|wind) + (1|block), family = poisson(link = log))

