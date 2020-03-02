library(ggplot2)
plants <- read.csv("week3exercise_plantdamage(1).csv")

# a. fitting 
mod1 <- lm(d.biomass~light*damage, data = plants)

# b. pulling out coefficients 
beta <- coef(mod1)

# c. predicting 
prediction <- predict(mod1, newdata = data.frame(light = "L", damage = 0.15)) # d.biomass = 5.632843 

# d. sumulating 
mu <- fitted(mod1) 
sigma <- sigma(mod1)
p_sim <- plants
nsim <- 1000 
beta_mat <- matrix(nrow = nsim, ncol = length(beta))
for (i in 1:nsim) {
  p_sim$y_sim <- rnorm(n = length(p_sim$d.biomass), mean = mu, sd = sigma)
  mod2 <- update(mod1, y_sim ~ ., data = p_sim)
  beta_sim <- coef(mod2)
  beta_mat[i,] <- beta_sim
}
beta_0_ci <- quantile(beta_mat[,1], c(0.025, 0.975))
beta_1_ci <- quantile(beta_mat[,2], c(0.025, 0.975))
beta_2_ci <- quantile(beta_mat[,3], c(0.025, 0.975))
beta_3_ci <- quantile(beta_mat[,4], c(0.025, 0.975))
pars <- data.frame(cbind(beta, rbind(beta_0_ci, beta_1_ci, beta_2_ci, beta_3_ci)))

# e. plotting 
qplot(data=pars, x = rownames(pars), y = pars$beta, ymin = pars$X2.5., ymax = pars$X97.5., geom='pointrange') +
  labs(y = "Estimate", x = "Parameter") +
  theme_bw()


