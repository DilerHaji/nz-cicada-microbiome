
    model{

    # process model
    for (i in 1:N){
pdis[i] ~ dbeta(alpha[i], beta[i])
  alpha[i] <- mu[i] * phi
  beta[i]  <- (1-mu[i]) * phi
  logit(mu[i]) <- 
  a + 
  (b + bj[library_block[i]])*copho[i] + 
  (c + bj[library_block[i]])*geodis[i] + 
  (d + bj[library_block[i]])*elev_diff[i] + 
  (e + bj[library_block[i]])*habitat_diff[i] +
  bj[library_block[i]]
}
    
    # latent factors
    for(j in N2) {
       bj[j] ~ dnorm(0, taubj)
    }

    # priors
    phi ~ dgamma(.001,.001)
    taubj ~ dgamma(.001,.001)
    a ~ dnorm(0,.001)
b ~ dnorm(0,.001)
c ~ dnorm(0,.001)
d ~ dnorm(0,.001)
e ~ dnorm(0,.001)
    
  }
  