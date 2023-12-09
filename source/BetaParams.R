#function to get beta parameters 
BetaParams <- function(mu, var) {
  alpha <- mu*((mu*(1-mu)/var-1))
  beta <- (1-mu)*((mu*(1-mu)/var-1))
  return(c(alpha = alpha, beta = beta))
}