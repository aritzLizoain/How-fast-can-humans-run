# Monte Carlo simulation algorithm
library(Rlab) # for rbern()
set.seed(42)

MC_simulation = function(B, years, beta11, beta12, beta21, beta22, last_record){
  
  t = years_to_t(years)
  n = length(t)
  simulation_matrix = matrix(nrow=B, ncol=n)
  
  for (j in 1:B){
    
    # 1. Randomly generate observations from h_t(e), where theta(t) is theta_hat_t
    theta_t_hat = theta_t(beta21, beta22, t)
    a = rexp(n, rate=1/theta_t_hat)
    # Same as:
    # a = numeric(n)
    # for (i in 1:n){
    #   a[i] = rexp(1, rate=1/theta_t_hat[1])
    # }
    # But with B=100,000 it goes from 2s to 11s
    
    # 2. Randomly generate observations from a Bernoulli distribution with mean p(t), where p(t) is p_hat_t
    p_t_hat = p_t(beta11, beta12, t)
    b = rbern(n, prob=p_t_hat)
    
    # 3. Random walk error: e_t = a_t*b_t
    e_t = a*b
    
    # 4. Random walk: X_t = X_{t-1} - e_t
    simulation_estimates = numeric(n)
    simulation_estimates[1] = last_record
    for (i in 2:n){
      simulation_estimates[i] = simulation_estimates[i-1] - e_t[i]
    }
    # Written in a more efficient way
    # e_t[1] = 0
    # simulation_estimates = cumsum(simulation_estimates-e_t)
    # With B=100,000 it goes from 1.7s to 1.85s
    
    # We can visualize different paths
    # plot(year, record, type='l', lwd=2, xlim=c(min(year), max(years)), ylim=c(9.4, max(record)), xlab=quote(bold('Year')), ylab=quote(bold('Record (s)')))
    # lines(years, simulation_estimates, col='red', lwd=2)
     
    # 5. Save simulation in matrix and repeat
    simulation_matrix[j,] = simulation_estimates
    
  }
  return(simulation_matrix)
}
