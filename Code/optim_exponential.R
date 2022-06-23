# LIKELIHOOD
optim_l_exponential = function(initial_values, year, record){
  
  # Define the likelihood
  l_exponential <- function(beta2,x) {
    n = length(x)
    beta21 = beta2[1]
    beta22 = beta2[2]
    likelihood = exp(-beta21*n) * exp(-beta22*sum(x)) * exp(-sum(et*exp(-beta21-beta22*x)))
    return(-likelihood)
  }
  
  # parametrization of time variable (explanation on page 10 of the paper)
  # year is only the years when the record is broken
  t = years_to_t(year)
  # et (amount by which a record is broken)
  et = -diff(record)
  # we don't know by how much the first record was broken
  # so we don't take into account the first record
  x = t[-1]
  
  # Maximize (minimize the negative)
  optimization = optim(par=initial_values, fn=l_exponential, x=x)
}

# -------------------------------------------------------------------------------------------

# LOG-LIKELIHOOD
optim_ll_exponential = function(initial_values, year, record){
  
  # Define the log-likelihood
  ll_exponential <- function(beta2,x) {
    n = length(x)
    beta21 = beta2[1]
    beta22 = beta2[2]
    log_likelihood = -beta21*n - beta22*sum(x) - sum(et*exp(-beta21-beta22*x))
    return(-log_likelihood)
  }
  
  # parametrization of time variable (explanation on page 10 of the paper)
  # year is only the years when the record is broken
  t = years_to_t(year)
  # et (amount by which a record is broken)
  et = -diff(record)
  # we don't know by how much the first record was broken
  # so we don't take into account the first record
  x = t[-1]
  
  # Maximize (minimize the negative)
  optimization = optim(par=initial_values, fn=ll_exponential, x=x)
}

# -------------------------------------------------------------------------------------------

# # GET THE PARAMETERS
# 
# # Data
# year = c(1884,1886,1892,1912,1921,1930,1932,1948,1958,1960,1964,1968,1983,1988,1991,1994,1996,1999,2002)
# record = c(11.49,11.44,11.04,10.84,10.64,10.54,10.38,10.34,10.29,10.24,10.06,9.95,9.93,9.92,9.86,9.85,9.84,9.79,9.78)
# 
# # Time parametrization function
# years_to_t = function(year){(year-min(year)+1)/100}
# 
# # Beta2 with MLE
# initial_values = c(10,10)
# optimization_likelihood = optim_l_exponential(initial_values, year, record)
# optimization_log = optim_ll_exponential(initial_values, year, record)
# beta2 = optimization_log$par
# beta21 = beta2[1] 
# beta22 = beta2[2]
# 
# beta21_paper = -1.02
# beta22_paper = -2.181
