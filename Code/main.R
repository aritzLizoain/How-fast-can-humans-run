# A Simple Random Walk Model for Predicting Track and Field World Records
# An implementation of: https://www.degruyter.com/document/doi/10.2202/1559-0410.1067/html

# ------------------------------------------------------------------------------
# DATA
# ------------------------------------------------------------------------------

year = c(1884,1886,1892,1912,1921,1930,1932,1948,1958,1960,1964,1968,1983,1988,1991,1994,1996,1999,2002)
record = c(11.49,11.44,11.04,10.84,10.64,10.54,10.38,10.34,10.29,10.24,10.06,9.95,9.93,9.92,9.86,9.85,9.84,9.79,9.78)
plot(year,record, main = '100m world record progression', type='l', xlim=c(1884,2002), ylim=c(9.5,11.5), lwd=2, xlab=quote(bold('Year')), ylab=quote(bold('Record (s)')))

# ------------------------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------------------------

# Time conversion
min_year = min(year)
years_to_t = function(year){(year-min_year+1)/100}

p_t = function(beta11, beta12, t){exp(beta11+beta12*t)/(1+exp(beta11+beta12*t))}
theta_t = function(beta21, beta22, t){exp(beta21+beta22*t)}

# MC simulation algorithm
setwd("C:/Users/Aritz/Desktop/1.2/(o) Computational Statistics/Project")
source("MC_simulation.R")

# MLE function
source("optim_exponential.R")

# ------------------------------------------------------------------------------
# GET PARAMETERS
# ------------------------------------------------------------------------------

# PAPER
beta11_paper = -2.676
beta12_paper = 1.447
beta21_paper = -1.02
beta22_paper = -2.181

# MINE

# Beta1 with logistic regression
time = min(year):max(year)
response = rep(0, length(time))
response[match(year, time)] = 1
model = glm(response ~ years_to_t(time),family=binomial(link='logit'))
# summary(model)
beta11 = summary(model)$coeff[1,1] 
beta12 = summary(model)$coeff[2,1]

# Beta2 with MLE
# initial values
initial_values = c(-1,-1)
# get MLE
optimization_likelihood = optim_l_exponential(initial_values, year, record)
optimization_log = optim_ll_exponential(initial_values, year, record)
beta2 = optimization_log$par
beta21 = beta2[1] 
beta22 = beta2[2] 

# ------------------------------------------------------------------------------
# MONTE CARLO SIMULATION - 50 years
# ------------------------------------------------------------------------------

# Simulation settings
future_years = 50
years = max(year):(max(year)+future_years)
n = length(years)
B = 10000 # 100,000 in the paper

simulation_matrix_paper = MC_simulation(B, years, beta11_paper, beta12_paper, beta21_paper, beta22_paper, last_record=9.78)
simulation_matrix = MC_simulation(B, years, beta11, beta12, beta21, beta22, last_record=9.78)

# Simulation means as future predictions
simulation_predictions_paper = colMeans((simulation_matrix_paper))
simulation_predictions = colMeans((simulation_matrix))

# ------------------------------------------------------------------------------
# CONFIDENCE INTERVALS
# ------------------------------------------------------------------------------

# Get 2.5th and the 97.5th percentiles of the simulated data (95% confidence intervals)
alpha = 0.05
# Paper
upper_limits_paper = numeric(n)
lower_limits_paper = numeric(n)

for (i in 1:n){
  year_predictions_paper = simulation_matrix_paper[,i]
  upper_limits_paper[i] = quantile(year_predictions_paper, 1-alpha/2)
  lower_limits_paper[i] = quantile(year_predictions_paper, alpha/2)
}
# Mine
upper_limits = numeric(n)
lower_limits = numeric(n)
for (i in 1:n){
  year_predictions = simulation_matrix[,i]
  upper_limits[i] = quantile(year_predictions, 1-alpha/2)
  lower_limits[i] = quantile(year_predictions, alpha/2)
}

# Plot records + simulation predictions + confidence intervals
plot(year, record, main = '100m world record prediction', type='l', lwd=2, xlim=c(min(year), max(years)), ylim=c(min(lower_limits)-0.1, max(record)+0.1), xlab=quote(bold('Year')), ylab=quote(bold('Record (s)')))
# Paper
lines(years, simulation_predictions_paper, col='gray', lwd=2)
lines(years, upper_limits_paper, col='gray', lwd=1.8, lty=2)
lines(years, lower_limits_paper, col='gray', lwd=1.8, lty=2)
# Mine
lines(years, simulation_predictions, col='red', lwd=2)
lines(years, upper_limits, col='red', lwd=1.8, lty=2)
lines(years, lower_limits, col='red', lwd=1.8, lty=2)

# ------------------------------------------------------------------------------
# HOW FAST CAN HUMANS RUN? - 400 years predictions
# ------------------------------------------------------------------------------

# Get predictions for the future 400 years
future_years_long = 400
years_long = max(year):(max(year)+future_years_long)
n = length(years_long)
B = 10000 # 100,000 in the paper

# Get simulations (with my parameters)
simulation_matrix_long = MC_simulation(B, years_long, beta11, beta12, beta21, beta22, last_record=9.78)
# Get mean
simulation_predictions_long = colMeans((simulation_matrix_long))
# Get confidence intervals
upper_limits_long = numeric(n)
lower_limits_long = numeric(n)
alpha = 0.05
for (i in 1:n){
  year_predictions_long = simulation_matrix_long[,i]
  upper_limits_long[i] = quantile(year_predictions_long, 1-alpha/2)
  lower_limits_long[i] = quantile(year_predictions_long, alpha/2)
}

# Get simulations (with paper parameters)
simulation_matrix_long_paper = MC_simulation(B, years_long, beta11_paper, beta12_paper, beta21_paper, beta22_paper, last_record=9.78)
# Get mean
simulation_predictions_long_paper = colMeans((simulation_matrix_long_paper))
# Get confidence intervals
upper_limits_long_paper = numeric(n)
lower_limits_long_paper = numeric(n)
for (i in 1:n){
  year_predictions_long_paper = simulation_matrix_long_paper[,i]
  upper_limits_long_paper[i] = quantile(year_predictions_long_paper, 1-alpha/2)
  lower_limits_long_paper[i] = quantile(year_predictions_long_paper, alpha/2)
}

# Plot
plot(year, record, type='l', main = '100m world record prediction', lwd=2, xlim=c(min(year), max(years_long)), ylim=c(min(lower_limits_long)-0.1, max(record)+0.1), xlab=quote(bold('Year')), ylab=quote(bold('Record (s)')))
lines(years_long, simulation_predictions_long_paper, col='gray', lwd=2)
lines(years_long, upper_limits_long_paper, col='gray', lwd=1.8, lty=2)
lines(years_long, lower_limits_long_paper, col='gray', lwd=1.8, lty=2)
lines(years_long, simulation_predictions_long, col='red', lwd=2)
lines(years_long, upper_limits_long, col='red', lwd=1.8, lty=2)
lines(years_long, lower_limits_long, col='red', lwd=1.8, lty=2)

# On the paper they get ultimate record of of 9.296
min(simulation_predictions_long)

# ------------------------------------------------------------------------------
# 6 NEW WORLD RECORDS - EVALUATING THE MODEL'S PREDICTIONS
# ------------------------------------------------------------------------------

# The paper is from 2007
# In the paper they claim: Note first that the prediction interval for year 2006 contains the current world record of 9.77 seconds.
# In fact all records are within the prediction intervals, all except the 2009 9.58 record by Usain Bolt
# This huge improvement was outside the predicted confidence intervals
# New records: page 33 http://iaaf-ebooks.s3.amazonaws.com/2015/Progression-of-IAAF-World-Records-2015/projet/IAAF-WRPB-2015.pdf 

# Same plot as before but changed xlim
plot(year, record, type='l', main = '100m world record prediction', lwd=2, xlim=c(min(year), 2150), ylim=c(min(lower_limits_long)-0.1, max(record)+0.1), xlab=quote(bold('Year')), ylab=quote(bold('Record (s)')))
lines(years_long, simulation_predictions_long_paper, col='gray', lwd=2)
lines(years_long, upper_limits_long_paper, col='gray', lwd=1.8, lty=2)
lines(years_long, lower_limits_long_paper, col='gray', lwd=1.8, lty=2)
lines(years_long, simulation_predictions_long, col='red', lwd=2)
lines(years_long, upper_limits_long, col='red', lwd=1.8, lty=2)
lines(years_long, lower_limits_long, col='red', lwd=1.8, lty=2)

year2 = c(2002,2005,2007,2008,2008,2009)
record2 = c(9.78,9.77,9.74,9.72,9.69,9.58)
lines(year2,record2, col='blue', lwd=3)

# In the paper they also say: the intervals can be used to estimate when the record will be broken next
# This will occur when the interval no longer contains the current record.
# If this is true, Bolt's record will not be broken until...
record_breaking_year = years_long[which(upper_limits_long<9.58)[1]]
record_breaking_year_paper = years_long[which(upper_limits_long_paper<9.58)[1]]

clip(1,record_breaking_year, -100, 9.58) # clip with my results because it happens later in time
abline(h=9.58, col='darkgray', lwd=2)
text(9.58, col='darkgray', x = 1890, y = 9.51, cex=1.5, font=2)
abline(v=record_breaking_year, col='darkgray', lwd=2)
text(eval(record_breaking_year), col='darkgray', x = record_breaking_year-16, y = 9.05, cex=1.5, font=2)
abline(v=record_breaking_year_paper, col='darkgray', lwd=1)
text(eval(record_breaking_year_paper), col='darkgray', x = record_breaking_year_paper-11, y = 9.04, cex=1, font=2)

# ------------------------------------------------------------------------------
# UPDATING THE MODEL PARAMETERS ACCORDING TO THE NEW WORLD RECORDS
# ------------------------------------------------------------------------------

# Beta1 with logistic regression
time_update = min(year):max(year2)
response_update = rep(0, length(time_update))
response_update[match(c(year,year2), time_update)] = 1
model_update = glm(response_update ~ years_to_t(time_update),family=binomial(link='logit'))
# summary(model_update)
beta11_update = summary(model_update)$coeff[1,1] 
beta12_update = summary(model_update)$coeff[2,1]

# Beta2 with MLE
# initial values same as before
# get MLE
# optimization_likelihood_update = optim_l_exponential(initial_values, c(year,year2), c(record,record2)) # doesn't work
optimization_log_update = optim_ll_exponential(initial_values, c(year,year2), c(record,record2))
beta2_update = optimization_log_update$par
beta21_update = beta2_update[1] 
beta22_update = beta2_update[2] 

# Compare paper vs. mine vs. updated coefficients
coefficients = cbind(c(beta11_paper, beta12_paper, beta21_paper, beta22_paper),c(beta11, beta12, beta21, beta22),c(beta11_update, beta12_update, beta21_update, beta22_update))
colnames(coefficients) = c("Paper","Mine","Updated")
rownames(coefficients) = c("Beta11","Beta12", "Beta21", "Beta22")
print(round(coefficients,3))

# ------------------------------------------------------------------------------
# MONTE CARLO SIMULATION 2 + CONFIDENCE INTERVALS 2
# ------------------------------------------------------------------------------

# Simulation settings
years2 = max(year2):(max(year)+future_years)
n2 = length(years2)
B = 10000 # 100,000 in the paper

simulation_matrix_update = MC_simulation(B, years2, beta11_update, beta12_update, beta21_update, beta22_update, last_record=9.58)
# Simulation means as future predictions
simulation_predictions_update = colMeans((simulation_matrix_update))

# Get 2.5th and the 97.5th percentiles of the simulated data (95% confidence intervals)
upper_limits_update = numeric(n2)
lower_limits_update = numeric(n2)
for (i in 1:n2){
  year_predictions_update = simulation_matrix_update[,i]
  upper_limits_update[i] = quantile(year_predictions_update, 1-alpha/2)
  lower_limits_update[i] = quantile(year_predictions_update, alpha/2)
}

# Plot records + simulation predictions + confidence intervals
plot(year, record, main = '100m world record prediction', type='l', lwd=2, xlim=c(min(year), max(years)), ylim=c(min(lower_limits_update)-0.1, max(record)+0.1), xlab=quote(bold('Year')), ylab=quote(bold('Record (s)')))
# Paper
lines(years, simulation_predictions_paper, col='gray', lwd=2)
lines(years, upper_limits_paper, col='gray', lwd=1.8, lty=2)
lines(years, lower_limits_paper, col='gray', lwd=1.8, lty=2)
# Mine
lines(years, simulation_predictions, col='red', lwd=2)
lines(years, upper_limits, col='red', lwd=1.8, lty=2)
lines(years, lower_limits, col='red', lwd=1.8, lty=2)
# Add new records
lines(year2, record2, lwd=2)
# Add updated model predictions
lines(years2, simulation_predictions_update, col='blue', lwd=2)
lines(years2, upper_limits_update, col='blue', lwd=1.8, lty=2)
lines(years2, lower_limits_update, col='blue', lwd=1.8, lty=2)

# With the new model the record should've been broken in 2018
years2[which(upper_limits_update<9.58)[1]]
# abline(h=9.58, col='gray')
# abline(v=2018, col='gray')

# New ultimate records
min(simulation_predictions_update) # 8.79

text(round(min(simulation_predictions),2), x=2300, y=round(min(simulation_predictions),2)+0.3, cex=1.5, col='red', font=2)
text(round(min(simulation_predictions_update),2), x=2300, y=round(min(simulation_predictions_update),2)-0.35, cex=1.5, col='blue', font=2)
text(round(min(simulation_predictions_paper),2), x=2300, y=round(min(simulation_predictions_paper),2), cex=1.5, col='darkgray', font=2)
