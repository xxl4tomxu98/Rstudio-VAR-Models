rm(list = ls())
graphics.off()
par(mar = c(1, 1, 1, 1)) # Set the margin on all sides to 2
library(devtools)
devtools::install_github("KevinKotze/tsm")
library(tsm)
library(vars)
library(mFilter)

# Download data
data <- read.table("http://www.jmulti.de/download/datasets/e1.dat", skip = 6, header = TRUE)

# Only use the first 76 observations so that there are 73 observations
# left for the estimated VAR(2) model after taking first differences.
data <- data[1:76, ]
invest <-  ts(data$invest, start = c(1960, 1), frequency = 4)
income <-  ts(data$income, start = c(1960, 1), frequency = 4)
cons <-  ts(data$cons, start = c(1960, 1), frequency = 4)

#As the persistence in all data is relatively high we need to check for
# a unit root, where we note that we can't reject null of unit root for
# invest, was able for income and consumption when using the Dickey-Fuller test.
adf.invest <- ur.df(invest, type = "trend", selectlags = "AIC")
summary(adf.invest)

adf.income <- ur.df(income, type = "trend", selectlags = "AIC")
summary(adf.income)

adf.cons <- ur.df(cons, type = "trend", selectlags = "AIC")
summary(adf.cons)

# persistency of data by auto-correlation function
invest.acf <- ac(invest, main = "investment")
income.acf <- ac(income, main = "income")
cons.acf <- ac(cons, main = "consumption")

# Convert to time series object
data <- ts(data, start = c(1960, 1), frequency = 4)

# Take logs and differences
diff_data <- diff(as.matrix(log(data)))

# Plot data
plot(diff_data, main = "Dataset E1 from Lütkepohl (2007)")

#information criteria to decide upon the number of lags to include
info.bv <- VARselect(diff_data, lag.max = 12, type = "const")
info.bv$selection

#This data is used to estimate a VAR(2) model with a constant term:
# yt=v+∑i=12Aiyt−i+ut,
# where ut∼N(0,Σ).
# Estimate model
model <- VAR(diff_data, p = 2, type = "const")

# Look at summary statistics, rejecting null hypothesis of non-correlations
summary(model)

#To consider the model fit we can perform some diagnostic tests on residuals
#of the model. To test for serial correlation we can apply a Portmanteau-test,
#can't reject null therefore pass if p > 0.05
bv.serial <- serial.test(model, lags.pt = 12, type = "PT.asymptotic")
bv.serial

plot(bv.serial, names = "invest")
plot(bv.serial, names = "income")
plot(bv.serial, names = "cons")

#To test for heteroscedasticity in the residuals we can perform a 
#multivariate ARCH Lagrange-Multiplier test. p-value that is greater than 5%
#would indicate the absence of heteroscedasticity
bv.arch <- arch.test(model, lags.multi = 12, multivariate.only = TRUE)
bv.arch

#distribution of the residuals, we could apply a normality test. null is rejected
bv.norm <- normality.test(model, multivariate.only = TRUE)
bv.norm

#structural break in the residuals we can apply a CUSUM test where there does
#not appear to be a break in the respective confidence intervals
bv.cusum <- stability(model, type = "OLS-CUSUM")
plot(bv.cusum)

#note that the null hypothesis of no Granger causality is able to reject for invest
bv.cause.invest <- causality(model, cause = "invest")
bv.cause.invest
#null hypothesis of no Granger causality is rejected for income and cons
bv.cause.income <- causality(model, cause = "income")
bv.cause.income
bv.cause.cons <- causality(model, cause = "cons")
bv.cause.cons


#Since variables in a VAR model depend on each other, individual 
#coefficient estimates only provide limited information on the reaction of 
#the system to a shock. In order to get a better picture of the model’s 
#dynamic behavior, forecast error impulse response (FEIR) function are used.
irf.cons <- irf(model, impulse = "income", response = "cons",
            n.ahead = 8, ortho = FALSE, runs = 1000)

plot(irf.cons, cex.main=1.25, cex.lab=1.0, cex.axis=0.75)

#to see shock effect of invest to consumption
irf.cons1 <- irf(model, impulse = "invest", response = "cons", 
               n.ahead = 8, boot = TRUE)
plot(irf.cons1, ylab = "consumption", main = "Shock from investment")

#to see shock effect of invest to income
irf.income <- irf(model, impulse = "invest", response = "income", 
                 n.ahead = 8, boot = TRUE)
plot(irf.income, ylab = "income", main = "Shock from investment")

#to see shock effect of cons to income
irf.income1 <- irf(model, impulse = "cons", response = "income", 
                  n.ahead = 8, boot = TRUE)
plot(irf.income1, ylab = "income", main = "Shock from consumption")

#to see shock effect of cons to invest
irf.invest <- irf(model, impulse = "cons", response = "invest", 
                  n.ahead = 8, boot = TRUE)
plot(irf.income, ylab = "investment", main = "Shock from consumption")


#Information on contemporaneous relations is rather contained in the 
#off-diagonal elements of the symmetric variance-covariance matrix Σ. 
#For the used data set it is estimated to be:
# Calculate summary statistics
model_summary <- summary(model)

# Obtain variance-covariance matrix
model_summary$covres

#Since the off-diagonal elements of the estimated variance-covariance matrix
#are not zero, we can assume that there is contemporaneous correlation 
#between the variables in the VAR model. This is confirmed by the correlation
#matrix, which corresponds to Σ:
model_summary$corres

t(chol(model_summary$covres))

#In R the irf function of the vars package can be used to optain OIRs by 
#setting the argument ortho = TRUE:
oir <- irf(model, impulse = "income", response = "cons",
           n.ahead = 8, ortho = TRUE, runs = 1000, seed = 12345)
#Note that the output of the Choleski decomposition is a lower triangular 
#matrix so that the variable in the first row will never be sensitive to a 
#contemporaneous shock of any other variable and the last variable in the 
#system will be sensitive to shocks of all other variables. Therefore, the
#results of an OIR might be sensitive to the order of the variables and it 
#is advised to estimate the above VAR model with different orders to see 
#how strongly the resulting OIRs are affected by that.
plot(oir)

#forecast error variance de-compositions (FEVD), consumption is influenced by
#shocks in income most strongly which is intuitive 
bv.vardec <- fevd(model, n.ahead = 10)
plot(bv.vardec)

#we are forecasting 8 steps ahead. We are also looking to make use of 95% 
#confidence intervals for the forecast
predictions <- predict(model, n.ahead = 8, ci = 0.95)
plot(predictions, names = "invest")
plot(predictions, names = "income")
plot(predictions, names = "cons")
fanchart(predictions, names = "invest")
fanchart(predictions, names = "income")
fanchart(predictions, names = "cons")




