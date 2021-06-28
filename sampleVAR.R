#VAR in R using Phillipine data
#Tom Xu
#Sources:Kevin Kotzr, vector auto-regression in R

#Load required packages for VAR

library(urca)
library(vars)
library(mfilter)
library(tseries)
library(forecast)
library(tidyverse)

#load dataset
okun = read.csv(file.choose())
head(okun)

#simple growth
ggplot(data = okun) + geom_point(mapping = aes(x = unem, y = real_gdp_growth))

# declare time series variables
gdp <- ts(okun$real_gdp_growth, start = c(1999,3), frequency = 4)
unem <- ts(okun$unem, start = c(1999,3), frequency = 4)

# plot the series
autoplot(cbind(gdp, unem))

#OLS
ols1 = lm(gdp ~ unem)
summary(ols1)

#determine the persistence of the model
acf(gdp, main = "ACF for real gdp growth")
pacf(gdp, main = "PACF for real gdp growth")

acf(unem, main = "ACF for Unemployment")
pacf(unem, main = "PACF for Unemployment")

#find the optimal lags
okun.bv <- cbind(gdp, unem)
colnames(okun.bv) <- cbind("GDP", "UNEMPLOYMENT")
lagselect <- VARselect(okun.bv, lag.max = 10, type = "const")
lagselect$selection

#building model VAR
ModelOkun1 <- VAR(okun.bv, p = 4, type = "const", season = NULL, exog = NULL)
summary(ModelOkun1)

#Diagnose VAR
#Serial correlation pass when p-value > 0.05
Serial1 <- serial.test(ModelOkun1, lags.pt = 12, type = "PT.asymptotic")
Serial1

#Heteroscedasticity pass when p-value > 0.05
Arch1 <- arch.test(ModelOkun1, lags.multi = 12, multivariate.only = TRUE)
Arch1

#Normal distribution of residuals pass if p-value > 0.05 
Norm1 <- normality.test(ModelOkun1, multivariate.only = TRUE)
Norm1

#Structural breaks in the residuals
Stability1 <- stability(ModelOkun1, type = "OLS-CUSUM")
plot(Stability1)

#Granger causality: if a signal X1 "Granger-causes" (or "G-causes") a signal X2, 
#then past values of X1 should contain information that helps predict X2 above
#and beyond the information contained in past values of X2 alone
GrangerGDP <- causality(ModelOkun1, cause = "GDP")
GrangerGDP

GrangerUnemployment <- causality(ModelOkun1, cause = "UNEMPLOYMENT")
GrangerUnemployment

#impulse response functions
GDPirf <- irf(ModelOkun1, impulse = "UNEMPLOYMENT", response = "GDP", n.ahead = 20, boot = TRUE)
plot(GDPirf, ylab = "GDP", main = "Shock from Unemployment")
