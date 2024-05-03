# This is the sample size calculation for BASELINE models #
# settings: train/test = 2/1 #
# R package: pmsampsize #
# Ref: Ensor, Joie, Emma C. Martin, and Richard D. Riley. 2021. Pmsampsize: Calculates the Minimum Sample Size Required for Developing a Multivariable Prediction Model. https://CRAN.R-project.org/package=pmsampsize. #

library(rms)
library(pmsampsize)

# function provided by Riley et al. 2021 on estimating the Cox-Snell R-squared from a reported ROC-AUC #
# Riley, RD, Van Calster, B, Collins, GS. A note on estimating the Cox-Snell R2 from a reported C statistic (AUROC) to inform sample size calculations for developing a prediction model with a binary outcome. Statistics in Medicine. 2021; 40: 859– 864. https://doi.org/10.1002/sim.8806 #
approximate_R2 <- function(auc, prev, n = 1000000){
  # define mu as a function of the C statistic
  mu <- sqrt(2) * qnorm(auc)
  # simulate large sample linear prediction based on two normals
  # for non-eventsN(0, 1), events and N(mu, 1)
  LP <- c(rnorm(prev*n, mean=0, sd=1), rnorm((1-prev)*n, mean=mu, sd=1))
  y <- c(rep(0, prev*n), rep(1, (1-prev)*n))
  # Fit a logistic regression with LP as covariate;
  # this is essentially a calibration model, and the intercept and
  # slope estimate will ensure the outcome proportion is accounted
  # for, without changing C statistic
  fit <- lrm(y~LP)
  max_R2 <- function(prev){
    1-(prev^prev*(1-prev)^(1-prev))^2
  }
  return(list(R2.nagelkerke = as.numeric(fit$stats['R2']),
              R2.coxsnell = as.numeric(fit$stats['R2']) * max_R2(prev)))
}

# outcome could be binary or time-to-event #
# steps:
# 1. Choose number of predictor variables/parameters that will be included in the model.
# 2. Choose the adjusted Cox-Snell R-squared (based on previously published, internally validated model in same setting/population) or the maximal apparent x-Snell R-squared (maximum value determined by population prevalence/risk).
# 3. Calculate sample size so that the estimated Van Houwellingen’s global shrinkage factor (S_VH) is greater than some chosen constant (commonly S_VH >= 0.90).
# 4. Calculate the sample size so that the absolute difference in the adjusted nagelkerke R-squared and the apparent nagelkerke R-squared is less than some constant (commonly <= 0.05).
# 5. Calculate the sample size required so that the absolute margin of error for the population outcome risk is less than a chosen constant (commonly <= 0.05).
# 6. The minimum sample size required for the model is the maximum sample size from steps 3-5.

######################################
# what's required for general input: #
######################################
# - rsquared: specifies the expected value of the (Cox-Snell) R-squared of the new model, which can be obtained from published studies based on the function above
# NOTE: If taking a value from a previous prediction model development study, users should input the model’s adjusted R-squared value, not the apparent R-squared value, as the latter is optimistic (biased).
# NOTE: If taking the R-squared value from an external validation of a previous model, the apparent R-squared can be used (as the validation data was not used for development, and so R-squared apparent is then unbiased).
# - parameters: specifies the number of candidate predictor parameters for potential inclusion in the new prediction model. 
# NOTE: This may be larger than the number of candidate predictors, as some categorical predictors require more parameters (degrees of freedom) to be estimated.
# - shrinkage: specifies the level of shrinkage desired at internal validation after developing the new model, which measures overfitting. 
# NOTE: This can range from 0 to 1, with higher values denoting less overfitting. A shrinkage = 0.9 is recommended (the default in pmsampsize), which indicates that the predictor effect (beta coefficients) in the model would need to be shrunk by 10% to adjust for overfitting.
##########################################
# what's in addition for binary outcome: #
##########################################
# - prevalence: specifies the overall outcome proportion (for a prognostic model) or overall prevalence (for a diagnostic model) expected within the model development dataset.
# NOTE: This should be derived based on previous studies in the same population.
# - cstatistic: specifies the C-statistic reported in an existing prediction model study to be used in conjunction with the expected prevalence to approximate the Cox-Snell R-squared using the approach of Riley et al. 2021. 
# NOTE: Ideally, this should be an optimism-adjusted C-statistic. The approximate Cox- Snell R-squared value is used as described above for the rsquared() option, and so is treated as a baseline for the expected performance of the new model.
# - seed: specifies the initial value of the random-number seed used by the random-number functions when simulating data to approximate the Cox-Snell R-squared based on reported C-statistic and expect prevalence as described by Riley et al. 2021
############################################
# what's in addition for survival outcome: #
############################################
# - rate: specifies the overall event rate in the population of interest, for example as obtained from a previous study, for the survival outcome of interest. 
# NOTE: rate must be given in time units used for meanfup and timepoint options.
# - timepoint: specifies the timepoint of interest for prediction. 
# NOTE: time units must be the same as given for meanfup option (e.g. years, months).
# - meanfup: specifies the average (mean) follow-up time anticipated for individuals in the model development dataset, for example as taken from a previous study in the population of interest. 
# NOTE: time units must be the same as given for timepoint option.


# there are overall 30862 catheter episodes in the playground dataset
# thus it is reasonable to assume 20575 catheter episodes (more or less) in the data_train_BASE 

#######################################
# Binary outcome (for logistic model) #
#######################################

# based on our systematic review, the reported internally validated C-index are: 0.72, 0.82, 0.72, 0.82, 0.71, 0.76, 0.75
# mean: 0.75

# According to the UZ Leuven EHRs and CLABSI calculation algorithm, the occurrence rate of CLABSI among all catheter episodes from 2012 to 2013 is around 3.14%
pmsampsize(type = "b", cstatistic=0.75, parameters = 21, prevalence = 0.0314)

# all these above use all admission days data

# If restricted to 7-day data only
# then prevalence of CLABSI is only 0.0131
pmsampsize(type = "b", cstatistic=0.75, parameters = 21, prevalence = 0.0131)

###############################################
# Survival outcome (for Cox prediction model) #
###############################################

# approximate_R2(0.75, 0.0314, n = 1000000)
# $R2.coxsnell 0.02670039

# we select a timepoint of interest for prediction using the newly developed model of 7 days

# if using the mean follow up (19 days) and event rate (0.0314) from our UZ Leuven playground dataset
pmsampsize(type = "s", csrsquared = 0.0267, parameters = 21, rate = 0.0314, timepoint = 7, meanfup = 19)


# approximate_R2(0.75, 0.0131, n = 1000000)
# $R2.coxsnell 0.01175451

# if using the event rate 0.0131
pmsampsize(type = "s", csrsquared = 0.0118, parameters = 21, rate = 0.0131, timepoint = 7, meanfup = 19)

