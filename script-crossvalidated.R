# Load libraries

library(RCurl)
library(lme4)
library(ggplot2)
library(tidyverse)
library(car)
library(MASS)
library(vcd)
library(DHARMa)
library(plyr)
library(fitdistrplus)
library(optimx)
library(dfoptim)

# Load data from Github

mydata = read.csv(text=getURL("https://raw.githubusercontent.com/samuel-walker/glmm-help/master/LTDB_stackexchange.csv"), header = TRUE)

# Overdispersion parameter estimation function from Ben Bolker: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor
# Comment additions from http://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html

overdisp_fun <- function(model) {
  ## number of variance parameters in an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  # The next two lines calculate the residual degrees of freedom
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  # extracts the Pearson residuals
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  # Generates a p-value. If less than 0.05, the data are overdispersed.
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# Set tract and decade as factors

mydata$TRTID10 <- factor (mydata$TRTID10)
mydata$decade <- factor (mydata$decade)

# Group mean center variables by decade

head(mydata)

scaled.mydata <- ddply(mydata, c("decade"), transform, P_NONWHT = scale(P_NONWHT), 
                                                       a_hinc = scale(a_hinc))

scaled.mydata <- na.omit(scaled.mydata)

head(scaled.mydata)

# Check missing values

table(is.na(scaled.mydata)) 
sapply(scaled.mydata, function(x) sum(is.na(x)))
table(rowSums(is.na(scaled.mydata)))
table(complete.cases(scaled.mydata))

# # Descriptives
# 
# # Histogram of dependent
# 
# p0 <- ggplot(data=mydata, aes(mydata$R_VAC))
# p0 + geom_histogram(binwidth=50)
# 
# # Boxplot of vacancy by decade; shows increasing proportion of outliers in later years
# 
# p1 <- ggplot(scaled.mydata, aes(decade, R_VAC, color=decade))
# p1 + geom_boxplot()
# 
# # Tracing tract vacancy by decade
# 
# p2 <- ggplot(scaled.mydata, aes(decade, R_VAC, group=TRTID10, color=decade))
# # A basic dot plot
# p2 + geom_point(size=1) +
#      geom_line()
# 
# # Percent vacant by decade
# 
# p3 <- ggplot(scaled.mydata, aes(decade, R_VAC/HU, group=TRTID10, color=decade))
# p3 + geom_point(size=1) +
#      geom_line() +
#      facet_grid(. ~ D_suburb)
# 
# # Percent vacant by nonwhite grouped by decade
# 
# p4 <- ggplot(mydata, aes(P_NONWHT, R_VAC/HU, color=decade))
# p4 + geom_point(size=1) +
#      facet_grid(. ~ D_suburb)
# 
# # Decade group mean vacancy by nonwhite
# 
# p5 <- ggplot(mydata, aes(P_NONWHT, R_VAC/HU, color=D_suburb))
# p5 + geom_point() + 
#      facet_grid(decade ~ .)
# 
# # Rootogram
# 
# fit <- goodfit(scaled.mydata$R_VAC) 
# summary(fit) 
# p6 <- rootogram(fit)
# 
# # Ord plot
# 
# p7 <- Ord_plot(scaled.mydata$R_VAC)
# 
# # Poisson-ness and negative binomial-ness plots, suggest Poisson.
# 
# p8 <- distplot(scaled.mydata$R_VAC, type="poisson")
# 
# p9 <- distplot(scaled.mydata$R_VAC, type="nbinomial")
# 
# # Check distributions
# 
# p10 <- plot(hist(mydata$R_VAC))
# 
# # Normal
# 
# p11 <- qqp(mydata$R_VAC, "norm")
# 
# # Negative binomial
# 
# nbinom <- fitdistr(mydata$R_VAC, "Negative Binomial")
# p12 <- qqp(mydata$R_VAC, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
# 
# # Poisson
# 
# poisson <- fitdistr(mydata$R_VAC, "Poisson")
# p13 <- qqp(mydata$R_VAC, "pois", poisson$estimate)
# 
# # Gamma
# 
# gamma <- gamma <- fitdistr(mydata$R_VAC[mydata$R_VAC!=0], densfun="gamma")
# p14 <- qqp(mydata$R_VAC, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
# 
# # Visual comparision of distribution fit using fitdistrplus
# 
# R_VAC <- mydata$R_VAC
# fitn <- fitdist(R_VAC, "norm")
# fitp <- fitdist(R_VAC, "pois")
# fitnb <- fitdist(R_VAC, "nbinom")
# 
# # Empirical distribution versus fitted distribution in CDF
# 
# p15 <- cdfcomp(list(fitn, fitp, fitnb), legendtext = c("Normal", "Poisson", "Negative Binomial"),
#         main = "Vacancy fits", xlab = "Vacant housing units", xlegend = "topright")
# 
# # Density plot with histogram
# 
# p16 <- denscomp(list(fitn, fitp, fitnb), legendtext = c("Normal", "Poisson", "Negative Binomial"),
#          main = "Vacancy fits", xlab = "Vacant housing units", xlim = c(0, 1000), xlegend = "topright")
# 
# # Probabilities of fitted versus empirical probabilities
# 
# p17 <- ppcomp(list(fitn, fitp, fitnb), legendtext = c("Normal", "Poisson", "Negative Binomial"),
#          main = "Vacancy fits", xlab = "Vacant housing units", xlegend = "topright", line01 = TRUE)
# 
# # Quantiles of theoretical distribution versus empirical quantiles
# 
# p18 <- qqcomp(list(fitn, fitp, fitnb), legendtext = c("Normal", "Poisson", "Negative Binomial"),
#        main = "Vacancy fits", xlab = "Vacant housing units", xlegend = "topright", line01 = TRUE)

# GLMMs

# Poisson

# Null model

null.p <- glmer(R_VAC ~ 1 + offset(HU_ln) + (1|decade/TRTID10), family="poisson", data=scaled.mydata)

summary(null.p)
overdisp_fun(null.p)

# Exponentiated coefficients

exp(summary(null.p)$coefficients[,"Estimate"])

# Model

m1.p <- glmer(R_VAC ~ decade + P_NONWHT + a_hinc +  P_NONWHT*a_hinc + offset(HU_ln) + (1|decade/TRTID10), family="poisson",
            data=scaled.mydata)

summary(m1.p)
overdisp_fun(m1.p)

# Exponentiated coefficients

exp(summary(m1.p)$coefficients[,"Estimate"])

out <- capture.output(summary(m1.p)+summary(null.p))

cat("My title", out, file="results.txt", sep="n", append=TRUE)

# Negative binomial

# Null model

null.nb <- glmer.nb(R_VAC ~ 1 + offset(HU_ln) + (1|decade/TRTID10), data=scaled.mydata)

summary(null.nb)
overdisp_fun(null.nb)

# Exponentiated coefficients

exp(summary(null.nb)$coefficients[,"Estimate"])

# Model

m1.nb <- glmer.nb(R_VAC ~ decade + P_NONWHT + a_hinc +  P_NONWHT*a_hinc + offset(HU_ln) + (1|decade/TRTID10), data=scaled.mydata)

summary(m1.nb)
overdisp_fun(m1.nb)

# Exponentiated coefficients

exp(summary(m1.nb)$coefficients[,"Estimate"])

# Diagnostics

# Compare fit by distribution

anova(null.p, m1.p)
anova(null.nb, m1.nb)
anova(m1.p, m1.nb)

# Model validation

# DHARMa to check residuals

# Poisson

# Generate simulated residuals

simulationOutput <- simulateResiduals(fittedModel = m1.p, n = 1000)

# Plot simulated residuals

plotSimulatedResiduals(simulationOutput = simulationOutput)

par(mfrow = c(1,3))
plotResiduals(scaled.mydata$P_NONWHT,  simulationOutput$scaledResiduals)
plotResiduals(scaled.mydata$a_hinc,  simulationOutput$scaledResiduals)
plotResiduals(scaled.mydata$P_NONWHT*scaled.mydata$a_hinc,  simulationOutput$scaledResiduals)

# Tests

# K-S test for uniformity of scaled residuals; significant = cannot reject non-uniformity (i.e. evidence of non-uniformity)

testUniformity(simulationOutput = simulationOutput)

# Overdispersion cannot do glmer.nb objects, so use Bolker code instead.

# Test for zero inflation (Observed versus expected); significant = cannot reject non-zero inflation (i.e. evidence of zero inflation). 
# Could also be caused by overdipersion, though.

testZeroInflation(simulationOutput)

# Negative binomial

# Generate simulated residuals

simulationOutput2 <- simulateResiduals(fittedModel = m1.nb, n = 1000)

# Plot simulated residuals

plotSimulatedResiduals(simulationOutput = simulationOutput2)

par(mfrow = c(1,3))
plotResiduals(scaled.mydata$P_NONWHT,  simulationOutput2$scaledResiduals)
plotResiduals(scaled.mydata$a_hinc,  simulationOutput2$scaledResiduals)
plotResiduals(scaled.mydata$P_NONWHT*scaled.mydata$a_hinc,  simulationOutput2$scaledResiduals)

# Tests

# K-S test for uniformity of scaled residuals; significant = cannot reject non-uniformity (i.e. evidence of non-uniformity)

testUniformity(simulationOutput = simulationOutput2)

# Overdispersion cannot do glmer.nb objects, so use Bolker code instead.

# Test for zero inflation (Observed versus expected); significant = cannot reject non-zero inflation (i.e. evidence of zero inflation). 
# Could also be caused by overdipersion, though.

testZeroInflation(simulationOutput2)

# Optimize model fit using allFit, see http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#troubleshooting

source(system.file("utils", "allFit.R", package="lme4"))
modelfit.all.m1.p <- allFit(m1.p)
ss1 <- summary(modelfit.all.m1.p)
ss1

source(system.file("utils", "allFit.R", package="lme4"))
modelfit.all.m1.nb <- allFit(m1.nb)
ss2 <- summary(modelfit.all.m1.nb)
ss2
