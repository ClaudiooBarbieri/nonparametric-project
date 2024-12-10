
covariate <- read.csv("../datasets/LUCAS_workon.csv")
LC_legend <- read.csv("../datasets/LC_description.csv")

data <- subset(covariate , select = c(BDsample_0, LC_0))

data$LC_0 <- as.factor(data$LC_0)
head(data)

table(data$LC_0)
g <- nlevels(data$LC_0)
n <- dim(data)[1]
plot(data$LC_0, data$BDsample_0, xlab='treat',col=rainbow(g),main='Original Data')

#test gaussianity
# Set up an 8-panel plotting layout (2 rows, 4 columns)
par(mfrow = c(2, 4)) # 2 rows and 4 columns of plots

# Loop through each category and create the Q-Q plots
for (i in c("A", "B", "C", "D", "E", "F", "G", "H")) {
  qqnorm(data$BDsample_0[data$LC_0 == i], 
         main = paste("LC_0 =", i)) # Title for each plot
  qqline(data$BDsample_0[data$LC_0 == i]) # Adds a reference line
}

# Reset the plotting layout to default (single plot per window)
par(mfrow = c(1, 1))
p_values=NULL
for (i in c("A", "B", "C", "D", "E", "F","H")) {
  p_values=c(p_values, shapiro.test(data$BDsample_0[data$LC_0 == i]))
}
p_values

# H0: tau1 = tau2 = tau3 = tau4 = tau5 = tau6 = tau7 = tau8 = 0
# the BDs belong to the same population

# H1: (H0)^c
# the BDs belong to several different population
# Parametric test:
fit <- aov(BDsample_0 ~ LC_0, data = data)
summary(fit)
# Permutation test:
# Test statistic: F stat
T0 <- summary(fit)[[1]][1,4]
T0

# what happens if we permute the data?
permutazione <- sample(1:n)
weight_perm <- data$BDsample_0[permutazione]
fit_perm <- aov(weight_perm ~ data$LC_0)
summary(fit_perm)

plot(data$LC_0, weight_perm, xlab='treat',col=rainbow(g),main='Permuted Data')


# CMC to estimate the p-value
B <- 1000 # Number of permutations
T_stat <- numeric(B) 

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  weight_perm <- data$BDsample_0[permutation]
  fit_perm <- aov(weight_perm ~ data$LC_0)
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}

layout(1)
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)

plot(ecdf(T_stat))
abline(v=T0,col=3,lwd=4)

# p-value
p_val <- sum(T_stat>=T0)/B
p_val
# we reject the null hypothesis
