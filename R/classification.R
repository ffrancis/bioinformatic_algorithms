### Classification
### Logistic Regression
library(MASS)
library(ISLR)

names(Smarket)
dim(Smarket)
summary(Smarket)

# produces a matrix that contains all of the pairwise correlations among the predictors in a data set
cor(Smarket[,-9])
attach(Smarket)
plot(Volume)

### use Logistio Regression to predict Direction using Lag1 through Lag5 and Volume
# The syntax generalized of the glm() function is similar to that of lm(), 
# except that we must pass in linear model the argument family=binomial in order to tell R 
# to run a logistic regression rather than some other type of generalized linear model.

glm.fit=glm(Direction~Lag1+Lag2+Lag3+Lag4+Lag5+Volume, data=Smarket, family=binomial)
summary(glm.fit)


names(Default)

### plot Balance vs Income
plot(Default$balance, Default$income, xlab="Balance", ylab="Income", 
     main="Balance vs Income", pch=3, cex.main=1.5,
     frame.plot=FALSE, col=ifelse(Default$default=='Yes', "blue","grey"))
legend(2300, 65500, pch=c(3,3), col=c("grey", "blue"), c("Default:No", "Default:Yes"),
            bty="o",  box.col="darkgreen", cex=.8)

### Box plot
boxplot(Default$income~Default$default, col=(c("grey","blue")))
boxplot(Default$balance~Default$default, col=(c("grey","blue")))








