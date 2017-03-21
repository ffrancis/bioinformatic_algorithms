### Classification
### Logistic Regression
library(MASS)
library(ISLR)

names(Smarket)
dim(Smarket)
summary(Smarket)
head(Smarket)
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


coef(glm.fit) # to assess the coeffients of the fitted model
summary(glm.fit)$coef # can also be used to the above purpose







names(Default)






### plot Balance vs Income
plot(Default$balance, Default$income, xlab="Balance", ylab="Income", 
     main="Balance vs Income", pch=3, cex.main=1.5,
     frame.plot=FALSE, col=ifelse(Default$default=='Yes', "blue","grey"))
legend(2300, 65500, pch=c(3,3), col=c("grey", "blue"), c("Default:No", "Default:Yes"),
            bty="o",  box.col="darkgreen", cex=.8)

### Box plot
par(mfrow=c(1,2))
boxplot(Default$income~Default$default, notch=TRUE, col=(c("grey","blue")),  xlab="Default", ylab="Income")
boxplot(Default$balance~Default$default, notch=TRUE, col=(c("grey","blue")),  xlab="Default", ylab="Balance")






##############################################
### plot the probability of Default vs Balance
##############################################

# creates a quartz window with title
quartz(title="Default vs Balance") 

# plot with body size on x-axis and survival (0 or 1) on y-axis
plot(Default$balance,Default$default,xlab="Balance",ylab="Probability of default")
# change Yes=1; No=0
default_boolean <- ifelse(Default$default=="Yes", 1, 0)

###run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
g=glm(default_boolean~Default$balance,family=binomial,Default) 
balance <- Default$balance
# draws a curve based on prediction from logistic regression model
curve(predict(g,data.frame(balance=x),type='resp'),add=TRUE) 



# sample data


# This set of codes will produce plots for logistic regression. Text that follows # sign is ignored by R when running commands, so you can just copy-and-paste these straight into your R console or R document.

#First, we'll create a fake dataset of 20 individuals of different body sizes:
bodysize=rnorm(20,30,2) # generates 20 values, with mean of 30 & s.d.=2
bodysize=sort(bodysize) # sorts these values in ascending order. 
survive=c(0,0,0,0,0,1,0,1,0,0,1,1,0,1,1,1,0,1,1,1) # assign 'survival' to these 20 individuals non-randomly... most mortality occurs at smaller body size
dat=as.data.frame(cbind(bodysize,survive))

#quartz(title="bodysize vs. survival") # creates a quartz window with title

plot(bodysize,survive,xlab="X",ylab="p(X)", col = 'red') # plot with body size on x-axis and survival (0 or 1) on y-axis
g=glm(survive~bodysize,family=binomial,dat) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm

curve(predict(g,data.frame(bodysize=x),type="resp"),add=TRUE, col='blue') # draws a curve based on prediction from logistic regression model

points(bodysize,fitted(g),pch=20, col='blue')


#################
### Logistic Regression based prediction using ISLR stock market data
#################




