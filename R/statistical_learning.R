?matrix

x=matrix(data=c(1,2,3,4), nrow=2, ncol=2)
x

x=matrix(data=c(1,2,3,4), nrow=2, ncol=2,byrow=TRUE)
x

sqrt(x) # only entry wise; not real sqrt of matrix x

x%*%x   # matrix multiplication
t(x)    # transpose


### for random num generation
# 1st specify the distribution it follows
?rnorm
x=rnorm(50)
y=x+rnorm(50,mean=50,sd=.1)
cor(x,y)

rnorm(10)
rnorm(10)

set.seed(53)
rnorm(10)
set.seed(53)
rnorm(10)


mean(y)
var(y)
sqrt(var(y))


x=rnorm(100)
y=rnorm(100)

x
y

plot(x,y)
plot(y~x)
?plot
plot(x,y,xlab ="xlabel", ylab = "ylabel", main="title")

### save fig in directory
plot(x,y,xlab ="xlabel", ylab = "ylabel", main="title")
pdf("Figure.pdf")
plot(x,y, col='green')
dev.off()



x=1:10
x

x=seq(-pi,pi,length=50)
y=x
f=outer(x,y,function(x,y)cos(y)/(1+x^2))
dim(f)
?outer
contour(x,y,f)
contour(x,y,f, nlevels=45,add=T)

fa=(f-t(f)/2)
contour(x,y,fa, nlevels=45,add=T)
image(x,y,fa)
?image


persp(x,y,fa)
persp(x,y,fa,theta=30,phi=20)
#####


### indexing data
A=matrix(1:16,4,4)
A
A[2,3]
A[c(1,3),c(2,4)]
A[1:3,2:4]

nrow(A)
ncol(A)
length(A)

?read.csv

Auto=read.table("Auto.data.txt")
### change working dir to folder where u have data
### Click: Session > Set working dir > choose working dir
### Alternatively use import function


#####################################################
### writing your own function in R
functionname = function(input1, input2, ... )
{
    statements
    return(outout) or output
}

### function to calculat the sum of a vector:
myfun1 = function(x) # one input, one output
{
    sum(x)
}

b = c(1,4,2)
myfun1(b)

### forloop ex in R
myfun3 = function(x) # one input, one output
{
    n = length(x)
    x0 = 0
    for(i in 1:n)
    {
        x0=x0+x[i]
    }
    x0
}

b = c(1,4,2)
myfun3(b)


### fucntion with 2 inputs
### matrix product

myfun4 = function(x,y)
{
    x%*%y
    
}

A = matrix(1:4,2,2)
B = matrix(1:9,2,3)
myfun4(A,B)

myfun4 = function(x,y) ### not functional!!!
{
    if (is.matrix(x) == FALSE | is.matrix(y) == FALSE)
    {
        return(Check your input matrix)
    }
    
    else if (ncol(x)!=nrow(y))
    {
        return(NULL)
    }
    
    else
    {
        x%*%y
    }
    
}
A = matrix(1:4,2,2)
B = matrix(1:9,3,3)
A = c(4,2,2)
myfun4(A,B)
#####################################################


#####################################################
### one input, two outputs
### return mean, sd

myfun6 = function(x)
{
    list(x, sum(x))
    
}

a = c(1,4,2)
myfun6(a)

z = myfun6(a)
class(z)

z[[1]]
z[[2]]


### naming each component instead of list indexing
myfun7 = function(x)
{
    list(INPUT=x, SUMMATION=sum(x))
    
}

a = c(1,4,2)
myfun7(a)

z = myfun7(a)
class(z)

z
length(z)

z[[1]]
z[[2]]  # or use
z$INPUT
z$SUMMATION

#####################################################

#####################################################
### INSTALLING R PACKAGES
?install.packages
### LOAD A PACKAGE
library(MASS)
library(ISLR)

?Boston
fix(Boston )
names(Boston )
attach (Boston )
lm.fit =lm(medv~lstat,data=Boston)
summary(lm.fit)
plot(medv~lstat)
lm.fit =lm(medv~lstat)
plot(lm.fit)

lm.fit =lm(medv~lstat,data=Boston)
plot(lstat ,medv)
abline (lm.fit)

par(mfrow =c(2,2))
plot(lm.fit)
par(mfrow =c(1,1))
par(mfrow =c(1,2))

plot(predict (lm.fit), residuals (lm.fit))
plot(predict (lm.fit), rstudent (lm.fit))
plot(hatvalues (lm.fit ))
which.max (hatvalues (lm.fit)) #leverage values: should be between 0 and 1...largest one here is still small.so its a relative result..



### residual standard error is an estimate of sigma..as long as its small ur model has a good fit
### Multple R-squared: amount of variation explained by the regression
### F-stat value sould be the same as lstat(slope) Standard error value

### information lm.fit contains
names(lm.fit)

lm.fit$coefficients
# or
coef(lm.fit)
### confidence interval  ; by default is 95% CI
?confint
confint(lm.fit)
confint(lm.fit, level=0.99)


### CI for regression function, after we fit lm
?predict
predict(lm.fit, data.frame(lstat=(c(5,10,15))), interval = "confidence")

## CI is for fixed things
## PREDICTION interval is for random things
predict(lm.fit, data.frame(lstat=(c(5,10,15))), interval = "prediction")


###########################################

### multiple linear regression fitting
lm.fit0 =lm(medv~lstat+age ,data=Boston )   ### instead of ~ (for simple linear), for multiple lin reg we use lstat+age
summary (lm.fit0)

# deg of freedon = n-3; here its 503

dim(Boston)
# n = 506; so df n-3 = 503

# RESIDUAL ST ERROR can be used to est sigma
# deg of freeding for s tat is diff: here 2 and 503; 2= no of paramters-1; 503 = n-(no of parameters)=n-3 = 506-3



#if f stat tells that teh reg function is not sig, then no point in looking at t test or anything else
# Since p value for age is small, it means that: In the presence of lstat, age is a useful predictor


lm.fit =lm(medv~.,data=Boston )   ### to include all predictors
summary (lm.fit)

### see f stats and df changes..but pvalue still small
### since crim (crime) has small p value, it means given all the other predictors crim is still significant...and
### crime would not be a good predictor to remove to make the model simple


### comparing two models:
lm.fit0 =lm(medv~lstat+age ,data=Boston )   ## only two predictors
summary (lm.fit0)

lm.fit =lm(medv~.,data=Boston )   ### to include all predictors
summary (lm.fit)


anova(lm.fit0, lm.fit)

### partial f test : additional inf explained by the whole model, divide by vaurance
### so it might be useful to use lm.fit (with more predictors) rather than lm.fit0

library(car)
## to detect collinearity
vif(lm.fit) # vif larger than 5 or 10 (more conservative) means the regression you have with this predictor will have R square more than 20% (when vif >5)
## > vif means ??? pros no collinearity, cons: vey hard to explain


### interaction terms

summary(lm(medv~lstat*age, data=Boston)) # better than the one given below; gives same results
summary(lm(medv~lstat+age+lstat:age, data=Boston))

# F-statistic says function is sig; interaction is sig if alpha < 0.05
# so should we remove age? since age does not seems sig, in the presence of others..but since the interaction term
# is a higher order term than linear term. here since interaction term is sig, we should include even if the liner term is not sig
# so dont remove age

# for 3 way interactions:
summary(lm(medv~lstat*age*rm, data=Boston))
# surprisingly all signific
# as long as 3 way int is sig, you keep all others even in individual level is not sig



### non linear transformation of predictors:
# if liner fit is bad: transform perdictor or transform responce
# when to transoform responce: when error terms are non-normal

# here to include quadratic term
lm.fit2 = lm(medv~lstat+I(lstat^2))
summary (lm.fit2) ## quadratic term is significant
anova(lm.fit, lm.fit2)
# t test for quadratic term is equivanet to this




plot(lm.fit2)
# residual vs fitted
# much betetr than linear fit
# normal qq still not satisfactory; so prob transofrm the response to make it as normally distributed as possible: it could also be cause by > sample size(then it may always look bad)
# if any points outside the top right corner red line, they are highly infulential


## if quadratif fit does not work, use cubic fit

# here all

## quadratic/cubic may introduce collinearity prob
# becasue if ur obs is within a small range say0-1
# to avoid this, r has poly function (orthognal polynomials) ### by cebtering the predictors; used to avoid collinearity problem
?poly
lm.fit33 = lm(medv~poly(lstat,3))

summary(lm.fit33)
### the coefficients might be different, but the significane test will give same results








###########################################
###########################################
###########################################
# CHAPTER 1
advertising <- read.csv("http://www-bcf.usc.edu/~gareth/ISL/Advertising.csv")

#head(advertising, n=3)

names(advertising)

library(ggplot2)
library(ggthemes)
library(gridExtra)
library(grid)

p1 <- ggplot(advertising, aes(TV, Sales))+ geom_point(shape=1,color="red") +
geom_rangeframe() + theme_tufte() + geom_point(alpha = 0.4) +
stat_smooth(method = "loess")
p2 <- ggplot(advertising, aes(Radio, Sales))+ geom_point(shape=1,color="red") +
geom_rangeframe() + theme_tufte() + geom_point(alpha = 0.4) +
stat_smooth(method = "loess")
p3 <- ggplot(advertising, aes(Newspaper, Sales))+ geom_point(shape=1,color="red") +
geom_rangeframe() + theme_tufte() + geom_point(alpha = 0.4) +
stat_smooth(method = "loess")

grid.arrange(p1, p2, p3, ncol=3)


### Ch2

income <- read.csv("http://www-bcf.usc.edu/~gareth/ISL/Income2.csv")
names(income)
#install.packages("scatterplot3d")
library(scatterplot3d)
attach(income)

scatterplot3d(Education, Seniority, Income , main="3D Scatterplot")


#install.packages("Rcmdr")
library(Rcmdr)
attach(mtcars)
scatter3d(wt, disp, mpg)

scatter3d(Education, Seniority, Income)


### prediction accuracy vs model interpretability

plot(1:10, xaxt = "n", xlab='Some Letters')
axis(1, at=1:10, labels=letters[1:10])

plot(1:2, xaxt = "n", xlab='Some Letters')
axis(1, at=1:2, labels=letters[1:2])


### clustering
library(cluster)
#install.packages("fpc")
library(fpc)

data(iris)
dat <- iris[, -5] # without known classification
clus <- kmeans(dat, centers=3)# Kmeans clustre analysis
clusplot(dat, clus$cluster, color=TRUE, shade=TRUE,
labels=2, lines=0)

### Latex equation, image, url, html generator: http://www.codecogs.com/latex/eqneditor.php
### Assessing model accuracy
### Measuring quality of fit
### MSE: MSE = \frac{1}{n} \sum _{i=1}^{n} (y_{i} - \hat{f}(x_{i}))^2


### simulate artificial data for logistic regression
# create data
set.seed(666)
x1 = rnorm(0, 100)
x2 = rnorm(100, 200)



x1 = runif(100000, min=0, max=100)
x2 = runif(100000, min=2, max=20)
plot(x1,x2)


z = 1 + 2*x1 + 3*x2        # linear combination with a bias
pr = 1/(1+exp(-z))         # pass through an inv-logit function
y = rbinom(1000,1,pr)      # bernoulli response variable
#now feed it to glm:
df = data.frame(y=y,x1=x1,x2=x2)
glm( y~x1+x2,data=df,family="binomial")






