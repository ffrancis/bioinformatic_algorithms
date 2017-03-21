########### STAT612: Advanced Regression ############

# Lecture Note 4: Cross-Validation and the Bootstrap

### Part 1: Cross-Validation

# The Validation Set Approach

library(ISLR)
set.seed(1)
train=sample(392,196) 
# randomly assign training and validation sets

# linear regression
lm.fit=lm(mpg~horsepower,data=Auto,subset=train) 
# subset=train, use training set to fit a model

attach(Auto)
mean((mpg-predict(lm.fit,Auto))[-train]^2)
# test MSE; index=-train: obs not in training set

# repeat for quadratic regression
lm.fit2=lm(mpg~poly(horsepower,2),data=Auto,subset=train)
mean((mpg-predict(lm.fit2,Auto))[-train]^2)

# repeat for cubic regression
lm.fit3=lm(mpg~poly(horsepower,3),data=Auto,subset=train)
mean((mpg-predict(lm.fit3,Auto))[-train]^2)

# set a different seed => different training and validation sets
set.seed(2)
train=sample(392,196)

# linear regression
lm.fit=lm(mpg~horsepower,subset=train) 
# already run attach(Auto)
mean((mpg-predict(lm.fit,Auto))[-train]^2)

# quadratic regression
lm.fit2=lm(mpg~poly(horsepower,2),data=Auto,subset=train)
mean((mpg-predict(lm.fit2,Auto))[-train]^2)

# cubic regression
lm.fit3=lm(mpg~poly(horsepower,3),data=Auto,subset=train)
mean((mpg-predict(lm.fit3,Auto))[-train]^2)

# Leave-One-Out Cross-Validation

# Artificial implementation of LOOCV 
attach(Auto)
n=nrow(Auto)
cv.auto=0
for(i in 1:n)
{
  lm.i.fit=lm(mpg~horsepower, data=Auto[-i,])
  cv.auto=cv.auto+(mpg - predict(lm.i.fit,Auto))[i]^2
}
cv.auto=cv.auto/n

# Automatic implementaion of LOOCV using cv.glm()
?glm
# compare glm and lm
glm.fit=glm(mpg~horsepower,data=Auto)
coef(glm.fit)
lm.fit=lm(mpg~horsepower,data=Auto)
coef(lm.fit)


library(boot)
glm.fit=glm(mpg~horsepower,data=Auto)
?cv.glm # only works for glm objects
cv.err=cv.glm(Auto,glm.fit)
cv.err$delta
cv.err$delta[1] # test MSE estimate

# compute LOOCV for polynomial regression with degrees=1-5.

# artificial

cv.error=rep(0,5)
n=nrow(Auto)
for(d in 1:5)
{
  cv.auto=0
  for(i in 1:n)
  {
    lm.i.fit=lm(mpg~poly(horsepower,d), data=Auto[-i,]) 
    # polynormial regression
    cv.auto=cv.auto+(mpg - predict(lm.i.fit,Auto))[i]^2
  }
  cv.error[d]=cv.auto/n
}



# automatic
cv.error=rep(0,5)
for (i in 1:5)
{
  glm.fit=glm(mpg~poly(horsepower,i),data=Auto)
  cv.error[i]=cv.glm(Auto,glm.fit)$delta[1]
}
cv.error
which.min(cv.error)

# K-Fold Cross-Validation

# artificial

# K=10 for linear regression

K=10

n=nrow(Auto)
foldsize=floor(n/K) # size of each fold

foldidx=list(idx=0)
foldidx=rep(foldidx,K)

idx=sample(n) # randomize indices

for(k in 1:(K-1)) # split to K folds
{
  foldidx[[k]]=idx[(1+(k-1)*foldsize):(k*foldsize)]
}
foldidx[[K]]=idx[(1+(K-1)*foldsize):n]

cv.kfold=0

### the remiaining ones are test set
### need to record test mse for each evaluation
### for each k, use all those, except case 1 as training set
for(k in 1:K)
{
  lm.k.fit=lm(mpg~horsepower, data=Auto[-foldidx[[k]],]) # foldidx[[k]] will return index for case set; so negative of it means use all rows except these indices
  cv.kfold=cv.kfold+sum((mpg - predict(lm.k.fit,Auto))[foldidx[[k]]]^2)/length(foldidx[[k]]) ### predictove fn to calc predicted value
}
cv.kfold=cv.kfold/K


# automatic for linear regression

?cv.glm # input option K
glm.fit=glm(mpg~horsepower,data=Auto) # 1st obtain glm object; use that obj as one input; then specify k=10
cv.10fold=cv.glm(Auto,glm.fit,K=10)$delta[1]


# automatic for polynomial regression for degree=1-10


set.seed(17)
# K=10
cv.error.10=rep(0,10) # allocate new memory
for (i in 1:10) # poly degree = 1-10
{
  glm.fit=glm(mpg~poly(horsepower,i),data=Auto) # glm fn to fit a polynomial regression
  cv.error.10[i]=cv.glm(Auto,glm.fit,K=10)$delta[1]
}
cv.error.10
plot(cv.error.10,type='l')
which.min(cv.error.10)  # to find which model has the lowest cv value


# we can see that linear model as smallest cv value...when deg val >=8 an inflation of cv value, since the model is more complicated

# K=5
set.seed(9)
cv.error.5=rep(0,10)
for (i in 1:10)
{
  glm.fit=glm(mpg~poly(horsepower,i),data=Auto)
  cv.error.5[i]=cv.glm(Auto,glm.fit,K=5)$delta[1]
}
cv.error.5
plot(cv.error.5,type='l')
which.min(cv.error.5)

### this is only half of the story, since we only decide which degree is the best...
### we still need to decide which coefficient to use
summary(lm(mpg~poly(horsepower,6),data=Auto))


### the best model may not include the predictor you are interested in 
### minimizing test mse and finding the stat model may not overlap; but each can provide some usefulness



#### Part 2: Permutation Tests
# widely used to compare the characteristics of two populations
# all below illustrations based on simulations


# normal population (pretent that we have 2 pop, both normally dist; same variance; may/mayno have same mean)

x.null=rnorm(100) # randomly generate some obs from norm dist; mean = 0; sd =1; 100 random numbers
y.null=rnorm(200) # independenty generate another vector of 200 random numbers, from the same distribution

?rnorm  # generates a random number from a normal disctribution
?dnorm # value of desity function
?pnorm # probability for a given point for a given distribution
?qnorm # you tell a probabilty and we ll get he left tail and right tail probability



qqnorm(x.null) # check normality
qqline(x.null)
### when the obs are from same normal pop, all points should be close to the diagonal line (as seen in the figure)

sd(x.null)
sd(y.null)

shapiro.test(x.null) # Shapiro-Wilk normality test; will give you a test stat and p value; p-value should be > 0.05 to proove normality

qqnorm(y.null)
qqline(y.null)

shapiro.test(y.null)

t.test(x.null, y.null, var.equal=TRUE)
### here p-vale is large...so null hypothesis is true; and 2 pop means are the same



x.alt=rnorm(100)
y.alt=rnorm(200,mean=0.5)
t.test(x.alt, y.alt, var.equal=TRUE)


# use permuation test
# permute a large number of times; for each permutations, we first combine 2 vectors and compute their index and then split into two vectors
# x and y can have different lengths
oneperm=function(x,y)
{
  n=length(x)
  # m=length(y)
  z=c(x,y)    # combine vectors x and y; which is vector z
  z.perm=sample(z) # permutation once; sample function to permute the vector itself here (before it was used to permute index, remeber??)
  # z.perm is permutation of z
  mean(z.perm[1:n]) - mean(z.perm[-(1:n)])  # another vector of same length as z, but elements arranged separtaely; then compare the mean to two subsets (of same length) from z 
  # test statistic: difference bewtween sample means
  
  # mean(z.perm[1:n]) - mean(z.perm[(n+1):(n+m)])
}


oneperm(x.null,y.null)
# we need to repeat this B times; and allocate the values of rep distribtion as given below

# null is true
B=999
null.dist=rep(NA,B)
for(b in 1:B)
{
  null.dist[b]=oneperm(x.null,y.null)
}  

null.obs=mean(x.null)-mean(y.null)

hist(null.dist) # hist of ref distribution
abline(v=null.obs, lwd=2, col="purple") # draw a line to show the observed d (test statistic)

mean(abs(null.dist) > abs(null.obs)) # p-value; since p-value is >0.05, we accept the null hypothesis    ;the 2 samples genetrated from 2 norm pop with normal means...we should not reject null hypothesis if our procedure is correct

# null is false; ie two pop have diff pop means
# use x.alt=rnorm(100); y.alt=rnorm(200,mean=0.5)
B=999
null.dist=rep(NA,B)
for(b in 1:B)
{
  null.dist[b]=oneperm(x.alt,y.alt)
}  

alt.obs=mean(x.alt)-mean(y.alt)

hist(null.dist)
abline(v=alt.obs, lwd=2, col="purple")

mean(abs(null.dist) > abs(alt.obs)) # p-value

# here p val is small; so reject null hypothesis; observed test stat is different from the mean(0)




# use R function replicate()
?replicate
B=999
null.dist=replicate(B, oneperm(x.null,y.null))

hist(null.dist)
abline(v=null.obs, lwd=2, col="purple")

mean(abs(null.dist) > abs(null.obs))


null.dist=replicate(B, oneperm(x.alt,y.alt))

hist(null.dist)
abline(v=alt.obs, lwd=2, col="purple")

mean(abs(null.dist) > abs(alt.obs)) # p-value


# Cauchy population: CANNOT use t.test (why?)

xc.null=rcauchy(100)
yc.null=rcauchy(200)

qqnorm(xc.null) # check normality
qqline(xc.null)

shapiro.test(xc.null) # Shapiro-Wilk normality test


oneperm.med=function(x,y)
{
  n=length(x)
  # m=length(y)
  z=c(x,y)
  z.perm=sample(z) # permutation once
  
  median(z.perm[1:n]) - median(z.perm[-(1:n)]) 
  # test statistic: difference bewtween sample means
  
  # mean(z.perm[1:n]) - mean(z.perm[(n+1):(n+m)])
}

B=999
null.dist=replicate(B, oneperm.med(xc.null,yc.null))

null.obs=median(xc.null)-median(yc.null)

hist(null.dist)
abline(v=null.obs, lwd=2, col="purple")

mean(abs(null.dist) > abs(null.obs)) # p-value



xc.alt=rcauchy(100)
yc.alt=rcauchy(200,location=0.5)

B=999
null.dist=replicate(B, oneperm.med(xc.alt,yc.alt))

null.obs=median(xc.alt)-median(yc.alt)

hist(null.dist)
abline(v=null.obs, lwd=2, col="purple")

mean(abs(null.dist) > abs(alt.obs)) # p-value



#### Part 3: The Bootstrap

library(boot)

alpha.fn=function(data,index){
  X=data$X[index]
  Y=data$Y[index]
  return((var(Y)-cov(X,Y))/(var(X)+var(Y)-2*cov(X,Y)))
}
alpha.fn(Portfolio,1:100)
set.seed(1)
alpha.fn(Portfolio,sample(100,100,replace=T))
?boot
boot(Portfolio,alpha.fn,R=1000)

# Estimating the Accuracy of a Linear Regression Model

boot.fn=function(data,index) # resample observations: linear regression
  return(coef(lm(mpg~horsepower,data=data,subset=index)))
boot.fn(Auto,1:392)
set.seed(1)
boot.fn(Auto,sample(392,392,replace=T))
boot.fn(Auto,sample(392,392,replace=T))

bootauto.linear=boot(Auto,boot.fn,1000)
names(bootauto.linear)
bootauto.linear$t0
bootauto.linear$t
dim(bootauto.linear$t) # bootstrap statistics

apply(bootauto.linear$t,2,mean) # bootstrap mean for each column
apply(bootauto.linear$t,2,mean)-bootauto.linear$t0 # bias
apply(bootauto.linear$t,2,sd) # bootstrap sd for each column, i.e., SE

apply(bootauto.linear$t,2,quantile, probs=c(0.025, 0.975)) # 95% CI

summary(lm(mpg~horsepower,data=Auto))$coef
boot.fn=function(data,index)# resample observations: quadratic regression
  coefficients(lm(mpg~horsepower+I(horsepower^2),data=data,subset=index))
set.seed(1)
boot(Auto,boot.fn,1000)
summary(lm(mpg~horsepower+I(horsepower^2),data=Auto))$coef


