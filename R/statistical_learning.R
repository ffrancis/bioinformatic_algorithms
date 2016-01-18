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






