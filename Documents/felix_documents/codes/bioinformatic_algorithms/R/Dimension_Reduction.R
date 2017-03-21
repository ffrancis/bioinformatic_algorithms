### PCA

library (ISLR)
attach (Carseats ) 

head (Carseats, 4)
sub_carseats <- head (Carseats, 100)
subset <- subset(Carseats,(Advertising/Population) >0.04  )
plot(subset$Population,subset$Advertising)
plot(Advertising~Population, data = subset)


x <- subset$Population
y <- subset$Advertising

plot(x, y, xlab ="Population", ylab ="Advertising spending", col = "orange", cex=.3)
abline(7, .033, col="dark green", lwd=3)
abline(55, -.229, col="dark blue", lwd=3, lty=5)



####################
### PCA example
####################
# Load data
data(iris)
head(iris, 3)

# log transform
log.ir <- log(iris[,1:4])
head(log.ir, 3)
ir.species <- iris[,5]

# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE.
ir.pca <- prcomp(log.ir, center=TRUE, scale. = TRUE)

print (ir.pca)


# plot method
plot(ir.pca, type = "l")

#summary method
summary(ir.pca)


#Predict PCs
predict(ir.pca, newdata = tail(log.ir,2))


#PCA biplot generation
#install.packages("devtools")
library(devtools)
install_github("ggbiplot", "vqv")
library(ggbiplot)
g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
              groups = ir.species, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
biplot(ir.pca)





















