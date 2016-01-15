# CHAPTER 1
advertising <- read.csv("http://www-bcf.usc.edu/~gareth/ISL/Advertising.csv")

head(advertising, n=3)

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
