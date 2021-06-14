library(here)
library(magrittr)
library(plotrix)

# import data
bh <- read.table(here("ddPCR_counts.txt"),
                 header = TRUE,
                 sep = "\t")


# - - - - bacteria - - - - #

# calculate BH ratios (values were corrected by subtracting NTC background signal)
bh$ratio <- bh$Bcorrected/bh$Hcorrected

# subset
bh <- bh[c(2:71,73),c(1:4,9)]

# drop redundant levels
bh <- droplevels(bh)

library(nlme)
# build model
mod <- gls(ratio ~ day * treatment, data = bh)

library(car)
Anova(mod)
#               Df   Chisq     Pr(>Chisq)    
# day            1   14.524    0.0001384 ***
# treatment      1   22.641    1.953e-06 ***
# day:treatment  1   21.505    3.529e-06 ***

# there was a sig dif by time & treatment, and a sig interaction.


# tests for data normality and homogeneity of variance
for (i in c(0,1,3,7,14,21)) {
  vvv <- subset(bh,bh$day==i) %>% leveneTest(ratio~treatment,.)
  print(i)
  print(vvv)
  uuu <- subset(bh,bh$day==i & bh$treatCode=="1")[,5] %>% shapiro.test(.) # control
  print(i)
  print(uuu)
  ttt <- subset(bh,bh$day==i & bh$treatCode=="2")[,5] %>% shapiro.test(.) # treated
  print(i)
  print(ttt)
}

# Mann-Whitney U-tests used due to data deviation from normality and homogeneity of variance
# default argument in the wilcox.test is 'paired = FALSE', which runs the Mann-Whitney U-test
for (i in c(0,1,3,7,14,21)) {
  yyy <- subset(bh,bh$day==i) %>% wilcox.test(ratio ~ treatment, data = .)
  print(i)
  print(yyy)
}

# day 0 vs day 1 treated
bh[c(37:42,43:48),] %>%  wilcox.test(ratio ~ day, data = .)
# p-value = 0.3095

# day 0 vs day 21 treated
bh[c(37:42,67:71),] %>%  wilcox.test(ratio ~ day, data = .)
# p-value = 0.0303

# day 0 vs day 1 control
bh[c(1:6,7:12),] %>%  wilcox.test(ratio ~ day, data = .)
# p-value = 0.8182

# day 0 vs day 21 control
bh[c(1:6,31:36),] %>%  wilcox.test(ratio ~ day, data = .)
# p-value = 0.04113


# - - - - artemia - - - - #

# import data
bh <- read.table(here("ddPCR_counts.txt"),
                 header = TRUE,
                 sep = "\t")

# calculate BH ratios (values were corrected by subtracting NTC background signal)
bh$ratio <- bh$Bcorrected/bh$Hcorrected

bh <- bh[74:79,c(2,9)]

leveneTest(ratio ~ treatment, bh)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df   F value  Pr(>F)
# group  1   1.2934   0.3189

bh[1:3,2] %>% shapiro.test(.)
# p = 0.3363

bh[4:6,2] %>% shapiro.test(.)
# p = 0.3772

# Student's t-test used because data normality and homogeneity of variance criteria were met
t.test(bh[1:3,2], bh[4:6,2], var.equal=TRUE)
# Two Sample t-test
# t = 3.3314, df = 4, p-value = 0.02907


rm(list = ls())
