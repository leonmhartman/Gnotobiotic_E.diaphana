library(here)
library(magrittr)

source(here("08_outlier_deletion.R"))

# normalise by rarefying
anemones <- rarefy_even_depth(anemones, sample.size = min(sort(sample_sums(anemones))), rngseed = 1)
# rarefying depth on the scaled reads = 48176 reads per sample

# get diversity values
div <- estimate_richness(physeq = anemones, measures = c("Observed"))
div <- cbind(div, estimate_richness(physeq = anemones, measures = c("Simpson")))
div <- cbind(div, estimate_richness(physeq = anemones, measures = c("Shannon")))

# add meta data
div$"treatcode" <- 1
div$treatcode[37:69] <- 2
div <- cbind(div, anemones@sam_data$Day)
names(div)[5] <- c("day")

# change row order so days are in chronological order within each condition
library(data.table)
setorder(div, cols = "day")
setorder(div, cols = "treatcode")

div$group <- c(rep(1:6, each = 6),rep(7:8, each = 5),rep(9:11, each = 6),rep(12, 5))
div$group <- factor(div$group)
div$uniq <- 1:69

# change treatcode to factor
div$treatcode <- factor(div$treatcode)

# - - - - - - - - - - - - - - - - - - - #

# observed

library(nlme)
modo <- gls(Observed ~ treatcode * day, data = div, weights=varIdent(form = ~1 | group))

library(car)
Anova(modo)
# Analysis of Deviance Table (Type II tests)
# 
# Response: Observed
#               Df  Chisq   Pr(>Chisq)    
# treatcode      1  62.216  3.077e-15 ***
# day            1  42.561  6.853e-11 ***
# treatcode:day  1  31.000  2.580e-08 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# tests for data normality and homogeneity of variance
for (i in c(0,1,3,7,14,21)) {
  vvv <- subset(div,div$day==i) %>% leveneTest(Observed ~ treatcode,.)
  print(i)
  print(vvv)
  uuu <- subset(div,div$day==i & div$treatcode=="1")[,1] %>% shapiro.test(.) # control
  print(i)
  print(uuu)
  ttt <- subset(div,div$day==i & div$treatcode=="2")[,1] %>% shapiro.test(.) # treated
  print(i)
  print(ttt)
}

# Mann-Whitney U-tests used due to data deviation from normality and homogeneity of variance
# default argument in the wilcox.test is 'paired = FALSE', which runs the Mann-Whitney U-test
for (i in c(0,1,3,7,14,21)) {
  yyy <- subset(div,div$day==i) %>% wilcox.test(Observed ~ treatcode, data = .)
  print(i)
  print(yyy)
}

# day 0 vs day 21 treated
div[c(37:41,65:69),] %>%  wilcox.test(Observed ~ day, data = .)
# p-value = 0.007937

# day 0 vs day 21 control
div[c(1:6,31:36),] %>%  wilcox.test(Observed ~ day, data = .)
# p-value = 1

# - - - - - - - - - - - - - - - - - - - #

# simpsons

modsi <- gls(Simpson ~ day * treatcode, data = div, weights=varIdent(form = ~1 | group))

Anova(modsi)
# Analysis of Deviance Table (Type II tests)
# 
# Response: Simpson
#               Df   Chisq     Pr(>Chisq)   
# day            1   5.9341    0.01485 * 
# treatcode      1   0.0128    0.90984   
# day:treatcode  1   7.8049    0.00521 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# check homogeneity of variance & normality
for (i in c(0,1,3,7,14,21)) {
  vvv <- subset(div,div$day==i) %>% leveneTest(Simpson ~ treatcode,.)
  print(i)
  print(vvv)
  uuu <- subset(div,div$day==i & div$treatcode=="1")[,2] %>% shapiro.test(.) # control
  print(i)
  print(uuu)
  ttt <- subset(div,div$day==i & div$treatcode=="2")[,2] %>% shapiro.test(.) # treated
  print(i)
  print(ttt)
}

# Mann-Whitney U-tests used due to data deviation from normality and homogeneity of variance
# default argument in the wilcox.test is 'paired = FALSE', which runs the Mann-Whitney U-test
for (i in c(0,1,3,7,14,21)) {
  yyy <- subset(div,div$day==i) %>% wilcox.test(Simpson ~ treatcode, data = .)
  print(i)
  print(yyy)
}

# day 0 vs day 21 treated
div[c(37:41,65:69),] %>%  wilcox.test(Simpson ~ day, data = .)
# p-value = 0.05556

# day 0 vs day 21 control
div[c(1:6,31:36),] %>%  wilcox.test(Simpson ~ day, data = .)
# p-value = 0.4848

# - - - - - - - - - - - - - - - - - - - #

# shannon

modsh <- gls(Shannon ~ treatcode * day, data = div, weights=varIdent(form = ~1 | group))

Anova(modsh)
# Analysis of Deviance Table (Type II tests)
# 
# Response: Shannon
#               Df   Chisq    Pr(>Chisq)    
# treatcode      1   2.1801   0.1398056    
# day            1   10.3143  0.0013200 ** 
# treatcode:day  1   13.3926  0.0002526 ***
#  ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# check homogeneity of variance & normality
for (i in c(0,1,3,7,14,21)) {
  vvv <- subset(div,div$day==i) %>% leveneTest(Shannon~treatcode,.)
  print(i)
  print(vvv)
  uuu <- subset(div,div$day==i & div$treatcode=="1")[,3] %>% shapiro.test(.) # control
  print(i)
  print(uuu)
  ttt <- subset(div,div$day==i & div$treatcode=="2")[,3] %>% shapiro.test(.) # treated
  print(i)
  print(ttt)
}

# Mann-Whitney U-tests used due to data deviation from normality and homogeneity of variance
# default argument in the wilcox.test is 'paired = FALSE', which runs the Mann-Whitney U-test
for (i in c(0,1,3,7,14,21)) {
  yyy <- subset(div,div$day==i) %>% wilcox.test(Shannon ~ treatcode, data = .)
  print(i)
  print(yyy)
}

# day 0 vs day 21 treated
div[c(37:41,65:69),] %>%  wilcox.test(Shannon ~ day, data = .)
# p-value = 0.007937

# day 0 vs day 21 control
div[c(1:6,31:36),] %>%  wilcox.test(Shannon ~ day, data = .)
# p-value = 0.4848

rm(list = ls())
