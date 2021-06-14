library(here)
library(mvabund)

source(here("08_outlier_deletion.R"))

# statistical comparison of bacterial communities in the Day 0 and Day 21 anemones
# comparisons:
# control Day 0 vs Day 21
# treated Day 0 vs Day 21
# control Day 0 vs treated Day 0
# control Day 21 vs treated Day 21

# subset Day 0 and Day 21 anemones and remove zero sum taxa
anemones <- subset_samples(anemones, SampleType == "Anemone" & Day == "0" | Day == "21")
anemones <- prune_taxa((taxa_sums(anemones) > 0), anemones)

# collapse data to Species
glomSpec <- tax_glom(anemones, 'Species', NArm = FALSE)

# get the OTU table & transpose it
otuTab <- t(glomSpec@otu_table)

# subset the metadata
mvaMeta <- data.frame(glomSpec@sam_data)[,2:3]

mvaMeta$Day <- factor(mvaMeta$Day)

# combine metadata & OTU table
combo <- cbind(mvaMeta, otuTab)

# # make mvabund object
# comboMva <- mvabund(combo[,3:565])
# is.mvabund(comboMva)

# - - - - - - - -

# control Day 0 vs Day 21

comboC0C21 <- subset(combo, Treatment == "Control")

comboC0C21mva <- mvabund(comboC0C21[,3:565])

is.mvabund(comboC0C21mva)

mod <- manyglm(comboC0C21mva ~ comboC0C21$Day, family="negative_binomial")

plot(mod)

anova(mod)

# Time elapsed: 2 hr 19 min 35 sec
# Analysis of Deviance Table
# 
# Model: comboC0C21mva ~ comboC0C21$Day
# 
# Multivariate test:
#                Res.Df   Df.diff     Dev    Pr(>Dev)   
# (Intercept)        11
# comboC0C21$Day     10         1    1126    0.003 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Arguments:
# Test statistics calculated assuming uncorrelated response (for faster computation) 
# P-value calculated using 999 iterations via PIT-trap resampling.

# - - - - - - - -

# treatment Day 0 vs Day 21

comboT0T21 <- subset(combo, Treatment == "Treatment")

comboT0T21mva <- mvabund(comboT0T21[,3:565])

is.mvabund(comboT0T21mva)

mod <- manyglm(comboT0T21mva ~ comboT0T21$Day, family="negative_binomial")

plot(mod)

anova(mod)

# Time elapsed: 0 hr 36 min 26 sec
# Analysis of Deviance Table
# 
# Model: comboT0T21mva ~ comboT0T21$Day
# 
# Multivariate test:
#                Res.Df  Df.diff   Dev     Pr(>Dev)   
# (Intercept)         9                         
# comboT0T21$Day      8        1   1428    0.006 **
#  ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Arguments:
# Test statistics calculated assuming uncorrelated response (for faster computation) 
# P-value calculated using 999 iterations via PIT-trap resampling.

# - - - - - - - -

# control Day 0 vs treated Day 0

comboC0T0 <- subset(combo, Day == "0")

comboC0T0mva <- mvabund(comboC0T0[,3:565])

is.mvabund(comboC0T0mva)

mod <- manyglm(comboC0T0mva ~ comboC0T0$Treatment, family="negative_binomial")

plot(mod)

anova(mod)

# Time elapsed: 0 hr 55 min 19 sec
# Analysis of Deviance Table
# 
# Model: comboC0T0mva ~ comboC0T0$Treatment
# 
# Multivariate test:
#                     Res.Df    Df.diff    Dev    Pr(>Dev)
# (Intercept)             10                     
# comboC0T0$Treatment      9          1    492    0.187
# Arguments:
# Test statistics calculated assuming uncorrelated response (for faster computation) 
# P-value calculated using 999 iterations via PIT-trap resampling.

# - - - - - - - -

# control Day 21 vs treated Day 21

comboC21T21 <- subset(combo, Day == "21")

comboC21T21mva <- mvabund(comboC21T21[,3:565])

is.mvabund(comboC21T21mva)

mod <- manyglm(comboC21T21mva ~ comboC21T21$Treatment, family="negative_binomial")

plot(mod)

anova(mod)

# Time elapsed: 1 hr 38 min 4 sec
# Analysis of Deviance Table
# 
# Model: comboC21T21mva ~ comboC21T21$Treatment
# 
# Multivariate test:
#                       Res.Df   Df.diff    Dev    Pr(>Dev)    
# (Intercept)               10                          
# comboC21T21$Treatment      9         1   1750    0.001 ***
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Arguments:
# Test statistics calculated assuming uncorrelated response (for faster computation) 
# P-value calculated using 999 iterations via PIT-trap resampling.

# - - - - - - - -

rm(list = ls())
