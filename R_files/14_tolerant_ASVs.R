library(here)
library(dplyr)
library(plotrix)

source(here("08_outlier_deletion.R"))

# subset
anemones <- subset_samples(anemones, SampleType != "Artemia" & Treatment == "Treatment")
anemones <- prune_taxa((taxa_sums(anemones) > 0), anemones)

# extract otu table
otuTab <- as(otu_table(anemones), Class = "matrix")

# transpose otu table
otuTab <- t(otuTab)

# convert otu table to data frame
otuTab <- as.data.frame(otuTab)

# add Day
otuTab$Day <- anemones@sam_data$Day

# calculate means
means <- round(aggregate(otuTab, by = otuTab[c("Day")], FUN = "mean"))

# calculate sems
sems <- round(aggregate(otuTab, by = otuTab[c("Day")], FUN = "std.error"))

# remove redundant columns
means <- means[,-(ncol(means))]
means <- means[,-1]
sems <- sems[,-(ncol(sems))]
sems <- sems[,-1]

# transpose
means <- t(means)
sems <- t(sems)

# add col names
colnames(means) <- c("Day0","Day1","Day3","Day7","Day14","Day21")
colnames(sems) <- c("Day0","Day1","Day3","Day7","Day14","Day21")

# make data frames
means <- data.frame(means)
sems <- data.frame(sems)


# identify tolerant ASVs

# retain only ASVs with absolute abundance ≥1000 at every timepoint
minVal <- 1000
tolerant <- which(means$Day0 >= minVal & means$Day1 >= minVal & means$Day3 >= minVal & means$Day7 >= minVal & means$Day14 >= minVal & means$Day21 >= minVal)
means <- means[tolerant,]
sems <- sems[tolerant,]
  
# get taxa names
vvv <- match(row.names(means), row.names(gnotoTax))
means <- cbind(means, gnotoTax[vvv,2:7])


# which ASVs were also present in the antibiotic-treated artemia?
artemia <- subset_samples(phy, SampleType == "Artemia" & Treatment == "Treatment")
artemia <- prune_taxa((taxa_sums(artemia) > 0), artemia)
intersect(row.names(means), rownames(artemia@tax_table))
# [1] "dc984a63d7d685f878587c30e0d8f18f" "ae8a6381f4e78d50e3fad0e2faad8ea0"
# [3] "11bb58be038678b5e2b1878c01140c05" "53818a706e38c1584f139d2f90fbd8df"
# [5] "107e04ddc4519998efd0c7f0d4a3c619" "acd191a5f307d0579b3126166f102435"



# did the abundance of the tolerant ASVs change significantly from
# Day 0 to Day 21 in the treated anemones?
# calculate using Mann-Whitney U test and add values to means data frame

counts <- otuTab[tolerant]
counts <- counts[c(1:5, 17:21),]
means$day0_21 <- 1:16

for (i in c(1:16)) {
  yyy <- wilcox.test(counts[1:5,i], counts[6:10,i])
  means$day0_21[i] <- yyy[[3]]
}


# write data to file
write.table(means, file='tolerant_ASVs.tsv', sep='\t')


# plots

# subset data for each ASV and create plot file

vibrio1 <- as.data.frame(t(means[1,1:6]))
colnames(vibrio1) <- "mean"
vibrio1$sem <- as.integer(sems[1,1:6])
vibrio1$Day <- c(0,1,3,7,14,21)

p1 <- ggplot(vibrio1, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (vibrio1)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


gamma2 <- as.data.frame(t(means[2,1:6]))
colnames(gamma2) <- "mean"
gamma2$sem <- as.integer(sems[2,1:6])
gamma2$Day <- c(0,1,3,7,14,21)

p2 <- ggplot(gamma2, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (gamma2)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


sphing3 <- as.data.frame(t(means[3,1:6]))
colnames(sphing3) <- "mean"
sphing3$sem <- as.integer(sems[3,1:6])
sphing3$Day <- c(0,1,3,7,14,21)

p3 <- ggplot(sphing3, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (sphing3)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


thala4 <- as.data.frame(t(means[4,1:6]))
colnames(thala4) <- "mean"
thala4$sem <- as.integer(sems[4,1:6])
thala4$Day <- c(0,1,3,7,14,21)

p4 <- ggplot(thala4, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (thala4)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


cox5 <- as.data.frame(t(means[5,1:6]))
colnames(cox5) <- "mean"
cox5$sem <- as.integer(sems[5,1:6])
cox5$Day <- c(0,1,3,7,14,21)

p5 <- ggplot(cox5, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (cox5)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


alter6 <- as.data.frame(t(means[6,1:6]))
colnames(alter6) <- "mean"
alter6$sem <- as.integer(sems[6,1:6])
alter6$Day <- c(0,1,3,7,14,21)

p6 <- ggplot(alter6, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (alter6)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


vibrio7 <- as.data.frame(t(means[7,1:6]))
colnames(vibrio7) <- "mean"
vibrio7$sem <- as.integer(sems[7,1:6])
vibrio7$Day <- c(0,1,3,7,14,21)

p7 <- ggplot(vibrio7, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (vibrio7)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


eryth8 <- as.data.frame(t(means[8,1:6]))
colnames(eryth8) <- "mean"
eryth8$sem <- as.integer(sems[8,1:6])
eryth8$Day <- c(0,1,3,7,14,21)

p8 <- ggplot(eryth8, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (eryth8)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


oligo9 <- as.data.frame(t(means[9,1:6]))
colnames(oligo9) <- "mean"
oligo9$sem <- as.integer(sems[9,1:6])
oligo9$Day <- c(0,1,3,7,14,21)

p9 <- ggplot(oligo9, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (oligo9)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


sm2d10 <- as.data.frame(t(means[10,1:6]))
colnames(sm2d10) <- "mean"
sm2d10$sem <- as.integer(sems[10,1:6])
sm2d10$Day <- c(0,1,3,7,14,21)

p10 <- ggplot(sm2d10, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (sm2d10)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


marino11 <- as.data.frame(t(means[11,1:6]))
colnames(marino11) <- "mean"
marino11$sem <- as.integer(sems[11,1:6])
marino11$Day <- c(0,1,3,7,14,21)

p11 <- ggplot(marino11, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (marino11)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


rhodo12 <- as.data.frame(t(means[12,1:6]))
colnames(rhodo12) <- "mean"
rhodo12$sem <- as.integer(sems[12,1:6])
rhodo12$Day <- c(0,1,3,7,14,21)

p12 <- ggplot(rhodo12, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (rhodo12)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


sapro13 <- as.data.frame(t(means[13,1:6]))
colnames(sapro13) <- "mean"
sapro13$sem <- as.integer(sems[13,1:6])
sapro13$Day <- c(0,1,3,7,14,21)

p13 <- ggplot(sapro13, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (sapro13)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


steno14 <- as.data.frame(t(means[14,1:6]))
colnames(steno14) <- "mean"
steno14$sem <- as.integer(sems[14,1:6])
steno14$Day <- c(0,1,3,7,14,21)

p14 <- ggplot(steno14, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (steno14)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


vibrio15 <- as.data.frame(t(means[15,1:6]))
colnames(vibrio15) <- "mean"
vibrio15$sem <- as.integer(sems[15,1:6])
vibrio15$Day <- c(0,1,3,7,14,21)

p15 <- ggplot(vibrio15, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (vibrio15)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


teras16 <- as.data.frame(t(means[16,1:6]))
colnames(teras16) <- "mean"
teras16$sem <- as.integer(sems[16,1:6])
teras16$Day <- c(0,1,3,7,14,21)

p16 <- ggplot(teras16, aes(x=Day, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0, color = "#FF6400") +
  geom_line(size=.7, color = "#FF6400") +
  geom_point(size=1.1, color = "#FF6400") +
  theme_bw() +
  labs(y = "16S/H × 10^3", x = "Day (teras16)") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  theme(text = element_text(family = "Arial Narrow",size = 8))


# print plots
library(gridExtra)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,
             nrow = 4)
# export to PDF

rm(list = ls())
