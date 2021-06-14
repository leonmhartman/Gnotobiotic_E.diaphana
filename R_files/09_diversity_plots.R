library(here)
library(gridExtra)

source(here("08_outlier_deletion.R"))

# normalise by rarefying
anemones <- rarefy_even_depth(anemones, sample.size = min(sort(sample_sums(anemones))), rngseed = 1)
# rarefying depth on the BH-adjusted ASV counts = 48176 counts per sample
# 11 ASVs were removed

# get diversity values
div <- estimate_richness(physeq = anemones, measures = c("Observed"))
div <- cbind(div, estimate_richness(physeq = anemones, measures = c("Simpson")))
div <- cbind(div, estimate_richness(physeq = anemones, measures = c("Shannon")))

# add meta data
div$"treatcode" <- 1
div$treatcode[37:69] <- 2
div <- cbind(div, anemones@sam_data$Day)
names(div)[5] <- c("day")

# calculate means
means <- aggregate(div, by = div[c("day", "treatcode")], FUN = "mean")
means[2] <- rep(c("control", "gnoto"), each = 6)

# calculate sems
library(plotrix)
sems <- aggregate(div, by = div[c("day", "treatcode")], FUN = "std.error")
sems[2] <- rep(c("control", "gnoto"), each = 6)

# split data for plotting
ob <- means[1:3]
ob <- cbind(ob,sems[3])
names(ob)[2:4] <- c("group","mean","sem")

sp <- means[,c(1:2,4)]
sp <- cbind(sp,sems[4])
names(sp)[2:4] <- c("group","mean","sem")

sh <- means[,c(1:2,5)]
sh <- cbind(sh,sems[5])
names(sh)[2:4] <- c("group","mean","sem")

# generate plots with standard error bars
obg <- ggplot(ob, aes(x=day, y=mean, group=group, color=group)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1, position=position_dodge(0.5)) +
  geom_line(position=position_dodge(0.5), size=1) +
  geom_point(position=position_dodge(0.5), size=3) +
  theme_bw() +
  theme(legend.position="none") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  scale_y_continuous(breaks = seq(100, 400, 50))

spg <- ggplot(sp, aes(x=day, y=mean, group=group, color=group)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1, position=position_dodge(0.5)) +
  geom_line(position=position_dodge(0.5), size=1) +
  geom_point(position=position_dodge(0.5), size=3) +
  theme_bw() +
  theme(legend.position="none") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  ylim(0.85,0.98)

shg <- ggplot(sh, aes(x=day, y=mean, group=group, color=group)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1, position=position_dodge(0.5)) +
  geom_line(position=position_dodge(0.5), size=1) +
  geom_point(position=position_dodge(0.5), size=3) +
  theme_bw() +
  theme(legend.position="none") +
  scale_x_continuous(breaks = c(0,1,3,7,14,21))

# Print plots on one page
grid.arrange(obg, spg, shg, nrow = 1)

rm(list = ls())
