library(here)
library(magrittr)
library(plotrix)
library(ggplot2)

# import data
bh <- read.table(here("ddPCR_counts.txt"),
                 header = TRUE,
                 sep = "\t")


# - - - - bacteria - - - - #

# calculate BH ratios (values corrected by subtracting NTC background signal)
bh$ratio <- bh$Bcorrected/bh$Hcorrected

# remove gt215
bh <- bh[-72,]

# calculate means
bact <- aggregate(bh, by = bh[c("day", "treatCode")], FUN = "mean") %>% .[,c(1,2,11)]

# calculate sems
bact$sem <- aggregate(bh[,c(3,4,9)], by = bh[,c(3,4,9)][c("day", "treatCode")], FUN = "std.error") %>% .[,5]

# add names for plotting
bact$treatCode <- c(rep("control",6),rep("treated",6))

# # plot of treated only
# ggplot(bact[7:12,], aes(x=day, y=ratio)) +
#   geom_errorbar(aes(ymin=ratio-sem, ymax=ratio+sem), width=.1, position=position_dodge(0.5)) +
#   geom_line(position=position_dodge(0.5), size=1) +
#   geom_point(position=position_dodge(0.5), size=3) +
#   theme_bw() +
#   theme(legend.position="none") +
#   scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
#   ylim(0,0.12)

# plot of controls & treated
ggplot(bact, aes(x=day, y=ratio, group=treatCode, color=treatCode)) + 
  geom_errorbar(aes(ymin=ratio-sem, ymax=ratio+sem), width=.1, position=position_dodge(0.5)) +
  geom_line(position=position_dodge(0.5), size=1) +
  geom_point(position=position_dodge(0.5), size=3) +
  theme_bw() +
  theme(legend.position=c(0.1,0.85)) +
  scale_x_continuous(breaks = c(0,1,3,7,14,21)) +
  scale_y_continuous(breaks = seq(0, 2, by = 0.1))


# - - - - artemia - - - - #

artemia <- bh[73:78,c(2,9)]

meanControl <- mean(artemia$ratio[1:3])
meanTreated <- mean(artemia$ratio[4:6])
means <- c(meanControl,meanTreated)
semsControl <- std.error(artemia$ratio[1:3])
semsTreated <- std.error(artemia$ratio[4:6])
sems <- c(semsControl,semsTreated)
treatCode <- c("control","treated")
art <- data.frame(treatCode,means,sems)

ggplot(art, aes(fill=treatCode, y=means, x=treatCode)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=means-(sems), ymax=means+(sems)), width=.1) +
  theme_bw() +
  theme(legend.position=c(0.9,0.9)) +
  ylim(0,0.12)

rm(list = ls())
