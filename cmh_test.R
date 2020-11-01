
### CMH test
###Sonja_10_2018

# Load libraries

library(poolSeq)
library(ggplot2)
library(reshape2)
library(data.table)


# load sync file
Sync_Ace_rising <- read.sync("BASE_EVOLVED_ALL_Portugal_ACE_region_indelRepeatMasked_mq20.sync", gen=c(rep(0, 5), rep(15, 5), rep(37,5), rep(59,5), rep(75,5), rep(89,5)), repl=c(rep(c(1:5), 6)), polarization = "rising")



## Make the CMH F89
snps <- alleles(Sync_Ace_rising)
af0 <- af(Sync_Ace_rising, repl = 1:5, gen = 0)
cov0 <- coverage(Sync_Ace_rising, repl = 1:5, gen = 0)
A0 <- af0*cov0
a0 <- cov0 - A0

aft <- af(Sync_Ace_rising, repl = 1:5, gen = 89)
covt <- coverage(Sync_Ace_rising, repl = 1:5, gen = 89)
At <- aft*covt
at <- covt - At


rm(af0, aft, cov0, covt)

cmh_89 <- cmh.test(t(A0), t(a0), t(At), t(at), min.cov = 15, max.cov = .95, min.cnt = 2, log = T)


# Plot CMH test


data <- data.frame(snps$pos, cmh_89)

plot(data, xlab = "3R chromosome", ylim= c(0,120))




