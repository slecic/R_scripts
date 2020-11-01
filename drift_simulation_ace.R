
## simulate drift and Plot resistant allele trajectories
## Sonja_11_2018

library(poolSeq)
library(ggplot2)
library(reshape2)
library(data.table)


Coverages = coverage(Sync_Ace_rising, repl = 1:5, gen = c(0,15,37,59,75,89))

rownames(Coverages) = snps$pos

Coverages2 = data.frame(Coverages)

traj <- wf.traj(p0=rep(0.4, times=1000), Ne=219, t=tp, s=0, h=0.5, haploid=FALSE)


all <- rbind(sample.alleles(traj[,1], size = as.matrix(sample(Coverages2["11986861" , grepl( "F0" , names( Coverages2 ) ) ], nrow(traj), replace = T)), mode = "coverage"),
             sample.alleles(traj[,2], size = as.matrix(sample(Coverages2["11986861" , grepl( "F15" , names( Coverages2 ) ) ], nrow(traj), replace = T)), mode = "coverage"),
             sample.alleles(traj[,3], size = as.matrix(sample(Coverages2["11986861" , grepl( "F37" , names( Coverages2 ) ) ], nrow(traj), replace = T)), mode = "coverage"),
             sample.alleles(traj[,4], size = as.matrix(sample(Coverages2["11986861" , grepl( "F59" , names( Coverages2 ) ) ], nrow(traj), replace = T)), mode = "coverage"),
             sample.alleles(traj[,5], size = as.matrix(sample(Coverages2["11986861" , grepl( "F75" , names( Coverages2 ) ) ], nrow(traj), replace = T)), mode = "coverage"),
             sample.alleles(traj[,6], size = as.matrix(sample(Coverages2["11986861" , grepl( "F89" , names( Coverages2 ) ) ], nrow(traj), replace = T)), mode = "coverage"))


traj95 <- apply(t(all), 2, FUN=function(x) quantile(x,probs=c(.025,.975)))
