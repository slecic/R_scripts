
## Selection coeeficients of the Ace region
## Sonja_Lecic_10_2018

library(poolSeq)
library(ggplot2)
library(reshape2)
library(data.table)


timepoints = c(0,15,37,59,75,89)

# get positions
temp <- alleles(Sync_Ace_rising)
my.pos <- temp[, 1:2]
rm(temp)

# get coverage
all.cov <- coverage(Sync_Ace_rising, chr = my.pos$chr, pos= my.pos$pos, repl = 1:5, gen = timepoints)
cov.id <- splitLocusID(row.names(all.cov), ".")
cov.names <- paste0("F", rep(timepoints, times=5), ".R", rep(1:5, each=6), ".cov")

library(doParallel)
registerDoParallel(20) #nb of cores
system.time({
  results <- foreach(l = 1:nrow(my.pos), .combine=rbind) %dopar% { 
    
    # extract trajectories
    my.traj <- af.traj(Sync_Ace_rising, my.pos$chr[l], my.pos$pos[l], repl = 1:5)
    
    # extract coverages
    my.cov <- all.cov[which(my.pos$chr[l] == cov.id$chr & my.pos$pos[l] == cov.id$pos),]
    my.cov <- matrix(my.cov[cov.names], ncol=6, byrow=TRUE)
    
    # estimate s    
    my.sel <- estimateSH(my.traj, t = timepoints, Ne = 219, h = .5, haploid = F,
                         N.ctraj = 0, simulate.p.value = T, N.pval = 1000, 
                         cov = my.cov, approximate = T, method = "LLS")
    
    return(c(s = my.sel$s, pval = my.sel$p.value))
  }
})

results_good <- cbind(as.data.frame(my.pos), as.data.frame(results))


### Plot results
options(scipen=5)

resuls_filtered = results_good[results_good$pval <= 0.01,]

plot(resuls_filtered$pos, abs(resuls_filtered$s), col = "grey60", ylab = "Selection coefficient", xlim=c(11950000, 12050000), ylim=c(0.02, 0.10), xlab = "Chromosome 3R position", frame=FALSE, xaxt="n", yaxt="n", pch=19, cex.lab = 1.2)

axis(1)
axis(2)

points(resuls_filtered$pos[resuls_filtered$pos >= "11971691" & resuls_filtered$pos <= "12007976"], resuls_filtered$s[resuls_filtered$pos >= "11971691" & resuls_filtered$pos <= "12007976"], col = "grey20", pch = 19)

points(resuls_filtered$pos[resuls_filtered$pos == "11987174" | resuls_filtered$pos == "11986861" | resuls_filtered$pos =="11987483"], 
       resuls_filtered$s[resuls_filtered$pos == "11987174" | resuls_filtered$pos == "11986861" | resuls_filtered$pos =="11987483"],
       col = "red", pch = 19)


