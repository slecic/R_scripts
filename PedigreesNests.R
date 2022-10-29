#!/usr/bin/env Rscript


# Author Sonja 04.04.2021
# This script plots pedigrees for any nest in the data set

# the code

# load kinship2 package
require(kinship2)
# make a function to plot any list of samples
samplefull <- read.csv("/Volumes/LaCie/PoolSeq_BlueTit/samplefull.csv")
samplefulluniqueall <- read.csv("/Users/sonjalecic/Documents/samplefulluniqueall.csv")
Asamplesf <- read.csv("/Volumes/LaCie/PhD/data/RRBSsamples/Asamples.csv")
Bsamplesf <- read.csv("/Volumes/LaCie/PhD/data/RRBSsamples/Bsamplesfinal.csv")
b2 <- subset(samplefull, new_objective_sample == "B2")
b1 <- subset(samplefull, new_objective_sample == "B1")
a2 <- subset(samplefull, new_objective_sample == "A2")
a1 <- subset(samplefull, new_objective_sample == "A1")

PlotNests <- function(samples, plotgrid, subsample1, subsample2, simbolsize, density1, density2) {
  
  myNest <- na.omit(unique(samples$focal_nest_box))
  plotgrid = plotgrid
  NestNumber <- c()
  for(i in myNest) {
    # fix parents
    fixpar <- fixParents(na.omit(samples[samples$focal_nest_box == i,]$ID), 
                         na.omit(samples[samples$focal_nest_box == i,]$genfather), 
                         na.omit(samples[samples$focal_nest_box == i,]$genmother), 
                         na.omit(samples[samples$focal_nest_box == i,]$sex), 
                         missid = 0)
    fixpar$affected <- 0
    fixpar$affected[match(samples[samples$epy == 1,]$ID, fixpar$id)] <- 1
    fathers <- subset(samples, ID %in% samples$genfather)
    aff2 <- data.frame(epy = fixpar$affected)
    aff2$WPsire <- 0
    aff2$EPWPsire <- 0
    aff2$EPsire <- 0
    fathers2 <- subset(subsample2, ID %in% samples$genfather)
    fathers1 <- subset(subsample1, ID %in% samples$genfather)
    aff2$WPsire[match(fathers2[fathers2$EPsire == 0,]$ID, fixpar$id)] <- 1
    aff2$WPsire[match(fathers1[fathers1$EPsire == 0,]$ID, fixpar$id)] <- 1
    aff2$EPWPsire[match(fathers2[fathers2$EPsire == 0.5,]$ID, fixpar$id)] <- 1
    aff2$EPWPsire[match(fathers1[fathers1$EPsire == 0.5,]$ID, fixpar$id)] <- 1
    aff2$EPsire[match(fathers2[fathers2$EPsire == 1,]$ID, fixpar$id)] <- 1
    aff2$EPsire[match(fathers1[fathers1$EPsire == 1,]$ID, fixpar$id)] <- 1
    
    nest_year <- unique(na.omit(samples[samples$focal_nest_box == i,]$focal_nest_year))
    
    NestNumber <- append(NestNumber, length(nest_year))
    
    secondClutch <- subset(samples, focal_nest_secondClutch == 1)$focal_nest_box
    
    ## compute the pedigree file
    ped <- with(fixpar, pedigree(id = id, dadid = dadid, 
                                 momid = momid, sex = sex, affected = as.matrix(aff2), missid = 0))
    
    # plot pedigree
    peds = plot(ped, id = ped$id, cex = 0.01, symbolsize = simbolsize, density = density1)
    title(main = paste("Nest box", i, "\n",
                       length(nest_year),"nest(s)"))
    print(peds)
    print(paste("Total number of nests is:", sum(NestNumber)))
  }
  pedigree.legend(ped, new = F, location = "bottomright", radius = 3, density = density2,
                  cex = 1.5)
  
}

# plot objective B nests
pdf("/Volumes/LaCie/PhD/figures/Bluetit/BNests.pdf",width=170,height = 120, pointsize=150)
PlotNests(samples = Bsamplesf, plotgrid = par(mfrow = c(3, 6)), subsample1 = b1, subsample2 = b2,
          simbolsize = 130, density1 = c(-1, 4, 4, 5), density2 = c(-1, 4, 4, 5))
dev.off()
# plot objective A nests
pdf("/Volumes/LaCie/PhD/figures/Bluetit/ANests.pdf",width=190,height = 120, pointsize=150)
PlotNests(samples = Asamplesf, plotgrid = par(mfrow = c(3, 7)), subsample1 = a1, subsample2 = a2,
          simbolsize = 130, density1 = c(-1, 4, 4, 5), density2 = c(-1, 4, 4, 5))
dev.off()
# plot all nests
pdf("/Volumes/LaCie/PhD/figures/Bluetit/NestsAllsamplefull.pdf",width=400,height = 600, pointsize=200)
PlotNests(samples = samplefull, plotgrid = par(mfrow = c(15, 10)), subsample1 = samplefull, 
          subsample2 = samplefull, simbolsize = 100, density1 = c(-1, 3, 3, 4), density2 = c(-1, 3, 3, 4))
dev.off()
