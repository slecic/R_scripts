#!/usr/bin/env Rscript

relatives <- readRDS("/Users/sonjalecic/Desktop/rereadyfortakeoff/relatives.RDS")

# make a function that takes the list of relatives and plots each group relatives as a pedigree plot
PlotRelatives <- function(samples, plotgrid, simbolsize, density1, density2) {
  samplefulluniqueall <- read.csv("/Users/sonjalecic/Documents/samplefulluniqueall.csv")
  plotgrid = plotgrid
  # select a list of relatives with unique elements
  rel <- relatives$relatives
  relunique <- unique(lapply(rel, sort))
  # because relatives may contain some individuals that are not 
  #present in our samples we have to subset the sample list
  uniquerels <- c()
  for(i in 1:length(relunique)) {
    relativsamp <- subset(samplefulluniqueall, ID %in% relunique[[i]])
    family <- c(unique(na.omit(relativsamp$ID)), unique(na.omit(relativsamp$genfather)), unique(na.omit(relativsamp$genmother)))
    famunique <- unique(family)
    uniquerels[[i]] <- famunique
  }
  # then, again select a list of relatives with unique groups of relatives
  reluni <- unique(lapply(uniquerels, sort))
  # subset the samples list and plot if the number of focal_nest_year is greater than 4. This give us a deeper pedigree, with 
  # more that 4 generations to check if 'true' transgenerational epigenetic inheritance beyond grandparents is present and not an inheritance pattern due
  # 'transgenerational genetic effects')
  for(i in 1:length(reluni)) {
    relativ <- subset(samplefulluniqueall, ID %in% reluni[[i]])
    if(length(unique(na.omit(relativ$focal_nest_year))) >= 6) {
      relativmulti = relativ[, c(1, 12, 13, 14, 15, 16)]
      fixrel <- fixParents(relativmulti$ID, relativmulti$genfather, relativmulti$genmother, relativmulti$sex)
      fixrel$affected <- 0
      fixrel$affected[match(relativmulti[relativmulti$epy == 1,]$ID, fixrel$id)] <- 1
      aff2 <- data.frame(epy = fixrel$affected)
      fathers <- subset(relativmulti, ID %in% relativmulti$genfather)
      aff2$WPsire <- 0
      aff2$EPWPsire <- 0
      aff2$EPsire <- 0
      aff2$WPsire[match(fathers[fathers$EPsire == 0,]$ID, fixrel$id)] <- 1
      aff2$EPWPsire[match(fathers[fathers$EPsire == 0.5,]$ID, fixrel$id)] <- 1
      aff2$EPsire[match(fathers[fathers$EPsire == 1,]$ID, fixrel$id)] <- 1
      relped <- with(fixrel, pedigree(id = id, dadid = dadid, 
                                      momid = momid, sex = sex, affected = as.matrix(aff2), missid = 0))
      relpeds = plot(relped, cex = 0.01, symbolsize = simbolsize, density = density1)
      title(main = paste("Relatives"))
      print(relpeds)
    }
    
  }
  pedigree.legend(relped, new = F, location = "bottomright", radius = 3, density = density2,
                  cex = 1.5)
}

pdf("/Volumes/LaCie/PhD/figures/Bluetit/relativesEpy.pdf",width=500,height = 450, pointsize=83)
PlotRelatives(samples = samplefulluniqueall, plotgrid = par(mfrow = c(15, 13)), simbolsize = 17, 
              density1 = c(-1, 20, 20, 20), density2 = c(-1, 5, 5, 5))
dev.off()



