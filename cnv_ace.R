
## Check for CNVs in the Ace region using CnMops package
## Sonja_02_2019

## Using R software

# Libraries to load
library(IRanges)
library(cn.mops)
library(data.table)



# Load BAM file list
BAMfiles <- c(
  list.files(pattern = "Dsim_Portugal*.*.bam.sorted$") 
) 

# List of chromosomes (we have only 3R in our bam files)

Chromosomes = c("3R")

# Do the CNV analysis 
CNVs = NULL
CNVRegions = NULL

for (i in 1:length(Chromosomes)){
  
  bamDataRanges <- getReadCountsFromBAM(BAMfiles, refSeqName= Chromosomes[i], WL = 1000)
  res <- cn.mops(bamDataRanges)
  result <- calcIntegerCopyNumbers(res)
  
  segm <- as.data.frame(segmentation(result))
  CNVs <- rbind(CNVs, as.data.frame(cnvs(result)))
  CNVRegions <- rbind(CNVRegions, as.data.frame(cnvr(result)))
}




# Create a list of haplotypes with numbers (ordered by name)
List = data.frame(BAMfiles)
List = sort(List$BAMfiles)
List$number = c(1:nrow(List))
# merge the CNVs and haplotype list
CNVs2 = merge(CNVs, List, by.x = "sampleName", by.y = "List")



## Variables to be changed for the plot
chrom = "3R"
min = 11960000
max = 12050000



# Plot CNVs along the 3R chromosome

plot(NULL, xlim=c(min,max), xlab = "Coordinates 3R chromosomes", ylab = "", main = "", frame=FALSE, xaxt="n", yaxt="n", ylim= c(0, 11))
axis(1)
for (i in 1:nrow(CNVs)) {
  lines(c(CNVs2[i,3], CNVs2[i,4]), c(CNVs2[i,10], CNVs2[i,10]), lwd = 5)
  text(CNVs2[i,3], CNVs2[i,10], label = paste("Haplo ", CNVs2[i,10]), pos = 2, cex = 0.8)
  
}

