
## Nucleotide diversity pi analysis (sliding window)
## Sonja_Lecic_08_2019


setwd("/Volumes/Temp2/Haplotypes/Base_populations/VCF_analysis/splited_VCFs/")

pi_susceptibles <- read.delim("susceptibles.windowed.pi", h = T)

pi_other_susceptibles <- read.delim("other_susceptibles.windowed.pi", h = T)

pi_resistants <- read.delim("resistants.windowed.pi", h = T)

pi_susceptibles2 <- subset(pi_susceptibles, pi_susceptibles$BIN_START > 10000000 & pi_susceptibles$BIN_START < 14000000)
pi_other_susceptibles2 <- subset(pi_other_susceptibles, pi_other_susceptibles$BIN_START > 10000000 & pi_other_susceptibles$BIN_START < 14000000)
pi_resistants2 <- subset(pi_resistants, pi_resistants$BIN_START > 10000000 & pi_resistants$BIN_START < 14000000)

View(pi_susceptibles)


matplot(pi_susceptibles2$BIN_START, pi_susceptibles2$PI, type = "l", col="grey60", xlim = c(11500000, 12500000), ylim=c(0.00001, 0.1),log= "y", lty = 1, xlab="3R chromosome", ylab= expression(paste(pi, " diversity")),  frame=FALSE, cex.lab=1.2)
lines(pi_resistants2$BIN_START, pi_resistants2$PI, col="brown3")
lines(pi_other_susceptibles2$BIN_START, pi_other_susceptibles2$PI, col="darkolivegreen3")
dev.off()

# Add Ace gene

# Extract GTF file

gtf = read.delim("/Volumes/Temp/Reference/Dsim/annotation/JMCE01.2-exon-CDS-utr.gtf", h = FALSE)
gtf$gene = gsub("gene_id ", "", gtf$V9)
gtf$gene = gsub("; transcript_.*", "", gtf$gene)
#gtf$gene = gsub("; gene_symbol.*", "", gtf$gene)


#extract min value and max
gtf_min = aggregate(gtf$V4, by = list(gtf$gene), min)
gtf_max = aggregate(gtf$V4, by = list(gtf$gene), max)
gtf2 = unique(data.frame(gtf$V1, gtf$gene))
gtf3 = merge(gtf_min, gtf_max, by = "Group.1")
gtf3 = merge(gtf2, gtf3, by.x="gtf.gene", by.y = "Group.1")
colnames(gtf3) = c("gene", "chromosome", "min", "max")

# Add genes on the plot
for (i in 1:nrow(gtf2)){
  rect(gtf3[gtf3$chromosome == "3R",][i,3], 0.05, gtf3[gtf3$chromosome == "3R",][i,4], 0.025, col = "grey80", border="transparent")
  #text(((gtf3[gtf3$chromosome == chrom,][i,3]+ gtf3[gtf3$chromosome == chrom,][i,4])/2), -4, , gtf3[gtf3$chromosome == chrom,][i,1], pos = 1)
}

text(expression(italic("Ace")), x = 11990000, y = 0.035)

legend(x=12200000, y=0.0005, legend=c("Resistant haplotypes", "Susceptible haplotypes: 4,14,16,69,73,75", "Susceptible haplotypes: 10,18,22,24,31,52,57"), col=c("brown3", "grey60", "darkolivegreen3"),lty=c(1,1,1), bty="n", cex=0.8)

# Add a legend
legend("bottomleft", 
       legend = c("Group 1", "Group 2"), 
       col = c(rgb(0.2,0.4,0.1,0.7), 
               rgb(0.8,0.4,0.1,0.7)), 
       pch = c(17,19), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))



