#### fisher test #####
### library testing
##Author:Sonja_Lecic_06_2018

Ultra100r2_FS10r1_fisher_test <- read.table("/Volumes/Temp1/allele_frequencies_libs/mpileup/Ultra100_r2_FS10_r1_masked.fet", header = F, sep = "\t")

Ultra100r2_FS10r1_fisher_testd <- data.frame(chr = Ultra100r2_FS10r1_fisher_test$V1, 
                                             wind = Ultra100r2_FS10r1_fisher_test$V2,
                                             pval = Ultra100r2_FS10r1_fisher_test$V6)
Ultra100r2_FS10r1_fisher_testds <- Ultra100r2_FS10r1_fisher_testd%>%separate(pval,into=c("pair","p.val"),sep="=")
print(summary(as.numeric(Ultra100r2_FS10r1_fisher_testds$p.val)))
hist(as.numeric(Ultra100r2_FS10r1_fisher_testds$p.val))
Ultra100r2_FS10r1_fisher_testds$p.val <- as.numeric(Ultra100r2_FS10r1_fisher_testds$p.val)

### Manhattan plot
chromosomes <- c("2L", "2R", "3L", "3R", "X")
chromosomes_length <- matrix(c(21106064, 18996374, 22263316, 26981190, 20637937), nrow = 1, dimnames = list(NULL, chromosomes))
off <- 0
pval_lim <- max(as.numeric(Ultra100r2_FS10r1_fisher_testds$p.val))
png("Ultra100r2_FS10r1.png", w = 850)
for(i in seq_along(chromosomes)){
  chr <- chromosomes[i]
  if(i%%2){col <- "gray"}else{col <- "black"}
  temp <- Ultra100r2_FS10r1_fisher_testds[which(Ultra100r2_FS10r1_fisher_testds$chr == chr), ]
  print(head(temp))
  if(chr == chromosomes[1]){
    plot(temp$wind, temp$p.val, xlim = c(0, sum(chromosomes_length[1, ])), ylim = c(0, pval_lim), col = col, pch=19)
  }else{
    points(temp$wind+off, temp$p.val, col = col, pch=19)
  }
  off <- off + chromosomes_length[1, chr] 
}
dev.off()
