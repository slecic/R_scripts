
## Simulate resistance allele trajectories based on different h
## Sonja_12_2018

library(poolSeq)
library(ggplot2)
library(reshape2)
library(data.table)


tp=c(0,15,37,59,75,89)
aftraj<-1-af.traj(sync = Sync_Ace_rising,chr=rep("3R", times=3),pos = c(11986861, 11987174, 11987483),repl=c(1:5))
aftraj[7,2] <- NA


aftraj1<-1-af.traj(sync = Sync_Ace_rising,chr=rep("3R"),pos = c(11986861),repl=c(1:5))
aftraj2<-1-af.traj(sync = Sync_Ace_rising,chr=rep("3R"),pos = c(11987174),repl=c(1:5))
aftraj3<-1-af.traj(sync = Sync_Ace_rising,chr=rep("3R"),pos = c(11987483),repl=c(1:5))

aftraj1 = aftraj1[-3,]
estimateSH(aftraj1, t=tp, Ne=219, h=NA, haploid=FALSE, method="NLS")
estimateSH(aftraj2, t=tp, Ne=219, h=NA, haploid=FALSE, method="NLS")
estimateSH(aftraj3, t=tp, Ne=219, h=NA, haploid=FALSE, method="NLS")


s.h05 <-estimateSH(aftraj, t=tp, Ne=219, h=0.5, haploid=FALSE, method="NLS")
s.h0 <-estimateSH(aftraj, t=tp, Ne=219, h=0, haploid=FALSE, method="NLS")
s.h1 <-estimateSH(aftraj, t=tp, Ne=219, h=1, haploid=FALSE, method="NLS")


# plot 
plot(1, type="n", xlim=range(tp), ylim=c(0, 0.6), xlab = "Generations", ylab = "Resistant allele frequencies")
for(i in c(1,2,3)) {
  lines(tp, colMeans(aftraj[c(i+0, i+3, i+6, i+9, i+12),], na.rm = TRUE))
}
lines(tp, wf.traj(p0=mean(aftraj[,1]), Ne=NA, t=tp, s=s.h05$s, h=s.h05$h.given), col="orange", lwd=4)
lines(tp, wf.traj(p0=mean(aftraj[,1]), Ne=NA, t=tp, s=s.h0$s, h=s.h0$h.given), col="green", lwd=4)
lines(tp, wf.traj(p0=mean(aftraj[,1]), Ne=NA, t=tp, s=s.h1$s, h=s.h1$h.given), col="red", lwd=4)

legend(x=50, y=0.55, legend=c("Dominant", "Co-dominant", "Recessive", "Mean trajectories of \n resistant mutations"), col=c("red", "orange", "green", "black"),lty=c(1,1,1), lwd = c(4, 4, 4, 1), bty="n", cex=0.9, y.intersp = 2)


