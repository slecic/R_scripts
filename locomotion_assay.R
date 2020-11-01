###Sopntaneous_locomotion_assay
###three_ace_haplotypes
###Sonja_Lecic_04_2019

###### CHECK THE FIRST SEVEN HOURS #####

#########################################################
###### RAMPING: 18C(12h), 28C(2h), 33C(3h), 35C(2h) #######
### 3 classes; 5 lines per class #######

monitor11 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor11.txt", sep="\t", header = F)
monitor21 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor21.txt", sep="\t", header = F)
monitor31 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor31.txt", sep="\t", header = F)
monitor41 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor41.txt", sep="\t", header = F)
monitor51 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor51.txt", sep="\t", header = F)
monitor12 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor12.txt", sep="\t", header = F)
monitor22 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor22.txt", sep="\t", header = F)
monitor32 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor32.txt", sep="\t", header = F)
monitor42 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor42.txt", sep="\t", header = F)
monitor52 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor52.txt", sep="\t", header = F)
monitor13 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor13.txt", sep="\t", header = F)
monitor23 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor23.txt", sep="\t", header = F)
monitor33 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor33.txt", sep="\t", header = F)
monitor43 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor43.txt", sep="\t", header = F)
monitor53 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping182833353days5lines20fliesperline/Monitor53.txt", sep="\t", header = F)



mon11=monitor11[10 : 240,]
mon21=monitor21[10 : 240,]
mon31=monitor31[10 : 240,]
mon41=monitor41[10 : 240,]
mon51=monitor51[10 : 240,]
mon12=monitor12[10 : 240,]
mon22=monitor22[10 : 240,]
mon32=monitor32[10 : 240,]
mon42=monitor42[10 : 240,]
mon52=monitor52[10 : 240,]
mon13=monitor13[10 : 240,]
mon23=monitor23[10 : 240,]
mon33=monitor33[10 : 240,]
mon43=monitor43[10 : 240,]
mon53=monitor53[10 : 240,]

##### Day 1 #####
mon111=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon11[,1])/20)){
  a=as.data.frame(t(colSums(mon11[seq(20*i-4,20*i,1),11:42]))) #this part sums up passes for 10 minutes
  b=as.data.frame(mon11[20*i,c(2,3,4,5,12)]) #this part sums up passes for 10 minutes
  cc=cbind(b,a)
  colnames(cc)=colnames(mon111)
  mon111=rbind(mon111,cc)
}

mon211=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon21[,1])/20)){
  a=as.data.frame(t(colSums(mon21[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon21[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon211)
  mon211=rbind(mon211,cc)
}

mon311=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon31[,1])/20)){
  a=as.data.frame(t(colSums(mon31[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon31[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon311)
  mon311=rbind(mon311,cc)
}

mon411=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon41[,1])/20)){
  a=as.data.frame(t(colSums(mon41[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon41[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon411)
  mon411=rbind(mon411,cc)
}

mon511=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon51[,1])/20)){
  a=as.data.frame(t(colSums(mon51[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon51[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon511)
  mon511=rbind(mon511,cc)
}

##### Day 2 #####
mon112=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon12[,1])/20)){
  a=as.data.frame(t(colSums(mon12[seq(20*i-4,20*i,1),11:42]))) #this part sums up passes for 10 minutes
  b=as.data.frame(mon12[20*i,c(2,3,4,5,12)]) #this part sums up passes for 10 minutes
  cc=cbind(b,a)
  colnames(cc)=colnames(mon112)
  mon112=rbind(mon112,cc)
}

mon212=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon22[,1])/20)){
  a=as.data.frame(t(colSums(mon22[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon22[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon212)
  mon212=rbind(mon212,cc)
}

mon312=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon32[,1])/20)){
  a=as.data.frame(t(colSums(mon32[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon32[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon312)
  mon312=rbind(mon312,cc)
}

mon412=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon42[,1])/20)){
  a=as.data.frame(t(colSums(mon42[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon42[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon412)
  mon412=rbind(mon412,cc)
}

mon512=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon52[,1])/20)){
  a=as.data.frame(t(colSums(mon52[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon52[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon512)
  mon512=rbind(mon512,cc)
}


###### Day 3 #####
mon113=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon13[,1])/20)){
  a=as.data.frame(t(colSums(mon13[seq(20*i-4,20*i,1),11:42]))) #this part sums up passes for 10 minutes
  b=as.data.frame(mon13[20*i,c(2,3,4,5,12)]) #this part sums up passes for 10 minutes
  cc=cbind(b,a)
  colnames(cc)=colnames(mon113)
  mon113=rbind(mon113,cc)
}

mon213=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon23[,1])/20)){
  a=as.data.frame(t(colSums(mon23[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon23[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon213)
  mon213=rbind(mon213,cc)
}

mon313=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon33[,1])/20)){
  a=as.data.frame(t(colSums(mon33[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon33[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon313)
  mon313=rbind(mon313,cc)
}

mon413=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon43[,1])/20)){
  a=as.data.frame(t(colSums(mon43[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon43[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon413)
  mon413=rbind(mon413,cc)
}

mon513=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon53[,1])/20)){
  a=as.data.frame(t(colSums(mon53[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon53[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon513)
  mon513=rbind(mon513,cc)
}

head(mon513)


positions=read.delim("/Users/sonjalecic/Desktop/Ace manuscript/locomotion_test2.txt",head=F,stringsAsFactors = F) # deswegen besser alle slots bef??llen #die excel tabelle kopieren und in ein txt file einf??gen
positions <- positions[, c(5, 4, 1, 2)]
colnames(positions) <- c("sampleID", "monitor", "row", "genotype")
positions$slot = c(rep(c(1:32), 5))
positions$sampleID = c(1:160)
positions <- positions[, c(1, 2, 5, 4, 3)]
positions = positions[, -5]


monall=cbind(mon111,mon211,mon311,mon411,mon511,
             mon112,mon212,mon312,mon412,mon512,
             mon113,mon213,mon313,mon413,mon513) 


monall=monall[,-c(grep("time", colnames(monall))[-1])] #removes extra time/light bla bla columns
monall=monall[,-c(grep("day", colnames(monall))[-1])]
monall=monall[,-c(grep("month", colnames(monall))[-1])]
monall=monall[,-c(grep("year", colnames(monall))[-1])]
monall=monall[,-c(grep("light", colnames(monall))[-1])]
colnames(monall)[-1:-5]=paste(positions$genotype,sep = "_")
head(monall)


#png("/Volumes/Temp1/shengkai/behavior_assay/spontaneous_locomotion.png",width = 24,height = 12,pointsize = 10,units = "cm",res=600)
#par(mfrow=c(1,2))
plot(NA,NA, ylim=c(0,7),xlim=c(0,length(monall[,1])), xaxt="n", ylab="Passes/10min",xlab="Time (hr)",type="n",axes=F, cex.lab=1.6)
#axis(1, 1:length(monall[,1]))
axis(1, at=c(seq(0,length(monall[,1]),6)+1),labels=0:2, cex.axis=1.6)
axis(2, cex.axis=1.6)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),6,6),substr(colnames(monall),7,7))=="Ib"],1,mean)),col="orange",lwd=3)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),6,6),substr(colnames(monall),7,7))=="II"],1,mean)),col="grey40",lwd=3)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),6,6),substr(colnames(monall),7,7))=="Ia"],1,mean)),col="red",lwd=3)
legend("topleft", c("Class Ia", "Class Ib", "Class II"), col=c("red", "orange","grey40"), lty=1,bty="n", lwd = 3, cex=1)
points(c(1,30,42,60,73),c(18,18,28,33,35)*6.9/35,col="black", type="l",lty=2, lwd = 2.5, cex=1)
text(c(1,28,40,58,73),(c(18,18,28,33,35)+1)*6.9/35,paste(c(18,18,28,33,35),"C"), cex=1.2, col="black")
text(c(15,33,49,65),(c(18,23,31,34)+1)*6.9/35,paste("Stage",1:4),cex=1, col="grey30")



plot(NA,NA, ylim=c(0,7),xlim=c(0,length(monall[,1])), xaxt="n", ylab="Passes/10min",xlab="Time (hr)",type="n",axes=F)
axis(1, at=c(seq(0,length(monall[,1]),6)+1),labels=0:12)
axis(2)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="69_1"],1,mean)),col="red",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="69_2"],1,mean)),col="red",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10))=="45"],1,mean)),col="red",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10))=="13"],1,mean)),col="red",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9))=="5"],1,mean)),col="red",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10))=="75"],1,mean)),col="orange",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10))=="11"],1,mean)),col="orange",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9))=="2"],1,mean)),col="orange",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="26_1"],1,mean)),col="orange",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="26_2"],1,mean)),col="orange",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9), substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="31_1"],1,mean)),col="grey",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9), substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="31_3"],1,mean)),col="grey",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9), substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="31_4"],1,mean)),col="grey",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9), substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="22_2"],1,mean)),col="grey",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9), substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="22_4"],1,mean)),col="grey",lwd=1)
legend("topleft", c("resistance class", "highly expressed", "lowly expressed"), col=c("red", "orange","grey40"), lty=1,bty="n")
points(c(1,30,42,60,73),c(18,18,28,33,35)*6.9/35,col="black", type="l",lty=2)
text(c(1,28,40,58,73),(c(18,18,28,33,35)+1)*6.9/35,paste(c(18,18,28,33,35),"C"),cex=1, col="black")


legend("topleft", c("11", "69", "22", "31"), col=c("#d95f0e","#fe9929", "#525252", "#969696"), lty=1,bty="n")
points(c(1,30,37,49,62),c(18,18,28,33,35)*6.9/35,col="red", type="l",lty=2)
text(c(1,30,37,49,62),(c(18,18,28,33,35)+1)*6.9/35,paste(c(18,18,28,33,35),"C"),cex=1, col="red")
#text(c(30,37,49,62),(c(18,18,28,33,35)+1)*6.9/35,paste("Stage",1:3),cex=1, col="black")






#########################################################
###### RAMPING: 18C(12h), 28C(2h), 33C(3h), 35C(2h) #######
### 3 classes; 5 lines per class #######

monitor11 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/first batch May 19th/Monitor11fem.txt", sep="\t", header = F)
monitor21 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/first batch May 19th/Monitor21fem.txt", sep="\t", header = F)
monitor31 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/first batch May 19th/Monitor31fem.txt", sep="\t", header = F)
monitor41 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/first batch May 19th/Monitor41fem.txt", sep="\t", header = F)
monitor51 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/first batch May 19th/Monitor51fem.txt", sep="\t", header = F)
monitor12 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/second batch May20th/Monitor1-4.txt", sep="\t", header = F)
monitor22 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/second batch May20th/Monitor2-4.txt", sep="\t", header = F)
monitor32 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/second batch May20th/Monitor3-4.txt", sep="\t", header = F)
monitor42 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/second batch May20th/Monitor4-4.txt", sep="\t", header = F)
monitor52 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/second batch May20th/Monitor5-4.txt", sep="\t", header = F)
monitor13 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/third batch May21st/Monitor1-5.txt", sep="\t", header = F)
monitor23 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/third batch May21st/Monitor2-5.txt", sep="\t", header = F)
monitor33 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/third batch May21st/Monitor3-5.txt", sep="\t", header = F)
monitor43 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/third batch May21st/Monitor4-5.txt", sep="\t", header = F)
monitor53 = read.table("/Users/sonjalecic/Desktop/Ace manuscript/ramping18_28_33_35_females/third batch May21st/Monitor5-5.txt", sep="\t", header = F)


mon11=monitor11[40 : 240,]
mon21=monitor21[40 : 240,]
mon31=monitor31[40 : 240,]
mon41=monitor41[40 : 240,]
mon51=monitor51[40 : 240,]
mon12=monitor12[40 : 240,]
mon22=monitor22[40 : 240,]
mon32=monitor32[40 : 240,]
mon42=monitor42[40 : 240,]
mon52=monitor52[40 : 240,]
mon13=monitor13[40 : 240,]
mon23=monitor23[40 : 240,]
mon33=monitor33[40 : 240,]
mon43=monitor43[40 : 240,]
mon53=monitor53[40 : 240,]


##### Day 1 #####
mon111=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon11[,1])/20)){
  a=as.data.frame(t(colSums(mon11[seq(20*i-4,20*i,1),11:42]))) #this part sums up passes for 10 minutes
  b=as.data.frame(mon11[20*i,c(2,3,4,5,12)]) #this part sums up passes for 10 minutes
  cc=cbind(b,a)
  colnames(cc)=colnames(mon111)
  mon111=rbind(mon111,cc)
}

mon211=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon21[,1])/20)){
  a=as.data.frame(t(colSums(mon21[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon21[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon211)
  mon211=rbind(mon211,cc)
}

mon311=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon31[,1])/20)){
  a=as.data.frame(t(colSums(mon31[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon31[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon311)
  mon311=rbind(mon311,cc)
}

mon411=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon41[,1])/20)){
  a=as.data.frame(t(colSums(mon41[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon41[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon411)
  mon411=rbind(mon411,cc)
}

mon511=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon51[,1])/20)){
  a=as.data.frame(t(colSums(mon51[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon51[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon511)
  mon511=rbind(mon511,cc)
}

##### Day 2 #####
mon112=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon12[,1])/20)){
  a=as.data.frame(t(colSums(mon12[seq(20*i-4,20*i,1),11:42]))) #this part sums up passes for 10 minutes
  b=as.data.frame(mon12[20*i,c(2,3,4,5,12)]) #this part sums up passes for 10 minutes
  cc=cbind(b,a)
  colnames(cc)=colnames(mon112)
  mon112=rbind(mon112,cc)
}

mon212=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon22[,1])/20)){
  a=as.data.frame(t(colSums(mon22[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon22[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon212)
  mon212=rbind(mon212,cc)
}

mon312=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon32[,1])/20)){
  a=as.data.frame(t(colSums(mon32[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon32[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon312)
  mon312=rbind(mon312,cc)
}

mon412=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon42[,1])/20)){
  a=as.data.frame(t(colSums(mon42[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon42[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon412)
  mon412=rbind(mon412,cc)
}

mon512=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon52[,1])/20)){
  a=as.data.frame(t(colSums(mon52[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon52[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon512)
  mon512=rbind(mon512,cc)
}


###### Day 3 #####
mon113=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon13[,1])/20)){
  a=as.data.frame(t(colSums(mon13[seq(20*i-4,20*i,1),11:42]))) #this part sums up passes for 10 minutes
  b=as.data.frame(mon13[20*i,c(2,3,4,5,12)]) #this part sums up passes for 10 minutes
  cc=cbind(b,a)
  colnames(cc)=colnames(mon113)
  mon113=rbind(mon113,cc)
}

mon213=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon23[,1])/20)){
  a=as.data.frame(t(colSums(mon23[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon23[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon213)
  mon213=rbind(mon213,cc)
}

mon313=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon33[,1])/20)){
  a=as.data.frame(t(colSums(mon33[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon33[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon313)
  mon313=rbind(mon313,cc)
}

mon413=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon43[,1])/20)){
  a=as.data.frame(t(colSums(mon43[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon43[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon413)
  mon413=rbind(mon413,cc)
}

mon513=data.frame(day=1, month="A",year=1,time="A", light=1,p01=0, p02=0, p03=0, p04=0, p05=0, p06=0, p07=0, p08=0, p09=0, p10=0, p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0, p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0, p31=0, p32=0)
for(i in 1:floor(length(mon53[,1])/20)){
  a=as.data.frame(t(colSums(mon53[seq(20*i-4,20*i,1),11:42])))
  b=as.data.frame(mon53[20*i,c(2,3,4,5,12)])
  cc=cbind(b,a)
  colnames(cc)=colnames(mon513)
  mon513=rbind(mon513,cc)
}

head(mon513)


positions=read.delim("/Users/sonjalecic/Desktop/Ace manuscript/locomotion_test2.txt",head=F,stringsAsFactors = F) # deswegen besser alle slots bef??llen #die excel tabelle kopieren und in ein txt file einf??gen
positions <- positions[, c(5, 4, 1, 2)]
colnames(positions) <- c("sampleID", "monitor", "row", "genotype")
positions$slot = c(rep(c(1:32), 5))
positions$sampleID = c(1:160)
positions <- positions[, c(1, 2, 5, 4, 3)]
positions = positions[, -5]


monall=cbind(mon111,mon211,mon311,mon411,mon511,
             mon112,mon212,mon312,mon412,mon512,
             mon113,mon213,mon313,mon413,mon513) 


monall=monall[,-c(grep("time", colnames(monall))[-1])] #removes extra time/light bla bla columns
monall=monall[,-c(grep("day", colnames(monall))[-1])]
monall=monall[,-c(grep("month", colnames(monall))[-1])]
monall=monall[,-c(grep("year", colnames(monall))[-1])]
monall=monall[,-c(grep("light", colnames(monall))[-1])]
colnames(monall)[-1:-5]=paste(positions$genotype,sep = "_")
head(monall)


#png("/Volumes/Temp1/shengkai/behavior_assay/spontaneous_locomotion.png",width = 24,height = 12,pointsize = 10,units = "cm",res=600)
#par(mfrow=c(1,2))
plot(NA,NA, ylim=c(0,7),xlim=c(0,length(monall[,1])), xaxt="n", ylab="Passes/10min",xlab="Time (hr)",type="n",axes=F, cex.lab=1.6)
#axis(1, 1:length(monall[,1]))
axis(1, at=c(seq(0,length(monall[,1]),6)+1),labels=0:1, cex.axis=1.6)
axis(2, cex.axis=1.6)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),6,6),substr(colnames(monall),7,7))=="Ib"],1,mean)),col="orange",lwd=3)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),6,6),substr(colnames(monall),7,7))=="II"],1,mean)),col="grey40",lwd=3)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),6,6),substr(colnames(monall),7,7))=="Ia"],1,mean)),col="red",lwd=3)
legend("topleft", c("Class Ia", "Class Ib", "Class II"), col=c("red", "orange","grey40"), lty=1,bty="n", lwd = 3, cex=1)
points(c(1,30,42,60,78),c(18,18,28,33,35)*6.9/35,col="black", type="l",lty=2, lwd = 2.5, cex=1)
text(c(1,28,40,58,78),(c(18,18,28,33,35)+1)*6.9/35,paste(c(18,18,28,33,35),"C"), cex=1.2, col="black")
text(c(15,33,49,65),(c(18,23,31,34)+1)*6.9/35,paste("Stage",1:4),cex=1, col="grey30")



plot(NA,NA, ylim=c(0,7),xlim=c(0,length(monall[,1])), xaxt="n", ylab="Passes/10min",xlab="Time (hr)",type="n",axes=F)
axis(1, at=c(seq(0,length(monall[,1]),6)+1),labels=0:12)
axis(2)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="69_1"],1,mean)),col="red",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="69_2"],1,mean)),col="red",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10))=="45"],1,mean)),col="red",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10))=="13"],1,mean)),col="red",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9))=="5"],1,mean)),col="red",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10))=="75"],1,mean)),col="orange",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10))=="11"],1,mean)),col="orange",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9))=="2"],1,mean)),col="orange",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="26_1"],1,mean)),col="orange",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9),substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="26_2"],1,mean)),col="orange",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9), substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="31_1"],1,mean)),col="grey",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9), substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="31_3"],1,mean)),col="grey",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9), substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="31_4"],1,mean)),col="grey",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9), substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="22_2"],1,mean)),col="grey",lwd=1)
lines(1:length(monall[,1]),(apply(monall[,paste0(substr(colnames(monall),9,9), substr(colnames(monall),10,10), substr(colnames(monall),11,11), substr(colnames(monall),12,12))=="22_4"],1,mean)),col="grey",lwd=1)
points(c(1,30,42,60,78),c(18,18,28,33,35)*6.9/35,col="black", type="l",lty=2, lwd = 2.5, cex=1)
text(c(1,28,40,58,78),(c(18,18,28,33,35)+1)*6.9/35,paste(c(18,18,28,33,35),"C"), cex=1.2, col="black")
text(c(15,33,49,65),(c(18,23,31,34)+1)*6.9/35,paste("Stage",1:4),cex=1, col="grey30")



legend("topleft", c("11", "69", "22", "31"), col=c("#d95f0e","#fe9929", "#525252", "#969696"), lty=1,bty="n")
points(c(1,30,37,49,62),c(18,18,28,33,35)*6.9/35,col="red", type="l",lty=2)
text(c(1,30,37,49,62),(c(18,18,28,33,35)+1)*6.9/35,paste(c(18,18,28,33,35),"C"),cex=1, col="red")
#text(c(30,37,49,62),(c(18,18,28,33,35)+1)*6.9/35,paste("Stage",1:3),cex=1, col="black")



