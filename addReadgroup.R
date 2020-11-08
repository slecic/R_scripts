require(data.table)
require(sys)
require(glue)
require(digest)

### script for adding read groups ###
##Author:Sonja_Lecic

setwd ("/home/vetlinux02/Sonja/projec/Ace")

# build a data table of file names, descriptions, replicate/generation/condition, barcode (if any)
vet = fread("Info/Popgen_Illumina_Data - 1-Vetmed.csv")
bgi = fread("Info/Popgen_Illumina_Data - 3-BGI+CSF+Fasteris.csv")
bar = fread("Info/Popgen_Illumina_Data - Barcoded samples.csv")

# start with original filenames
fn = fread("10RepGen/original-filenames.txt",header=FALSE)$V1

### 
t1 = vet[sub("(.bam|.txt)","",`File name`)%in%fn,.(`Unique brief description`,`File name`,`ID in FASTQ header`,`Flowcell`,
                                                   `Total fragment size (incl. 120 bp adapters/primers)`)]
t2 = bgi[sub("(.bam|.txt)","",`File name`)%in%fn,.(`Unique brief description`,`File name`,`ID in FASTQ header`,
                                                   `Total fragment size (incl. 120 bp adapters/ primers)`,
                                                   `Run folder name or Submission ID`)]
t3 = bar[`Submission ID`%in%t2$`Run folder name or Submission ID`,.(`Brief sample description`,`Barcode sequence (1)`,`Submission ID`)]

# Pool 336 also contains samples we do not need
t3 = t3[!grepl("(D. sim|OregonR|Kahlenberg)",`Brief sample description`)]

# now need to fuse these temporary data tables
# filename, fastq ID, and read length t2 -> t3
t3$`File name` = t3$`ID in FASTQ header` = t3$`Total fragment size (incl. 120 bp adapters/primers)` = as.character(NA)
for (id in t2$`Run folder name or Submission ID`)
{
  t3[`Submission ID`==id,"File name"] = t2[`Run folder name or Submission ID`==id,"File name"]
  t3[`Submission ID`==id,"ID in FASTQ header"] = t2[`Run folder name or Submission ID`==id,"ID in FASTQ header"]
  t3[`Submission ID`==id,"Total fragment size (incl. 120 bp adapters/primers)"] = 
    t2[`Run folder name or Submission ID`==id,"Total fragment size (incl. 120 bp adapters/ primers)"]
}

t4 = merge(t3,t1,
           by.x=c("Brief sample description","File name","ID in FASTQ header","Total fragment size (incl. 120 bp adapters/primers)"),
           by.y=c("Unique brief description","File name","ID in FASTQ header","Total fragment size (incl. 120 bp adapters/primers)"),
           all=TRUE)
setnames(t4,"Barcode sequence (1)","Barcode")

# now extract experiment info from description
x= ".*Portugal (base|HOT|Hot|COLD) .*"
t4$`Regime` = factor(toupper(sub(x,"\\1",t4$`Brief sample description`)),levels=c("BASE","HOT","COLD"))
x = ".*F([0-9]+) ?.*"
t4$`Generation` = as.integer(ifelse(grepl(x,t4$`Brief sample description`),sub(x,"\\1",t4$`Brief sample description`),"0"))
x = ".* r\\. ?([0-9]*).*"
t4$`Replicate` = as.integer(sub(x,"\\1",t4$`Brief sample description`))

# insert size
t4$`Insert size` = as.integer(t4$`Total fragment size (incl. 120 bp adapters/primers)`)-120

# drop unused columns
t4[,`Total fragment size (incl. 120 bp adapters/primers)`:=NULL]
t4[,`Submission ID`:=NULL]
# reorder again
setcolorder(t4,c("Regime","Replicate","Generation","Barcode","Insert size","Brief sample description","ID in FASTQ header","Flowcell","File name"))
setorderv(t4,cols=c("Regime","Replicate","Generation"))
t4$`File name` = sub("(.bam|.txt)","",t4$`File name`)
t4$`File name` = sub("_1_sequence","_sequence",t4$`File name`)
# fix flowcell entries
t4[`Flowcell`=="n.a.","Flowcell"] = NA
t4[is.na(`Flowcell`),"Flowcell"] = ifelse(t4$`ID in FASTQ header`!="",
                                          sub("ILLUMINA-","",lapply(strsplit(t4$`ID in FASTQ header`,":"),"[",1)),
                                          lapply(strsplit(t4$`File name`,"_"),"[",3))[is.na(t4$Flowcell)]
# find lane
t4$`Lane` = unlist(ifelse(t4$`ID in FASTQ header`!="",
                   lapply(strsplit(t4$`ID in FASTQ header`,":"),"[",2),
                   sub("L","",lapply(strsplit(t4$`File name`,"_"),"[",4))))

# fix _Loaner ext.
t4$`File name` = sub("_Loaner","",t4$`File name`)

# construct readgroup entries
t4$RG_PU = paste(t4$Flowcell,t4$Lane,t4$Barcode,sep=".")
t4$RG_PL = "ILLUMINA"
t4$RG_PI = t4$`Insert size`
t4$RG_SM = paste("Portugal_",t4$Regime,"_r",t4$Replicate,".F",t4$Generation,sep="")
t4$RG_LB = sub("(.*)[a-d]$","\\1",lapply(lapply(strsplit(sub("s_._","",t4$`File name`),"_"),"[",1:2),collapse,sep="-"))
t4$RG_DS = t4$`Brief sample description`
t4$RG_ID = apply(t4[,.(RG_SM,RG_LB,RG_PL,RG_PU)],1,digest,algo="xxhash32")

# math cold and new file names (generated on cmd line from sam @PG headers)
id = fread("10RepGen/file-match.txt")
id$Replicate = as.integer(sub(".*_r([0-9]+)_.*","\\1",id$`New name`))
id$Generation = as.integer(sub(".*_F([0-9]+)_.*","\\1",id$`New name`))
id$`Sample ID` = as.integer(sub(".*__([0-9]+).bam","\\1",id$`New name`))

t4 = merge(t4,id,by=c("Replicate","Generation","File name"),all=T)

running = list()
for (f in dir("12MapNovo/",pattern=".bam"))
{
  print(paste("Processing file",f))
  t_ = t4[`New name`==f,.(RG_ID,RG_PU,RG_PL,RG_PI,RG_SM,RG_LB,RG_DS)]
  cmd = "picard"
  args = with(t_,c("AddOrReplaceReadGroups",
                   paste("I=12MapNovo",f,sep="/"),
                   paste("O=98Temp",f,sep="/"),
                   paste("ID",RG_ID,sep="="),
                   paste("PU",RG_PU,sep="="),
                   paste("PL",RG_PL,sep="="),
                   paste("PI",RG_PI,sep="="),
                   paste("SM",RG_SM,sep="="),
                   paste("LB",RG_LB,sep="="),
                   paste("DS",RG_DS,sep="="),
                   paste("PG","picard-tools",sep="=")))
  std_out = "12MapNovo/add-readgroup.out"
  std_err = "12MapNovo/add-readgroup.err"
  #std_err = paste("12MapNovo/",sub(".bam","",f),".err",sep="")
  # Wait until a slot for running a picard instance becomes available
  repeat
  {
    running = Filter(function(pid) is.na(exec_status(pid,wait=FALSE)),running)
    if (length(running)<12)
      break
    Sys.sleep(10)
  }
  print(paste("Executing",cmd,collapse(args,sep=" "),"in the background ..."))
  running[[length(running)+1]] = exec_background(cmd,args,std_out=std_out,std_err=std_err)
}
