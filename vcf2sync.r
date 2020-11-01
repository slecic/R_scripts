#!/usr/bin/env Rscript

require(sys)
require(data.table)

vcffile="08CallSNPsNovo/markers.vcf.gz"
syncfile="output.sync"

cmd="bcftools"
args=c("query","-f","%CHROM\t%POS\t%REF\t%ALT[\t%RO:%AO]\n",vcffile)
exec_wait(cmd,args,std_out="foo.txt",std_err="bcftools.err")

dt=fread("foo.txt")
names(dt)[1:4]=c("CHROM","POS","REF","ALT")
dt=melt(dt,id.vars=names(dt)[1:4],variable.name="POP",value.name="COUNTS")
dt[,c("COUNT1","COUNT2"):=tstrsplit(COUNTS,':')]
dt[,COUNTS:=NULL]
dt[,c("COUNT1","COUNT2"):=.(as.integer(COUNT1),as.integer(COUNT2))]
bases=c('A','C','G','T','N')
dt=na.omit(dt[REF%in%bases & ALT%in%bases])
mkcstring=function(r,a,c1,c2)
{
  cvec=c(A=0L,C=0L,G=0L,T=0L,N=0L,DEL=0L)
  cvec[c(r,a)]=c(c1,c2)
  return(paste(cvec,collapse=':'))
}
dt[,cstring:=mkcstring(REF,ALT,COUNT1,COUNT2),
   by=.(CHROM,POS,POP)]
fwrite(na.omit(dcast(dt,CHROM+POS+REF~POP,value.var="cstring")),
       file=syncfile,
       sep='\t',
       col.names=FALSE)
