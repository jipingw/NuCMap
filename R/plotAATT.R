plotAATT=function(seqname="default",genfile,center)
{
centerfile=read.table(center,header=T)
if(sum(seqname=="default")==0){
seqname=seqname
}
else if(sum(seqname=="default")==1){
seqname=levels(centerfile[,1])
}

number=length(seqname)
len=numeric(number)
center1=numeric()
combine=numeric()

for(i in 1:number)
{
 LgF=nchar(genfile[i]); gFd=as.character(genfile[i]); zstart=0; ind=0
 readLength=.Fortran("readLength", LgF, gFd, z=as.integer(zstart), ind=as.integer(ind))
 len[i]=readLength$z
 ind=readLength$ind
  if(ind==1){
    print("Stop! Error: the genome file should be a fasta file, i.e. the first row should be '>sequence name'")
  } 
  if(ind==2){
    print("Stop! Error: the length of each line in fasta file should be less than 400bp!")
 }
 wstart=numeric(as.integer(len[i]))
 readGenome=.Fortran("readGenome", LgF, gFd, z=as.integer(len[i]), w=as.integer(wstart))
 w=readGenome$w
 combine=c(combine,as.integer(w))
 pos=centerfile[centerfile[,1]==seqname[i],]
 if(i==1){
 center1=c(center1, pos[,2])
 }
 else{
 center1=c(center1,pos[,2]+sum(len[1:(i-1)]))
 }
}
 n=sum(len)
 center1=center1[center1>73]
 center1=center1[center1<n-73]
 freq=numeric(146)

 result=.C("AATT", seq=as.integer(combine), len1=as.integer(n), len2=as.integer(length(center1)), pos=as.integer(center1), aatt=as.double(freq))
 freq=result$aatt
 wcp=freq[1:146]+freq[146:1]
 plot(c(-73:72),wcp/2/length(center1),type="l",col="blue",xlab="distance to nucleosome center", ylab="frequency",main="AA/TT/AT/TA frequency")
}
