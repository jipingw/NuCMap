trainTEMP4=function(seqname="default",genfile,watsonfile,crickfile,center){
cutw=read.table(watsonfile)
cutc=read.table(crickfile)
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
template=numeric()
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
 tem=numeric(as.integer(len[i]))
 result=.C("Template", size=as.integer(len[i]), seq=as.integer(w), temp=as.integer(tem))
 temp=result$temp
 template=c(template,as.integer(temp))
 pos=centerfile[centerfile[,1]==seqname[i],]
 if(i==1){
 center1=c(center1, pos[,2])
 }
 else{
 center1=c(center1,pos[,2]+sum(len[1:(i-1)]))
 }
}

n=sum(len)
watson=numeric(n)
crick=numeric(n)

indw=cutw[cutw[,1]==seqname[1],2]
indc=cutc[cutc[,1]==seqname[1],2]
clvw=cutw[cutw[,1]==seqname[1],3]
clvc=cutc[cutc[,1]==seqname[1],3]

if (number>1){
for(i in 2:number)
{
 indw=c(indw,cutw[cutw[,1]==seqname[i],2]+sum(len[1:(i-1)]))
 indc=c(indc,cutc[cutc[,1]==seqname[i],2]+sum(len[1:(i-1)]))
 clvw=c(clvw,cutw[cutw[,1]==seqname[i],3])
 clvc=c(clvc,cutc[cutc[,1]==seqname[i],3])
} 
}

watson[indw]=clvw
crick[indc]=clvc

te1=center1[template[center1]==1]
te2=center1[template[center1]==2]
te3=center1[template[center1]==3]
te4=center1[template[center1]==4]
te1=te1[te1>20]
te1=te1[te1<n-21]
te2=te2[te2>20]
te2=te2[te2<n-21]
te3=te3[te3>20]
te3=te3[te3<n-21]
te4=te4[te4>20]
te4=te4[te4<n-21]

freq1=numeric(41)
freq2=numeric(41)
freq3=numeric(41)
freq4=numeric(41)
for(i in -20:20){
freq1[i+21]=mean(c(watson[te1+i],crick[te1-i]))
}

for(i in -20:20){
freq2[i+21]=mean(c(watson[te2+i],crick[te3-i]))
}

for(i in -20:20){
freq3[i+21]=mean(c(watson[te3+i],crick[te2-i]))
}

for(i in -20:20){
freq4[i+21]=mean(c(watson[te4+i],crick[te4-i]))
}

par(mfrow=c(2,2))
plot(c(-20:20),freq1,type="l",col="red",xlab="distance from center",ylab="average cut frequency",main="A(-3)T(+3)")
points(c(-20:20),freq1,pch=19,col="red")
plot(c(-20:20),freq2,type="l",col="red",xlab="distance from center",ylab="average cut frequency",main="A(-3)-T(+3)")
points(c(-20:20),freq2,pch=19,col="red")
plot(c(-20:20),freq3,type="l",col="red",xlab="distance from center",ylab="average cut frequency",main="-A(-3)T(+3)")
points(c(-20:20),freq3,pch=19,col="red")
plot(c(-20:20),freq4,type="l",col="red",xlab="distance from center",ylab="average cut frequency",main="-A(-3)-T(+3)")
points(c(-20:20),freq4,pch=19,col="red")

cluster=c(-2,-1,0,1,4,5,6,7)

lambda=matrix(0,nrow=4,ncol=8)

for(i in 1:8){
lambda[1,i]=freq1[21+cluster[i]]
lambda[2,i]=freq2[21+cluster[i]]
lambda[3,i]=freq3[21+cluster[i]]
lambda[4,i]=freq4[21+cluster[i]]
}

 FilePath=getwd()
 write.table(signif(lambda,digits=5),paste("~/trainTEMP4result",".txt",sep=""),row.names=c("+A+T","+A-T","-A+T","-A-T"),col.names=c("Pos(-2)","Pos(-1)","Pos(0)","Pos(+1)","Pos(+4)","Pos(+5)","Pos(+6)","Pos(+7)"))
 cat(paste("Output: '"), FilePath, "/", "trainTEMP4result",".txt'",sep="")
}
