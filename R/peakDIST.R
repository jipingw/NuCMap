peakDIST=function(seqname="default",watsonfile,crickfile){
watson=read.table(watsonfile)
crick=read.table(crickfile)
if(sum(seqname=="default")==0){
seqname=seqname
}
else if(sum(seqname=="default")==1){
seqname=levels(watson[,1])
}
number=length(seqname)

len=numeric(number)

for(i in 1:number)
{
 maxw=max(watson[watson[,1]==seqname[i],2])
 maxc=max(crick[crick[,1]==seqname[i],2])
 len[i]=max(maxw,maxc)
}

n=sum(len)
cutw=numeric(n)
cutc=numeric(n)

indw=watson[watson[,1]==seqname[1],2]
indc=crick[crick[,1]==seqname[1],2]
clvw=watson[watson[,1]==seqname[1],3]
clvc=crick[crick[,1]==seqname[1],3]

if (number>1){
for(i in 2:number)
{
 indw=c(indw,watson[watson[,1]==seqname[i],2]+sum(len[1:(i-1)]))
 indc=c(indc,crick[crick[,1]==seqname[i],2]+sum(len[1:(i-1)]))
 clvw=c(clvw,watson[watson[,1]==seqname[i],3])
 clvc=c(clvc,crick[crick[,1]==seqname[i],3])
} 
}

cutw[indw]=clvw
cutc[indc]=clvc


K1=cbind(cutw,c(1:n))
K3=cbind(cutc,c(1:n))
K2=K1[order(K1[, 1],decreasing=T), ]
K4=K3[order(K3[, 1],decreasing=T), ]
seqw=c(1:n)
result1=.C("getpeak",size=as.integer(n),ID=as.integer(seqw),position=as.integer(K2[,2]),package="NuCMap")

seqc=c(1:n)
result2=.C("getpeak",size=as.integer(n),ID=as.integer(seqc),position=as.integer(K4[,2]),package="NuCMap")

seqw=result1$ID
seqc=result2$ID

seqw=seqw[seqw>0]
seqc=seqc[seqc>0]

pos=numeric()
for(i in 1:length(seqc)){
dist=abs(seqc[i]-seqw)
pos[i]=which.min(dist)
}

mindist=seqc-seqw[pos]
freq=numeric(81)
for(i in 1:81){
freq[i]=sum(mindist==(i-41))
}

plot(c((-40):40),freq,type="h",col="red",lwd=3,main="Crick-Watson cleavage peak-peak distance",xlab="Crick-Watson",ylab="Frequency")

}






