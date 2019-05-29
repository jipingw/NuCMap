trainTEMP1=function(seqname="default",watsonfile,crickfile){
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

n=sum(len)+7
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


Pold=c(0.125,0.25,0.125,0,0,0.125,0.25,0.125)
Pnew=numeric(8)
score=numeric(n)
ss=1;

while(ss>10^(-2)){
Pnew=Pold
result1=.C("getscore",cut1=as.integer(cutw),cut2=as.integer(cutc),size=as.integer(n),p=as.double(Pnew),scores=as.double(score),package="NuCMap")
score=result1$scores
K1=cbind(score,c(1:n))
K2=K1[order(K1[, 1],decreasing=T), ]
seq=c(1:n)
result2=.C("getpeak",size=as.integer(n),ID=as.integer(seq),position=as.integer(K2[,2]),package="NuCMap")
seq=result2$ID
seq=seq[seq>2]
Pold=c(sum(cutw[seq-2],cutc[seq+2]),sum(cutw[seq-1],cutc[seq+1]),sum(cutw[seq],cutc[seq]),sum(cutw[seq+1],cutc[seq-1]),sum(cutw[seq+4],cutc[seq-4]),sum(cutw[seq+5],cutc[seq-5]),sum(cutw[seq+6],cutc[seq-6]),sum(cutw[seq+7],cutc[seq-7]))/length(seq)/2
ss=sum((Pold-Pnew)^2)
}
 FilePath=getwd()
 write.table(signif(t(Pold),digits=5),paste("~/trainTEMP1result",".txt",sep=""),row.names="Freq",col.names=c("Pos(-2)","Pos(-1)","Pos(0)","Pos(+1)","Pos(+4)","Pos(+5)","Pos(+6)","Pos(+7)"))
 cat(paste("Output: '"), FilePath, "/", "trainTEMP1result",".txt'",sep="")
}
