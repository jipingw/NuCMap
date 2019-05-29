callRED=function(estresults,seqname="default",threshold="default")
{
 estres=as.character(estresults)
 data=read.table(estres,header=T)
 if(sum(seqname=="default")==0){
 seqname=seqname
 }
 else if(sum(seqname=="default")==1){
 seqname=levels(data[,1])
 }
 number=length(seqname)
 output=data.frame()
 for(i in 1:number){
 NCPr=data[data[,1]==seqname[i],]
 n=NCPr[nrow(NCPr),2]
 seq=c(1:n)
 ratio=numeric(n)
 NCP=numeric(n)
 ratio[NCPr[,2]]=NCPr[,4]
 NCP[NCPr[,2]]=NCPr[,3]
 K1=cbind(ratio,NCP,c(1:n))
 K2=K1[order(K1[, 1],decreasing=T), ]
 result=.C("NonRedun",size=as.integer(n),ID=as.integer(seq),position=as.integer(K2[,3]),package="NuCMap")
 Red=result$ID
 Red=Red[Red>0]
 cutoff=as.numeric(quantile(ratio[Red],0.1))
 if(threshold=="default"){
 cutoff=cutoff}
 else{
 cutoff=threshold
 }
 Redund=K1[K1[,1]>cutoff,3]
 summary=cbind(rep(seqname[i],length(Redund)),Redund,NCP[Redund],ratio[Redund])
 output=rbind(output,summary)
}
 FilePath=getwd()
 write.table(output, paste("~/REDcenters.txt",sep=""),row.names=F,col.names=c("Chromosome","Position","NCPscore","Ratio"),quote=F)
 cat(paste("Output: '"), FilePath, "/", "REDcenters.txt'",sep="")
}
