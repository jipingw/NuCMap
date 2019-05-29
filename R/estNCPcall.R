estNCPcall=function(seqname="default",genfile,watsonfile,crickfile)
{
 watson=read.table(watsonfile)
 if(sum(seqname=="default")==0){
 seqname=seqname
 }
 else if(sum(seqname=="default")==1){
 seqname=levels(watson[,1])
 }
 number=length(seqname)
 output=data.frame()

 for(i in 1:number){
  cN=as.character(seqname[i]); LcN=nchar(cN); LcN=as.integer(LcN); gF=as.character(genfile[i]); LgF=nchar(gF); LgF=as.integer(LgF); wF=as.character(watsonfile); LwF=nchar(wF);
  LwF=as.integer(LwF); cF=as.character(crickfile); LcF=nchar(cF); LcF=as.integer(LcF); ind=0

  tempdata=pattern[2:5,]
  
    results=.Fortran("EM4",LcN,cN,LgF,gF,LwF,wF,LcF,cF,pattern=as.matrix(tempdata),ind=as.integer(ind),index=as.integer(i),PACKAGE = "NuCMap")
    ind=results$ind
    FilePath=getwd()
    
 data=read.table(paste(FilePath,"/","NCPscore.ratio_4temp.txt",sep=""),header=T)
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
 NonRed=result$ID
 NonRed=NonRed[NonRed>0]
 cutoff=as.numeric(quantile(ratio[NonRed],0.1))
 NonRed=NonRed[ratio[NonRed]>cutoff]
 summary=cbind(rep(seqname[i],length(NonRed)),NonRed,NCP[NonRed],ratio[NonRed])
 output=rbind(output,summary)
} 
 write.table(output,paste("~/UNIQUEcenters.txt",sep=""),row.names=F,col.names=c("Chromosome","Position","NCPscore","Ratio"),quote=F)
 cat(paste("Output: '"), FilePath, "/", "NCPscore.ratio_4temp.txt'","\n",sep="")
 cat(paste("Output: '"), FilePath, "/", "UNIQUEcenters.txt'",sep="")
}
