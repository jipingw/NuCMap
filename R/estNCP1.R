estNCP1=function(seqname="default",genfile,watsonfile,crickfile,temp1="default")
{
watson=read.table(watsonfile)  
if(sum(seqname=="default")==0){
seqname=seqname
}
else if(sum(seqname=="default")==1){
seqname=levels(watson[,1])
}
  number=length(seqname)
  for(i in 1:number){
  cN=as.character(seqname[i]); LcN=nchar(cN); LcN=as.integer(LcN); gF=as.character(genfile[i]); LgF=nchar(gF); LgF=as.integer(LgF); wF=as.character(watsonfile); LwF=nchar(wF);
  LwF=as.integer(LwF); cF=as.character(crickfile); LcF=nchar(cF); LcF=as.integer(LcF);ind=0
  if(temp1=="default")
  {
    tempdata=pattern[1,]
  }
  else{
  temp1=as.character(temp1)
  tempdata=read.table(temp1)
}
  results=.Fortran("EM1",LcN,cN,LgF,gF,LwF,wF,LcF,cF,pattern=as.matrix(tempdata),ind=as.integer(ind),index=as.integer(i),PACKAGE = "NuCMap")
  ind=results$ind    

  if(ind==1){
    print("Stop! Error: the genome file should be a fasta file, i.e. the first row should be '>sequence name'")
  } else if(ind==2){
    print("Stop! Error: the length of each line in fasta file should be less than 400bp!")
  } else if(ind==3){
    print("The EM algorithm does not converge within 5000 repetitions!")
  }

}
if(ind==0){
 FilePath=getwd()
 cat(paste("Output: '"), FilePath, "/", "NCPscore.ratio_1temp.txt'",sep="") 
}
}
