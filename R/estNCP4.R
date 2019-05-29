estNCP4=function(seqname="default",genfile,watsonfile,crickfile,temp4="default")
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
  LwF=as.integer(LwF); cF=as.character(crickfile); LcF=nchar(cF); LcF=as.integer(LcF); ind=0
  if(temp4=="default")
  {
    tempdata=pattern[2:5,]
  }
  else{
  temp4=as.character(temp4)
  tempdata=read.table(temp4)
}
    results=.Fortran("EM4",LcN,cN,LgF,gF,LwF,wF,LcF,cF,pattern=as.matrix(tempdata),ind=as.integer(ind),index=as.integer(i),PACKAGE = "NuCMap")
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
      cat(paste("Output: '"), FilePath, "/", "NCPscore.ratio_4temp.txt'",sep="")
    }
}
