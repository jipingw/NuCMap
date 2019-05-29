calOccup=function(estresults,genfile,rednufile,seqname="default")
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
 len=numeric(number)
 center=read.table(rednufile,header=T)
 output=data.frame()

 for(i in 1:number){
 NCPr=data[data[,1]==seqname[i],]
 genf=as.character(genfile[i])
 LgF=nchar(genf); gFd=as.character(genf); zstart=0; ind=0
 readLength=.Fortran("readLength", LgF, gFd, z=as.integer(zstart), ind=as.integer(ind))
 len[i]=readLength$z
 ind=readLength$ind
  if(ind==1){
    print("Stop! Error: the genome file should be a fasta file, i.e. the first row should be '>sequence name'")
  } 
  if(ind==2){
    print("Stop! Error: the length of each line in fasta file should be less than 400bp!")
 }
 leng=as.integer(len[i])
 cNCP=numeric(leng)
 crNCP=numeric(leng)
 pos=center[center[,1]==seqname[i],]
 cNCP[NCPr[,2]]=NCPr[,5]
 crNCP[pos[,2]]=cNCP[pos[,2]]
 crNCP=c(rep(0,73),crNCP,rep(0,73)) 

 occup=function(m){
 return(sum(crNCP[(m-73):(m+73)]))
 }
 occupancy=sapply(c(74:(leng+73)),occup)
 summary=cbind(rep(seqname[i],leng),c(1:leng),occupancy)
 output=rbind(output,summary)
}

 FilePath=getwd()
 write.table(output, paste("~/NuOccupancy.txt",sep=""),row.names=F,col.names=c("Chromosome","Position","Occupancy"),quote=F)
 cat(paste("Output: '"), FilePath, "/", "NuOccupancy.txt'",sep="")
}
