plotCUTS=function(seqname,watsonfile,crickfile,startpos,endpos)
{
  wF=as.character(watsonfile)
  cF=as.character(crickfile)
  startpos=as.integer(startpos)
  endpos=as.integer(endpos)
  watson=read.table(wF)
  crick=read.table(cF)
  watson=watson[watson[,1]==seqname,]
  crick=crick[crick[,1]==seqname,]
  len=max(max(watson[,2]),max(crick[,2])) 
  if(endpos>len){
  print("Stop! Error: no cleavage in this region or the ending position is out of range!")
  }
  else{
  watson_in=watson[watson[,2]>=startpos,]
  watson_in=watson_in[watson_in[,2]<=endpos,c(2,3)]
  crick_in=crick[crick[,2]>=startpos,]
  crick_in=crick_in[crick_in[,2]<=endpos,c(2,3)]
  plot(NA,NA,xlim=c(startpos,endpos),ylim=c(-max(crick_in[,2]),max(watson_in[,2])),main="Watson and Crick cleavage frequency",xlab="position",ylab="cuts for watson(+) and crick(-)")
  lines(watson_in[,1],watson_in[,2],type="h",col="blue")
  lines(crick_in[,1],-crick_in[,2],type="h",col="red") 
  } 
}
