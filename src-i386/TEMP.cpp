extern "C" {
void Template(int *size, int *seq, int *temp)
{
  int b=*size;   
  for(int i=3;i<b-3;i++)
  {
    if((seq[i-3]==1) && (seq[i+3]==4))
    {temp[i]=1;}
    else if((seq[i-3]==1) && (seq[i+3]!=4))
    {temp[i]=2;}
    else if((seq[i-3]!=1) && (seq[i+3]==4))
    {temp[i]=3;}
    else if((seq[i-3]!=1) && (seq[i+3]!=4))
    {temp[i]=4;}        
  }       
  temp[0]=0;
  temp[1]=0;
  temp[2]=0;
  temp[b-3]=0;
  temp[b-2]=0;
  temp[b-1]=0;
}

}














