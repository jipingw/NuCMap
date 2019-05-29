extern "C"{
void AATT(int *seq, int *len1, int *len2, int *pos, double *aatt)
{
 double dinu[146]={};

 for(int j=146;j>=1;j--)
 {
    for(int i=0;i<*len2;i++)
     {
           if(seq[pos[i]+72-j]==1)
            {
             if((seq[pos[i]+73-j]==1) || (seq[pos[i]+73-j]==4))                         
             {dinu[146-j]=dinu[146-j]+1;}
            }
            else if(seq[pos[i]+72-j]==4)
            {
             if((seq[pos[i]+73-j]==1) || (seq[pos[i]+73-j]==4))             
             {dinu[146-j]=dinu[146-j]+1;}
            }
	  }
      }


for(int i=0;i<146;i++)
{
   aatt[i]=dinu[i];
}

}
}
