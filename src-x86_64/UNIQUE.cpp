extern "C" {
void NonRedun(int *size, int *ID, int *position)
{
  for(int i=0;i<*size;i++)
  {
    if(ID[position[i]-1]!=0)
    { 
      for(int j=1;j<=107;j++)
      {
	if(position[i]>j)
	{ID[position[i]-j-1]=0;}
        if(position[i]+j<=*size)
	{ID[position[i]+j-1]=0;}	
      }	      
    }
  }	  
}

}
















