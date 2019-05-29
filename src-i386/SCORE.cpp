extern "C" {
void getscore(int *cut1, int *cut2, int *size, double *p, double *scores)
{   	
  int b=*size;	
  for(int i=100;i<(b-7);i++){
  scores[i]=p[0]*(cut1[i-2]+cut2[i+2])+p[1]*(cut1[i-1]+cut2[i+1])+p[2]*(cut1[i]+cut2[i])+p[3]*(cut1[i+1]+cut2[i-1])
          +p[4]*(cut1[i+4]+cut2[i-4])+p[5]*(cut1[i+5]+cut2[i-5])+p[6]*(cut1[i+6]+cut2[i-6])+p[7]*(cut1[i+7]+cut2[i-7]);             
   }
}

}

















