subroutine readGenome(LgF,gF,z,w)

  implicit none
  integer  i,j,k,l,m,n,z,LgF; integer w(z)
  character*500 tpc; character*1 a(z); character(len=LgF) gF

  open(1,file=gF,action='read')
  read(1,*) tpc; read(1,*) tpc; k=len_trim(tpc)
  backspace(1)
  do i=1,z/k; read(1,'(500a1)') a((i-1)*k+1:i*k); end do
  if(z/k*k<z) read(1,'(500a1)') a(z/k*k+1:z)
  close(1); w=0
  do i=1,z
     if(a(i)=='A'.or.a(i)=='a') then
        w(i)=1
     else if(a(i)=='C'.or.a(i)=='c') then
        w(i)=2
     else if(a(i)=='G'.or.a(i)=='g') then
        w(i)=3
     else if(a(i)=='T'.or.a(i)=='t') then
        w(i)=4
     end if
  end do

end subroutine readGenome
