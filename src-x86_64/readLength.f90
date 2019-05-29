subroutine readLength(LgF,gF,z,ind)

  implicit none
  integer  i,j,k,l,m,n,z,ind,LgF
  character*500 tpc; character(len=LgF) gF

  open(1,file=gF,action='read'); ind=0
  read(1,*) tpc; if(tpc(1:1)/='>') then; ind=1; return; end if
  read(1,*,iostat=l) tpc; k=len_trim(tpc); n=1; m=k
  do; read(1,*,iostat=l) tpc; if(l/=0) exit; n=n+1
      m=len_trim(tpc)
  end do
  if(n==1.and.k>400) then; ind=2; return; end if
  close(1)
  z=k*(n-1)+m

end subroutine readLength
