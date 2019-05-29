
subroutine EM1(LcN,cN,LgF,gF,LwF,wF,LcF,cF,pattern,ind,index)

  implicit none
  integer  i,j,k,l,m,n,z,rd,wd,ind,tep,step,index
  integer  LcN,LgF,LwF,LcF,tp(0:500),id(0:500)
  real*8   pattern(8),lam(-2:7),temp,olik,lik,sw,sc,s1,s2,tmp,thd1,thd2,av,sd
  character*500 tpc; character*60 tepc
  character(len=LcN) cN; character(len=LgF) gF; character(len=LwF) wF; character(len=LcF) cF

  integer,allocatable:: cutW(:),cutC(:)
  real*8,allocatable:: O(:),Q(:),r(:),lamW(:),lamC(:),tw(:),tc(:)
  ind=0

  open(1,file=gF,action='read')
  read(1,*) tpc; if(tpc(1:1)/='>') then; ind=1; return; end if
  read(1,*,iostat=l) tpc; k=len_trim(tpc); n=1; m=k
  do; read(1,*,iostat=l) tpc; if(l/=0) exit; n=n+1
      m=len_trim(tpc)
  end do
  if(n==1.and.k>400) then; ind=2; return; end if
  close(1); z=k*(n-1)+m
  allocate(cutW(z),cutC(z),r(z),lamW(z),lamC(z))

  !---------- input cutW and cutC

  open(1,file=wF,action='read'); cutW=0
  do; read(1,*,iostat=l) tepc,i,k; if(l/=0) exit
      if(trim(tepc)==cN) cutW(i)=k
  end do
  close(1)
  open(2,file=cF,action='read'); cutC=0
  do; read(2,*,iostat=l) tepc,i,k; if(l/=0) exit
      if(trim(tepc)==cN) cutC(i)=k
  end do
  close(2)

  !---------- adjust pattern

  pattern=pattern/sum(pattern)
  lam=0.0; lam(-2:1)=pattern(1:4); lam(4:7)=pattern(5:8)
  temp=(sum(cutW)+sum(cutC))/2.0d0/z; temp=temp/lam(6)
  lam=temp*lam

  rd=73; wd=2*rd+1
  !---------- calculate cut ratio (Crick/Watson)

  m=sum(cutW(1:wd)); n=sum(cutC(1:wd))
  if(m==0.and.n==0) then
     r(1:rd+1)=1.0
  else if(m==0) then
     r(1:rd+1)=min(n*1.0d0,1000.0)
  else if(n==0) then
     r(1:rd+1)=max(0.001,1.0d0/m)
  else
     r(1:rd+1)=min(max(0.001,n*1.0d0/m),1000.0)
  end if
  do i=rd+2,z-rd
     m=m-cutW(i-1-rd)+cutW(i+rd); n=n-cutC(i-1-rd)+cutC(i+rd)
     if(m==0.and.n==0) then
        r(i)=1.0
     else if(m==0) then
        r(i)=min(n*1.0d0,1000.0)
     else if(n==0) then
        r(i)=max(0.001,1.0d0/m)
     else
        r(i)=min(max(0.001,n*1.0d0/m),1000.0)
     end if
  end do
  r(z-rd+1:z)=r(z-rd)
  av=0.0; do i=1,z; av=av+log10(r(i)); end do; av=av/z
  sd=0.0; do i=1,z; sd=sd+(log10(r(i))-av)**2; end do; sd=sqrt(sd/(z-1.0))
  thd1=10.0**(av-2.0*sd); thd2=10.0**(av+2.0*sd)

  !---------- calculate noise on Watson

  tp(1:wd)=cutW(1:wd); do i=1,wd; id(i)=i; end do
  do i=1,wd-1; do j=wd,i+1,-1
     if(tp(j-1)>tp(j)) then
        tep=tp(j-1); tp(j-1)=tp(j); tp(j)=tep
        tep=id(j-1); id(j-1)=id(j); id(j)=tep
     end if
  end do; end do
  lamW(1:rd+1)=sum(tp(1:rd))*1.0d0/rd

  do i=rd+2,z-rd
     do j=1,wd; if(cutW(i+rd)<=tp(j)) exit; end do
     do k=0,j-2; tp(k)=tp(k+1); id(k)=id(k+1); end do; tp(j-1)=cutW(i+rd); id(j-1)=i+rd
     do k=0,wd; if(id(k)==i-1-rd) exit; end do
     do j=k,1,-1; tp(j)=tp(j-1); id(j)=id(j-1); end do
     lamW(i)=sum(tp(1:rd))*1.0d0/rd
  end do

  lamW(z-rd+1:z)=lamW(z-rd)

  !---------- calculate noise on Crick

  tp(1:wd)=cutC(1:wd); do i=1,wd; id(i)=i; end do
  do i=1,wd-1; do j=wd,i+1,-1
     if(tp(j-1)>tp(j)) then
        tep=tp(j-1); tp(j-1)=tp(j); tp(j)=tep
        tep=id(j-1); id(j-1)=id(j); id(j)=tep
     end if
  end do; end do
  lamC(1:rd+1)=sum(tp(1:rd))*1.0d0/rd

  do i=rd+2,z-rd
     do j=1,wd; if(cutC(i+rd)<=tp(j)) exit; end do
     do k=0,j-2; tp(k)=tp(k+1); id(k)=id(k+1); end do; tp(j-1)=cutC(i+rd); id(j-1)=i+rd
     do k=0,wd; if(id(k)==i-1-rd) exit; end do
     do j=k,1,-1; tp(j)=tp(j-1); id(j)=id(j-1); end do
     lamC(i)=sum(tp(1:rd))*1.0d0/rd
  end do

  lamC(z-rd+1:z)=lamC(z-rd)

  allocate(tw(z),tc(z),O(z),Q(z))

  !---------- EM start

  id(1:8)=(/-2,-1,0,1,4,5,6,7/)

  lik=0.0; step=0; O(1:z)=0.05 !<< initial K
  !!!!
  do !
  !!!!
  olik=lik

  lik=0.0
  do i=1,z
     tw(i)=lamW(i); tc(i)=lamC(i)
     do k=1,8; j=id(k)
        if(1<=i-j.and.i-j<=z) tw(i)=tw(i)+O(i-j)*lam(j)
        if(1<=i+j.and.i+j<=z) tc(i)=tc(i)+O(i+j)*lam(j)*r(i+j)
     end do
     lik=lik+cutW(i)*log(tw(i))-tw(i)+cutC(i)*log(tc(i))-tc(i)
  end do

  if(step>=500.and.abs(lik-olik)<abs(olik)*1.0d-9) then !<<<<<
     exit
  else if(step>=5000) then
     if(abs(lik-olik)<abs(olik)*1.0d-8) then !<<<<<
        exit
     else
        ind=3; return
     end if
  end if

  Q=O
  do i=1,z
     sw=0.0; sc=0.0; s1=0.0; s2=0.0
     do k=1,8; j=id(k)
        if(1<=i+j.and.i+j<=z) then
           sw=sw+cutW(i+j)*Q(i)*lam(j)/tw(i+j); s1=s1+lam(j)
        end if
        if(1<=i-j.and.i-j<=z) then
           sc=sc+cutC(i-j)*Q(i)*r(i)*lam(j)/tc(i-j); s2=s2+r(i)*lam(j)
        end if
     end do
     O(i)=(sw+sc)/(s1+s2); O(i)=max(1.0d-6,O(i)) !#
  end do

  step=step+1
  !!!!!!!!
  end do !
  !!!!!!!!

  if(index==1) then; tepc='Chromosome'
     open(1,file='NCPscore.ratio_1temp.txt',status='replace')
     write(1,'(a,3x,a8,2x,a8,6x,a5,5x,a9)') tepc(1:LcN),'Position','NCPscore','Ratio','cNCPscore'
  else
     open(1,file='NCPscore.ratio_1temp.txt',position='append')
  end if
  do i=1,z
     O(i)=O(i)*(1.0+r(i))/2.0; temp=O(i)/max(0.5,(lamW(i)+lamC(i))/2.0)
     tmp=O(i)
     if(r(i)>thd2) then
        tmp=tmp*2.0/(1.0+r(i))*r(i)
     else if(r(i)<thd1) then
        tmp=tmp*2.0/(1.0+r(i))
     end if
     if(temp>=1.0d-3) write(1,'(a,1x,i10,1x,f9.3,1x,f10.5,1x,f9.2)') cN,i,O(i),temp,tmp
  end do
  close(1)

  deallocate(cutW,cutC,O,Q,r,lamW,lamC,tw,tc)
end subroutine EM1
