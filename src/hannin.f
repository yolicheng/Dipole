      subroutine hannin(a,np, iha)
      implicit real*8(a-h,o-z)
      dimension a(np)
      integer, parameter :: dp = kind(1.d0)
      
      real(kind=dp), parameter :: pi2 = 8.d0 *  datan(1.d0)
      
!      print*,'np, iha', np, iha
!      print*; read*
      
!      common /contro/pi2,tol,lecpas,irepas,iha,iat,tofre,nrec1,nrec2
!     1,tolra,ftolra,nminpa,ilat
     
      if(iha.eq.0)return
      coe=1
      if(iha.eq.2)coe=2.d0/3
      if(iha.eq.3)coe=.4d0
      del=pi2/np
      c0=1
      c1=dcos(del)
      z=2*c1
      a(1)=0
      a(2)=a(2)*(1-c1)**iha*coe
      do 1 i=3,np
         c2=z*c1-c0
         c0=c1
         c1=c2
         a(i)=a(i)*(1-c1)**iha*coe
1     continue
      return
      end
