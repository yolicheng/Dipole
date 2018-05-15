! Modified to free format,  2016-09-23 11:53:49 

!  This routine is an implementation of a runge-kutta-fehlberg
!  method of orders 7 and 8. using a total of 13 steps (and
!  evaluations of the vectorfield) it computes two different
!  estimations of the next point. the difference between both
!  estimations (with local errors of order 8 and 9) is computed
!  and the l1 norm is obtained. this norm is divided by n (the
!  number of equations). the number obtained in this way is required
!  to be less than a given tolerance e1 times (1+0.01*dd) where dd
!  is the l1 norm of the point computed to order 8. if this
!  requirement is satisfied the order 8 estimation is taken as the
!  next point. if not, a suitable value of the step h is obtained
!  and the computation is started again.
!  in any case, when the next point is computed, a prediction of
!  the step h, to be used in the next call of the routine, is
!  done.
!
!  Input data:
!
!       x  current value of the independent variable.
!    y(i) i=1,n  the current value of the dependent variable.
!       n  the dimension of the dependent variable.
!       h  the time step to be used.
!     hmin  the minimum allowed value for the absolute value of h.
!    hmax  the maximum allowed value for the absolute value of h.
!      e1  a tolerance.
!   deriv  the name of the routine computing the vector field (to
!          be declared external in the calling program).
!
!  Output data:
!
!       x  the next value of the independent variable.
!     y(i) i=1,n  the estimated next value for the dependent
!          variable.
!       h  the time step to be used in the next call of this
!          routine.
!
!  Auxiliary parameters:
!
!  r,b,f   a matrix of dimensione

!  Routines used: deriv

!  Revised to f90 free format by yu  

!********************************************************************
subroutine gr_rk78_f90(x,y,n, h,hmin,hmax,e, r,b,f,  deriv)

implicit none
integer, parameter :: dp = kind(1.d0) 
  
integer, intent(in) :: n 
real(kind=dp), intent(in) ::  hmin, hmax, e
real(kind=dp), intent(inout) :: x, y(n), h
external deriv
  
! local varaibles
integer :: ii, j, j1, jk, k, l, iter  
real(kind=dp) ::  a,bet,alfa(13), r(13,n),b(n),f(n), &
                  beta(79),c(11),cp(13),d,dd,e3,e4, fac
            
data ii/0/
save ii, alfa, beta, c, cp  
 
print*, '---- start gr_rk78_f90 ----- '         
if ( ii == 0) then 

! assignment of the constants alpha and beta, only for the first time  

  ii=1
  alfa(1)=0.d0
  alfa(2)=2.d0/27.d0
  alfa(3)=1.d0/9.d0
  alfa(4)=1.d0/6.d0
  alfa(5)=5.d0/12.d0
  alfa(6)=.5d0
  alfa(7)=5.d0/6.d0
  alfa(8)=1.d0/6.d0
  alfa(9)=2.d0/3.d0
  alfa(10)=1.d0/3.d0
  alfa(11)=1.d0
  alfa(12)=0.d0
  alfa(13)=1.d0
  
  beta(1)=0.d0
  beta(2)=2.d0/27.d0
  beta(3)=1.d0/36.d0
  beta(4)=1.d0/12.d0
  beta(5)=1.d0/24.d0
  beta(6)=0.d0
  beta(7)=1.d0/8.d0
  beta(8)=5.d0/12.d0
  beta(9)=0.d0
  beta(10)=-25.d0/16.d0
  beta(11)=-beta(10)
  beta(12)=.5d-1
  beta(13)=0.d0
  beta(14)=0.d0
  beta(15)=.25d0
  beta(16)=.2d0
  beta(17)=-25.d0/108.d0
  beta(18)=0.d0
  beta(19)=0.d0
  beta(20)=125.d0/108.d0
  beta(21)=-65.d0/27.d0
  beta(22)=2.d0*beta(20)
  beta(23)=31.d0/300.d0
  beta(24)=0.d0
  beta(25)=0.d0
  beta(26)=0.d0
  beta(27)=61.d0/225.d0
  beta(28)=-2.d0/9.d0
  beta(29)=13.d0/900.d0
  beta(30)=2.d0
  beta(31)=0.d0
  beta(32)=0.d0
  beta(33)=-53.d0/6.d0
  beta(34)=704.d0/45.d0
  beta(35)=-107.d0/9.d0
  beta(36)=67.d0/90.d0
  beta(37)=3.d0
  beta(38)=-91.d0/108.d0
  beta(39)=0.d0
  beta(40)=0.d0
  beta(41)=23.d0/108.d0
  beta(42)=-976.d0/135.d0
  beta(43)=311.d0/54.d0
  beta(44)=-19.d0/60.d0
  beta(45)=17.d0/6.d0
  beta(46)=-1.d0/12.d0
  beta(47)=2383.d0/4100.d0
  beta(48)=0.d0
  beta(49)=0.d0
  beta(50)=-341.d0/164.d0
  beta(51)=4496.d0/1025.d0
  beta(52)=-301.d0/82.d0
  beta(53)=2133.d0/4100.d0
  beta(54)=45.d0/82.d0
  beta(55)=45.d0/164.d0
  beta(56)=18.d0/41.d0
  beta(57)=3.d0/205.d0
  beta(58)=0.d0
  beta(59)=0.d0
  beta(60)=0.d0
  beta(61)=0.d0
  beta(62)=-6.d0/41.d0
  beta(63)=-3.d0/205.d0
  beta(64)=-3.d0/41.d0
  beta(65)=-beta(64)
  beta(66)=-beta(62)
  beta(67)=0.d0
  beta(68)=-1777.d0/4100.d0
  beta(69)=0.d0
  beta(70)=0.d0
  beta(71)=beta(50)
  beta(72)=beta(51)
  beta(73)=-289.d0/82.d0
  beta(74)=2193.d0/4100.d0
  beta(75)=51.d0/82.d0
  beta(76)=33.d0/164.d0
  beta(77)=12.d0/41.d0
  beta(78)=0.d0
  beta(79)=1.d0
   
  c(1)=41.d0/840.d0
  c(2)=0.d0
  c(3)=0.d0
  c(4)=0.d0
  c(5)=0.d0
  c(6)=34.d0/105.d0
  c(7)=9.d0/35.d0
  c(8)=c(7)
  c(9)=9.d0/280.d0
  
  c(10)=c(9)
  c(11)=c(1)
  cp(1)=0.d0
  cp(2)=0.d0
  cp(3)=0.d0
  cp(4)=0.d0
  cp(5)=0.d0
  cp(6)=c(6)
  cp(7)=c(7)
  cp(8)=c(8)
  cp(9)=c(9)
  cp(10)=c(10)
  cp(11)=0.d0
  cp(12)=c(1)
  cp(13)=c(1)
    
else ! (ii .ne. 1)
  iter = 0  
  do 
    jk = 1
    iter = iter + 1 
     
    do j = 1, 13 ! 3
     
      do l = 1,n !6
        b(l) = y(l)
      enddo 

      a = x + alfa(j)*h
         
      if(j .ne. 1) then  
        
        j1 = j - 1

        do k = 1, j1, 1
          jk = jk + 1
          bet = beta(jk) * h
          do l = 1, n
            b(l) = b(l) + bet * r(k,l)
          enddo 
        enddo 
           
      endif

      call deriv (a,b,n,f)
      do   l = 1, n
        r(j,l) = f(l)
      enddo 
      
    enddo  ! 3     

    d  = 0
    dd = 0
         
    do l = 1, n !--1
      
      b(l)=y(l)
      f(l)=y(l)
      
      do k=1,11 !--5 
        bet = h * r(k,l)
        b(l) = b(l) + bet*c(k)
        f(l) = f(l) + bet*cp(k)
      enddo  ! -5
      
      f(l) = f(l) + h *( cp(12)*r(12,l) + cp(13)*r(13,l) )
      d  = d + dabs( f(l) - b(l) )
      dd = dd + dabs( f(l) )
    enddo !--1

    d   = d/n 
    fac = 1.d0 + dd*1.d-2
    e3  = e*fac
    
    ! --- ck
    print*, 'f90: jk, e3,d,h', jk, e3,d,h;  read* 
            
    if (dabs(h) .lt. hmin  .or.  d .lt. e3)  exit !go to 7
    h = h * 0.9d0 * (e3/d)**0.125d0
    if( dabs(h) .lt. hmin )  h = hmin * h / dabs(h)
    
  enddo 
       
  x = x + h ! -- 7 
  if(d.lt.e3)  d = dmax1(d, e3/256)
  h=h*0.9d0*(e3/d)**0.125d0
         
  if(dabs(h).gt.hmax)  h=hmax*h/dabs(h)
  if(dabs(h).lt.hmin)  h=hmin*h/dabs(h)
   
  do l=1,n
    y(l)=f(l)
  enddo 
       
  b(1)=d
endif 
  
return
end subroutine gr_rk78_f90  
         
