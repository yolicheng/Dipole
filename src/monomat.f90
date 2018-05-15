subroutine monomat(yi, n, tp, mmat, deriv, gr_cj)
! compute the Monodramy matrix of the p.o. of period tp
! 	Input Varialbles
!   yi	initial state(x0,y0,z0,vx0,vy0,vz0)
!   n      dimension of the phase space
!   tp 	period of the P.O.      

! 	Output Variables
!   mmat	Monodramy matrix, which is also the state transition matrix after one period

!  Routine used: dfcr, fctn, champ, adams, gr_cjlf

 
! 20150408 --- the vertical lyapunov checked  
!     imax is set to 1 to get the first intersection, pay attention to
!     the computation of g for new PO

! Version 0.1     : 20150525
! Finally revised :  2017-05-21 12:13:57 
!***********************************************************************
use dp_mod
implicit none
  
integer, intent(in)        :: n
real(kind=dp), intent(in)  :: yi(n), tp 
real(kind=dp), intent(out) :: mmat(n,n)

external ::  deriv, gr_cj ! here should be gr_lf (the vector field for lorentz force problem)

integer :: nvar, i , debug 
real(kind=dp) :: y(n*(n+1)), t, h, hmin, hmax, e, r(13,n*(n+1)), b(n*(n+1)), f(n*(n+1)), &

                 cj, cj2, dcj

debug = 0
nvar = n*(n+1)
if(debug == 1) print*, 'n, nvar', n, nvar

h    = 1.d-3

! this is the same... should we use the value from the module pomod? no...
hmin = 1.d-6
hmax = 1.d-1
e    = 1.d-13

t    = 0.d0 ! start time is 0
! initialize vector field and variational matrix(identity)
y(1:n)          = yi
y(n+1:nvar)     = 0.d0
y(n+1:nvar:n+1) = 1.d0

if(debug == 1)  then 
  write(*,'(6d20.10)') y
  print*; read*
endif 

! check if the energy is conversative 
call gr_cj(yi, cj)
 
do while (dabs(t+h) .le. dabs(tp))

!SUBROUTINE GR_RK78 (X,Y,N,H,HMI,HMAX,E1,R,B,F,DERIV)
  call gr_rk78(t,y, nvar,h,hmin,hmax,e,r,b,f, deriv)
  if(debug == 1)   print*, t, y(1:n)
enddo  


if (dabs(t-tp) > 1.d-9 ) h = tp - t 

! h for sure will be smaller than the previous h, so 1 step is enough
call gr_rk78(t,y, nvar,h,hmin,hmax,e,r,b,f, deriv )

! check if the energy integral is conservative
call gr_cj(y(1:n),cj2) 
dcj = cj2 - cj 


mmat = reshape(y(n+1:nvar), (/n,n/)) 
  
! check the final state, should be the same with the initial state 
if(debug == 1) then
  print*,'check  if the energy integral is conservative -- monomat'  
  print*, 'dcj, cj, cj2', dcj, cj, cj2
  print*,'check  the final state -- monomat'  
  print*, 'yf, y0', y(1:6), yi(1:6)
   
  do i = 1, n
    write(*,'(6d20.10)') mmat(i,:)
  enddo  
endif 

return
end subroutine monomat
	 
	 
