	subroutine dfcr(y0,imax,e2,id, xfi,yf,g, deriv)

C++++++++++++++++++++++++++  
c
c**********************************************************************
c  Using modified Newton method to refine the p.o, provided with the initial conditon.
c    with 3 free variables instead of 2 variables( inderterminant equations), least-square solution

c  	Input Variables:
c  y0(6)	initial condition(x,y=0,z, vx=0, vy, vz=0)
c  imax 	the times of the crossing through the section(y=0)
c  e2		error control for vx,vz to be zero(1.d-10)

c  	Output Variables:
c  y0, yf	the refined initial and final conditon
c  xfi		half period or one period, depends on the number of intsec-imax
c  g(2,3)	JACOBIAN MATRIX OF F(*), 2 by 3,  3 variables and 2 equations
c  
c     the equations: f1= vx = 0, f2 = vz = 0 symmetric wrt xz plane : halo and planar lyapunov orbit, 1st crossing in T/2
c                    f1= z = 0,  f2 = vx = 0 symmetric wrt x-aixs : vertical lyapunov orbit, first crossing in T/4
  
c The point is get the derivative of the  target variables with respect to the control variables


c  20141215
c     integrate to the first intersection with y=0, then, use it as the inital conditon.


c  20150408 -Finally caution!
c    deal with vt_lypt orbit individually, different from halo and planar lypt

c  if the initial condition changes to x,vy,vz,then the g should 
c    be different, the dfcr routine still can't be share for three kinds of orbits  
c    between 10 and 12 are the special approach to compute g for vt_lypt

c  function used: poinc, deriv, fctn, deltx
 
c***********************************************************************
 	implicit real*8 (a-h,o-z)
    	
c both the state and the variational equations are desired, so n=42
	integer, parameter:: n=42 
	dimension y0(6),dx(3)
 	dimension f(2),g(2,3),yi(n),yf(n),pf(n)
   
   
!  	print*, 'start dfcr, y0!',y0  
c because we need the jacobi matrix with fctn subroutine, and poinc is 
c called in it, so no need to call poinc before fctn
   	iter = 0 
   	
10   	if (id .ne. 3)  goto 12  ! planar lypt or halo

! vertical lyapunov orbit, use the t/4 as target point
	
! vt_lyap cannot be refined by setting the intersection to 2, so we choose to do the 
!   the refinement also at the first intersection of the y=0 plane, 
!   but the target variables are z=0,vx=0, thus the target point is actually at T/4 
  
  
!  	write(*,*) 'looking for the first intersection for vtly'--- need to be modified!!
! the point is to compute the derivative G
 	yi = 0.d0
 	yi(1:6) = y0  
 	
 	!yi(1:5:2) = y0  
	yi(7:42:7) = 1.d0
 	call poinc(yi, 1, xfi, yf ) 
 	call gr_lf(0.D0, yf, 42, pf) ! at this moment, keep the name gr_lf, but replace with general deriv later
 	
c Instead of using ftcn subroutine to get f,g as for lyap and halo orbits
c User-define f and g seperately. 
 
!     write(*,*) 'found the first intersection for vtly'
      	
	dv2 = yf(3)*yf(3) + yf(4)*yf(4)
	if (dsqrt(dv2) .lt. e2) return

    	f = yf(3:4)  ! f1=z = 0, f2=vx = 0 for vertical lyapunov
     	f1 = -pf(3)/pf(2)  !vz/vyf
	f2 = -pf(4)/pf(2)  !axf/vyf

! there should be a subroutine to compute the G as required
! todo ---  20151223
	
! G is the Jacobi matrix of f(f1=z, f2=vx) with respect to x0,z0,vy0
! the variational eqs are stored by columns, be careful with the index

! instead of using one-dimensional index, use index of 2-dimensional matrix is better 
 
!******************** start here tomorrow ****************************** 
	g(1,1) = yf(9)  + f1*yf(8)  ! d z/d x  = phi(3,1) - vz/vfy*phi(2,1)  
	g(1,2) = yf(21) + f1*yf(20) ! d z/d z  = phi(3,3) - vz/vfy*phi(2,3) 
	g(1,3) = yf(33) + f1*yf(32) ! d z/d vy = phi(3,5) - vz/vfy*phi(2,5) 
	
	g(2,1) = yf(10) + f2*yf(8)  ! d vx/d x  = phi(4,1) - axf/vfy*phi(2,1)
	g(2,2) = yf(22) + f2*yf(20) ! d vx/d z  = phi(4,3) - axf/vfy*phi(2,3)
	g(2,3) = yf(34) + f2*yf(32) ! d vx/d vy = phi(4,5) - axf/vfy*phi(2,3)
      
	goto 15

	
! common part for planar lyapunov and halo, just to get G to pass to DELTX
!	subroutine fctn(x,ind,f,g,z,y,xfi,imax, deriv )
 
12	call fctn(y0, 0, f, g, yf, pf, xfi, imax, deriv )  

	dv2 = yf(4)*yf(4) + yf(6)*yf(6)
	if (dsqrt(dv2) .lt. e2)  return
      
15	ipo = 0
 	call deltx(f,g, ipo,dx)
 	
	y0(1:5:2) = y0(1:5:2) + dx
	iter = iter+1
	if(iter.lt. 10 .and. ipo .eq. 0) goto 10
	
!for Newton method, if the iteration is great than 5 or 6, then there is something wrong, provided that the initial guess is a good one, which is 
! close to the solution
	
	if (iter .gt. 10) then 
	  write(*,*) 'Maximum iteration exceeded!'
	else if (ipo .eq. 1) then
	  write(*,*) 'dx less than the precision !'
	endif
	
  	return
  	end 
