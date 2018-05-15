  	subroutine poinc_z0(yi,imax, tmax, ind, z0, xfi,yf, ispc, deriv) 
  	
!c Determination of the imax-th passage of the orbit through the poincare section 
!!  defined by subroutine sect (here, z = z0) 
!c   only one intersection is record here!
!     
!c 	input parameters:
!c yi(*)  initial point 
!c imax	 the number of time for intersecting the section

!c  	output parameters:
!c xfi 	 time spent by the orbit to go from yi to yf
!c yf(*)  first cut of the orbit which passes by yi with the surface 
!c	 of section

!c hminim	 minimus step used in the integration of the p.o.
!c
!c function used: gr_rk78, sect

!c to make this more general, use deriv to represent to subroutine to 
!c compute the vector field

!c **********************************************************************  
!   
 
	implicit real*8(a-h,o-z)
	dimension yi(42),yf(42), dg(42), r(13,42),b(42),f(42)  
	
 	external deriv ! here should be gr_lf (the vector field for lorentz force problem)
 
  	neq = 42
 	ispc = 1
 	
! 	print*,'poinc_z0, x0', yi(1:6)
!	read(*,*) aaa

 	call sect(yi,ind, z0, g,dg)
! 	print*, 'ck sect, ind,z0, g', ind, z0, g
 	
! 	read(*,*) aaa

	if(dabs(g) .lt. 1e-9) g=0.d0  

	h    = 1.d-2
	hmin = 1.d-5
	hmax = 1.d-0
	e    = 1.d-13 ! there is something wrong here, the error in energy is supposed to be less than 1.d-13

    	x    = 0.d0 ! start time
    	iint = 0  !
	
1    	gi = g
!SUBROUTINE GR_RK78 (X,Y,N,H,HMI,HMAX,E1,R,B,F,DERIV)
	call gr_rk78(x,yi,neq,h,hmin,hmax,e,r,b,f, deriv)
 	call sect(yi, ind, z0, g,dg)
 	
! 	write(*,'(3e20.12)') x, gi, g !yi(1:6)
  	if(gi*g .ge. 0 )  goto 1
  	
!  	print*,  'crossing the poincare section', gi, g
!  	read(*,*) aaa
!  	
 	if (x > tmax ) then 
 	  ispc = 0
 	  return
 	endif

	iint = iint+1 ! number of time of intersection with y=0 plane
	if(iint .lt. imax) goto 1
    
! use newton method, with yi as initial guess, to find the intersecion
! with df*dg, we get vy, the derivative of y respect to t, and dt

! 20160115, update for newton method, if there is a sign change in g
!           use the previous point, (x-dh, ypre) as the new initial guess
	
	
2	call deriv(x, yi, 42, f)

	dy = 0
	do i=1,42
	  dy = dy + f(i)*dg(i)
	enddo
	
	dh = - g/dy
	
!	print*, 't, dh', x, dh

! here is the newton method process
	call gr_rk78(x,yi,neq,dh,hmin,hmax,e,r,b,f,deriv)
 	call sect(yi, ind, z0, g,dg)
 	
! 	print*, 'newton, t,g', x,g
 	if(dabs(g) .gt. 1.d-10)   goto 2
! ********************************************************** 	
 	
! we get y=0, x is tp, yi is the intersection	
 	yf = yi
 	xfi = x
c  	write(*,*)'Poinc finished, g, xfi,yf:',g, xfi, yf(1:6) 
 	return
 	end
