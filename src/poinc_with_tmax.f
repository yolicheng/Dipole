  	subroutine poinc(yi,imax, tmax, xfi,yf, ispc, deriv) 
c
C++++++++++++++++++++++++++  
c Determination of the imax-th passage of the orbit through the poincare section defined by subroutine sect(here,y=0) 
c   only one intersection is record here!
     
c 	input parameters:
c yi(*)  initial point 
c imax	 the number of time for intersecting the section

c  	output parameters:
c xfi 	 time spent by the orbit to go from yi to yf
c yf(*)  first cut of the orbit which passes by yi with the surface 
c	 of section

c hminim	 minimus step used in the integration of the p.o.
c
c function used: gr_rk78, sect

c to make this more general, use deriv to represent to subroutine to 
c compute the vector field

c **********************************************************************  
   
 !	use lfmod, only : beta, cs, sgn ! we don't need to use it in poinc
 !        GR_LF is the only subroutine that uses these three parameters
 
	implicit real*8(a-h,o-z)
	dimension yi(42),yf(42), ypre(42), dg(42), r(13,42),b(42),f(42) !, cj, cj2, dcj
	
 	external deriv ! here should be gr_lf (the vector field for lorentz force problem)
	
!	print*, 'the first poinc, Initial state & variational matrix'
!	print*, '(the transpose(stored by 1 dimensioanal array))!'
!	print*, yi(1:6)
	
!	do i = 1, 6
!	  write(*,'(6f10.6)') yi(6*i+1: 6*i+6) 
!	enddo 
	
	!read(*,*) aaa
  	neq = 42
 	ispc = 1
 	
	call sect(yi,g,dg)
	if(dabs(g) .lt. 1e-9) g=0.d0  

	h    = 1.d-2
	hmin = 1.d-10
	hmax = 1.d-1
	e    = 1.d-13 ! there is something wrong here, the error in energy is supposed to be less than 1.d-13

    	x    = 0.d0 ! start time
    	iint = 0  !
 
 
	!call gr_cjlf(yi(1:6), cj)
	
1    	gi = g
	xpre = x
	ypre = yi
!SUBROUTINE GR_RK78 (X,Y,N,H,HMI,HMAX,E1,R,B,F,DERIV)
	call gr_cjlf(yi(1:6), cj)
	call gr_rk78(x,yi,neq,h,hmin,hmax,e,r,b,f, deriv)
	
!	print*, 't,y', x, yi(1:6)
	write(*,'(7f10.6)') x, yi(1:6)
	
!	if(x > tmax) then ! no intersecion during the allowable integration time
!	   ispc = 0
!	   print*, 'exceed maximum time, t, tmax', x, tmax
!	   read(*,*) aaa
!	   return 
!	endif
	
! check if the energy integral is conservative
	call gr_cjlf(yi(1:6),cj2) 
	dcj = cj2 - cj 

!	print*,'check  if the energy integral is conservative  --  poinc'  
! 20151225 -- chiristmas!
! there is some problem here, please check the routine really carefully
! because the energy is not conservative!!!!  the error is of order 1.d-5 every step	
	!print*, 'dcj, cj, cj2', dcj, cj, cj2
	!read(*,*) aaa

 	call sect(yi,g,dg)
  	if(gi*g .ge. 0 .or. x < 3.d-1)  goto 1 	
!  	if(gi*g .ge. 0 )  goto 1
  	
  	print*,  'crossing the poincare section', gi, g
!  	read(*,*) aaa
  	
	iint = iint+1 ! number of time of intersection with y=0 plane
	if(iint .lt. imax) goto 1
    
! use newton method, with yi as initial guess, to find the intersecion
! with df*dg, we get vy, the derivative of y respect to t, and dt

! 20160115, update for newton method, if there is a sign change in g
!           use the previous point, (x-dh, ypre) as the new initial guess
	
	x  = xpre
	yi = ypre
	g = gi
	
2	call deriv(x, yi, 42, f)

	dy = 0
	do i=1,42
	  dy = dy + f(i)*dg(i)
	enddo
	
	dh = - g/dy
	
	print*, 't, dh', x, dh

! here is the newton method process
        hmin = dabs(dh)
	call gr_rk78(x,yi,neq,dh,hmin,hmax,e,r,b,f,deriv)
 	call sect(yi,g,dg)
 	
 	print*, 'newton, t,g', x,g
 	if(dabs(g) .gt. 1.d-10)   goto 2
! ********************************************************** 	
 	
! we get y=0, x is tp, yi is the intersection	
 	yf = yi
 	xfi = x
c  	write(*,*)'Poinc finished, g, xfi,yf:',g, xfi, yf(1:6) 
 	return
 	end
