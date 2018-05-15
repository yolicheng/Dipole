c  with position and velocity
c  calculate 6 orbital elements of the 2-body problem

c  Input Variables
c 	sv  	state vector, 1*6, (x,y,z,vx,vy,vz), 
c 	oe  	orbital elements (h,i,omega,e,w,theta)
c 	xmu 	is the mass prarameter corresponding to mu*(M*m). Depending 
c 		on the units of the units of the position and velocity. It
c 		has a value of other. If the coordinates are canonical then
c 		xmu = 1 
c  Where,
c    	hm 	specific angular momentum 
c 	xi 	inclination
c 	omega	right acension(RA) of the ascending node 
c  	em	eccentricity
c 	w 	argument of perigee
c 	theta   ture anomaly
c 
c  An alternative is, replace h and w to  a and Mean Anomaly M
 
	subroutine sv2oe(sv,xmu, hm,xi,omega,em,w,theta)
	implicit real*8(a-h, o-z)

!with variables starts out of this domin, put a x before the name
	dimension sv(6),h(3),xn(3),e(3)
	pi = 3.141592653589793D0
		
!      DATA PI/3.141592653589793D0/, PI2/6.283185307179586/
!      SAVE PI,PI2
!      RETRO=.FALSE.
      
c	here, for 2 body problem, the mu should be the moon w.p.t the sun
	xmu = 1 ! this should be update to the real value 
	x=sv(1); y=sv(2); z=sv(3); vx=sv(4); vy=sv(5); vz=sv(6) 
c 	h: angular momentum
	rm = dsqrt(x*x + y*y + z*z )
	vm = dsqrt(vx*vx + vy*vy + vz*vz )
	vr = ( x*vx + y*vy + z*vz )/ rm
	h(1) = y*vz - z*vy !hx
	h(2) = z*vx - x*vz !hy
	h(3) = x*vy - y*vx !hz
	hm = dsqrt(h(1)**2 + h(2)**2 + h(3)**2)
	
	if( dabs(hm) .lt. 1.d-10) then
	  write(*,*) 'problem: orbit degenerated into a line'
	  stop
	endif
	
c 	i: inclination
	xi = dacos(h(3)/hm)  
	
c 	Omega: right acension(RA) of the ascending node 
	xn(1) = -h(2) !nx
	xn(2) =  h(1) !ny
	xn(3) =  0.d0  !nz
	xnm = dsqrt(xn(1)**2 + xn(2)**2 + xn(3)**2)
	omega = dacos(xn(1)/xnm)
 
	if(xn(2) .lt. 0)  omega = 2*pi - omega
	
c 	e: eccentricity
c 	to reduce the cost of computation, use aux=v*v-xmu/rm, compute once
	aux = vm*vm-xmu/rm 
!	e(1) = ( aux*x - rm*vr*vx )/xmu !ex
!	e(2) = ( aux*y - rm*vr*vy )/xmu !ex
!	e(3) = ( aux*z - rm*vr*vz )/xmu !ex
!	write(*,*) 'e1', e 
	e  = ( aux*sv(1:3) - rm*vr*sv(4:6) )/xmu !e 
	em = dsqrt( e(1)**2 + e(2)**2 + e(3)**2 )
	
c	w : argument of perigee
	w = dacos( ( xn(1)*e(1) + xn(2)*e(2) + xn(3)*e(3) ) /xnm/em )
	if(e(3) .lt. 0.d0 )  w = 2*pi - w
	
c 	theta : ture anomaly
	theta = dacos( ( ( e(1)*x + e(2)*y + e(3)*z ) /em/rm ) )
	if(vr .lt. 0.d0)  theta = 2*pi - theta
	
c	p = hm*hm
c 	oe = (/hm,xi,omega,em,w,theta/) !return the scalar
	
	return
	end		
	
	
