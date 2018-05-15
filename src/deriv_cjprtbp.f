subroutine deriv_cjprtbp(x, dcj)
! ******************************************************
! computation of the differential of Jacobi CSonstant w.r.t. the all the coordinates 
! in the planar rtbp.

! in this rtbp the big primary is at (xmu,0,0) with
! mass 1-xmu, and the small one at (xmu-1,0,0) with
! mass xmu.

!******************************************************
use dp_mod 
use rtbpconst_mod, only : mu 
implicit none 

real(kind=dp), intent(in)       ::  x(4)
real(kind=dp), intent(out)      ::  dcj(4)

!     local  parameters
real(kind=dp) ::  y1, y12, y22, r1, r13, r15, r2, r23, r25, p1, p2, q & ! velocity + acceleration

y1  = x(1) - mu
y12 = y1*y1

y22  = x(2)**2

r1  = y12 + y22   ! the square of the distance to the big primary
r13 = r1 * dsqrt(r1)
r15 = r13 * r1


r2  = (y1+1.d0)**2 + y2 !+ y3
r23 = r2 * dsqrt(r2)
r25 = r23*r2

p1 = (1.d0-mu) / r13
p2 = mu/r23
q  = -1.d0*(p1 + p2)

! this is from gr_cjrtbp.f --- Provided by Gerard
!  OMEGA = (X1+X2)*.5+(1-MU)/R1+MU/R2+ .5*MU*(1-MU)
!  CJA = 2*OMEGA- Vx**2 - Vy**2 - Vz**2 
        
! the differential of d CJ / d {x, y} = 2 * d Omega / d {x, y}
dcj(1) = 2 * ( x(1) - 2y1 * p1 - (y1 + 1.d0) * p2 ) ! dcjdx 
dcj(2) = 2 * ( x(2) * (1.d0 + q) )                  ! d cj / d y

!  d CJ / d {vx, vy}  = - 2 {vx, vy}
dcj(3) = -2.d0 * x(3) ! d cj / d vx
dcj(4) = -2.d0 * x(4) ! d cj / d vy

return
end

