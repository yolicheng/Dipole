!****************************************************************************  
!
!  subroutine gr_rtbp 
!      
!  vectorfield of the spatial rtbp and its variational equations if desired.
!
!  in this rtbp the big primary is located at (mu,0,0) with mass 1-mu and
!  the small one at (mu-1,0,0) with mass mu.
!
!  input parameters:   
!
!    t     rtbp time
!    x(*)  rtbp-variational coordinates (x(1),x(2),...,x(42))
!    n     number of equations (6 or 42 according if we want just the
!  vectorfield or moreover the variational flow).
!    mu   rtbp mass parameter through common rtbpc. 
!
!  output parameters:
!
!    f(*)  vectorfield. the first 6 components correspond to the rtbp
!  equations, and the remaining ones to the variational equations
!  stored by columns. 
!
!  note: if n equals 6 the variational equations are skipped 
!
!  revised by yu to free format and use the jacobi matrix and matmul--- 20160401  
!***************************************************************************** 
subroutine gr_rtbp(t,x,n,f)
implicit none 

integer, parameter :: dp=kind(1.d0)

! 	input and output declaration 
integer, intent(in) :: n ! demension of x
real(kind=dp), intent(in)  :: t, x(n)
real(kind=dp), intent(out) :: f(n)
 

!	local  parameters
real(kind=dp) ::  y1, y12, y2, y3, r1, r1a, r15, r2, r2a, r25, p1, p2, q, & ! velocity + acceleration
		  rr1, rr2, qp, qpp, a(6,6), phi(6,6), dphi(6,6) ! variational matrix


!--------test------------
real(kind=dp) :: mu ! for test only 
integer :: i ! 
! include this subroutine in the specified systemt-based module to get the mass 
! parameter mu 

! used to check pomod 
mu =  0.1215058561000000d-01 
!------------------



! velocity
f(1) = x(4)
f(2) = x(5)
f(3) = x(6)

! acceleration
! omega = 1/2(x**2 + y**2) + (1-mu)/r1 + mu / r2
! r1**2 = (x-mu)**2 + y**2 + z**2
! r2**2 = (x-mu+1)**2 + y**2 + z**2

! f4 = ax =  2vy + d omega / d x (donated as omega_x)
! f5 = ay = -2vx + d omega / d y (donated as omega_y)
! f6 = az =  d omega / d z       (donated as omega_z)

! omega_x = x - 1/r1**3 * (1-mu) * (x-mu) - 1/r2**3 * mu * (x-mu+1) 
! omega_y = y - 1/r1**3 * (1-mu) * y      - 1/r2**3 * mu * y
! omega_z =    -1/r1**3 * (1-mu) * z      - 1/r2**3 * mu * z

     
y1  = x(1) - mu
y12 = y1*y1

y2  = x(2)*x(2)
y3  = x(3)*x(3)

r1 = y12 + y2 + y3 ! the square of the distance to the big primary
r1a = r1 * dsqrt(r1)
r15 = r1a * r1
r1 = r1a

r2 = (y1+1.d0)**2 + y2 + y3
r2a = r2 * dsqrt(r2)
r25 = r2a*r2
r2  = r2a

p1 = (1.d0-mu) / r1
p2 = mu/r2
q  = -1.d0*(p1 + p2)

f(4) = 2.d0 * x(5) + x(1) - y1 * p1 -(y1 + 1.d0) * p2
f(5) =-2.d0 * x(4) + x(2) * (1.d0 + q)
f(6) = x(3) * q

if (n.eq.6) return


! --------------------------- variational matrix ------------------------------
phi =  reshape(x(7:42), (/6,6/))  ! do the reshape to work with matrix

!!----------------------------------
!print*, 'check phi'
!do i = 1, 6
!  write(*,  '(6f20.14)')  phi(i,:)
!enddo 
!print*; read*
!!----------------------------------

! jacobi matrix a
! a = |    0         i ( 3 x 3 )|         l =  | 0   -2    0 |
!     | omega_xx         l      |              |-2    0    0 |
!                                              | 0    0    0 |

a  = 0.d0 
a(1,4) = 1.d0 
a(2,5) = 1.d0
a(3,6) = 1.d0

a(4,5) =  2.d0 
a(5,4) = -2.d0

!  comments:  
!     d ( 1/r1) / dx =  1/r1 *(x-mu) 
!     d ( 1/r1) / dy =  1/r1 * y    
!     d ( 1/r1) / dz =  1/r1 * z 

!     d ( 1/r1**3) / dx =  -3/r1**5 *(x-mu) 
!     d ( 1/r1**3) / dy =  -3/r1**5 * y    
!     d ( 1/r1**3) / dz =  -3/r1**5 * z 

! derivative of the variational matrix  
! -- omega_x      
! omega_xx = 1 -  1/r1**3 * (1-mu) +  3/r1**5 * (1-mu) * (x-mu)**2 
!              -  1/r2**3 *  mu    +  3/r2**5 *  mu    * (x-mu+1)**2 

! omega_xy =   3/r1**5 * (1-mu) * (x-mu)   * y 
!	     + 3/r1**5 *  mu    * (x-mu+1) * y
 
! omega_xz =   3/r1**5 * (1-mu) * (x-mu)   * z 
! 	     + 3/r2**5 *  mu    * (x-mu+1) * z

! -- omega_y
! omega_yx = omega_xy   

! omega_yy = 1 -  1/r1**3 * (1-mu) +  3/r1**5 * (1-mu) * y**2 
!              -  1/r2**3 *  mu    +  3/r2**5 *  mu    * y**2 
 
! omega_yz =   3/r1**5 * (1-mu) * y * z 
!   	     + 3/r2**5 *  mu    * y * z 

! -- omega_z
! omega_zx =  omega_xz
! omega_zy =  omega_yz

! omega_zz =  -  1/r1**3 * (1-mu) +  3/r1**5 * (1-mu) * z**2 
!             -  1/r2**3 *  mu    +  3/r2**5 *  mu    * z**2 

rr1 = 3.d0*(1.d0-mu) / r15
rr2 = 3.d0*mu / r25

qp  = rr1 * y1 + rr2 * (y1+1.d0)
qpp = (rr1+rr2) * x(2)

a(4,1) =  1.d0 + q + rr1 * y12 + rr2 * (y1+1.d0)**2 ! omega_xx
a(4,2) = qp * x(2) ! omega_xy
a(4,3) = qp * x(3) ! omega_xz

a(5,1) = a(4,2) ! omega_yx = omega_xy
a(5,2) = 1.d0 + q + (rr1+rr2)*y2 ! omega_yy
a(5,3) = qpp*x(3) ! omega_yz

a(6,1) = a(4,3) ! omega_zx = omega_xz
a(6,2) = a(5,3)! omega_zy = omega_yz
a(6,3) = q + (rr1+rr2)*y3 ! omega_zz

! d phi = a * phi
dphi = matmul(a, phi)  
f(7:42) = reshape( dphi, (/36/) )
 
return
end

