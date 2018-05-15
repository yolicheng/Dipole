!****************************************************************************  
!  vectorfield of the Planar RTBP and its variational equations if desired.
!
!  in this rtbp the big primary is located at (mu,0,0) with mass 1-mu and
!  the small one at (mu-1,0,0) with mass mu.
!
!      Input Parameters:   
!
!    t     rtbp time
!    x(*)  rtbp-variational coordinates (x(1),x(2),...,x(42))
!    n     number of equations (6 or 42 according if we want just the
!  vectorfield or moreover the variational flow).


!     Output Parameters:
!
!    f(*)  vectorfield. the first 6 components correspond to the rtbp
!  equations, and the remaining ones to the variational equations
!  stored by columns. 
!
!     Module-based Parameter
!  mu       rtbp mass parameter through common Module massrat_rtbp_mod
!  ndim     the dimension of the phase space, 4: planar problem; 6: spatial problem

!  note: if n equals 6(4 for planar case) the variational equations are skipped 
!
!  revised by Yu to free format and use the jacobi matrix and matmul--- 20160401
  
!***************************************************************************** 
subroutine gr_prtbp(t, x, n, f)
use  rtbpconst_mod  ! ndim and mu 
implicit none 

integer, parameter  ::  dp = kind(1.d0)

! 	input and output declaration 
integer, intent(in) :: n ! demension of x
real(kind=dp), intent(in)  :: t, x(n)
real(kind=dp), intent(out) :: f(n)
 

!     local  parameters
real(kind=dp) ::  y1, y12, y2, y3, r1, r13, r15, r2, r23, r25, p1, p2, q, & ! velocity + acceleration
                rr1, rr2, qp, qpp,  a(ndim, ndim), phi(ndim, ndim), dphi(ndim, ndim), & ! variational matrix
                f2(n), a2(ndim, ndim) ! ck 
                
!real(kind=dp), allocatable ::   a(:,:), phi(:,:), dphi(:,:), a2(:,:) ! variational matrix


!--------test------------
integer :: i,  ndimhf, indvx, indvy, debug 
debug = 0
        
! for the dimension of the Pos + Vel, and the index for x and y components
ndimhf = ndim / 2
indvx = 1 + ndimhf
indvy = 2 + ndimhf

if(debug == 1) then
  print*, 'mu=', mu, 'ndim=', ndim 
  print*, 'Dimension of the state: ',    ndim
  print*, 'check the index of vx and vy: ', indvx, indvy
  read*
endif 
                
! vx - vy 
f(1) = x( indvx )
f(2) = x( indvy )

!-- vz(ndim=6)
if(ndim == 6) f(3) = x(6)
      
                
!f(1) = x(4)
!f(2) = x(5)
!f(3) = x(6)

! acceleration
! Omega = 1/2(x**2 + y**2) + (1-mu)/r1 + mu / r2
! r1**2 = (x-mu)**2 + y**2 + z**2
! r2**2 = (x-mu+1)**2 + y**2 + z**2

! f4 = ax =  2vy + d Omega / d x (denoted as Omega_x)
! f5 = ay = -2vx + d Omega / d y (denoted as Omega_y)
! f6 = az =  d Omega / d z       (denoted as Omega_z)

! Omega_x = x - 1/r1**3 * (1-mu) * (x-mu) - 1/r2**3 * mu * (x-mu+1) 
! Omega_y = y - 1/r1**3 * (1-mu) * y      - 1/r2**3 * mu * y
! Omega_z =    -1/r1**3 * (1-mu) * z      - 1/r2**3 * mu * z

     
y1  = x(1) - mu
y12 = y1*y1

y2  = x(2)**2

y3 = 0.d0 
if(ndim == 6) y3 = x(3)**2

!y3  = x(3)*x(3)

r1  = y12 + y2 + y3 ! the square of the distance to the big primary
r13 = r1 * dsqrt(r1)
r15 = r13 * r1


r2  = (y1+1.d0)**2 + y2 + y3
r23 = r2 * dsqrt(r2)
r25 = r23*r2

p1 = (1.d0-mu) / r13
p2 = mu/r23
q  = -1.d0*(p1 + p2)

! ax - ay
f(indvx) =  2.d0 * x(indvy) + x(1) - y1 * p1 -(y1 + 1.d0) * p2
f(indvy) = -2.d0 * x(indvx) + x(2) * (1.d0 + q)

!-- az(ndim=6)
if(ndim == 6)  f(6) = x(3)*q


if(debug == 1) then  ! --ckd
  print*, 'Check the vector field: '
  print*, 'x0= ', x
  print*, 'vf=', f 
  print*, 'vf by element computation: ', f2
  read*
endif 

if (n .eq. ndim) return


! --------------------------- variational matrix ------------------------------
phi =  reshape( x(ndim+1 : n), (/ndim, ndim/))  ! do the reshape to work with matrix

!phi =  reshape( x(7:42), (/6,6/))  ! do the reshape to work with matrix

!!----------------------------------
!print*, 'check phi'
!do i = 1, 6
!  write(*,  '(6f20.14)')  phi(i,:)
!enddo 
!print*; read*
!!----------------------------------

! jacobi matrix a
! a = |    0         i ( 3 x 3 )|         l =  | 0   -2    0 |
!     | Omega_xx         l      |              |-2    0    0 |
!                                              | 0    0    0 |

a  = 0.d0 
do i = 1, ndimhf, 1
  a(i, i+ndimhf) = 1.d0 
end do

a(indvx, indvy) =  2.d0 
a(indvy, indvx) = -2.d0

!  comments:  
!     d ( 1/r1) / dx =  1/r1 *(x-mu) 
!     d ( 1/r1) / dy =  1/r1 * y    
!     d ( 1/r1) / dz =  1/r1 * z 

!     d ( 1/r1**3) / dx =  -3/r1**5 *(x-mu) 
!     d ( 1/r1**3) / dy =  -3/r1**5 * y    
!     d ( 1/r1**3) / dz =  -3/r1**5 * z 

! derivative of the variational matrix  
! -- d Omega_x   /  d {x, y, z}   
! Omega_xx = 1 -  1/r1**3 * (1-mu) +  3/r1**5 * (1-mu) * (x-mu)**2 
!              -  1/r2**3 *  mu    +  3/r2**5 *  mu    * (x-mu+1)**2 

! Omega_xy =    3/r1**5 * (1-mu) * (x-mu)   * y 
!	       +  3/r1**5 *  mu    * (x-mu+1) * y
 
! Omega_xz =   3/r1**5 * (1-mu) * (x-mu)   * z 
! 	     +   3/r2**5 *  mu    * (x-mu+1) * z

! -- d Omega_y  /  d {x, y, z}  
! Omega_yx = Omega_xy   

! Omega_yy = 1 -  1/r1**3 * (1-mu) +  3/r1**5 * (1-mu) * y**2 
!              -  1/r2**3 *  mu    +  3/r2**5 *  mu    * y**2 
 
! Omega_yz =   3/r1**5 * (1-mu) * y * z 
!   	     + 3/r2**5 *  mu    * y * z 

! -- d Omega_z /  d {x, y, z}  
! Omega_zx =  Omega_xz
! Omega_zy =  Omega_yz

! Omega_zz =  -  1/r1**3 * (1-mu) +  3/r1**5 * (1-mu) * z**2 
!             -  1/r2**3 *  mu    +  3/r2**5 *  mu    * z**2 

rr1 = 3.d0 / r15 * (1.d0-mu) 
rr2 = 3.d0 / r25 * mu 

qp  = rr1 * y1 + rr2 * (y1+1.d0) 
qpp = (rr1+rr2) * x(2)

! d vx / d (x, y) 
a(indvx,1) =  1.d0 + q + rr1 * y12 + rr2 * (y1+1.d0)**2 ! Omega_xx
a(indvy,2) = qp * x(2) ! Omega_xy

! d vy / d (x, y, z) 
a(indvy,1) = a(indvx, 2) ! Omega_yx = Omega_xy
a(indvy,2) = 1.d0 + q + (rr1+rr2)*y2 ! Omega_yy


! Components related to z-components: d v_{x,y} / d {z}, and d vz / d {x, y, z}
if(ndim == 6) then 
  a(4,3) = qp * x(3)   ! Omega_xz
  a(5,3) = qpp * x(3)  ! Omega_yz
  
  a(6,1) = a(4,3) ! Omega_zx = Omega_xz
  a(6,2) = a(5,3) ! Omega_zy = Omega_yz
  a(6,3) = q + (rr1+rr2)*y3 ! Omega_zz
endif

!print*,'Jocabi matrix A:'
!do i = 1, ndim, 1
!  print*, a(i,:)
!end do
!print*; read*  !ck 

! d phi = A * phi ! DO not forget to multiply the Jocabi matrix A 
dphi = matmul(a,  phi)  

!print*, 'ndim = ', ndim, 'n=', n !ckd 

f(ndim+1 : n) = reshape( dphi, (/ndim * ndim/) )

! check phi  -- ck
!print*, 'Phi:'
!do i = 1, ndim, 1
!  print*, dphi(i,:)
!end do
!print*; read*  !ck 

 
return
end

