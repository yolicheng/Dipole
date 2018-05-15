subroutine gr_lf(t, y, neq, f)
!  computation of the vector field +  variational equations to be integrated for lorentz force problem 
!  include this model subroutine into the respective module, so we can directly use the system-based paramters
!  without adding to the input parameter list,  keep the form deriv(t, y, n, f) to be called in gr_rk78

!  Input parameters:
!     x  	time
!     y(*)      point (y(1),y(2),....y(42))
!     neq 	neqnumber of equations (in this case 42), if it equals 6 the variational equation are skipped

!  Output parameters:
!     f(*)       vector field
!                in y(*) and f(*) the first 6 components correspond to the position and velocity
!                the others to the variational equations
!
!  lf_mod Variable 
!    cs, sgn, beta, dp, n !  inside the module, no need for use declaration

! ------------------Discard this option to delcare dummy argument and optional attributes 
!  as in gr_lf case, keep to first four of the same dummy argument, declare the rest with optional attributes

!  but then the subroutine containing optional arguments need to be explicit interface 
!  Arguments could be OPTIONAL, requiring the caller to somehow indicate that an argument was omitted. Even more fun, a function could return an array or a character value whose length was dependent on some function of the arguments. All of these things, and more, required the ability to describe in detail the procedure and its arguments so that the compiler could call the procedure the right way. This information also helped correctness as the compiler could now check to make sure you were passing the correct arguments. 
!  This set of information is called an "explicit interface".


! Subroutine use:  dflrtz, dflrtz_g3,  swap2gr,  inv_swap2gr
! Finaly revised by Yu -- 20160306
! -----------------------------------------------

implicit none 
 
! Input and Output
integer, intent(in) :: neq  ! dimension of the state
real(kind=dp), intent(in)  :: t, y(neq)  
real(kind=dp), intent(out) :: f(neq)

! Local Variable
integer ::  debug 
real(kind=dp) :: fg(6), dlf(6,6), yg(6), x,x1,x2, dx,dx1,dx2, &
                 a,abt,acst, b, c, cbt, ccst, d, dbt, dcst, & ! dlfrtz
                 r2, r, r5,  tmp, &  ! mat(6,6),
                 dphi(6,6), phi(6,6) ! variatioal matrix
                 
!equivalence (fcopy(7), dlf) ! obsolete usage, could be replaced by reshape..

! be careful, the array is stored in columnwise order 

debug = 1

! dx, dy, dz
f(1)=y(4)
f(2)=y(5)
f(3)=y(6)

! swap to general form
call swap2gr(y(1:6), yg)

fg(1:3) = yg(4:6) ! velocity in general form
 
x   = yg(1)
x1  = yg(2)
x2  = yg(3)
dx  = yg(4)
dx1 = yg(5)
dx2 = yg(6)


r2 = x*x + x1*x1 + x2*x2
r = dsqrt(r2)
r5 = r**5 
tmp = sgn / r5! the commom first term


! fl_i = sgn*1/R^ 5 * 3x_i*a
!a = beta*(x_i+2 * dx_i+1 - x_i+1 * dx_i+2) + x_i+1^ 2 + x_i+2^ 2
 
abt = x2 * dx1 - x1 * dx2 ! coefficient of beta
acst = x1 * x1 + x2 * x2  ! constant part
a = abt * beta + acst
fg(4) = tmp * 3 * x * a

!print*, 'check ax, a,abt,acst, x, r,r5, tmp', fg(4), a,abt,acst, x, r,r5, tmp
! **************** fl_i+1 ******************************
! fl_i+1 = sgn*1/R^ 5 * c
! b = x_i+1^ 2 + x_i+2^ 2 - 2*x_i^ 2
! c = -beta * b * dx_i+2 - 3*beta* x_i+2 * x_i * dx_i + x_i+1 * b

b = x1 * x1 + x2 * x2 - 2 * x * x

 cbt = -b * dx2 - 3 * x2 * x * dx ! coefficient of beta
 ccst = x1 * b   ! constant part
 c = cbt * beta + ccst
 
 fg(5) = tmp * c
 
 
 ! **************** fl_i+2 ******************************

! fl_i+2 = sgn*1/R^ 5 * d
! b = x_i+1^ 2 + x_i+2 ^ 2 - 2*x_i^ 2
! d = beta * b * dx_i+1 + 3*beta* x_i+1 * x_i * dx_i + x_i+2 * b

!b = x1 * x1 + x2 * x2 - 2 * x * x ! same as the one in fl_i+1

dbt = b * dx1 + 3 * x1 * x * dx ! coefficient of beta
dcst = x2 * b   ! constant part
d = dbt * beta + dcst
 
fg(6) = tmp * d 
 
!print*, 'test inv_swap2gr'
!print*, 'fg, general ', fg(1:6)

call inv_swap2gr(fg, f(1:6))

!print*, 'f, real ', f(1:6)


! the previous just got the lorentz force, follows is the acceleration
f(4) = f(4) + 2 * y(5) + 3 * y(1) ! fx + 2vy + 3x
f(5) = f(5) - 2 * y(4) ! fy - 2vx
f(6) = f(6) - y(3) ! fz - z


!print*, 'initial state', y(1:6)
!print*, 'f, real ', f(1:6)
!print*

!read(*,*) aaa
 
if(neq.eq.6)  return
! ------------------- Position + Velocity -------------------------------- 


! ----------------Variational Matrix --------------------------------------
! this part need to be check really carefully!
! how stupid can u be!!! - 20160317 -- big bug here! how can you make such a stupid mistake.

!              \dot phi = Dlf * phi
  
! so here, we need to do a matrix multiplication.
phi = reshape( y(7:42), (/6, 6/) )

! To compute Jacobi Matrix of the vector field, only y(1:n) is used, where n is the dimension of state vector, specified in lf_mod 
call dflrtz(y(1:n), dlf) 

dphi = matmul(dlf, phi)

f(7:42) =  reshape(dphi, (/36/)) ! checked, discard the use of equivalence

!if (debug == 1) then ! checked f(7:42), fine! ! ckd 2017-02-20 10:49:51 
!  print*, 'check the parameter of the system, beta, cs, sgn'
!  print*, beta, cs, sgn
!  
!  
!  print*, 'check the variational matrix'
!  do i = 1, 6
!    write(*,'(6e18.8)') dlf(i,:) 
!  enddo
!  print*
!  
!  print*, 'check the STM'
!  do i = 1, 6
!    write(*,'(6e18.8)') phi(i,:) 
!  enddo
!  print*
!  
!!print*, 'check last 36 components of f'
!  do i = 1, 6
!    write(*,'(6e18.8)') f(6*i+1: 6*i+6) 
!  enddo

!  print*, 'check the reshape of the variational matrix' !ckd
!  mat =  reshape(f(7:42), (/6,6/))
!  do i = 1, 6
!    write(*,'(6e18.8)') mat(i,:) 
!  enddo
!  read*
!endif 


return
end

