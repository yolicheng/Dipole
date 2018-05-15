subroutine deriv(t,y,neq,f)
!     computation of the vector field +  variational equations to be integrated
!     depend on the parameter beta, for lorentz force problem

!     gr_rtbp: (rtbp if mu.ne.0,and hill's problem if mu.eq.0  ) 
! 
!  input parameters:
!     x  time
!     y(*)       point (y(1),y(2),....y(42))
!     neqnumber of equations (in this case 42), if it equals 6 the variational equati
! care skipped

!  output parameters:
!     f(*)       vector field
!     in y(*) and f(*) the first 6 components correspond to the position and velocity
!      the others to the variational equations
!
! since this vector field computation is supposed to be called by gr_rk78,
! in which the deriv has 4 arguments, 3 inputs and 1 output

!  as in gr_lf case, keep to first four of the same dummy argument, 
!  declare the rest with optional attributes

!  but then the subroutine containing optional arguments need to be explicit interface 
!  Arguments could be OPTIONAL, requiring the caller to somehow indicate that an argument was omitted. Even more fun, a function could return an array or a character value whose length was dependent on some function of the arguments. All of these things, and more, required the ability to describe in detail the procedure and its arguments so that the compiler could call the procedure the right way. This information also helped correctness as the compiler could now check to make sure you were passing the correct arguments. 
!  This set of information is called an "explicit interface".

! the best way is to put the subroutine in a module, which needs more effort, give up at this moment 
! because gr_lf, and gr_rtbp are never called by the same problem together, a simpler option is :

! Finaly solution, 20151223 
! use module lf_mod to pass the variables: beta, cs, sgn, so the subroutine deriv to compute the vector field 
! in gr_rk78 need not to be modified 

!    just modify the call of deriv in gr_rk78, to make it suitable for gr_lf

! Subroutine used: swap2gr, inv_swap2gr, dflrtz
!-------------------------------------------------------------

use lf_mod, only : cs, sgn, beta
implicit none 

integer, parameter:: dp = kind(1.d0)

integer, intent(in) :: neq !, cs, sgn
real(kind=dp), intent(in)  :: t, y(neq) !, beta
real(kind=dp), intent(out) :: f(42)

integer :: i
real(kind=dp) :: fg(6), dlf(6,6), yg(6), x,x1,x2, dx,dx1,dx2, a,abt,acst, &
     		 b, c, cbt, ccst, d, dbt, dcst, &
    		 r2, r, r5, fcopy(42), mat(6,6), tmp, aaa
    		 
equivalence (fcopy(7), dlf)

! be careful, the array is stored in columnwise order 

print*, 't,y', t, y(1:6)

 print*, 'In deriv- tmax, hmin,hmax,e,ind,sec',tmax, hmin,hmax,e, ind,sec
      print*, 'tar1, tar2, ctr1, ctr2, ctr3 ', tar1, tar2, ctr1, ctr2, ctr3 
      read(*,*)

! dx, dy, dz
f(1)=y(4)
f(2)=y(5)
f(3)=y(6)

! swap to general form
call swap2gr(y(1:6), cs, yg)

fg(1:3) = yg(4:6) ! velocity in general form
 
 
x = yg(1)
x1 = yg(2)
x2 = yg(3)
dx = yg(4)
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

call inv_swap2gr(fg, cs, f(1:6))

!print*, 'f, real ', f(1:6)


! the previous just got the lorentz force, follows is the acceleration
f(4) = f(4) + 2 * y(5) + 3 * y(1) ! fx + 2vy + 3x
f(5) = f(5) - 2 * y(4) ! fy - 2vx
f(6) = f(6) - y(3) ! fz - z


!print*, 'initial state', y(1:6)
!print*, 'f, real ', f(1:6)
!print*

!read(*,*) aaa
 
if(neq.eq.6) return

! the variational matrix, call subroutine dflrtz 
call dflrtz( y, beta, cs, sgn, dlf)

!print*, 'check the parameter of the system, beta, cs, sgn'
!print*, beta, cs, sgn
!print*, 'check the variatioal matrix'
!do i = 1, 6
!  write(*,'(6f8.4)') dlf(i,:) 
!enddo



!print*, 'check the reshape of the variatioal matrix'
!mat =  reshape(fcopy(7:42), (/6,6/))
!do i = 1, 6
!  write(*,'(6f8.4)') mat(i,:) 
!enddo

!read(*,*) aaa

f(7:42) = fcopy(7:42) 

!print*, 'check last 36 components of f'
!do i = 1, 6
!  write(*,'(6f8.4)') f(6*i+1: 6*i+6) 
!enddo

return
end

