subroutine dflrtz_g3(x0,sgn, dlf)
! The general form of the differential of the lorentz force, as a function
! of parameter beta, only consider the last 3 row,(derirative wrt dvx,dvy,dvz)
! each one consists of two row, the first is the coefficient of beta, the second
! is constant coefficient 

!    Input
!  x0	the state of point of which to compute
!  sgn	the sign of q/m 

!	Output
!  dlf  the last three row of dlf, each row actually is split into 2 rows, 1-coefficient of beta, 2-constant part 

implicit none
integer, parameter :: dp = kind(1.d0) , n = 6

integer, intent(in) ::  sgn
real(kind=dp), intent(in) :: x0(n)
real(kind=dp), intent(out) :: dlf(n,n)

! local variables
real(kind=dp) :: x, x1, x2, dx, dx1, dx2, r2, r, r5, &
                 abt, acst, tmp, tmp0, tmp1, tmp2, &
                 b, cbt, ccst, dbt, dcst


x = x0(1)
x1 = x0(2)
x2 = x0(3)
dx = x0(4)
dx1 = x0(5)
dx2 = x0(6)

dlf = 0.d0 ! initialize as zero matrix

r2 = x*x + x1*x1 + x2*x2
r = dsqrt(r2)
r5 = r**5 

! fl_i = sgn*1/R⁵*3x_i*a
!a = beta*(x_i+2 * dx_i+1 - x_i+1 * dx_i+2) + x²_i+1 + x²_i+2 
abt = x2 * dx1 - x1 * dx2 ! coefficient of beta
acst = x1 * x1 + x2 * x2  ! constant part

! df_i/dx_i, df_i/dx_i+1, df_i/dx_i+2, df_i/d dx_i, df_i/d dx_i+1, df_i/d dx_i+2


! 1 column df_i/dx_i
! 1 row, coefficient of beta; 2 row, constant part
tmp = sgn * 3/r5 * (1 - 5*x*x /r2 ) ! independent on beta
dlf(1,1) =  abt * tmp  !beta
dlf(2,1) =  acst * tmp !cst


! 2 column df_i/dx_i+1
! 1 row, coefficient of beta; 2 row, constant part

tmp1  = sgn * 3 * x / r5
tmp2 = 5*x1 / r2 ! independent on beta

dlf(1,2) = tmp1 * ( -dx2 - tmp2 * abt ) !beta 
dlf(2,2) = tmp1 * ( 2*x1 - tmp2 * acst ) !cst

! 3 column df_i/dx_i+2
! 1 row, coefficient of beta; 2 row, constant part

!tmp1  = sgn * 3 * x / r5 ! the same as previous one
tmp2 = 5*x2 / r2 ! independent on beta

dlf(1,3) = tmp1 * ( dx1 - tmp2 * abt ) !beta 
dlf(2,3) = tmp1 * ( 2*x2 - tmp2 * acst ) !cst


! 4 column df_i/d dx_i
! 1 row, coefficient of beta; 2 row, constant part
! 0 
! 0 


! 5 column df_i/d dx_i+1
! 1 row, coefficient of beta; 2 row, constant part

dlf(1,5) = sgn * 3 * x * x2 / r5 !beta 
! 0 !cst


! 6 column df_i/d dx_i+2
! 1 row, coefficient of beta; 2 row, constant part

dlf(1,6) =  sgn * (-3) * x * x1 / r5 !beta 
! 0 !cst



! **************** fl_i+1 ******************************
! fl_i+1 = sgn*1/R⁵ * c
! b = x²_i+1 + x²_i+2 - 2*x²_i
! c = -beta * b * dx_i+2 - 3*beta* x_i+2 * x_i * dx_i + x_i+1 * b

b = x1 * x1 + x2 * x2 - 2 * x * x

 cbt = -b * dx2 - 3 * x2 * x * dx ! coefficient of beta
 ccst = x1 * b   ! constant part

! df_i+1/dx_i, df_i+1/dx_i+1, df_i+1/dx_i+2, df_i+1/d dx_i, df_i+1/d dx_i+1, df_i+1/d dx_i+2

tmp = sgn / r5 ! the commom first term


! 1 column df_i+1/dx_i
! 3 row, coefficient of beta; 4 row, constant part
tmp0 =  5*x / r2  ! independent on beta
dlf(3,1) = tmp *( 4 * x * dx2 - 3 * x2 * dx - tmp0 * cbt )  !beta
dlf(4,1) = tmp *( -4 * x * x1 - tmp0 * ccst ) !cst


! 2 column df_i+1/dx_i+1
! 3 row, coefficient of beta; 4 row, constant part
tmp1 =  5*x1 / r2  ! independent on beta
dlf(3,2) = tmp *( -2 * x1 * dx2  - tmp1 * cbt )  !beta
dlf(4,2) = tmp *( 3 * x1 * x1 + x2 * x2 - 2 * x * x - tmp1 * ccst ) !cst


! 3 column df_i+1/dx_i+1
! 3 row, coefficient of beta; 4 row, constant part
tmp2 =  5*x2 / r2  ! independent on beta
dlf(3,3) = tmp *( -2 * x2 * dx2 - 3 * x * dx  - tmp2 * cbt )  !beta
dlf(4,3) = tmp *( 2 * x1 * x2 - tmp2 * ccst ) !cst


! 4 column df_i+1 / d dx_i
! 3 row, coefficient of beta; 4 row, constant part
dlf(3,4) = sgn * (-3) * x2 * x / r5 !beta 
! 0 !cst

! 5 column df_i+1 / d dx_i+1
! 3 row, coefficient of beta; 4 row, constant part
! 0 
! 0 


! 6 column df_i+1 / d dx_i+2
! 3 row, coefficient of beta; 4 row, constant part

dlf(3,6) =  sgn * (- b) / r5 !beta 
! 0 !cst



! **************** fl_i+2 ******************************

! fl_i+2 = sgn*1/R^5 * d
! b = x²_i+1 + x²_i+2 - 2*x²_i
! d = beta * b * dx_i+1 + 3*beta* x_i+1 * x_i * dx_i + x_i+2 * b

!b = x1 * x1 + x2 * x2 - 2 * x * x ! same as the one in fl_î+1

 dbt = b * dx1 + 3 * x1 * x * dx ! coefficient of beta
 dcst = x2 * b   ! constant part

!tmp = sgn / r5 ! the commom first term, same with f_i+1

! df_i+2/dx_i, df_i+2/dx_i+1, df_i+2/dx_i+2, df_i+2/d dx_i, df_i+2/d dx_i+1, df_i+2/d dx_i+2


! 1 column df_i+2/dx_i
! 5 row, coefficient of beta; 6 row, constant part
!tmp0 =  5*x / r2  ! independent on beta
dlf(5,1) = tmp *( -4 * x * dx1 + 3 * x1 * dx - tmp0 * dbt )  !beta
dlf(6,1) = tmp *( -4 * x * x2 - tmp0 * dcst ) !cst

!print*, 'check dlf(6,1)', dlf(6,1)
!read(*,*)

! 2 column df_i+2/dx_i+1
! 5 row, coefficient of beta; 6 row, constant part
!tmp1 =  5*x1 / r2  ! independent on beta
dlf(5,2) = tmp *( 2 * x1 * dx1  +  3 * x * dx - tmp1 * dbt )  !beta
dlf(6,2) = tmp *( 2 * x1 * x2 - tmp1 * dcst ) !cst


! 3 column df_i+2/dx_i+1
! 5 row, coefficient of beta; 6 row, constant part
!tmp2 =  5*x2 / r2  ! same as in fl_i+1, independent on beta
dlf(5,3) = tmp *( 2 * x2 * dx1   - tmp2 * dbt )  !beta
dlf(6,3) = tmp *( 3 * x2 * x2 + x1 * x1 - 2 * x * x - tmp2 * dcst ) !cst



! 4 column df_i+2/d dx_i
! 5 row, coefficient of beta; 6 row, constant part
dlf(5,4) = sgn * 3 * x1 * x / r5 !beta 
! 0 !cst


! 5 column df_i+2/d dx_i+1
! 5 row, coefficient of beta; 6 row, constant part

dlf(5,5) =  sgn * b / r5 !beta 
! 0 !cst

! 6 column df_i+2/d dx_i+2
! 5 row, coefficient of beta; 6 row, constant part
! 0 
! 0 

 
return

end subroutine dflrtz_g3


























