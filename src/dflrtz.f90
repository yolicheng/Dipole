subroutine dflrtz(x0, dlf)
! The differential of the lorentz force, as a function of parameter beta
! 
! the general form of the differential of the lorentz force of only the last 3 row,(derivative wrt dvx,dvy,dvz) is 
! computed in subroutine dflrtz.f90

! 2016-08-01 21:21:59  update the index of case by 1:normal-2:radial-3:tangential 

! 201603027 --checked for all the components by finite difference, can believe this routine  

! Cautions!
! to use the general form, we have to deal with the order of the compoent of x0, and re-swith after the call of dflrtz 

! rule of     x_i, x_i+1, x_i+2
! x -> y  -> z  -> x -> y  
! case 1: N =  [0 0 1]  z  -> x  -> y  
! case 2: N =  [1 0 0]  x  -> y  -> z  !! the original case
! case 3: N =  [0 1 0]  y  -> z  -> x 
! 

!    Input
!  x0	the state of point of which to compute
!  sgn	the sign of q/m 

!	Output
!  dlf  the last three row of dlf, each row actually is split into 2 rows, 1-coefficient of beta, 2-constant part 

!	Modular public variables 
! cs, sgn, beta

! function used :   dflrtz_g3,  swap2gr, inv_swap2gr
! Finaly revised by Yu -- 20160307
!-------------------------------------------------

implicit none

real(kind=dp), intent(in) :: x0(n)  
real(kind=dp), intent(out) :: dlf(n,n)

! local variables
integer :: i, debug 
real(kind=dp) :: dlf3(n,n), dlfsub(3,6), xg(n), dlfsubcopy(3,6) 
                 
  debug = 0 
    
  if(debug == 1) then 
    print*, 'Debug mode for dflrtz.f90'
    print*; read*
  endif 
  
  ! initialize dlf
  dlf = 0.d0 
  dlf(1,4) = 1.d0
  dlf(2,5) = 1.d0
  dlf(3,6) = 1.d0

  ! swap the component of x0, to produce the appropriate general order to pass to dflrtz_g3
  ! x,y,z,dx,dy,dz -> y,z,x,dy,dz,dx
  if (debug == 1) then
    print*, 'dflrtz: original state'
    print*, 'case = ', cs,  'beta=', beta
    write(*,'(6f12.8)') x0
    read*
  endif 


! ***********  swap the state to general form ************
  call swap2gr(x0, xg)

  if (debug == 1) then
    print*, 'General form'
    write(*,'(6f12.8)') xg
    read* 
  endif

  ! compute the last three rows of differential of the vector field
  call dflrtz_g3(xg, dlf3) 

  ! this dlf3 is the one that seperated by 2 rows, one is the coefficient of beta, the other is the constant part.
  if (debug == 1) then
    print*, 'the last three rows of dfl, 3-by-6'
    do i = 1, 6
      write(*,'(6f12.8)') dlf3(i,:) 
    enddo
    print*
    read(*,*)
    print*, beta; read*
  endif 

  dlfsub(1,:) = dlf3(1,:) * beta + dlf3(2,:)
  dlfsub(2,:) = dlf3(3,:) * beta + dlf3(4,:)
  dlfsub(3,:) = dlf3(5,:) * beta + dlf3(6,:)

  dlfsubcopy = dlfsub

  if (debug == 1) then
    print*, 'original dlfsub, 3-by-6, beta', beta
    do i = 1, 3
     write(*,'(6f12.8)') dlfsubcopy(i,:) 
    enddo
    print*
  endif


! ***********  swap the jacobi matrix back to real form (x,y,z)************ 
!   -- tangential case 3, N = [0 1 0]
  if(cs == 3) then ! (y,z,x) --> (x, y, z)
  ! row switching transformation
  ! (R1, R2, R3) <- (R3, R1, R2)

  dlfsub = dlfsub((/3,1,2/),:)

  if (debug == 1) then 
    print*, 'Direct row swap:  (R1, R2, R3) <- (R3, R1, R2)'
    do i = 1, 3
      write(*,'(6f8.4)') dlfsub(i,:) 
    enddo
    print*
  endif

! reverse of column switching transformation 
! (C1, C2, C3, C4, C5, C6) <- (C3, C1, C2, C6, C4, C5)

! ---- ckd! do not need an intermediate variable, use array to address the index of array

  dlfsub = dlfsub(:, (/3,1,2,6,4,5/)) ! this seems the best way to do the swap

  if (debug == 1) then
    print*, ' (C1, C2, C3, C4, C5, C6) <- (C3, C1, C2, C6, C4, C5)'
    do i = 1, 3
     write(*,'(6f8.4)') dlfsub(i,:) 
    enddo
    print*
    endif

! -- Normal case 1, N = [ 0 0 1]
  elseif( cs == 1) then !
    ! (R1, R2, R3) <- (R2, R3, R1)
    dlfsub = dlfsub((/2,3,1/),:)

    if (debug == 1) then 
      print*, ' (R1, R2, R3) <- (R2, R3, R1)'
      do i = 1, 3
        write(*,'(6f8.4)') dlfsub(i,:) 
      enddo
      print*
    endif
 
! column switching transformation
!(C1, C2, C3, C4, C5, C6) <- (C2, C3, C1, C5, C6, C4).
    dlfsub = dlfsub(:, (/2,3,1,5,6,4/)) 

    if (debug == 1) then
      print*, '(C1, C2, C3, C4, C5, C6) <- (C2, C3, C1, C5, C6, C4).'
      do i = 1, 3
       write(*,'(6f8.4)') dlfsub(i,:) 
      enddo
      print*
    endif

  endif


  dlf((/4,5,6/), :) = dlfsub

! with the normal order, x -> y  -> z
! add the constant matrix 
!
  dlf(4,1) = dlf(4,1) + 3.d0
  dlf(4,5) = dlf(4,5) + 2.d0
  dlf(5,4) = dlf(5,4) - 2.d0
  dlf(6,3) = dlf(6,3) - 1.d0

  return 
end subroutine dflrtz



subroutine dflrtz_g3(x0, dlf)
! The general form of the differential of the lorentz force, as a function
! of parameter beta, only consider the last 3 row,(derirative wrt dvx,dvy,dvz)
! each one consists of two row, the first is the coefficient of beta, the second
! is constant coefficient 

!	Input
!  x0	the state of point of which to compute
!!  sgn	the sign of q/m 

!	Output
!  dlf  the last three row of dlf, each row actually is split into 2 rows, 1-coefficient of beta, 2-constant part 

! 	Module Global Variable 
! dp, n, sgn   variable precision, and the demension of the system

! Subroutine use: None
! Finaly revised by Yu -- 20160306
! -----------------------------------------------

implicit none

!integer, intent(in) ::  sgn
real(kind=dp), intent(in) :: x0(n)
real(kind=dp), intent(out) :: dlf(n,n)

! local variables
real(kind=dp) :: x, x1, x2, dx, dx1, dx2, r2, r, r5, &
                 abt, acst, tmp, tmp0, tmp1, tmp2, &
                 b, cbt, ccst, dbt, dcst

! check the variatioal matrix 
  if (debug == 1) then 
    print*, 'dflrtz_g3 -- check the variatioal matrix,  dp, n, sgn'
    print*,  dp, n, sgn
    print*, 'x0', x0
    read*
  endif 

  ! position
  x = x0(1)
  x1 = x0(2)
  x2 = x0(3)

  ! velocity
  dx = x0(4)
  dx1 = x0(5)
  dx2 = x0(6)

  dlf = 0.d0 ! initialize as zero matrix

  r2 = x*x + x1*x1 + x2*x2
  r = dsqrt(r2)
  r5 = r**5 

  ! fl_i = sgn * 1/ R⁵ * 3x_i * a
  !a = beta * (x_i+2 * dx_i+1 - x_i+1 * dx_i+2) + x²_i+1 + x²_i+2 
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

  dlf(3,6) =  - sgn * b / r5 !beta 
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


  ! 4 column df_i+2/d dx_i  -- dvx  -- to check
  ! 5 row, coefficient of beta; 6 row, constant part
  dlf(5,4) = sgn * 3 * x1 * x / r5 !beta 
  ! 0 !cst


  ! 5 column df_i+2/d dx_i+1  -- dvy -- to check   
  ! 5 row, coefficient of beta; 6 row, constant part

  !b = x1 * x1 + x2 * x2 - 2 * x * x 

  dlf(5,5) =  sgn * b / r5 !beta 
  ! 0 !cst

  ! 6 column df_i+2/d dx_i+2
  ! 5 row, coefficient of beta; 6 row, constant part
  ! 0 
  ! 0 

 
  return

end subroutine dflrtz_g3

subroutine swap2gr(x0,  xg)
! ***********  swap the state to general form ************

! rule of     x_i, x_i+1, x_i+2
! x -> y  -> z  -> x -> y   
! case 1: N =  [0 0 1]  z  -> x  -> y 
! case 2: N =  [1 0 0]  x  -> y  -> z  !! this is the initial one 
! case 3: N =  [0 1 0]  y  -> z  -> x 


!	Input
!  x0	the state of point of which to compute
!  sgn	the sign of q/m 

! 	Output 
!  xg   the general form of the state

! 	Module Global Variable 
! dp, n, cs   variable precision, and the demension of the system

! Subroutine use: None
! Finaly revised by Yu -- 20160801 -- update the index of case 
! -----------------------------------------------

implicit none

real(kind=dp), intent(in) ::  x0(n)
real(kind=dp), intent(out) :: xg(n)

  if(cs == 2) then 
    ! -- radial case 

    ! copy the initial one
    xg = x0 
    !  write(*,'(6f8.4)') xg
  
  elseif(cs == 3) then 
    ! -- tangential case

    !(C1, C2, C3, C4, C5, C6) → (C2, C3, C1, C5, C6, C4).
    !  print*, '(C1, C2, C3, C4, C5, C6) → (C2, C3, C1, C5, C6, C4)'
    xg = x0((/2,3,1,5,6,4/) )

  elseif( cs == 1) then 
    ! -- normal case 

    ! (C1, C2, C3, C4, C5, C6) → (C3, C1, C2, C6, C4, C5).
    !  print*, '(C1, C2, C3, C4, C5, C6) → (C3, C1, C2, C6, C4, C5)'
    xg = x0((/3,1,2,6,4,5/) )

  endif


end subroutine swap2gr


subroutine inv_swap2gr(x0,  xg)
! ***********  the inverse of the swap of the state to general form ************
! rule of     x_i, x_i+1, x_i+2
! x -> y  -> z  -> x -> y 
! case 1: N =  [0 0 1]  z  -> x  -> y   
! case 2: N =  [1 0 0]  x  -> y  -> z  !! the original case
! case 3: N =  [0 1 0]  y  -> z  -> x 


! 	Input
!  x0	the state of point of which to compute
!  sgn	the sign of q/m 

! 	Output 
!  xg   the general form of the state

! 	Module Global Variable 
! dp, n, cs   variable precision, and the demension of the system

! Subroutine use: None
! Finaly revised by Yu -- 20160801  -- update the index of case 
! -----------------------------------------------

implicit none

real(kind=dp), intent(in) ::  x0(n)
real(kind=dp), intent(out) :: xg(n)


  if(cs == 2) then 
    ! -- radial case
    ! copy the initial one
   
    xg = x0 
    !  write(*,'(6f8.4)') xg !ckd!
  
  elseif(cs == 3) then 
    ! -- tangential case 

    !(C1, C2, C3, C4, C5, C6) → (C2, C3, C1, C5, C6, C4).
    !  print*, '(C1, C2, C3, C4, C5, C6) → (C2, C3, C1, C5, C6, C4)'
    xg((/2,3,1,5,6,4/) ) = x0

  elseif( cs == 1) then 
    ! -- normal case 

    ! (C1, C2, C3, C4, C5, C6) → (C3, C1, C2, C6, C4, C5).
    !  print*, '(C1, C2, C3, C4, C5, C6) → (C3, C1, C2, C6, C4, C5)'
    xg((/3,1,2,6,4,5/) ) = x0
    
  endif


end subroutine inv_swap2gr










