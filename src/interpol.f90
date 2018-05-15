!***********************************************************************
!     ****   interpol   ****
! Linear interpolation for points of dimension np-by-ndim, according to 
! the parameter arg, in order to obtain npnew points equally spaced by arg within [0, 2pi] 

!   with the interpolated values for arg:  arg = k*(prd/N),   k=0..N-1, prd = 2pi

! **NOTE** the preprocess of two more points only fit this problem... 
!          if there is possible use for more general linear interpolation, we have to modify this routine

! since the angles for the first and the last point are not exactly 0 and 2pi, to avoid extrapolation, we add two more angles 
! one is added before the first point as  arg(npoinc) - 2*pi,  
! the other is added after the last point as arg(1)   + 2*pi 
! Such that we only deal with the interpolation 

! **NOTE** if pt has only one column, we have to preprocess the data into a two dimensional array with ndim=1
!       Input Variables 
!  pt         dimension np-by-ndim, number of original points 
!  np,ndim    the dimension of the state of the point, pt 
!  arg        dimension np, the standard according which to do the interpolation 
!  npnew      number of points to interpolate
!  fptnew     file tag to save all the new interpolated points.
  

!       Output Variables 
!  ptnew      dimension npnew-by-dimension, the interpolated points        

!  Routine Used: None

!  Finally Revised by Yu -- 20160811
!***********************************************************************
subroutine interpol( pt, np, ndim, arg, npnew, ptnew)
use dp_mod
use pi_mod 
implicit none

! Input  and Output Declaration   
integer, intent(in)      ::  np, ndim, npnew !, fptnew 
real(kind=dp), intent(in)      :: pt(np, ndim), arg(np)  
real(kind=dp), intent(out)     :: ptnew(npnew, ndim)
 
! Local Variable
integer :: i,  ipt 
real(kind=dp)  :: ptaux(np+2, ndim), argaux(np+2), arg_temp, darg, delt_arg 

!real(kind=dp)  :: pt222(ndim) !ck 
  
  ! We assume arg(1) >= 0.d0, and arg(np) <= 2pi 
  ! add two more points to make sure the values of argument cover the full interval [0, 2pi]
  
  ! -- deal with the first point, if arg(1) > 0. to obtain arg=0,  we need  a point less than 0, 
  !    take the last point, and substract the angle by pi2 
  if( arg(1) > 0.d0 ) then 
    argaux(1)  = arg(np) - pi2
    ptaux(1,:) = pt(np,:)
  else 
    argaux(1) = arg(1) 
    ptaux(1,:) = pt(1,:) 
  endif 
  
  argaux(2:np+1)   = arg 
  ptaux(2:np+1, :) = pt 
  
  ! deal with the last point, if it equals 2pi, assign as the last point,
  ! otherwise take the first point with the argument be the value after 1 period 
   if( arg(np) < pi2 ) then 
    argaux(np+2)  = arg(1) + pi2
    ptaux(np+2,:) = pt(1,:)
  
  else 
    argaux(np+2)  = arg(np) 
    ptaux(np+2,:) = pt(np,:) 
  endif 
  
  ! use linear interpolation, so find the two closest values of arg for the interpolated one 
  ipt = 1 ! initialize the index of the first point to use  for interpolation 
  
  darg = pi2 / npnew
  arg_temp = -darg  
  
  do i = 1, npnew, 1
!    arg = (i-1)*(prd/N),     k = 1_N, prd = 2pi
    arg_temp = arg_temp + darg  ! just a small trick to replace the multiplication by addition 
    
    do while (argaux(ipt+1) < arg_temp )
      ipt = ipt + 1
    end do 
    
    ! check the valus of arg(ipt) and the target value arg_temp
!    print*, 'arg_temp = ', arg_temp, 'ipt=', ipt, 'argaux(ipt:ipt+1) = ', argaux(ipt:ipt+1)
!    read*

    ! use ptaux(ipt:ipt+1) and argaux(ipt:ipt+1) for interpolation 
    ! ptnew(i)_{arg_temp} =   ptaux(ipt) + ( ptaux(ipt+1) -  ptaux(ipt) ) /  ( argaux(ipt+1) -  argaux(ipt) )  *  ( arg_temp -  argaux(ipt) )
    
    if( dabs(argaux(ipt+1) - arg_temp ) < 1.d-6 )  then 
      ptnew(i, :) = ptaux(ipt+1, :) 
    else   
      delt_arg = ( arg_temp -  argaux(ipt) ) / ( argaux(ipt+1) -  argaux(ipt) ) ! the commom part for all the components
      ptnew(i, :)  =  ptaux(ipt, :) + ( ptaux(ipt+1, :) -  ptaux(ipt, :) ) * delt_arg 
    endif 
    
!    print*, 'Orginal two points + interpolated point  - arg - pv(6)' !ckd
!    write(*, '(7f20.14)')  argaux(ipt),   ptaux(ipt, :) 
!    write(*, '(7f20.14)')  argaux(ipt+1), ptaux(ipt+1, :) 
!    write(*, '(7f20.14)')  arg_temp, ptnew(i,:)
!    read* 
!    
    ! -- save all the points to file ob_interpol.dat, better to do this in the routine interpol  
    ! TODO--- check the energy of the new points, then we have to specify everything from the very begining.....
!    call 
!    write(fptnew, '(7f20.14)') arg_temp, ptnew(i, :)
    
  end do   
   

  return  
end subroutine interpol

