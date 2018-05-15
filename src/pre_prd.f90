! Predict along the a single parameter:  period (isarc=2) or energy(isarc=3)

!  We observe that the curve of the control variables ie, (x, z, vy, tp) w.r.t energy or period are almost linear 
!  and for the asymmetric periodic orbit, it is hard to do continuation along the characteristic curve,
!  in this case we swith to the one-parameter continuation, there is no much difference for period or energy...  

!  ********************* Idea *******************************
!  Add one more constraint function :  
!     T = T0 for along the period ;   CJ = CJ0 for along the energy
!  Now we have the same number of unknowns and equations, but the Jacobi Matrix is 

!  because the control variables is not explicitly dependent on the period 
!  use the difference bewteen the previous two p.o.s as the derivative
!  --- disard this approach : for the second p.o.s, we have to use arc-length continuation(a very small step if necessary), 
!      and then use the finite difference to  compute the prediction for the new p.o. 

! 	Input  Varaible 
!  curv 	the control variables in the initial conditions of the previous p.o.s, for general case, the full state vector 
!           	for generality, use as curv(5,nctr+1), the last column is the parameter to follow 

!  nds		if incrds = -1, the number of current available rows in curv 
!		else, it has nds+1 rows .... 

!  incrds 	flag to decide how to change the current arc step, 1: 2*ds, 0: remain, -1: ds/2


! 	Input  and Output Varaible 
!  yctr		the initial (and the updated) guess  of the prediction for the next p.o. 
! 

! ******************************************************************************

subroutine pre_prd( ipo, curv, vf, g, nds, ds, dir, incrds, csangle, yctr, con_stop)

implicit none 
integer, parameter :: dp = kind(1.d0)

integer, intent(out)   :: con_stop
integer, intent(inout) :: ipo, nds, dir, incrds

real(kind=dp), intent(in)    ::  g(ntar, nctr)
real(kind=dp), intent(out)    :: yctr(nctr) 
real(kind=dp), intent(inout)  :: vf(nctr),ds, csangle, curv(5, nctr+1)
  
! Local Variables 
integer :: i, nrow 
real(kind=dp) :: x1(nctr+1), x2(nctr+1), dx(nctr), dtp, der(nctr),  &
                 curv2ds(3, nctr), csangle2ds 

print*, 'Compute the prediction of new p.o. along the period!'
    
!initialize the flags for control of stepsize to be zeros
con_stop = 0 ! too small stepsize required
 

if (incrds == -1 ) then     
  print* ,'Need decrease arc step!'; read*
  ipo = ipo-1
! the case when we need to decrease the stepsize 
  if (ds/2.d0 .le. ds_min) then 
    print*, 'Cannot continue with ds_min! Terminate the continuation process!'; read*
    con_stop = 1
    return
  else               
    ds =  ds / 2.d0
  endif      
  
  
else 

  nds = nds + 1 ! update the counter

! before doubling the step, check the possible new 3 points with stepsize of 2*ds
  if ( incrds == 1) then 
    print*;  print* ,'Need increase arc step!'; read*
    print*, 'nds >=5 to double ds! nds=',nds ; read* !ck
    
    ! I think there is no point in constraint the angle in these two approach, because it is not along the characteristic curve.
    ! anyway, put it here, for possible use in the futurn 
         
!    if( anglectr == 1) then 
    curv2ds(1:3, :) = curv(1:5:2, 1 : nctr)
    call curv_angle(curv2ds, nctr, csangle2ds)
    
    print*, 'check angle for possible doubling of ds, csangle=', csangle2ds
    read* 
    
    if ( dabs(csangle2ds) < csangle_min .and. anglectr == 1)  then ! if no control in angle, double the step
      incrds = 0 
      
    else 
      print*, 'Double ds!!'; read*
      ds = ds * 2.d0
      csangle = csangle2ds
    endif 
  
  endif  ! incrds == 1
  
endif  ! incrds = -1 

! compute the finite difference between the last two available rows in curv 
nrow  = min0(nds, 5)

x1 = curv(nrow-1, :)
x2 = curv(nrow, :)

print*, 'check the last 2 points'
print*, x1
print*, x2
print*;  read*

! instead of the tangent vector along the characteristic curve, try do continuation along the energy or period 

! the difference in control variables between two consecutive points 
dx =  x2(1:nctr) - x1(1:nctr) 

! the difference  in the continuation parameter 
dtp = x2(1+nctr) - x1(1+nctr)

! check difference 
print*, 'dx', dx 
print*, 'dtp =', dtp 

der = dx / dtp 

der = der / dnrm2(nctr, der, 1)

! compare the difference and the derivative 
print*, 'By difference: der=', der  
read*

if( isarc == 2) then 
  print*, 'By vector field: vf=', vf/dnrm2(nctr, vf, 1) 
  read*
endif 


print*, 'dir=', dir, 'ds=', ds  
read* 

! check 
print*, 'check if yctr = the last available row in curv'
print*, yctr
print*, curv(nrow,:)
read* 


  
if (nds > 2 .and. debug == 1) then 
  print*;  print*, 'Check the angle after updating cham!'
  print*, 'cos<v1,v2> =', csangle, 'only make sense for nds>2, nds=', nds
  print*; read* 
endif

yctr = yctr + dir * ds * der 

print*, 'check updated yctr'
print*, yctr 
read*
    
return
end subroutine pre_prd


