! this subroutine control the angle between 3 consecutive points , to produce a smooth characteristic curve 

! use interpolartion to obtain an intermediate point, the nds is updated to be 2, with the previous point + new interpolated point 

! 20160420  -- disard the interpolation approach, 
! Instead use the option that decrease the stepsize, and starts over, this the mose efficient way! 
    
    
! 20160419 
!  set debug==0, to  check this routine ..., after the check, set back to 1


! 	Input  Varaible
!  nds 		number of available vectors 
!  a 		the current vector, updated as nds-th row depending on the its angle with the previous vector cham(min0(nds,4)-1,:)
! 

! 	Input and Output Varaible
!  cham 	VECTOR FIELD ON THE CHARACTERISTIC CURVE AT THE LAST (up to 4)
!		POINTS
!  curv 	coordinate of the last (up to 3) points on the characteristic curve 


! 	 Output Varaible
!  con_stop  	flag to show if the continuation process needs to be terminated due to the small stepsize required (<ds_min)

! 	Module-based Varaible
!  anglectr, ntar, nctr, csangle_min

! Routine used: none
! Finally revised by Yu -20160419 
! ---------------------------------------------------------------------------------


subroutine angle_interpol(nds, ds, a,  csangle, cham, curv, con_stop ) 

implicit none 
integer, parameter :: dp = kind(1.d0)

! 	Input and Output Decalaration 
integer, intent(inout)    :: nds
integer, intent(out)      :: con_stop

real(kind=dp), intent(in)    :: a(nctr)  
real(kind=dp), intent(inout) :: ds, csangle, cham(4, nctr), curv(3, nctr) 
  
  
! 	Local Variables 
integer :: i, j  
  
real(kind=dp) ::   angle, v1(nctr), v2(nctr), x1(nctr), x2(nctr), x3(nctr),  vm(nctr), xm(nctr) ! 


! set to 0 by default  
con_stop = 0 

! nds >= 2, at least 2 available row in cham, and 3 computed p.o.s

x1 = curv(1, : )
x2 = curv(2, : )
x3 = curv(3, : )

v1 = x2 - x1 
v2 = x3 - x2

! normalize 
v1 = v1 / dnrm2(nctr, v1, 1)
v2 = v2 / dnrm2(nctr, v2, 1)
! check if csangle is the one  passed in 
print*, 'cos<v1,v2>=', csangle 
 csangle = dot_product( v1, v2 ) ! v1, v2 are both unit vector, so do not need to devide the norm 
print*, 'check if the value remains? cos<v1,v2>=', csangle 
read* 


! intermediate point need to be added to make sure the characteristic curve is smooth 

if (debug == 0) then !ck --todo  debug==1 
  
! the start, check the value of nds, should be >=2, but the available rows in curv should be 3 
  print*, 'check  before angle_interpol: nds>=2, and curv(with 3 available rows)!'
  print*, 'nds =', nds 
   
  print*, 'curv'
  do i = 1, 3
    print*, curv(i,:)
  enddo  
  print* 
  
  print*, 'cham'
  do i = 1, min0(nds-1, 4)
    print*, cham(i,:)
  enddo 
    
      
  print*, 'v1, v2'
  print*, v1
  print*, v2
    
  print*, 'check the norm = 1?', dnrm2(nctr, v1, 1), dnrm2(nctr, v2, 1)
  read* ! if == 1, discard this part 
endif 
! ------------------------- ck --------------------------------------

! check if csangle is the one  passed in 
print*, 'cos<v1,v2>=', csangle 
 csangle = dot_product( v1, v2 ) ! v1, v2 are both unit vector, so do not need to devide the norm 
print*, 'check if the value remains? cos<v1,v2>=', csangle 
read* 

! ---------------   use  interpolartion to obtain ds ----------------      
! it seems not a good option, try again.....

do 
  print*, 'csangle= ', csangle ;    read*
      
! it is possible that   c < 0, and c < -0.995, in this case, we need to reverse, otherwise, vm = (v1+v2)/2 will give an almost rectangle one.
  if (csangle < 0.d0) then
    print*, 'cos<v1,v2> < 0, donot know how to handle this!!!' !!
    
    ! 20160420  -- stop here, it is hard to use interpolation to obtain a good guess... because if <v1,v2> < 90 degree, then we can not guarantee the linear interpolartion locate in the middle between x2 and x3! so this option is quite inefficient .... 
    ! disard the interpolation approach, and decrease the stepsize, and start over, this the mose efficient way! 
    
    read* 
    
    v2 = -v2  
    print*,'Reverse v2!'
    read*
  endif 
	 
! just simple linear interpolation, and the option of  (vm = xm - x1 ) is a bad choice
  xm = ( x1 + x2) / 2.d0

! evaluate the angle between vm and v1, to see if we need to decrease the arc step 
  vm = ( v1 + v2) / 2.d0  ! simple linear interpolation 
  vm = vm/ dnrm2(nctr, vm, 1)  ! normalize 
        
  angle = dot_product( v1, vm) /  dnrm2(nctr, v1, 1) / dnrm2(nctr, vm, 1)

  x2 = xm 
  v2 = vm 
        
! the only case we need to modify the modulus of ds, otherwise, keep the modulus of ds 
  if ( dabs( angle ) .ge. csangle_min ) exit 
       
enddo 
      

! a big bug here-- to do , stuck at 27-th p.o. and get a guess for a smooth characteristic curve       
! Consider the available bound for arc step,  if ds < ds_min, fail of contuation.... 
! From Carles' book: The successive reduction of the arc step can laed to a too small step, showing some problem in the computation. Then one should stop and look for a detailed analysis of what happens.  Maybe it's time to terminate the continuation
     
ds = dnrm2(nctr, x2 - x1, 1)
      
if (ds < ds_min) then 
    ! stop the continuation, not possible to go further. Return with flag con_stop = 1, and deal it in the callee.
  con_stop = 1 
  print*, 'ds=', ds, '< ds_min', ds_min, 'Terminate the continuation!'
  return 
endif 
        
nds = 2 ! the number of available p.o. becomes 2 
 cham(1, :) = v1;   cham(2, :) = v2
 curv(1, :) = x1;   curv(2, :) = x2
    
if (debug == 0) then 
  print*, 'Check after the adjustment of ds'      
  print*, 'angle, ds', angle, ds
  read* 
endif  
    
    
print* , 'Finish angle_interpol!, nds=', nds, 'ds=', ds
return       
  
end  subroutine angle_interpol 
