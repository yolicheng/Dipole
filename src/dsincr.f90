subroutine dsincr( ds, nds, cham, cham5, curv )
! This routine deal with the case when the step size need to be doubled with niter < 3, and ds*2 < ds_max 
! if ds*2 > ds_max, keep the current ds without increasing 

! Varaibles need to be updated: cham, cham5, curv 

!  All the reassignments have been  checked carefully. -20160419     

! 	Input and Output Varaibles
!  ds 		the current and updated prediction of the stepsize for the next p.o. 
!  nds 		number of available p.o.s in cham5, with the same fixed step size 
!  cham		the vector filed of the latest(up to 4) points, used by Adams Multistep predictor
!               to predict the vector field of next p.o. 
!  cham5 	all the previous vector filed with fixed step size, used to double the step size 
!  curv 	the points(up to 3) on the characteristic curve, also with the fixed step size, corresponding to cham  

! 	Module-based Parameters
!  nctr, debug, ds_max 

implicit none 
integer, parameter :: dp = kind(1.d0)

integer, intent(inout) :: nds 
real(kind=dp), intent(inout) :: ds, cham(4, nctr), cham5(5, nctr), curv(3, nctr)

! Local Parameters 
integer :: i

! Here, the rows of curv and cham are corresponding to each other
! if nds = 3, we keep the available vectors as the last 2 rows in curv, because only the last 2 rows will be used, if necessary (in champ)
! if nds = 1, keep only available the first row in curv
 
! check cham, cham5, curv after the dscrease  ! -- ckd 
if (debug == 0) then 
  print*, 'Before dsincr, check!'
  print*, 'ds=',ds, 'nds=', nds

  print*, 'cham'
  do i = 1, 4
    print*, cham(i, :)
  enddo
  print*          

  print*, 'cham5'
  do i = 1, 5
    print*, cham5(i, :)
  enddo
  print*   

  print*, 'curv'
  do i = 1, 3
    print*, curv(i, :)
  enddo
  print*   
endif 

!   if ds*2 > ds_max, we cannot use the previous np:np-4: 2 points...., not fixed step size 

if ( ds * 2.d0 <= ds_max) then 
          
  print*, 'ds needs to be increased!'
 
  if(debug == 1) read*  !ck
        
  ! we have the previous 5 points to provide us 2 points with stepsize as 2*ds 
  ds = ds * 2.d0
  
  print*, 'nds>=5', nds 
  if(nds < 5) read*
  
  cham(1:3,:) = cham5(1:5: 2, : ) 
  
  curv(2,:) = curv(1,:) ! curv(2:3,:) = curv(1:3:2, :) ! the first row is of no meaning  
  nds = 3
  cham5(1:3,:) = cham(1:3, : )

         
else 
  
  ! this part will never be called, because the condition for tha clal of dsincr is   if( nds .ge. 4 .and. 2.d0*ds .le. ds_max ) 

  print*, 'Will never reach this part,  Error!!!! ds*2.d0 > ds_max'
  read* 

! if ds*2.d0 > ds_max, special operation is needed to guarantee ds not to exceed ds_max
! after the increase, we only have 1 available p.o. 

! if ds == ds_max already, keep everything unchanged 

  if(dabs(ds-ds_max) .gt. 1.d-14) then
    cham(1,:) = cham5(nds,: ) 
    curv(1,:) = curv(3, :)
    nds = 1
    
    cham5(1, :) = cham(1, : )
  endif 
  
  ds = ds_max  

endif 
  
! check cham, cham5, curv after the dscrease  ! -- ckd  
if (debug == 0) then 
  print*, 'After dsincr, check!'
  print*, 'ds=',ds, 'nds=', nds

  print*, 'cham'
  do i = 1, 4
    print*, cham(i, :)
  enddo
  print*          

  print*, 'cham5'
  do i = 1, nds
    print*, cham5(i, :)
  enddo
  print*   

  print*, 'curv'
  do i = 1, 3
    print*, curv(i, :)
  enddo
  print*
  
  print*, 'After the increase!, nds = ', nds, 'ds=', ds
  read* 
          
endif 

 
return      
end subroutine dsincr

      
       
