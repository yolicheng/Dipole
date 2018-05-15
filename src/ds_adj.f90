subroutine ds_adj( ds_org, niter, ds)  
! To adjust the stepsize ds for the continuation, based on the number of iteration of Newton Method
! But make sure, there are always 4 available vectors in cham, with the same step size

! bound [ds_min, ds_max] for ipo > 1
! for each value of ds, do at least 4 orbits, in order to provide cham(4, ctr) for adams
! if niter > 6,  ds = ds / 2;   
! if niter < 3,  ds = ds * 2.


! 	Input 
  ds_org	the previous stepsize 
  niter 	number of iterations of Newton Method to refine the previous p.o. 
     
  implicit none 
   
  ! Input and Output declaration 
  real(kind=dp), intent(in) :: ds_org


  
  ! Local varaibles
  real(kind=dp) :: ds_max, ds_min ! bound for stepsize of the continuation 
    
  ! bound vaules of stepsize ds
  ds_max = 1.d-1
  ds_min = 5.d-4  
  
  if(niter .gt. 6) then 
     
        
        ! detect the case, where ds=ds_min, but still fails to reach convergent after 6 iterates of Newton Method
        ! regard this case as failure of the continuation and stop
        
        if( ispo == 0) then 
          if (ds .eq. ds_min) then  
            return
          else 
          
           ! decrease the stepsize by half of the original one, and start again with the old initial guess, without updating the counter ipo
            ds = dmax1(ds_min,  ds / 2 )
            ! check y0_org 
            print*, 'Start again with y0_org:', y0_org
            read*
            y0 = y0_org
            cycle
            
          endif
        endif   
          
      elseif (niter < 3) then 
        ds = dmin1(ds_max, ds*2) ! increase the step size, without exceeding ds_max
      endif 
