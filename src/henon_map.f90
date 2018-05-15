!***********************************************************************
!     ****   henon_map   ****
! This is the standard Henon map, which is conservative, used to check the computation of invariant curve
! Todo: use this map the final decide, what is the necessary input and output of a map 
!       for Poincare map, which is special case, use module to pass the relevent parameters 

p¯=p+Ksinxx¯=x+p¯(1)

!       Input Variables 
!  a            

!       Output Variables 
!  a            

!  Routine Used:
!     

!  Finally Revised by Yu -- 20160
!***********************************************************************
!subroutine poinc_prtbp( pvin, ind, p0, h, n, ndim, ind_fun, dir, imax, tmax, &
!                        tf, pvf, dp_dx)

subroutine henon_map( pvin, n, dir, imax, pvf, dp_dx)

use dp_mod
implicit none
 
! Input  and Output Declaration   
real(kind=dp), intent(in)      ::    
real(kind=dp), intent(out)     ::  pvf(n), dp_dx(n, n)   
 
! Local Variable
real(kind=dp)  ::
  
    
  source

  return  
end subroutine henon_map

!***********************************************************************
!     ****   poinc_prtbp   ****
!  Compute the 2d poincare map and its differential for PRTBP, given initial point 
!  in this map, and compute its image and differential, for purpose to refine

!  since we have already specify the model to be PRTBP, so we put the called subroutines the real name.

!       Input Variables 
!  pv0            the input n-dimension state on the n-d map   
!  ind, p0        pv(ind) = p0 is the poincare section  
!  h              the fixed energy (Jacobi constant, in this case)
!  n              dimension of the Poincare map
!  ndim           the dimension of the full state of the problem is always:  ndim = n+2   
!  ind_fun        the index of the free components in the state, expect pos(ind) and vel(ind) compoments


!  parameters for the Poincare map:   
!    dir          direction to cross the poincare section 
!    imax         number of time of the intersection, 
!    tmax         maximum time between two consecutive crossing 
     

!       Output Variables 
!  tf            return time of Poincare map 
!  pvf           dimension ndim, the final full state of the next image of pvin            
! dpdx           the differential of Poincare map w.r.t. the state pvin 
 
!  Routine Used:
!     

!  Finally Revised by Yu -- 20160
!***********************************************************************
subroutine poinc_prtbp( pvin, ind, p0, h, n, ndim, ind_fun, dir, imax, tmax, &
                        tf, pvf, dp_dx)
