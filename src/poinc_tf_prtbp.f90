!***********************************************************************
!     ****   poinc_prtbp   ****

!  this is more for the computation of return time, add one more output tf 
!  and remove the computation of differential, do not need at this moment

!  Compute the 2d poincare map, given initial point on this map,  compute its first image 

!  since we have already specify the model to be PRTBP, so we put the called subroutines the real name.

!       Input Variables 
!  pv0            the input 2-dimension state on the 2-d map 

!       Output Variables 
!  tf            return time of Poincare map 
!  pvf           dimension ndim(4), the final full state of the next image of pvin            


!  **NOTE**: the followings are parameters for the Poincare map, and specified through module poinc_mod,
!            before the call, we have to call init_poinc      

! put all the Poincare map related parameters in the module poinc_mod, and use init_poinc for the assignment 

!     --- Advantage of 'include' approach 
!         if we include this subroutine into the module, we can share all the private parameters 

!     --- Disadvantage of 'include approach '   
!         otherwise, we have to declare them public, which will conflict the ones declared in the main routine 
!         we need to put different names for all of them 
 
!    -- I prefer the standalone approach, since we only want to put really related parameters and routines in the moudle 

!       we do not need to clean everytime we make some modification in this routine....      

  
!  ind, p0      pv(ind) = p0 is the poincare section  
!  h            the fixed energy (Jacobi constant, in this case)
!  n            dimension of the Poincare map
!  ndim         the dimension of the full state of the problem is always:  ndim = n+2   
!  ind_fun      the index of the free components in the state, expect pos(ind) and vel(ind) compoments

!  dir          direction to cross the poincare section 
!  imax         number of time of the intersection, 
!  tmax         maximum time between two consecutive crossing 
     
!  Routine Used:
!     

!  Finally Revised by Yu -- 20160905
!***********************************************************************
subroutine poinc_tf_prtbp(pvin, tf, pvf)

use dp_mod
use poinc_mod 
implicit none
integer, parameter :: ndim  = 4 ! dimension of phase space in PRTBP
integer, parameter :: nall = ndim*(ndim+1) ! PV + Variational matrix

!  poinc_mod module-based parameters 
!  integer, intent(in)     ::  n, ndim, ind_fun(n), ind, dir, imax
!   p0, h, tmax  

! Input  and Output Declaration   

real(kind=dp), intent(in)      ::  pvin(n)
real(kind=dp), intent(out)     ::  tf, pvf(ndim) 
 
external :: gr_prtbp, gr_cjprtbp, cj2v_prtbp, dvind_dx_prtbp  

!external :: gamm, deriv_gamm                    
 
! Local Variable
integer        :: npvar, ind_para(ndim), ispc, i, debug 
   
real(kind=dp)  :: fun0(n), pv0(ndim), cj, dgamm_dx(ndim, n)                ! gamm:  R2->R4 

!real(kind=dp), allocatable ::  pci(:),  pcf(:)    ! for this purpose, fixed size ndim(4) is enough 

real(kind=dp)  ::  pci(ndim), pcf(ndim), hminim ! poinc 


! --- start from pv0, compute the first return to the Poincare section  --- 

  
  debug = 0  
  
  ! copy the input, just be sure to keep pvin untouched 
  fun0 = pvin ! just to 
  
  ! Map \gamma: R2-> R4  (y, vy)_0 --> (x,y,vx,vy)_0
  ! subroutine gamm( pvin, n, ind_fun, h0, ind, p0, cj2v, pv)
  call gamm(fun0, n, ind_fun, h0, ind, p0, pv0, cj2v_prtbp)
  
  ! initialize the variational matrix
  pci = pv0 
  
    ! ------------- debug ---------------  
    if (debug == 1) then   
      print*, 'input of poinc_prtbp: ', pvin; print*;  read*

      ! -- check also the energy ! --ckd 
      call gr_cjprtbp(pv0, cj)
      print*, 'Energy: ', cj , 'Prescirbed value h0: ', h0 ;  read* 

      print*, 'initial state before poinc'
      write(*, '(4f18.8)') pci 
      read*    
    endif  
    ! ------------- debug ---------------  


  ! subroutine poinc(sti, ndim, n, 1, tf,stf, hminim, ispc,  deriv, gr_cj) 
  ! without plotting the orbit, ispl=0
  call poinc(pci, ndim, ndim, 1, tf, pcf, hminim, ispc, gr_prtbp, gr_cjprtbp)
  
  if (ispc == 0) then 
    print*, 'Fail to return to the Poincare section'
    read*; return
  endif
  pvf = pcf
 
  
    if(debug == 1) then 
      print*, 'final state after poinc! tf = ', tf   
      print*, pcf
    endif 
  

  return  
end subroutine poinc_tf_prtbp

