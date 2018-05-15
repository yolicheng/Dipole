!***********************************************************************
!     ****   TimeMap_lf   ****
!  Compute the 6d Time T map and its differential for lf, given initial point 
!  in this map, and compute its image and differential, for purpose to refine

!  We need to pass the period for TimeMap

!  since we have already specify the model to be lf, so we put the called subroutines the real name.

!       Input Variables 
!  pv0            the input n-dimension state on the n-d map   
! isdp            flag for the computation of STM or not 

!       Output Variables 
!  pvf           dimension n, the final full state of the next image of pvin            
! dpdx           the differential of Poincare map w.r.t. the state pvin 
 
!       Auxillary Varaibles 
!  n=6           dimension of the Time T map, equal to the dimension of the full phase space

!  Routine Used:
!     

!  Finally Revised by Yu -- 20160
!***********************************************************************
subroutine TimeMap_lf( pvin, pvf, isdp, dp_dx)

! subroutine std_map( ptin, pvf,  isdp,  dpdx)

! TODO --- put this routine in the TimeMap_mod, since we have one parameter tp to pass,
!          for the moment, we assign directly.... later use module... to specify

!          or use this as a standalone subroutine which uses the module poinc_mod
!          whihc option is more convenient? 

!     --- Advantage of 'include' approach 
!         if we include this subroutine into the module, we can share all the private parameters 

!     --- Disadvantage of 'include approach '   
!         otherwise, we have to declare them public, which will conflict the ones declared in the main routine 
!         we need to put different names for all of them 
 
!    -- I prefer the standalone approach, since we only want to put really related parameters and routines in the moudle 

!       we do not need to clean everytime we make some modification in this routine....      

use dp_mod
use lf_mod
use curve_time_mod, only : m, tp, h0 

implicit none

integer, parameter :: n  = 6 ! dimension of phase space in lf

! Input  and Output Declaration   
integer, intent(in)     :: isdp 

real(kind=dp), intent(in)      ::  pvin(n)
real(kind=dp), intent(out)     ::  pvf(n), dp_dx(n, n)  
 
! Local Variable
integer        :: npvar,  i, debug 
real(kind=dp), allocatable ::  yi(:),  yf(:) ! plob 

    
real(kind=dp)  :: cj, stm(n,n), dt, tdir,  &     ! plob
                  
                  ! the rest is to check dp_dx with center difference  
                  step, pvmh(n), cjmh, pvph(n), cjph,  & ! check STM
                  pvfmh(n), pvfph(n), stm_diff(n, n), det  
                  

! --- start from pvin, compute the first return to time T map --- 

  ! if we need the differential matrix, isdp=1, we have npvar = 20
  npvar = n + isdp*n*n
    
  dt   = tp / m 
  tdir = dsign(1.d0, dt)
    
  allocate( yi(npvar) )
  allocate( yf(npvar) )
  
  debug = 0 ! -- dp / d x checked!!!
  step  = 1.d-4
  
  ! Time T Map:  take the isotrophical map, record the point every tp unit time 

  ! initialize the state vector and the variational matrix
  yi      = 0.d0
  yi(1:n) = pvin 
  if(isdp == 1)  yi(n+1 :npvar: n+1) = 1.d0 
  
    ! ------------- debug ---------------  
    if (debug == 1) then  
 
      print*, 'input of TimeMap_lf: ', pvin; print*;  read*
  
      ! -- check also the energy ! --ckd 
      call gr_cjlf(pvin, cj)
      print*, 'Energy: ', cj, 'Prescirbed value h0: ', h0 ;  read* 
  
      print*, 'initial state befor plob '
      write(*, '(6e18.8)')  yi(1:n) 
      read*    
    endif  
    ! ------------- debug ---------------  
  
  call plob(yi, 0.d0, dt, n, npvar, tdir , 0, 0, gr_lf, gr_cjlf,  yf) 
  
  
  pvf = yf(1:n) ! the image of pvin under the Time T map 
  
  ! if we do not need to compute the differential of Poincare map, up to this point is enough
  if (isdp == 0 ) then 
     return 
  else  
     stm = reshape(yf(n+1 : npvar), (/n, n/) )
     dp_dx = stm 
  endif 
  
    ! ------------- debug ---------------  
    if (debug == 1) then  !--ckd

      write(*, '(6e18.8)') yf(1:n); print*
      call gr_cjlf(pvf, cj)
      print*, 'pvf = ', pvf(1:n), 'cj =', cj, 'dcj', dabs(cj-h0)  !ckd 
      read* 
  
    ! differential of Time T map w.r.t. the initial state 
  
    ! ---------------- check d P / d x  by finite difference --------------
    ! ( P(x+h) - P(x-h) / 2h 
  
      stm_diff = 0.d0  

      do i = 1, n
        pvmh = pvin; pvph = pvin
        
        pvmh(i) = pvin(i) - step 
        pvph(i) = pvin(i) + step 
      
        call plob(pvmh, 0.d0, dt, n, n, tdir, 0, 0, gr_lf, gr_cjlf, pvfmh) 
        call plob(pvph, 0.d0, dt, n, n, tdir, 0, 0, gr_lf, gr_cjlf, pvfph) 
    
        stm_diff(:,  i) =  ( pvfph - pvfmh ) / 2.d0 / step 
       
        call gr_cjlf(pvfmh, cjmh)
        call gr_cjlf(pvfph, cjph)
        print*, 'pvfmh = ', pvfmh, 'cjmh =', cjmh  !ckd
        print*, 'pvfph = ', pvfph, 'cjph =', cjph !ckd
        read*
      end do
  
      print*, 'check STM by center difference:  ' 
      do i =  1, n, 1
        write(*, '(10e22.12)') stm_diff(i, :)
      end do
     
!      subroutine detmat( a, n, det)
      call detmat(stm_diff, n, det)
      print*, 'det (STM_diff)  =', det 
      print*; read* 
  
   
      print*, 'STM by variational matrix: ' 
      do i = 1, n, 1
        write(*, '(10e22.12)') stm(i,:)
      end do
      
      print*; 
      print*, 'dif in dp_dx between two methods: maxnorm', maxval(dabs(stm_diff-stm))
      do i = 1, n, 1
        write(*, '(10e22.12)')  stm_diff(i,:)-stm(i,:)
      end do
      
      call detmat(stm, n, det)
      print*, 'det (STM)  =', det 
      print*; read* 

    endif 
   ! ------------- debug ---------------  

  return  
  
end subroutine TimeMap_lf

