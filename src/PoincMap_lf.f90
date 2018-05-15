!***********************************************************************
!     ****   PoincMap_lf   ****
!  Compute the 2d poincare map and its differential for Lorentz force problem, given initial point 
!  in this map, and compute its image and differential, for purpose to refine

!  The related parameters are passed by module 'lf_mod', claim the use of this module in the beginning 

!  since we have already specify the model to be lf_mod, so we put the called subroutines the real name.

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
 
!     

!  Finally Revised by Yu -- 20170519
!***********************************************************************
subroutine PoincMap_lf( pvin, pvf, isdp, dp_dx)

! subroutine std_map( ptin, ptf,  isdp,  dpdx)

! TODO: make a general routine for Poincare map, with the same input and output as std_map
!       in such a way it could be called by the refine_curve_poinc module 
!       so the other parameter should be specified by either module or the function itself. 

! TODO --- must not forget  !!!
! put all the Poincare map related parameters in the module poinc_mod, and use init_poinc 
! for the assignment 

! TODO --- put this routine in the module poinc_mod? 
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
use poinc_mod  
use lf_mod 

implicit none

integer, parameter :: ndim  = 6 ! dimension of phase space in LF
                      
integer, parameter :: nall = ndim*(ndim+1) ! PV + Variational matrix

!  poinc_mod module-based parameters 
!integer, intent(in)     ::  n, ndim, ind_fun(n), ind, dir, imax
!, p0, h, tmax  debug,


! Input  and Output Declaration   
integer, intent(in)     :: isdp 

real(kind=dp), intent(in)      ::  pvin(n)
real(kind=dp), intent(out)     ::  pvf(n), dp_dx(n, n)  
 
!external :: gr_lf , gr_cjlf, cj2v_lf, dvind_dx_lf 
  
external :: gamm, deriv_gamm                    
 
!  Local Variable
integer        :: npvar, ind_para(ndim), ispc, i, debug, isv2 
   
real(kind=dp)  :: fun0(n), pv0(ndim), cj, dgamm_dx(ndim, n)               ! gamm:  R4 -> R6

real(kind=dp), allocatable ::  pci(:),  pcf(:)  

real(kind=dp)  :: tf, stm(ndim,ndim), hminim, pv_pc(ndim), vf_pc(ndim), & ! poinc
                  dpi_dx(n, ndim),  &                                     ! auxillary differential
                  
                  ! the rest is to check dp_dx with center difference  
                  step, funmh(n), pvmh(ndim), cjmh, &   
                  funph(n), pvph(ndim), cjph, dgamm_dx_diff(ndim, n), &   ! check deriv_gamm
                  
                  pciph(nall), tfph, pcfph(nall), &                       ! ck d Pi / d x
                  pcimh(nall), tfmh, pcfmh(nall), dp_dx_diff(n, n), det 

! --- start from pv0, compute the first return to the Poincare section  --- 

  ! if we need the differential matrix, isdp=1, we have npvar = 20
  npvar = ndim + isdp*ndim*ndim  
  allocate( pci(npvar) )
  allocate( pcf(npvar) )
  
  
  debug = 0 ! -- dp / d x checked!!!
  step  = 1.d-3
  
  fun0 = pvin 
  
!  print*, 'Check PoincMap_lf: n, ind_fun, h0', n,  ind_fun, h0; print*; read*
  
  ! Map \gamma: R4-> R6  (y, vy)_0 --> (x,y,vx,vy)_0
  !subroutine gamm( pvin, n, ind_fun, h0, ind, p0, pv, isv, cj2v)
  call gamm(fun0, n, ind_fun, h0, ind, p0, pv0, isv, cj2v_lf) 
  if(isv == 0) return 
  
  if(isdp == 1) call deriv_gamm( pv0, ndim, ind, dgamm_dx, dvind_dx_lf) 
  
  ! initialize the variational matrix
  pci = 0.d0
  pci(1:ndim) = pv0 
  
  
  if(isdp == 1)  pci(ndim+1 :npvar: ndim+1) = 1.d0 

  ! ------------- debug dgamm_dx  ---------------  
  if (debug == 1) then   !  -------- ckd 2017-05-21 16:49:25 
!  !  set debug=2 after check, and then reset back to 1 after all the check
 
    print*, 'input of PoincMap_lf: ', pvin; print*;  read*
  
    print*, 'check pv0 by gamm: ', pv0; read* ! --ckd
  
    ! -- check also the energy ! --ckd 
    call gr_cjlf(pv0, cj)
    print*, 'Energy: ', cj , 'Prescirbed value h0: ', h0 ;  read* 
  
    ! -- check dgamm_dx   ! --ckd
    print*,'dgamm_dx'
    do i = 1, ndim, 1
      print*, dgamm_dx(i, :)  
    end do
    print*;  read*  
  
    ! ----- check dgamm_dx by finite difference: 
    !   ( \gamm(y+h) - \gamm(y-h) ) / (2h)
    do i = 1, n 
      funmh = fun0
      funmh(i) = fun0(i) - step 
      call gamm(funmh, n, ind_fun, H0, ind, p0, pvmh, isv2, cj2v_lf)
    
      funph = fun0
      funph(i) = fun0(i) + step 
      call gamm(funph, n, ind_fun, H0, ind, p0, pvph, isv2, cj2v_lf)
    
      call gr_cjlf(pvmh, cjmh)
      call gr_cjlf(pvph, cjph)
    
      print*, 'pvmh:', pvmh,  'cj=', cjmh
      print*, 'pvph:', pvph,  'cj=', cjph; read*
    
      dgamm_dx_diff(:, i) = (pvph - pvmh) / 2 / step 
    enddo 
  
    print*, 'check dgamm_dx by center difference:  ' 
    do i =  1, ndim, 1
      print*, dgamm_dx_diff(i, :)
    end do
    print*; read* 
   
    print*, 'initial state befor poinc '
    write(*, '(6f18.8)') pci(1:ndim) 
    read*    
  endif  
  ! ------------- debug ---------------  


  ! subroutine poinc(sti, ndim, nvar, tdir, tf, stf, hminim, ispc, deriv, gr_cj) 
  ! without plotting the orbit, ispl=0
  
  call poinc(pci, ndim, npvar, 1, tf, pcf, hminim, ispc, gr_lf, gr_cjlf)
!  print*, 'poinc from PoincMap_lf:', pcf(1:ndim);  !ckd
!  print*; read*
  
!  call poinc(pvf, ndim, ndim, 1, tpc, pvpc, hminim, ispc, gr_lf, gr_cjlf) 
    
  if (ispc == 0) then 
    print*, 'Fail to return to the Poincare section'
    read*; return
  endif

  if(debug == 1) print*, 'final state after poinc! tf = ', tf  
  
  pvf = pcf(ind_fun) ! the image of pv_in under the Poincare map 
  
  ! if we do not need to compute the differential of Poincare map, up to this point is enough
  if (isdp == 0 ) return 
  
  pv_pc = pcf(1:ndim) ! the final state 
  call gr_lf(0.d0, pv_pc,  ndim, vf_pc)  
  
  
    ! ------------- debug ---------------  
    if (debug == 1) then  !--ckd 

      write(*, '(6f18.8)') pcf(1:ndim); print*
    
      print*, 'pv_pc = ', pv_pc(1:ndim)  !ckd 
      print*, 'vf_pc = ', vf_pc(1:ndim) 
      read* 
    endif     
    ! ------------- debug ---------------  
    
    
  stm = reshape( pcf(ndim+1:npvar), (/ndim,ndim/) ) ! the STM 
  
  ! differential of Poincare map w.r.t. the initial state 
  ! in order to make it more general, compute all the components of dp_dx0
  
  !  subroutine diffpoinc( phi, f,  ndim, ind, nr, nc, para, fun, dpdx)
  ind_para = (/(i, i = 1, ndim)/)  ! implied do loop
  call diffpoinc( stm, vf_pc,  ndim, ind, n, ndim, ind_para, ind_fun, dpi_dx)
 

    ! ------------- debug ---------------  
    if (debug == 1)  then  !--- ckd    2017-05-21 16:56:14 
   
      print*, 'check d Pi / d x by diffpoinc:  '  
      do i =  1, n, 1
        print*, dpi_dx(i, :)
      end do
      print*; read* 
  
    ! ---------------- check d P / d x  by finite difference --------------
    ! ( P(x+h) - P(x-h) / h 
  
      dp_dx_diff = 0.d0  

      do i = 1, n
      funmh = fun0
      funmh(i) = fun0(i) - step 
      call gamm(funmh, n, ind_fun, H0, ind, p0, pvmh, isv2,  cj2v_lf)
      
      funph = fun0
      funph(i) = fun0(i) + step 
      call gamm(funph, n, ind_fun, H0, ind, p0, pvph, isv2, cj2v_lf)
    
      pciph = 0.d0
      pciph(1:ndim) = pvph 
      pciph(ndim+1 :npvar: ndim+1) = 1.d0 
    
      pcimh = pciph;  pcimh(1:ndim) = pvmh
    
      call poinc(pcimh, ndim, npvar, 1, tfmh, pcfmh, hminim, ispc, gr_lf, gr_cjlf)
    
      call poinc(pciph, ndim, npvar, 1, tfph, pcfph, hminim, ispc, gr_lf, gr_cjlf)
      dp_dx_diff(:,  i) =  ( pcfph(ind_fun) - pcfmh(ind_fun) ) / 2 / step 
      
      print*, 'final state after poinc, tf = ', tfph
      
      call gr_cjlf(pcfph(1:ndim), cj) 
      print*, 'Energy: ', cj ; print*; read* 
      
      write(*, '(6f18.8)') pcfph(1:ndim); !print*
      
      read*
    end do
  
      print*, 'check d P / d x by center difference:  ' 
      ! TODO: something is wrong here. 
      do i =  1, n, 1
        write(*, '(10f22.14)') dp_dx_diff(i, :)
      end do
      
      call detmat(dp_dx_diff, n, det)
      
      print*, 'det (D P / D X)  =', det 
      print*; read* 
    end if 
    ! ------------- debug ---------------  
  
   
  !    d Pi / d \hat P (x,y,vx, vy)_f  *  d P(x,y,vx, vy)_f / d (x,y,vx, vy)_0  
  ! =  d Pi / dxf  * dp_dx , dimension  n - by - ndim 
  !   where  d Pi / dx = [0 1 0 0; 0 0 0 1]  
   
  !  so it could be computed by seleciting the corresponding rows  of dp_dx 
  
  ! dimension n-by-n, n-by-ndim X ndim-by-n
  dp_dx = matmul(dpi_dx, dgamm_dx)
  
    ! ------------- debug ---------------  
    if (debug == 1) then  !--ckd, put condition to 2
    
      print*, 'dp_dx, dimension:', n, ' X ', n
      do i = 1, n, 1
        write(*, '(10f22.14)') dp_dx(i, :)
      end do
  
      call detmat(dp_dx, n, det)
      print*, 'det (D P / D X)  = 1 ?', det 
      print*; read*  
  
      print*, 'dp_dx - dp_dx_diff' 
      do i = 1, n, 1
        write(*, '(10f22.14)') dp_dx(i, :) - dp_dx_diff(i, :)
      end do
      print*; read*
      
    endif 
   ! ------------- debug ---------------  
      

  return  
end subroutine PoincMap_lf

