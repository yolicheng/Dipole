!  **** My curve Module*******
!  Compute and refine the invariant curve obtained by Poincare map method,  with the energy and rotation number are both fixed
!  Take PRTBP for example, choosing the Poincare section to be x = x_L4 

!  the idea is different from the Time t map, 
!  so do a module seperately, although they have some similar part, but do not mix the two methods


!  we have a 2d map that maps R^2(y, vy) to R^2 (y, vy), take a good look at Alex's note 
!  and vx is a function of vx = f(x0, y, vy, h0)
!  so vx is not free, the differential of the Poincare map should be w.r.t. the two unknowns (y, vy) 
!  taking the  vx = f(x0, y, vy, h0) into consideration 

!  It's suggested to write the equations in a compact way, and use matrix to express the Poincare map  
!  such that we will not miss anything

! Assume we have \hat P: (x,y,vx, vy)_0 --> (x,y,vx, vy)_f, actually the flow \phi
! introduce a map \gamma:  (y, vy)_0 --> (x,y,vx, vy)_0 
! A projection Pi (the inverse of gamma): (x,y,vx, vy)_f -->  (y, vy)_f 


! Poincare map is P:  (y, vy)_0 --> (y, vy)_f, can be expressed as a series of compositions:
  
!   Pi(   \hat P  ( \gamma( (y, vy)_0 ) )  )   )

!  so d P / d (y, vy)_0  = d Pi / d \hat P *  d \hat P / d \gamma *  d \gamma / d (y, vy)_0
! where the first term: 
!   d Pi / d \hat P (x,y,vx, vy)_f  = d Pi /  [0 1 0 0; 0 0 0 1]   ! the second and fourth rows of STM

! the second term 
! d \hat P / d \gamma = d (x,y,vx, vy)_f / d (x,y,vx, vy)_0 = STM 

! the third term --- ** NOTE ** be careful with this one, vx= f(x0, y, vy, h0) 
! d \gamma / d (y, vy)_0 = d (x,y,vx, vy)_0 / d (y, vy)_0   
!                        = [0    0; 1    0; d vx/ d y     d vx / d vy;  0    1 ]

! following several steps : 

!  -- 1. take a Poincare map, the section is specified by pv(ind) = p0, for example, x = x_L4 

!  -- 2. In order to get equispaced points to do Fourier analysis  
!        -- 2.1  compute the rotation number  within [0,1] by Carles' note, multiply by 2pi to obtain rho \in [0, 2pi]
!        -- 2.2  sort all the points by the value \rho, and then do linear interpolation within [0, 2pi] for y and vy
!                the last component vx is obtained by the prescribed energy vx = f(H0, x0, y, vy)
!        -- 2.3  using the general Fourier.f routine to compute the Fourier coefficients, which provide the initial seed for refinement
!                here, we need to check the approximated Fourier representation and the interpolated points 

!  -- 3. do the refinement using Newton Method by imposing the invariance equation 
!        ** NOTE **
!        3.1  there is only one indetermination in Poincare map method, the phase shift indetermination, 
!             it can be avoided by fixing C1_1 = 0, this is done by remove C1_1 from the unknows, although in the system equations
!             it is expressed by adding another equation  C1_1 = 0. 
!             since if we keep C1_1 as one unknown, then the Jacobi Matrix will be singular, which makes it difficult for Newton Method 
!        3.2  by fixing the energy H0 and the rotation number \rho, we fix a particular invariant curve 
!        3.3  the unknowns are only the Fourier coefficients: C0, CK, SK 

!  so, we have the number of unknowns : n * (2*nf+1) - 1 
!      where n is the dimension of the phase space to be evaluated,  n = ndim -2
!      since we fix h0 and x0, decrease dimension by 2
!      nf is the number of Fourier modes. 

!  the number of equations 
!      n * (2*nf+1)  

! ** Note ** 
!    1. we have 2 more equations than unknows, but Josep Maria claims the kernel dimension is zero... 
!    2. use the general name for vector field and energy (Jacobi costant)
!    3. subroutines from libraries, need to declare external attribute **
!         external  ::  deltx, trigrec_c, trigrec_s   ! -- libgr.a 

module curve_poinc_mod

! TODO: what if dp and pi, pi2 conflict with the ones in the main routine?
use dp_mod ! dp and pi, pi2 
use pi_mod 

! dimension declaration and termination control of Newton method
!use tori_mod, only: nf, nf2, nf2p1, ndim, nitmax, tol, tol_err 

implicit none
save 

! subroutines from libraries, declare external attribute, 
! this is only necessary for module to call external subroutines or functions 
external  ::  deltx,  trigrec_c, trigrec_s   ! -- libgr.a 

real(kind=dp), private ::  tmax, p0, &  ! poinc
                           tol, tol_err ! for tori_mod TODO
                           
! --- For declaration of array ---- 
integer, private   ::  nc_dg, nr_dg, nc_dgall, nr_dgall, nc_dgxi1, & ! dimension
                       ind, dir, imax,  & ! poinc 
                       nf, nf2, nf2p1, ndim, nitmax , & ! for tori_mod TODO
                       debug ! for debug 

integer :: ind_vel, n                        
integer, allocatable  :: ind_para(:), ind_fun(:) ! poinc                        

! --- For save -- 
! to write the output of refinment process
integer, private :: fwork, fout 


contains
  include  'refine_curve_prtbp.f90'
  include  'gamm.f90'       ! map from (y, vy) to (x0, y, vx, vy)
  include  'deriv_gamm.f90' ! differential of gamm w.r.t.  (y, vy)
  
!***********************************************************************
!     ****   init_curve_poinc   ****
!  the dimension of related arrays 
!  From tori_mod:  nf, nf2, nf2p1, ndim 
!  for Poincare map: ind, p0, tmax 
!***********************************************************************
subroutine init_curve_poinc(nf0, ndim0, nitmax0, tol0, tol_err0, ind0, p00, dir0, imax0, tmax0 ) 

implicit none
integer, intent(in)          ::  ind0, dir0, imax0, &
                                 nf0, ndim0, nitmax0  
                                 
real(kind=dp), intent(in)    ::  p00, tmax0, &
                                 tol0, tol_err0   

integer :: i, k
  print*, 'Debug (1) or not?'
  read(*,*)  debug
!  debug = 0  
  
  
! for the tori_mod... TODO: put this part into tori_mod later
  ndim  = ndim0
  n     = ndim - 2 ! the dimension of Poincare map 
  nf    = nf0
  nf2   = 2*nf  
  nf2p1 = nf2 + 1
  
  
  ! termination control for Newton method 
  tol     = tol0
  tol_err = tol_err0
  nitmax  = nitmax0 
  
!  print*, 'check the initialization of curve_poinc_mod module'
!  print*, 'nf, ndim, tol, tol_err, nitmax', nf, ndim, tol, tol_err, nitmax
!  !read*  
  
  
! ---- commonly used value for dimension declaration  
  ! number of free components of the phase 
  ! since one component is specified by Poincare map, we only need to compute the rest ones 
    
  ! number of columns and rows in dg, which is of dimension ndim * (nf2p1)+1 - by -  ndim * (nf2p1)-1
  nr_dg  = n * nf2p1  
  nc_dg  = nr_dg  
  
!  dicard this approach at this moment  
!  nc_dgm1 = nr_dg -1       ! fix c1_x = 0 to avoid the phase shift indetermination
  
  ! dimension of dgall, the number of row is the same as nr_dg 
  ! number of the column is 1 less, since we need to fix c1_(ind)
  
  ! dimension of dg for one value of xi 
!  nr_dgxi1 = nf2p1      ! withou the energy equation
!  nc_dgxi1 = ndim*nf2p1 ! number of unknowns, all the coefficients for (T,y,vx,vy), each is of dimension nf2p1 

  ! For Poincare map 
  ind = ind0
  
  ! the velocity in the direction of pv(ind) 
  ind_vel = ind + ndim / 2  
  
  p0   = p00
  tmax = tmax0  
  dir  = dir0
  imax = imax0 
  
  ! assign the index of the component in  free parameters and functions 
  allocate(ind_para(ndim)) 
  allocate(ind_fun(n)) 
  
!  ind_para =  (/ (i, i = 1, ndim) / ! implied do loop

  k = 1
  do i = 1, ndim, 1
    ind_para(i) = i 
    if(i == ind .or. i == ind_vel) cycle 
    ind_fun(k) = i
    k = k+1
  end do
  
!  for save - TODO: for the moment, not used
   fwork = 6
!  fwork = 222 
!   open(fwork, file='curve_poinc.work')
  
if(debug == 1) then 
  print*, 'check ind, p0, tmax, dir, imax',  ind, p0, tmax, dir, imax
  print*, 'finish init_curve_poinc!'
  read*  
endif    
  
  return  
end subroutine init_curve_poinc

  

!*********************************************************************** 
!  This routine computes the coefficient matrix and error for the linear equations to be solved
!  by Newton Method, as input of deltx_lnx 
!                  dg *  dX = - ferr
!                 
!  For the nf values of xi_i, i = 0, ..., 2*nf 
!  
!  1.  the whole diffrential dg: D F/ d   A(C_0, C's, S's)  
!  2.  the error: ferr,  to be refined by Newton method 

! where the matrix A include components T,y,vx,vy
! the Fourier coefficients of the return time T also needs to be refined
! to save memory, we store them in the column ind (the Poincare section) 


!  -- the Linear equation F to slove for Newton Method : 
!     F:  \varphi(xi_i + rho) - P (xi_i) = 0 ],  i = 1,..., Nf   ! the invariance condition  
!          Note that the invariance conditon is for all the  ndim components  
!          the one chosen to be Poincare section is replace by the return time T 
!          For instantce, if ind =1, we have (T, y, vx, vy)
  
!  the unknowns are written in column vector, of dimension n *(2*Nf+1) - 1:
!  ---  X = [A_y, ... , A_vy]  if ind=1... 
  
! the component of dg = d F / d A(C_0, C's, S's)     
!        for each xi, compute by subroutine gdg_inv_xi1
!        do not forget:  -d P( \varphi(xi_i) ) / d A  
  
! -- 2. TO make everything clear, we compute the full matrix to aviod confusion
!       Do the substraction of corresponding submatrix in the routine refine_curve_poinc

!       Input Variables 
!  h            the presribed vaule of energy 
! c0, ck, sk   the initial unkowns to be refined, in total is of dimension  n * (2*nf+1)


!       Output Variables 
! ferr_all    the error in F to be refined, dimension:   n*(2*nf+1) 
! dg_all       the differential of F w.r.t. the unkowns, 
!             dimension: n*(2*nf+1) - by -   n*(2*nf+1) -1       
! 
!        Module-Based Varaibles
! ndim, nf, nf2, nf2p1,  pi2 
! nr_dg = n * ( 2*nf+1 )  , nc_dg = n * ( 2*nf+1 ) - 1


!  routine Used: -- TODO 
!     None

!  Finally Revised by Yu -- 20160731
! ---- TODO : check 20160826
!------------------------------------------------------------------------
subroutine  gdg_curve_poinc( h, rho, c0, ck, sk, ferr_all, dg_all, deriv, gr_cj, cj2v, dvind_dx )
implicit none

! Input  and Output Declaration   
real(kind=dp), intent(in)    ::  h, rho,  c0(n), ck(n, nf), sk(n, nf) 
real(kind=dp), intent(out)   ::  ferr_all(nr_dg),  dg_all(nr_dg, nc_dg) 
 
! Local Variable
integer        :: i, k, irow  
real(kind=dp)  :: dxi, xi, dg_xi1(n, nc_dg),  ferr_xi1(n) 

external :: deriv, gr_cj,  cj2v, dvind_dx 
  
  ! Initialize as zeros, then we only need to deal with nonzero components 
  dg_all = 0.d0 
  
  ! ---  The Jacobi Matrix of the general function F ----
  !    xi_0   = A_0 ([y, vy])
!  ----- deal with argument xi_i one by one  ---------   

! time step size for discretisize the parameter space, dxi = 2*pi/(1+2*nf), 
! and we need  to evaluate xi_i = i * dxi, where i = 0, ..., 2*nf 

  ! step size for xi 
  dxi   = pi2 / nf2p1
  
  do i = 1, nf2p1, 1
  
    ! the current xi to be evaluated
    xi = (i-1) * dxi
    
    ! the starting row for each block of xi_i, i = 0, ..., 2*nf 
    irow = 1 + (i-1) * n
    
!    print*, 'i=', i, 'xi=', xi, 'starting row:', irow 
        
    !  subroutine gdg_inv_xi1( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, deriv, gr_cj, cj2v, dvind_dx)
    call gdg_inv_xi1( xi, h, rho, c0, ck, sk, dg_xi1, ferr_xi1, deriv, gr_cj, cj2v, dvind_dx)
    
    print*, 'ferr_xi1 = ', ferr_xi1 
    
!    print*, 'dg_xi1 = '
!    do k =  1, n, 1
!      write(*, '(12f18.8)') dg_xi1(k,:)
!    end do  
!    print*; read*    
    
    ! the sub-block in dg associated with one value of xi, ndim rows starting from irow
    dg_all(irow : irow+n-1,  :) = dg_xi1
    
    ! error in F (the invariance equations)
    ferr_all(irow : irow+n-1 ) =  -ferr_xi1 
    
  end do 
  
  print*, 'finish dg computation!'
  return 
end subroutine gdg_curve_poinc


!***********************************************************************
!     ****   gdg_inv_xi1  ****
!  for a certain value of xi, the Jacobi Matrix of the invariance equation 
!  we fix the energy level here, so h is eliminated from the unknowns 

!  **NOTE** 
!    \varphi(xi + rho) - P ( \varphi(xi) ) = 0,  i = 1,..., Nf, for components (y,vx,vy)

!-- where P_T is the return time , 
!   the C_T and S_T are computed in the first iteration   
!   so  the refinement in the first iteration only to improve the precision of the Fourier anay


!  the unknowns are written as a column vector, of dimension ndim *(2*Nf+1):
!  ---  X = [  A_T,  A_y, ... , A_vy ] (ind=1 as example) 

!   we denote A = [ C0, C's, S's] as the Fourier coefficients, then A is of dimenson ndim * (2*Nf+1), 
!   A = [A_x, A_y, A_z, A_vx, A_vy, A_vz], and ind-th column is replaced by the derivative w.r.t. the return time T in Poincare map method 
    
!  the component of dg =
!   -- d F  / d A_(y, vy) ---  ( as a block of dimension 2*Nf + 1) 
!      we also apply the chain rule to compute  d P / d A_(y, vy) 

! Finally, we have: 
!    d \varphi(xi+rho) / d A_(y, vy) - [ d P / d \varphi(xi) ] * d \varphi(xi) / d A_(y, vy)
!            
! where   \varphi is a truncated Fourier series evaluated at xi 
!         \varphi (xi+rho) is the same Fourier series evaluated at xi+rho 
  
!       Input Variables  (initial guess)
!  xi           the argument for the invariant curve  \varphi 
!  rho          rotation number:  2pi / omega2 * omega1
!  c0, ck, sk   dimension ndim-by-nf,  the unknowns C_k and S_k (k=1, ..., nf), for each component
!               (x,y,z,vx,vy,vz), the values are different. 

                    
!       Output Variables 
!  dg_xi1      dimension: n(4 or 2), n * ( 2*nf+1 ) 
!              the differential (Jacobi Matrix) of invariance equations w.r.t. A_{y, vy}
 
! ferr_xi1     dimesion n, the error to be refined by Netwon method 
!              F   = \varphi(xi_i + rho) - P (xi_i)  

!     Module based Varaibles
! n, ndim, nf      dimension of the phase space and the model, number of fourier modes
! for poinc:  ind, p0, dir, imax

!  routine Used:
!   trigrec_c, trigrec_s 

!  Finally Revised by Yu -- 20160826
!***********************************************************************
subroutine gdg_inv_xi1( xi, h, rho, c0, ck, sk, dg_xi1, ferr_xi1, deriv, gr_cj, cj2v, dvind_dx )
implicit none

! Input  and Output Declaration   
!integer, intent(in)            ::  nf 
real(kind=dp), intent(in)      ::  xi, h, rho, c0(n), ck(n, nf), sk(n, nf) 
                                   
real(kind=dp), intent(out)     ::  dg_xi1(n, n*nf2p1), ferr_xi1(n) 
 
! Local Variable
integer        :: i, j ,k, npvar, jcol_st, jcol_end,  ispc, isvy   

real(kind=dp)  :: cn(nf2), sn(nf2), cn_rho(nf2), sn_rho(nf2), &
                  fun0(n), pv0(ndim), vf0(ndim),  pv0_rho(n), cj,  &   ! the initial state 
                  pci( ndim*(ndim+1) ),  pcf( ndim*(ndim+1) ), pv_pc(ndim), stm(ndim,ndim), & !  poinc
                  dp_dx(n, n), dpi_dx(n, ndim),  dv_dx(n), dgamm_dx(ndim, n), & ! poinc
                  dp_da(n, n*nf2p1),  hminim, tf,  &  ! poinc , dp_da
                  dvphi_da(n, n*nf2p1), dvphi_rho_da(n, n*nf2p1), vf_pc(ndim), &  ! for dg_xi1
                  step, funmh(n), pvmh(ndim), cjmh, pvph(ndim), cjph, dgamm_dx_diff(ndim, n), & ! check deriv_gamm
                  pciph( ndim*(ndim+1) ), tfph, pcfph( ndim*(ndim+1) ),  dp_dx_diff(n, n), & ! ck d Pi / d x
                  c0ph(n), ckph(n,nf), skph(n,nf), funph(n), dvphi_da_diff(n, n*nf2p1), arg, & ! check dvphi_da_diff
                  dp_da_diff(n, n*nf2p1), &   ! check dp_da 
                  funph_rho(n), funmh_rho(n), dg_da_diff(n, n*nf2p1), & ! check dg_da the final check!!!!!
                  
                  ! (x-h)
                  c0mh(n), ckmh(n,nf), skmh(n,nf), pcimh( ndim*(ndim+1) ), tfmh, pcfmh( ndim*(ndim+1) ), det  

!  External subroutines                  
external :: deriv, gr_cj, cj2v, dvind_dx 

  
  step =   1.d-3  ! to check with center difference
  
  dg_xi1 = 0.d0 ! initialization
  
 ! the trigometric function at k* (xi+rho) for k = 1_2Nf
  call trigrec_c( xi+rho, nf2, cn_rho)  
  call trigrec_s( xi+rho, nf2, sn_rho)
 
  ! ---- The Jacobi Matrix of the general function F (the ndim coordinates) w.r.t. A---
  
  ! Evaluate the invariant curve \varphi(xi), with the initial guess of coefficients (c0, ck, sk)
  call varphi( xi, c0, ck, sk, fun0) 
  
  ! -- check \varphi, compare fun0 with the output of subroutine four_seri   -- 
  
  ! start from pv0, compute the first return to the Poincare section 
  ! we need to compute the variational matrix, so the size of initial state is  n = ndim*(ndim+1)
  
  ! initialize also the variational matrix
  npvar = ndim*(ndim+1)
  
  ! the inverse of map PI, \gamma: (y, vy)_0 --> (x,y,vx,vy)_0
  !  subroutine gamma_prtbp( pvin, n, ind_fun, h0, ind, p0, cj2v, pv)
  call gamm(fun0, n, ind_fun, h, ind, p0, pv0, cj2v)

if( debug == 1) then   
  print*, 'check pv0 by gamm: ', pv0; read* ! --ckd
  ! check also the energy ! --ckd 
  call gr_cj(pv0, cj)
  print*, 'Energy: ', cj , 'Prescirbed value h0: ', h 
  read*    !ck
endif 
  
!  subroutine deriv_gamm( pv, ndim, ind, dgamm_dx, dvind_dx) 
  call deriv_gamm( pv0, ndim, ind, dgamm_dx, dvind_dx) 

if(debug == 2) then   ! TODO -------- ckd, bingo! 
  ! check dgamm_dx   ! --ckd
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
    call gamm(funmh, n, ind_fun, h, ind, p0, pvmh, cj2v)
    
    funph = fun0
    funph(i) = fun0(i) + step 
    call gamm(funph, n, ind_fun, h, ind, p0, pvph, cj2v)
    
    call gr_cj(pvmh, cjmh)
    call gr_cj(pvph, cjph)
    print*, 'pvmh:', pvmh,  'cj=', cjmh
    print*, 'pvph:', pvph,  'cj=', cjph; read*
    
    dgamm_dx_diff(:, i) = (pvph - pvmh) / 2 / step 
  enddo 
  
  print*, 'check dgamm_dx by center difference:  ' 
  do i =  1, ndim, 1
    print*, dgamm_dx_diff(i, :)
  end do
  print*; read* 
  ! -----------ckd ----------------
endif 

  pci = 0.d0
  pci(1:ndim) = pv0 
  pci(ndim+1 : npvar : ndim+1) = 1.d0 

if(debug == 1) then  
  print*, 'initial state befor poinc '
  write(*, '(4f18.8)') pci 
  read*    
endif  


!subroutine poinc(sti, ndim, n, ind, p0, dir, imax, tmax, ispl, fob, tf,stf, hminim, ispc,  deriv, gr_cj) 
  ! without plotting the orbit, ispl=0
  call poinc(pci, ndim, npvar, ind, p0,  dir, imax, tmax, 0, 6, tf, pcf, hminim, ispc, deriv, gr_cj)
  
  print*, 'final state after poinc '
  print*, 'tf = ', tf 
  if(debug == 1) write(*, '(4f18.8)') pcf; print*
  !!read*    
  
  pv_pc = pcf(1:ndim) ! the final state 
  call deriv(0.d0, pv_pc,  ndim, vf_pc)  
  
!  if(debug == 1) then  !--ckd
!    print*, 'pv_pc = ', pv_pc  !ckd 
!    print*, 'vf_pc = ', vf_pc 
!    read* 
!  endif     
  
  stm = reshape( pcf(ndim+1:npvar), (/ndim,ndim/) ) ! the STM 
  
  ! differential of Poincare map w.r.t. the initial state 
  ! in order to make it more general, compute all the components of dp_dx0
  
  !  subroutine diffpoinc( phi, f,  ndim, ind, nr, nc, para, fun, dpdx)
  call diffpoinc( stm, vf_pc,  ndim, ind, n, ndim, ind_para, ind_fun, dpi_dx)
 
if(debug == 2)  then  !--- ckd, put condition to 2
  print*, 'check d Pi / d x by diffpoinc:  '  
  do i =  1, n, 1
    print*, dpi_dx(i, :)
  end do
  print*; read* 
endif 
  
  ! TODO: since A0(ind) is fixed, we should assign the ind-th column of dp_dx to zeros 
  if(ispc == 0) then 
      print*, 'Poinc failed!'  
      !!read*     
  endif 
  
  ! ---------------- check d P / d x  by finite difference --------------
  ! ( P(x+h) - P(x-h) / h 
  
if(debug == 2)  then   ! ckd, put to debug == 2
dp_dx_diff = 0.d0  

  do i = 1, n
    funmh = fun0
    funmh(i) = fun0(i) - step 
    call gamm(funmh, n, ind_fun, h, ind, p0, pvmh, cj2v)
    
    funph = fun0
    funph(i) = fun0(i) + step 
    call gamm(funph, n, ind_fun, h, ind, p0, pvph, cj2v)
    
    pciph = 0.d0
    pciph(1:ndim) = pvph 
    pciph(ndim+1 : npvar : ndim+1) = 1.d0 
    
    pcimh = pciph;  pcimh(1:ndim) = pvmh
    
    call poinc(pcimh, ndim, npvar, ind, p0,  dir, imax, tmax, 0, 6, tfmh, pcfmh, hminim, ispc, deriv, gr_cj)
    
    call poinc(pciph, ndim, npvar, ind, p0,  dir, imax, tmax, 0, 6, tfph, pcfph, hminim, ispc, deriv, gr_cj)
    
    dp_dx_diff(:,  i) =  ( pcfph(ind_fun) - pcfmh(ind_fun) ) / 2 / step 
    print*, 'final state after poinc, tf = ', tfph 
    write(*, '(4f18.8)') pcfph(1:ndim); !print*
    read*
  end do
  
  print*, 'check d P / d x by center difference:  ' 
  ! TODO: something is wrong here. 
  do i =  1, n, 1
    write(*, '(10f22.14)') dp_dx_diff(i, :)
  end do
  det = dp_dx_diff(1, 1)*dp_dx_diff(2, 2) - dp_dx_diff(1, 2)*dp_dx_diff(2, 1)
  print*, 'det (D P / D A)  =', det 
  print*; read* 
endif
  
!   check d Pi / d x by diffpoinc:  step = 1.d-6 
!   0.0000000000000000       0.58731372614356747      -0.18720114152169121       0.86947872443892493     
!   0.0000000000000000       -9.4576699639183026E-002 -0.29319804621519063        1.4137373025188946     

!!   check d Pi / d x by difference:  
!   0.0000000000000000       0.58799422220090491      -0.18679309887748019       0.86937764565142572     
!   0.0000000000000000       -9.2865562474273022E-002 -0.29250986260508788        1.4138317805043044 
  !------------------------------ ckd -----------------------------------------
   
  !    d Pi / d \hat P (x,y,vx, vy)_f  *  d P(x,y,vx, vy)_f / d (x,y,vx, vy)_0  
  ! =  d Pi / dxf  * dp_dx , dimension  n - by - ndim 
  !   where  d Pi / dx = [0 1 0 0; 0 0 0 1]  
   
  !  so it could be computed by seleciting the corresponding rows  of dp_dx 
  
  ! dimension n-by-n, n-by-ndim X ndim-by-n
  dp_dx = matmul(dpi_dx, dgamm_dx)
  
if(debug == 2) then  !--ckd, put condition to 2
  print*, 'dp_dx, dimension:', n, ' X ', n
  do i = 1, n, 1
    write(*, '(10f22.14)') dp_dx(i, :)
  end do
  
  det = dp_dx(1, 1)*dp_dx(2, 2) - dp_dx(1, 2)*dp_dx(2, 1)
  print*, 'det (D P / D A)  =', det 
  print*; read*  
  
  print*, 'dp_dx - dp_dx_diff' 
  do i = 1, n, 1
    write(*, '(10f22.14)') dp_dx(i, :) - dp_dx_diff(i, :)
  end do
  print*; read*
    
endif 
           
  ! \varphi( xi+rho ) = pv0_rho
  call varphi( xi+rho, c0, ck, sk, pv0_rho) ! the truncated Fourier series 
  
  ! the error to be corrected by Newton Method:  \varphi(xi+rho) - P( \varphi(xi) ) = 0
  ! for the components (y, vx, vy)
  ! first deal with the columns that are related to A_T for the phase components 

  ferr_xi1 = pv0_rho - pv_pc(ind_fun)  
 
 !  --- compute d \varphi / dA  ---------
  !  d \varphi(xi) / d A 
  call dvphi_da_cal( xi, dvphi_da)
  

if(debug == 1) then ! -- ckd -- put condition to 2
  ! print dvphi_da on screen to compare with the one computed by center difference 
  print*, 'pv0 =', pv0
  print*, 'pv0_rho = ', pv0_rho 
  print*, 'ferr_xi1 = ', ferr_xi1  
  read*
 
  ! ----------- check d \varphi / dA -----------------------------
  !  ( \varphi(A_i + h) -  \varphi(A_i + h) ) / h 
  dvphi_da_diff = 0.d0
  
  do i = 1, n
    do k = 1, nf2p1, 1
      c0ph = c0; ckph = ck; skph = sk 
      c0mh = c0; ckmh = ck; skmh = sk
      
      ! c0, ck, sk 
      if(k == 1 ) then 
        c0ph(i) = c0(i) + step
        c0mh(i) = c0(i) - step
        
      elseif( k <= nf+1) then 
        j = k -1
        ckph(i, j) = ck(i, j) + step
        ckmh(i, j) = ck(i, j) - step
      else 
        j = k - nf - 1 
        skph(i, j) = sk(i, j) + step
        skmh(i, j) = sk(i, j) - step
      endif 
    
      call varphi( xi, c0ph, ckph, skph, funph) 
      call varphi( xi, c0mh, ckmh, skmh, funmh) 
      
      dvphi_da_diff(:, (i-1)*nf2p1 + k) = (funph - funmh) / step / 2
      
      ! to compute d G / d A by finite difference, we compute \varphi(xi+rho) as well 
      call varphi( xi+rho, c0ph, ckph, skph, funph_rho)
      call varphi( xi+rho, c0mh, ckmh, skmh, funmh_rho)
      
      ! -------- check d P / d A ----------------
      pciph = 0.d0; pcimh = 0.d0
      
      ! -- x+h
!      subroutine gamm( pvin, n, ind_fun, h0, ind, p0, pv, cj2v)
      call gamm(funph, n, ind_fun, h, ind, p0, pvph, cj2v)
      pciph(1:ndim) = pvph 
      pciph(ndim+1 : npvar : ndim+1) = 1.d0 
    
      call poinc( pciph, ndim, npvar, ind, p0,  dir, imax, tmax, 0, 6,  &
                  tfph, pcfph, hminim, ispc, deriv, gr_cj)
      
      ! -- x- h
      call gamm(funmh, n, ind_fun, h, ind, p0, pvmh, cj2v)
      pcimh(1:ndim) = pvmh 
      pcimh(ndim+1 : npvar : ndim+1) = 1.d0 
      
      call poinc( pcimh, ndim, npvar, ind, p0,  dir, imax, tmax, 0, 6, & 
                  tfmh, pcfmh, hminim, ispc, deriv, gr_cj)
      
      ! --  check d P / d A by center difference----------------
      dp_da_diff(:,  (i-1)*nf2p1 + k ) =  ( pcfph(ind_fun) - pcfmh(ind_fun) ) / step  / 2
      
      ! -- check d G / d A by center difference ----
      dg_da_diff(:, (i-1)*nf2p1 + k ) = (funph_rho-funmh_rho  - pcfph(ind_fun)+pcfmh(ind_fun) ) / step / 2
 
    enddo   
  end do
  
  ! -- d vphi / d A by routine  -- ckd, put into comment
!  print*, 'check dvphi_da, dimension', size(dvphi_da), 'should be', n, n*nf2p1
!  
!  do i = 1, n !ck
!    write(*, '(10f16.8)') dvphi_da(i, :)
!  enddo   
!  read*  
!  
!  ! --- check d vphi / d A  ---
!  print*, 'check d vaphi / d A by difference:  '  ! -- ckd 
!  do i =  1, n, 1
!    write(*, '(10f16.8)'), dvphi_da_diff(i, :)
!  end do
!  print*; read*
!  
!  print*, 'check ddvphi_da_diff - dvphi_da:  '  ! -- ckd 
!  do i =  1, n, 1
!    write(*, '(10f16.8)'), dvphi_da_diff(i, :)-dvphi_da(i, :)
!  end do
!  print*; read*
  
  
!   check dvphi_da by dvphi_da_cal 
!   1.00000000     -0.80901699      0.30901699      0.58778525     -0.95105652      
!   0.00000000      0.00000000      0.00000000      0.00000000      0.00000000
!   
!   0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      
!   1.00000000     -0.80901699      0.30901699      0.58778525     -0.95105652

! check d vaphi / d A by difference:  
!   1.00000000     -0.80901699      0.30901699      0.58778525     -0.95105652      
!   0.00000000      0.00000000      0.00000000      0.00000000      0.00000000
!  
!   0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      
!   1.00000000     -0.80901699      0.30901699      0.58778525     -0.95105652
  ! ----------------------------------------------------------
  
  
  ! --- check d P  / d A ----
  print*, 'check d P  / d A by difference:  ' 
  do i =  1, n, 1
    write(*, '(12f18.8)')  dp_da_diff(i, :)
  end do
  print*; read* 

! step = 1.d-5   
! check d P  / d A by difference:  
!        0.10964325        0.10964325        0.10964325        0.00000000        0.00000000        
!        1.12030086        1.12030086        1.12030086        0.00000000        0.00000000
!       -0.72837197       -0.72837197       -0.72837197        0.00000000        0.00000000        
!        1.69524679        1.69524679        1.69524679        0.00000000        0.00000000


! dp_da = dp_dx X  dvphi_da
!        0.10845144        0.10845144        0.10845144        0.00000000        0.00000000        
!        1.11915660        1.11915660        1.11915660        0.00000000        0.00000000
!       -0.72934034       -0.72934034       -0.72934034        0.00000000        0.00000000        
!        1.69434312        1.69434312        1.69434312        0.00000000        0.00000000
 
endif ! if debug == 1   
 
  ! d P / dA  =  d pv_pc{y,vx,vy} / d \varphi(xi) * d \varphi(xi) / d A
  !           = Differential of poincare map  * dvphi_da
  
  ! dimension:  n-by-n  X  n-by-n*nf2p1
  dp_da = matmul(dp_dx, dvphi_da) 
      
  !  d \varphi( xi+rho ) / d A 
  call dvphi_da_cal( xi+rho, dvphi_rho_da)
  
  dg_xi1 = dvphi_rho_da - dp_da 
  
  
! --- check d G  / d A ----
if(debug == 1) then 
 
  print*, 'dp_da = dp_dx X  dvphi_da' 
  do i = 1, n, 1
    write(*, '(12f18.8)')  dp_da(i, :)
  enddo 
  print*;  read* 
  
  print*, 'dp_da_diff - dp_da' 
  do i = 1, n, 1
    write(*, '(12f18.8)')  dp_da_diff(i,:)-dp_da(i, :)
  enddo 
  print*;  read*
  read*
  
  print*, 'check dvphi_rho_da, dimension: ', size(dvphi_rho_da) !ck
  do i = 1, n !ck
    write(*, '(12f18.8)') dvphi_da(i, :)
  enddo   
  print*; read*  
 
 
  print*, 'check d G  / d A by difference:  ' 
  do i =  1, n, 1
    write(*, '(12f18.12)')  dg_da_diff(i, :)
  end do
  print*; read* 
 
  print*, 'd G / d A = dg_dx1 = dvphi_rho_da - dp_da ' 
  do i = 1, n, 1
    write(*, '(12f18.12)')  dg_xi1(i, :)
  enddo 
  print*; read* 
  
  print*, 'dg_da_diff - dg_da' 
  do i = 1, n, 1
    write(*, '(12f18.8)')  dg_da_diff(i,:) - dg_xi1(i, :)
  enddo 
  print*;  read*
  read*
   
endif  
 
!  print*, 'd G / d A = dg_dx1 = dvphi_rho_da - dp_da ' 
!  do i = 1, n, 1
!    write(*, '(12f18.12)')  dg_xi1(i, :)
!  enddo 
!  print*; read* 
  
  return  
end subroutine gdg_inv_xi1 


!***********************************************************************
!     ****   varphi   ****
!  Compute the truncated Fourier series of the invariant curve \varphi(xi) 
!  given the coefficients and initial phases

!       Input Variables 
!  xi           the argument for the invariant curve \varphi 
!  c0, ck, sk   dimension 6-by-nf,  the unknowns C_k and S_k (k=1, ..., nf), for each component
!               (x,y,z,vx,vy,vz), the values are different. 

!       Output Variables 
!  pv           the state vector: position + velocity            

!      Module-based Varaibles
!  n, nf  
! 
!  routine Used:
!     trigrec_c, trigrec_s

! TODO: probably we will not need to compute pv using a standalone subroutine 
!       because at the point when pv is needed, we have already computed  cn, sn 

!  Finally Revised by Yu -- 20160730
!***********************************************************************
subroutine varphi( xi, c0, ck, sk, pv)
implicit none

! Input  and Output Declaration 
real(kind=dp), intent(in)    :: xi, c0(n), ck(n, nf), sk(n, nf)   
real(kind=dp), intent(out)   :: pv(n) 
 
! Local Variable
integer :: i, k 
real(kind=dp)  :: cn(nf), sn(nf)
  
 ! compute the trigometric function at k*xi 
  call trigrec_c( xi, nf, cn)  
  call trigrec_s( xi, nf, sn)
  
  do i = 1, n, 1
    pv(i) = c0(i)
    do k = 1, nf, 1
       pv(i) = pv(i) + ck(i,k) * cn(k) + sk(i,k) * sn(k)
    end do
  end do  

  return  
end subroutine varphi


!***********************************************************************
!     ****   dvphi_da   ****
!   for a certain value of xi, the Jacobi Matrix of the invariant curve \varphi 
!   w.r.t. the  Fourier coefficients A.  
!   
!   so dvphi_da is a matrix with ndim same block in diagonal line 
!   we store by Fourier modes, k = 0, 1, ..., nf 

!   each block is a row vector of dimension 1 + 2Nf: 
!   1 ,  cos(xi)...  cos(nf*xi),  sin(xi) ... sin(nf*x) 

!   we denote A = [ C0, C's, S's] as the Fourier coefficients, then A is of dimenson n * (2*Nf+1), write as components A = [A_x, A_y, A_z, A_vx, A_vy, A_vz] 
  
!  where \varphi as a truncated Fourier series.
!        \varphi (xi) = C_0 + Sigma_ k from 1 to Nf of { ( C_k cos(k*xi) + S_k cos(k*xi)) }
  
!       Input Variables 
!  xi       the argument for the invariant curve \varphi 
                    
!       Output Variables 
!  dgda      the differential (Jacobi Matrix) of \varphi w.r.t. the Fourier coefficients A            
 
!       Module-based Varaibles
! nf, nf2, nf2p1

!  routine Used:
!   trigrec_c, trigrec_s 

! TODO: check this routine, to make sure every component of the differential is the right
!       one thing is use finite difference   

!  Finally Revised by Yu -- 20160729
!***********************************************************************
subroutine dvphi_da_cal( xi, dgda)
implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration   
real(kind=dp), intent(in)      ::  xi  
real(kind=dp), intent(out)     ::  dgda(n, n*nf2p1) 
 
! Local Variable
integer        :: i, j, k, jcol
real(kind=dp)  :: cn(nf), sn(nf)
  
!  if(debug == 2) then  !ckd, put condition to 2
!    print*;  print*, '-------- check dvphi/da ----------'
!    print*, 'xi = ', xi
!  endif   

  dgda = 0.d0 ! initialization
  
 ! compute the trigometric function  xi 
  call trigrec_c( xi, nf, cn)  
  call trigrec_s( xi, nf, sn)
  
  ! The Jacobi Matrix of the general function F (the n coordinates)  w.r.t. the  Foureier coefficients A 
  !  each is of dimension n X nf2p1 and are of the same block in the diagonal line, the rest components are zeros
  
  do i = 1, n, 1
  
    ! the starting  index of i-th block  for i-th component
    jcol = (i-1)*nf2p1 + 1  
     
    ! d \varphi / d C_0 = 1
    dgda( i,  jcol) = 1.d0       
     
    ! the 2*nf components are cos(k*xi) for k =1, ...,  nf, for jcol+1 : jcol+nf columns 
    ! sin(k*xi) for k =1, ...,  nf, for jcol+nf_1 : jcol_+2nf columns 
     
    dgda( i,  jcol+1 : jcol+nf)       =  cn  ! cos( k*xi ), k = 1_nf
    dgda( i,  jcol+nf+1 : jcol+2*nf)  =  sn  ! sin( k*xi ), k = 1_nf
     
    if(debug == 2) then  ! ckd -- put condition to 2
      print*, 
      print*, 'cn = '
      write(*, '(12f18.8)') cn
     
      print*, 'sn='
      write(*, '(12f18.8)') sn 
    endif   
     
  end do
 
  return  
end subroutine dvphi_da_cal

end module curve_poinc_mod
