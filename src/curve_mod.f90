!  **** ** The common routines for both method *******
!            curve_time_mod and curve_poinc_mod 
!  including 
!  -- 1. subroutine init_curve(n0, nitmax0, tol0, tol_err0, opt0) 
!  -- 2. subroutine init_nf(nf0) 
!  -- 3. subroutine gdg_da_xi1( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, gr_map)
!  -- 4. subroutine upd_four(isupd, nf0) 
!  -- 5. subroutine alloc_arr_2d(n0, nf0, ck0) 
!  -- 6. subroutine varphi( xi, n, nf, c0, ck, sk, pv) 
!  -- 7. subroutine dvphi_da_cal( xi, dgda)
 
!  -- 8. subroutine rot_fcs( ck, sk, ind_cs1, is_ck, n, nf)
!        rotate ann angle \xi_0 to avoid phase-shift indetermination 
!        to set one component c1 or s1 to be zero 
 
!  -- 9. subroutine fix_index(ck, sk, n, nf, ind_cs1, is_ck)
!        choose the index of c1 or s1 to be set to zero 
   
 
!  ** Note ** 
!    1. we have 2 more equations than unknows, but Josep Maria claims the kernel dimension is zero... 
!    2. use the general name for vector field and energy (Jacobi costant)
!    3. subroutines from libraries, need to declare external attribute **
!         external  ::  deltx, trigrec_c, trigrec_s    ! -- libgr.a 

!**********************  init_opt ************************************
! assignment of the approach by the value of opt 
subroutine init_opt( opt0)
  implicit none
  integer, intent(in)    ::  opt0 

  opt   = opt0 
  return
end subroutine init_opt

!**********************  init_cont ************************************
! the upper and lower bound for ds and nf 
subroutine init_cont(nf_min0, nf_max0, ds_min0, ds_max0)
!call init_cont(nf_min, nf_max, ds_min, ds_max)
  implicit none
  integer, intent(in)        :: nf_min0, nf_max0
  real(kind=dp), intent(out) :: ds_min0, ds_max0 
  
  nf_min = nf_min0;  nf_max = nf_max0
  ds_min = ds_min0;  ds_max = ds_max0
  return
end subroutine  init_cont
 
!**********************  init_curve ************************************
! The initialization for control  parameters of Newton method, and which approach to take 
! in fact, we will always choose opt = 3, where we remove one unknowns to get rid of the phase shift unknowns  
 
subroutine init_curve(n0, m0, nitmax0, tol0, tol_err0 ) 

implicit none
integer, intent(in)          :: n0, m0, nitmax0 
real(kind=dp), intent(in)    :: tol0, tol_err0 
logical :: isopen  

  debug = 0
!  print*, 'Debug (1) or not for the module curve_mod?'
!  read(*,*)  debug
  
  ! the dimension of the map, for Poincare map method n = ndim-2,
  ! for Time T map method, n = ndim  
  n     = n0 
  m     = m0 ! multiple shooting 
  nm    = n*m 
  
  ! termination control for Newton method 
  nitmax  = nitmax0 ! maximum time of iterations allowed
  tol     = tol0
  tol_err = tol_err0
    
!  for save - TODO: for the moment, not used
   fwork = 6
!  fwork = 222 
!  open(fwork, file='curve_poinc.work')
 
  fcurve_approx = 888
  inquire(file='curve_approx.dat', opened=isopen)
  if(isopen) then 
    print*, 'curve_approx.dat already opened!'; print*
  else 
     open(fcurve_approx, file='curve_approx.dat', status='replace',access='append')
  endif    
  
end subroutine init_curve 

!***********************************************************************
!     ****   init_nf   ****
!  the dimension of related arrays:  nf, nf2, nf2p1 
!  after the refinement of each curve, the value of nf will probably be updated, 
!  so as all the related parameters. 

!  --- we switch to Time-T map method 
!  10 -- fix T and rho, let h be free, unknonws: c0, ck, sk


! -- opt, the strategy for Poincare map method 
!  1 -- add rho as unknonw 
!     Unknowns: rho, c0, ck, sk, and remove c1_k or s1_k  
!    dimension: (n * 2nfp1 + 1) X  ( n*2nfp1 )

!  2 -- add energy as unknonw 
!     Unknowns: h0, c0, ck, sk, and remove c1_k or s1_k  
!    dimension: (n * 2nfp1 + 1) X  ( n*2nfp1 )

!  3 --- the current used one  for Poincare map method 
!     Unknowns: c0, ck, sk, and remove c1_k or s1_k, 
!    dimension: n * 2nfp1 X  (n*2nfp1 - 1)

!***********************************************************************
subroutine init_nf(nf0) 

implicit none
integer, intent(in)  :: nf0
  
  nf    = nf0
  nf2   = 2*nf  
  nf2p1 = nf2 + 1

! ---- commonly used value for dimension declaration  -----------  
    
  ! number of columns and rows in dg,  which is the differential of the linear equation 
  !  F:  \varphi(xi+rho) - P( \varphi(xi) ) = 0 
  !      dimension:   n * nf2p1  X  n*nf2p1 
  na1      = n*nf2p1
  na_all   = m*na1 
  
  nr_dgall = na_all + 1   ! inv eq.s + energy eq.s
  nc_dgall = na_all + 3   ! Fourie coefs + three unknonws rho, tp, h0 
  
!  print*, 'nr_dgall, nc_dgall',  nr_dgall, nc_dgall; print*; read*
   
  
  ! if we want to add more equations or remove certain unknowns, change the value 
  ! of nr_dg and nc_dg. 
  
  ! **NOTE** we always take opt=3 
  
  ! opt = 1, 2, 3 are for Poincare map method 
  ! opt = 11, 12 are for Time Map method 
  if( opt == 1) then ! 
 !  TODO :keep constant 2 unknown,  --- perhaps we can remove this one ...
 ! that is non zero c1,or s1, to make it simple, we take the one with ind_fxd, two much  
   
    nr_dg = na_all   
    nc_dg = na_all - 2 

! ---- avoid phase shift indetermination -------- 
 ! Assume \varphi(xi) gives us an invariant curve, \varphi(xi+ alpha) will also do 
 ! we fix C1(or S1)_J = 0, see more details about the idea in subroutine rot_fcs
 ! for this purpose, we have two approaches, 
 !  --1, add one more equation, C1(S1)_J = 0
 !  --2, remove the unknown C1(S1)_J after the call of rot_fcs to get C1(S1)_J = 0
 !  
 ! I prefer the second approach, since the first one will produce small correction on C1(S1)_J
 ! and we have to call rot_fcs in every iteration .
  
  
  elseif( opt == 2) then 
  ! only fix one value in non-zero  C1, S1 , avoid phase shift and falling back to p.o. at the same time 
  
    nr_dg = na_all  
    nc_dg = na_all - 1
  
  elseif( opt == 3) then  ! -- only use this one  
    ! -- approach 3 : remove one unknows to make  c1_1  = 0 to avoid phase shift indetermination 
    ! so we have one less unknonws, and the same number of equtaion
       
    nr_dg = na_all
    nc_dg = na_all - 1
  
  ! ---------------- Time T map approach -------------------------------  
  elseif(opt == 10 .or. opt == 20) then 
    ! only Foureier coefs + invariance eq. 
    ! fix c0(s0)_k to avoid curve indetermination
    ! rotate to set c1(s1)_k to zero to avoid phase shift indetermination 
    ! fix c1(s1)_k to avoid falling back to p.o. (all ck,sk with k>=1 is zeros)
    ! -- we have to add one more unknonw, here, we take rho
    !    but if we also let h free, we will fall back to p.o. with rho=0 
    ! -- add one more unknonw TP, and fix h by adding the energy equation. 
    
    nr_dg = na_all + 1
    nc_dg = na_all - 1 
    
    ! length of the new guess for the next torus used in dg_pre_time 
    ! but for convenience, we return dx_pre as a full matrix including all the unknonws except h0
    ! with the c0_k to be removed  as zero 
     
    len_dx_pre = na_all + 1
    call alloc_arr_1d(len_dx_pre+1, dx_pre)
    
!    print*, 'size(dx_pre) = len_dx_pre+1', size(dx_pre), len_dx_pre !ckd 
    
    ! index of c0, ck, sk used, 3 to be eliminated  
    len_ind_dg     = na_all - 3 
    
    ! two more than len_ind_dg, deal with the indeterminations in a different way w.r.t. the correction
    ! 1 removed to avoid curve indetermination 
    len_ind_dg_pre = na_all - 1   
 
  elseif(opt == 11 .or. opt == 21) then 
    ! JM's method for the continuation, and the system is rank-deficient, we use QR fracterization with column pivot
    ! opt = 11 : fix the energy, keep rho and tp as unknonws
    !            one more column than opt=10, since we didn't need to consider falling back to p.o. 
    nr_dg = na_all + 1
    nc_dg = na_all  
    
    ! length of the new guess for the next torus, all c0,ck,sk + rho, tp  but with c0_k=0    
    len_dx_pre = na_all + 1
    call alloc_arr_1d(len_dx_pre+1, dx_pre)
    
    ! index of c0, ck, sk used, 2 to be eliminated  c0_k and c1_k, no falling back to p.o. is necessary 
    len_ind_dg     = na_all - 2 
    
    ! two more than len_ind_dg, deal with the indeterminations in a different way w.r.t. the correction
    ! 1 removed to avoid curve indetermination, h0 removed 
    len_ind_dg_pre = na_all - 1
    
    ! -- ckd 
!    print*, 'size(dx_pre) = len_dx_pre+1 = len_ind_dg_pre+2?? ', size(dx_pre), len_dx_pre+1, len_ind_dg_pre+3
        
     
  elseif(opt == 100) then 
    ! Alex's suggestion, continue along tp, first fix tp and rho, free h 
    nr_dg = na_all 
    nc_dg = na_all - 2 
    len_ind_dg = nc_dg  ! index of c0, ck, sk used 
     
  endif 
    
  ! --- allocate memory for ind_dg and ind_dg_pre ---
  call alloc_arr_1d0(len_ind_dg,     ind_dg)
  call alloc_arr_1d0(len_ind_dg_pre, ind_dg_pre)
  
  ! --------------------------------------------------

!  print*, 'number of rows and columns in DG, Dg_ore ' 
!  print*, 'na_all, nr_dg, nc_dg, len_ind_dg_pre:', na_all, nr_dg, nc_dg, len_ind_dg_pre
!  print*, 'size(ind_dg), nc_dg', size(ind_dg), nc_dg; 
!  print*; read*

  return  
end subroutine init_nf


!***********************************************************************
!     ****   gdg_da_xi1  ****
!  for a certain value of xi, the Jacobi Matrix of the invariance equation w.r.t. c0,ck,sk
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
!!!  im           index of the segment in multiple shooting 
!  xi           the argument for the invariant curve  \varphi 
!  rho          rotation number:  2pi / omega2 * omega1
!  c0, ck, sk   dimension ndim-by-nf,  the unknowns C_k and S_k (k=1, ..., nf), for each component
!               (x,y,z,vx,vy,vz), the values are different. 

                    
!       Output Variables 
!  dg_xi1      dimension: n(4 or 6 ), m * n * ( 2*nf+1 ) 
!              the differential (Jacobi Matrix) of invariance equations w.r.t. A_{x,y,z,vx,vy,vz}
 
! ferr_xi1     dimesion n, the error to be refined by Netwon method 
!              \varphi_j+1 (xi_i + rho)    - P ( \varphi_j (xi_i) ), j = 0, m-2     
!              F   = \varphi_0(xi_i + rho) - P ( \varphi_m-1 (xi_i) )  

!     Module based Varaibles
! n, nf, c0, ck, sk

!  routine Used:
!   trigrec_c, trigrec_s, varphi, dvphi_da_cal

!  Finally Revised by Yu -- 20170612 -- ckd using Richardson extrapolation 
!***********************************************************************
subroutine gdg_da_xi1( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, gr_map)
implicit none

! Input  and Output Declaration   
!integer, intent(in) :: im 
real(kind=dp), intent(in)      ::  xi, rho, c0(nm), ck(nm, nf), sk(nm, nf) 
real(kind=dp), intent(out)     ::  dg_xi1(nm, na_all), ferr_xi1(nm) 

! External subroutines                  
external ::  gr_map 


! Local Variable
integer        :: i, k,  debug, im, row_st, col_st    
real(kind=dp)  :: fun0(n), fun1(n), pv0_rho(n),             &  ! the initial state 
                  pvf(n), dp_dx(n, n), dp_da(n, na1),       &  ! map 
                  dvphi_da(n, na1), dvphi_rho_da(n, na1)   ! for dg_xi1
                  

! -- debug 
real(kind=dp) ::  step,  dvphi_da_h(n,na1), dvphi_da_h2(n,na1), dvphi_da_diff(n,na1),  & 
                  dp_da_h(n,na1), dp_da_h2(n,na1),  dp_da_diff(n, na1),    & 
                  dg_da_h(nm,na_all), dg_da_h2(nm,na_all),  dg_da_diff(nm,na_all)     

! dvphi_da_diff(n,na1), dp_da_diff(nm,na_all), dg_da_diff(nm,na_all) 

  debug = 0      ! -- ckd! 2017-06-20 with m = 2
  step  = 1.d-4 
  
  dg_xi1 = 0.d0 ! initialization
  
  
   ! -- im == 0 --   \varphi_0( xi+rho ) = pv0_rho, and dvphi_rho_da^0 
   call varphi( xi+rho, n, 0, nf, c0, ck, sk, pv0_rho) ! the truncated Fourier series 
  
   !  d \varphi( xi+rho ) / d A^0    
   call dvphi_da_cal( xi+rho, dvphi_rho_da)
    
  ! ---- The Jacobi Matrix of the general function F (the ndim coordinates) w.r.t. A---
  do im = 0, m-1, 1 
    ! for im-th block 
    row_st = n*im 
    col_st = na1*im
    
    ! \varphi_0(xi), the first segment, and the rotation...  varphi_0(xi+rho)
  
    ! Evaluate  \varphi_im (xi) for im-th segment 
    call varphi( xi, n, im, nf, c0, ck, sk, fun0)
    
    ! use general routines to compute the map and its differential 
    !  subroutine std_map( ptin,  ptf, isdp,  dp)  isdp = 1, compute the differential             
    
    call gr_map( fun0,  pvf, 1, dp_dx)    
  
    if(debug == 1) then  ! -- ckd 
      print*, 'Finish gr_map! dp_dx, dimension:', n, ' X ', n
      do i = 1, n, 1
        write(fwork, '(10f22.14)') dp_dx(i,:)
      end do
    endif 
    
    !  --- compute d \varphi / dA  ---------
    !  d \varphi(xi) / d A 
    call dvphi_da_cal( xi, dvphi_da)   
    
    ! dvphi_da, dvphi_rho_da are only dependent on xi, and are independent on C0,CK,SK 
  
    ! d P / dA  =  d P / d \varphi   *   d \varphi  / d A
    !           = Differential of the map  * dvphi_da
  
    ! dimension:  n-by-n  X  n-by-n*nf2p1
   
    dp_da = matmul(dp_dx, dvphi_da)   ! for one xi, and \varphi_im 
     
    
    if( im < m-1 ) then  ! m has to be greater than 1
      ! -- m > 1 .and.  im < m-1
      ! Evaluate \varphi_im+1 (xi)
      call varphi( xi, n, im+1, nf, c0, ck, sk, fun1) 
    
      ! the error to be corrected by Newton Method:  \varphi_im+1 (xi) - P( \varphi_im (xi) ) = 0
      ferr_xi1(row_st+1 : row_st+n)       = fun1 - pvf
     
      dg_xi1(row_st+1 : row_st+n, col_st+1     : col_st+na1)     =  -dp_da    ! dg_da^im
      dg_xi1(row_st+1 : row_st+n, col_st+na1+1 : col_st+2*na1)   =  dvphi_da  ! dg_da^(im+1)
       
        
    elseif(m == 1) then 
      
      !  (m == 1)  ! \varphi_0(\xi+rho) -  P( \varphi_0 (xi) )
      dg_xi1(row_st+1 : row_st+n, 1 : na1)  = dvphi_rho_da - dp_da   
   
   
   ! dg_da to debug 
    else    
      ! (im == m-1) the last segment with m > 1 
      
      ! dg_da^0, 1:na1 columns  the last segment  \varphi_0(\xi+rho) -\phi_t ( \varphi_m-1 (xi) )
      !   if m = 1, both items are assoicated to A^0 
      !   if m > 1, this item only include  \varphi_0(\xi+rho)
      dg_xi1(row_st+1 : row_st+n, 1 : na1)                =  dvphi_rho_da
       
      dg_xi1(row_st+1 : row_st+n, col_st+1 : col_st+na1)  = -dp_da 
        
    endif 
    
    ! the last segment, invariance equation
    ferr_xi1(nm-n+1 : nm)       = pv0_rho - pvf
 
  enddo 
  
  
    !---------------------------- debug ----------------------------
    if(debug == 1) then   ! -- ckd! 2017-06-12 20:41:55 
      ! print dvphi_da on screen to compare with the one computed by center difference 
      print*, 'pv0      = ', fun0
      print*, 'pv0_rho  = ', pv0_rho 
      print*, 'ferr_xi1 = ', ferr_xi1  
      read*
      
      row_st = 30
!     subroutine dg_da_debug(step, gr_map, dvphi_da_diff, dp_da_diff, dg_da_diff)
! dimension:  dvphi_da_diff(n,na1), dp_da_diff(nm,na_all), dg_da_diff(nm,na_all) 
      call dvphi_da_debug(xi,  m-1,  step/2.d0, gr_map, dvphi_da_h,  dp_da_h)
      call dvphi_da_debug(xi,  m-1,  step,      gr_map, dvphi_da_h2, dp_da_h2)
      
      
      call dg_da_debug(xi, step/2.d0,   gr_map,  dg_da_h )
      call dg_da_debug(xi, step,        gr_map,  dg_da_h2)
      
      dp_da_diff    = (4.d0*dp_da_h - dp_da_h2) / 3.d0
      dg_da_diff    = (4.d0*dg_da_h - dg_da_h2) / 3.d0      
      dvphi_da_diff = (4.d0*dvphi_da_h - dvphi_da_h2) / 3.d0
      
      
      ! -- d vphi / d A by routine  -- ckd, put into comment
      print*, 'check dvphi_da, dimension', size(dvphi_da), 'should be', n, n*nf2p1
  
      do i = 1, nm !ck
        write(fwork, '(10f16.8)') dvphi_da(i, 1:row_st)
      enddo   
      read*  
  
      ! --- check d vphi / d A  ---
      print*, 'check d vaphi / d A by centre difference:  '  ! -- ckd 
      do i =  1, n, 1
        write(fwork, '(10f16.8)'), dvphi_da_diff(i, 1:row_st)
      end do
      print*; read*
  
      print*, 'check dvphi_da_diff - dvphi_da:, maxnorm', maxval(dabs(dvphi_da_diff-dvphi_da))  ! -- ckd 
      do i =  1, n, 1
        write(fwork, '(10f16.8)'), dvphi_da_diff(i, 1:row_st)-dvphi_da(i, 1:row_st)
      end do
      print*; read*
  
  
      print*, 'check d P  / d A by difference:  ' 
      do i =  1, n, 1
        write(fwork, '(10f18.8)')  dp_da_diff(i, 1:row_st)
      end do
      print*; read* 
 
      print*, 'dp_da = dp_dx X  dvphi_da' 
      do i = 1, n, 1
        write(fwork, '(10f18.8)')  dp_da(i, :)
      enddo 
      print*;  read* 
   
      print*, 'dp_da_diff - dp_da, maxnorm: ', maxval(dabs(dp_da-dp_da_diff)) 
      do i = 1, n, 1
        write(fwork, '(10f18.8)')  dp_da_diff(i,1:row_st)-dp_da(i, 1:row_st)
      enddo 
      print*;  read*
      
      ! --- check d G / D A -----
      print*, 'dg_da by centre difference' 
      do i = 1, nm, 1
        write(fwork, '(10f18.8)')  dg_da_diff(i, 1:m*row_st)
      enddo 
      print*;  read*
      
      print*,'shape dg_da_diff, dg_xi1:', shape(dg_da_diff), shape(dg_xi1)
      
      print*, 'd G / d A = dg_dx1' 
      do i = 1, nm, 1
        write(fwork, '(10f18.12)')  dg_xi1(i, 1:m*row_st)
      enddo 
      print*; read*
      
      print*, 'diff in d G / d A, maxnorm:', maxval(dabs(dg_xi1-dg_da_diff)), & 
      maxloc(dabs(dg_xi1-dg_da_diff)) !, dg_xi1(maxloc)-dg_da_diff(maxloc)
      
      do i = 1, nm, 1
        write(fwork, '(10f18.12)') dg_xi1(i,:)-dg_da_diff(i,:)  ! dg_xi1(i, 1:m*row_st)-dg_da_diff(i, 1:m*row_st)
      enddo 
      print*; read*
      
    endif ! if debug == 1   
! ---------- debug dp_da by centre difference -------------------------
  
  return  
end subroutine gdg_da_xi1   

!***********************************************************************
!     ****   upd_four   ****
!  Decide the approriate value of nf, and declare the Fourier coefficients CK, SK, 
!  and update the array size (if necessary) for Poincare map method. 

! The criteria is: 
!  the maximum norm of the last 1/4 Fourier Coef should be of one order 
!  of magnitude less than the tolrence, < tol / 1.d1, otherwise, we double nf=2nf 
!  and set the new ones as  zeros  
! 
!  ****** Input  Variables ******
!  isupd           returned by check_tail, othervalues than 0,-1,1: remain;   1: double; -1: half  
!                  if we want to initialize ck,sk, isupd = 0 

 
!  ****** Input and Output Variables ******
!  nf0              the old value of number of harmonics, to be update  is isupd=1
!                   ** NOTE **  Remerber to update nf0! since it is not a global value 


!  Module-base Varaibles:
!  private: n, tol 
!  public: c0, ck, sk 

!  Finally Revised by Yu -- 20170520
!***********************************************************************
subroutine upd_four(isupd, nf0) 

implicit none
! Input  and Output Declaration  
integer, intent(in)          :: isupd 
integer, intent(inout)       :: nf0
 
! Local Variable
integer        ::  nfmn, i, nf_copy,  debug 
real(kind=dp)  ::  ck_copy(nm,nf0), sk_copy(nm,nf0), c0_copy(nm)

  debug = 0
  nf_copy = nf0    ! in case that nf0 is passed into as nf  
  nf      = nf_copy
  
  
  if(isupd == 0) then 
!    print*, 'The initialization for C0, CK, SK.'; 
!    print*, 'upd_four:, nm, nf', nm, nf; print*; read*
   
    call init_nf(nf)
    
    call alloc_arr_1d(nm, c0)
    call alloc_arr_2d(nm, nf, ck); 
    call alloc_arr_2d(nm, nf, sk); 
    
    c0 = 0.d0; ck = 0.d0; sk = 0.d0
    return
  endif

  if(abs(isupd) /= 1)  return  
   
  ck_copy = ck; sk_copy = sk  
  
  if(isupd == 1) then   
     nf = 2 * nf_copy
  elseif(isupd == -1) then 
     nf = nf_copy / 2
  endif 
  
  
  call alloc_arr_2d(nm, nf, ck); call alloc_arr_2d(nm, nf, sk)
  ck = 0.d0; sk = 0.d0 
         
  nfmn = min0(nf, nf_copy)
  
  ! assignment for public module variables
  ck(:, 1:nfmn) = ck_copy(:, 1: nfmn)
  sk(:, 1:nfmn) = sk_copy(:, 1: nfmn)

  
  if(debug == 1) then 
    write(*, *) '# Updated c0, ck, sk : n, nf = ', n, nf 
    do i = 1, n, 1
      write(*, '(6e24.14)')   c0(i);       ! write(*,*) 
      write(*, '(10e24.14)')  ck(i, 1:10); ! write(*,*) 
      write(*, '(10e24.14)')  sk(i, 1:10) 
     write(*,*); ! write(8,*) ! two blank line2 to seperate components 
    end do
    print*; read*
  
!   commonly used value for dimension declaration 
    print*, 'Update nf for refine_curve_poinc: ', nf; read*
  endif 
  
  call init_nf(nf)
  
  nf0 = nf 
  return  
end subroutine upd_four


!***********************************************************************
!     ****   check_tail   ****
!  Check the tail of ck, sk to decide  if we need to update the number of Foureier modes nf 

! The criteria is: 
!  the maximum norm of the last 1/4 Fourier Coef should be of one order 
!  of magnitude less than the tolrence, < tol / 1.d1, otherwise, we double nf=2nf 
!  and set the new ones as zeros  
! 
!  ****** Output  Variables ******
!  isincr            10: keep the old nf, or any value except 0,-1,1 
!                    1: double it by 2*nf 
!                   -1: decrease it by half 

!  Module-base Varaibles:
!  private: n, tol 
!  public: c0, ck, sk 

!  Finally Revised by Yu -- 20170621
!***********************************************************************
subroutine check_tail(isincr) 

implicit none
integer, intent(out)       :: isincr
 
! Local Variable
integer        :: inf, nfmn, i, nf_copy, debug 
real(kind=dp)  :: ckm, skm, tolb10, ck_copy(nm,nf), sk_copy(nm,nf), c0_copy(nm)
real(kind=dp)  :: dlange 


  isincr =  10  ! by default keep the number of nf
  
  debug = 0
  nf_copy = nf 

  ck_copy = ck; sk_copy = sk  ! ; c0_copy = c0 
  
  ! -- check the Maximum norm of csf and sif, and choose an appropriate value for nf 
  ! FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
  inf = nf_copy / 4
  ckm = dlange('m', n, inf, ck(:, nf_copy-inf+1 : nf_copy), n, 0.d0)
  skm = dlange('m', n, inf, sk(:, nf_copy-inf+1 : nf_copy), n, 0.d0)
  
  print*; print*, 'maximum norm of CK and SK for the last ', inf, 'Fourier harmonics'
  print*,  ckm, skm 
  
  print*; ! read*

  tolb10 = tol / 1.d0 
  if(dmax1(ckm, skm)  < tolb10 ) then
    ! check the maximum norm of the first 1/4-1/2 Four Coefs
    !  if it's still within the tolerance, decrease nf by a half.
    ckm = dlange('m', n, inf, ck(:, inf+1 : 2*inf), n, 0.d0)
    skm = dlange('m', n, inf, sk(:, inf+1 : 2*inf), n, 0.d0)
    
    if(dmax1(ckm, skm)  < tolb10 .and. nf_copy >= nf_min*2.d0)  isincr = -1
    print*, inf+1, ' ---', 2*inf, 'Coefs:', ckm, skm; print*
  
  else 
    ! double the number of harmonics
!    nf = 2 * nf_copy
    isincr = 1
  endif 
  
  return  
end subroutine check_tail



!*******************   alloc_arr_1d0   *********************************
!  allocate 1d array of integer type 
!***********************************************************************
subroutine alloc_arr_1d0( nf0, arr0)  

implicit none
! Input  and Output Declaration  
integer, intent(in)       :: nf0 
integer, allocatable, intent(inout) :: arr0(:)  

  if (allocated(arr0) .eqv. .false.) then
    allocate(arr0(nf0)); ! allocate(sif(nf0))
  else 
    if(size(arr0) /= nf0) then 
      deallocate(arr0);    ! deallocate(sif)
      allocate(arr0(nf0)); ! allocate(sif(nf0)) 
    endif  
  endif
    
  return
end subroutine alloc_arr_1d0

!******************* alloc_arr_1d ***************************************
!  allocate 1d real array  
!***********************************************************************
subroutine alloc_arr_1d( nf0, csf)  

implicit none
! Input  and Output Declaration  
integer, intent(in)       :: nf0 
real(kind=dp), allocatable, intent(inout) :: csf(:)  

  if (allocated(csf) .eqv. .false.) then
    allocate(csf(nf0));  
  else 
    if(size(csf) /= nf0) then 
      deallocate(csf);    
      allocate(csf(nf0));  
    endif  
  endif
    
  return
end subroutine alloc_arr_1d


!***********************************************************************
!     ****   allococate 2D array  of dimention n-by-nf0
!   use 0 to avoid conflicts  with the variables from the module 
!***********************************************************************
subroutine alloc_arr_2d(nr, nc, arr) 

implicit none


integer, intent(in)       ::  nr,  nc 
real(kind=dp), allocatable, intent(inout) :: arr(:,:)  
  
  if (allocated(arr) .eqv. .false.) then
    allocate(arr(nr, nc)); 
  
  else 
    deallocate(arr);     
    allocate(arr(nr, nc));  
  endif 
    
  return
end subroutine alloc_arr_2d


!***********************************************************************
!     ****   varphi   ****
!  Compute the truncated Fourier series of the invariant curve \varphi(xi) 
!  given the coefficients and initial phases

!       Input Variables 
!  im           the subindex of the second map term, im = 0--m-1
!  xi           the argument for the invariant curve \varphi 
!  c0, ck, sk   dimension 6-by-nf,  the unknowns C_k and S_k (k=1, ..., nf), for each component
!               (x,y,z,vx,vy,vz), the values are different. 

!       Output Variables 
!  pv           the state vector: position + velocity            

!  routine Used:
!     trigrec_c, trigrec_s

! TODO: probably we will not need to compute pv using a standalone subroutine 
!       because at the point when pv is needed, we have already computed  cn, sn 

!  Finally Revised by Yu -- 20160730
!***********************************************************************
subroutine varphi( xi, n, im, nf, c0, ck, sk, pv)  
implicit none

! Input  and Output Declaration 
integer, intent(in)          ::  n,im,  nf
real(kind=dp), intent(in)    :: xi, c0(nm), ck(nm, nf), sk(nm, nf)   
real(kind=dp), intent(out)   :: pv(n) 
 
! Local Variable
integer :: i, k, row_st  
real(kind=dp)  :: cn(nf), sn(nf) 
  
  row_st = im*n
  
!  print*, 'm, im, row_st: ', m, im, row_st;  print*; read*

 
 ! compute the trigometric function at k*xi 
  call trigrec_c( xi, nf, cn)  
  call trigrec_s( xi, nf, sn)
  
  do i = 1, n, 1
    pv(i) = c0(row_st + i)
    do k = 1, nf, 1
       pv(i) = pv(i) + ck(row_st+i, k) * cn(k) + sk(row_st+i, k) * sn(k)
    end do
  end do  

  return  
end subroutine varphi


!***********************************************************************
!     ****   dvphi_da_cal   ****
!   for a certain value of xi, the Jacobi Matrix of the invariant curve \varphi 
!   w.r.t. the  Fourier coefficients A.  
!   

!   so dvphi_da is a matrix with ndim same block in diagonal line 
!   we store by Fourier modes, k = 0, 1, ..., nf 

!   ** NOTE ** 
!  -- we put the centre difference as check within each routine for integrety 

!   each block is a row vector of dimension 1 + 2Nf: 
!   1 ,  cos(xi)...  cos(nf*xi),  sin(xi) ... sin(nf*x) 

!   we denote A = [ C0, C's, S's] as the Fourier coefficients, then A is of dimenson n * (2*Nf+1), write as components A = [A_x, A_y, A_z, A_vx, A_vy, A_vz] 
  
!  where \varphi as a truncated Fourier series.
!        \varphi (xi) = C_0 + Sigma_ k from 1 to Nf of { ( C_k cos(k*xi) + S_k cos(k*xi)) }
  
!       Input Variables 
!  xi       the argument for the invariant curve \varphi 
                    
!       Output Variables 
!  dvphi_da      the differential (Jacobi Matrix) of \varphi w.r.t. the Fourier coefficients A            
 
!       Module-based Varaibles
!  n, nf, nf2, nf2p1

!  routine Used:
!   trigrec_c, trigrec_s 

! TODO: check this routine, to make sure every component of the differential is the right
!       one thing is use finite difference   

!  Finally Revised by Yu -- 20160729
!***********************************************************************
subroutine dvphi_da_cal(xi,  dvphi_da)  
implicit none

!  Input and Output Declaration   
real(kind=dp), intent(in)      ::  xi  
real(kind=dp), intent(out)     ::  dvphi_da(n, na1) 
 
!  Local Variable
integer        :: i, jcol, debug 
real(kind=dp)  :: cn(nf), sn(nf) 

                  
  debug  = 0
  
  dvphi_da = 0.d0 ! initialization
  
 ! compute the trigometric function  xi 
  call trigrec_c( xi, nf, cn)  
  call trigrec_s( xi, nf, sn)
  
  ! The Jacobi Matrix of the general function F (the n coordinates)  w.r.t. the  Foureier coefficients A 
  ! each is of dimension n X nf2p1 and are of the same block in the diagonal line, the rest components are zeros

   
  do i = 1, n, 1
  
    ! the starting  index of i-th block  for i-th component
    jcol = (i-1)*nf2p1 + 1   
     
    ! d \varphi / d C_0 = 1
    dvphi_da( i,  jcol) = 1.d0       
     
    ! the 2*nf components are cos(k*xi) for k =1, ...,  nf, for jcol+1 : jcol+nf columns 
    ! sin(k*xi) for k =1, ...,  nf, for jcol+nf_1 : jcol_+2nf columns 
     
    dvphi_da( i,  jcol+1    : jcol+nf  )   =  cn  ! cos( k*xi ), k = 1_nf
    dvphi_da( i,  jcol+nf+1 : jcol+2*nf)   =  sn  ! sin( k*xi ), k = 1_nf
     
  end do

  return  
end subroutine dvphi_da_cal



!***********************************************************************
!     ****   dvphi_drho_xi1   ****
!   for a certain value of xi, the Jacobi Matrix of the invariant curve \varphi 
!   w.r.t. the rotation number   
!   
!   so dvphi_drho is a column array with  n * nf2p1
!   Each n corresponds to one value of xi,  
!     dvphi_drho = Sigma_ k from 1 to Nf of { ( -C_k *k sin(k*(xi_rho) ) + S_k * K cos( k*(xi+rho) ) ) }

!    since    \varphi (xi+rho) = C_0 + Sigma_ k from 1 to Nf of { ( C_k cos(k*xi+k*rho) + S_k cos(k*xi+k*rho)) }
  
!       Input Variables 
!  xi       the argument for the invariant curve \varphi 
                    
!       Output Variables 
!  dgdrho      the differential (Jacobi Matrix) of \varphi w.r.t. the Fourier coefficients A            
 
!       Module-based Varaibles
!  n, nf, nf2, nf2p1

!  routine Used:
!   trigrec_c, trigrec_s 


!  Finally Revised by Yu -- 20170527
!***********************************************************************
subroutine dvphi_drho_xi1( xi, rho, dgdrho) ! only \varhi_0(xi+rho) use this one, A^0
implicit none

! Input  and Output Declaration   
real(kind=dp), intent(in)      ::  xi, rho 
real(kind=dp), intent(out)     ::  dgdrho(n) 
 
! Local Variable
integer        :: i,   k
real(kind=dp)  :: cn(nf), sn(nf), xi_rho, temp 
  

 
  dgdrho = 0.d0     ! initialization
  xi_rho = xi + rho 
  
 ! compute the trigometric function  xi 
  call trigrec_c( xi_rho, nf, cn)  
  call trigrec_s( xi_rho, nf, sn)
  
  ! The Jacobi Matrix of the general function F (the n coordinates)  w.r.t. the  Foureier coefficients A 
  !  each is of dimension n X nf2p1 and are of the same block in the diagonal line, the rest components are zeros
  
  do i = 1, n, 1
  
    temp = 0.d0  
    ! dvphi_drho = Sigma_ k from 1 to Nf of { ( -C_k *k sin(k*(xi_rho) ) + S_k * K cos( k*(xi+rho) ) ) }
    do k = 1, nf, 1 
      temp = temp - ck(i, k)*k*sn(k) + sk(i, k)*k*cn(k)
    enddo  
    
    dgdrho(i) = temp 
     
  end do
  
  return  
end subroutine dvphi_drho_xi1


!***********************************************************************
!     ****   rot_fcs   ****

!  This routine is to avoid the indetermination of the phase shift in  the 
!  computation of invariant curve. 

!  The phase shift indetermination is: given \varphi(xi) as a solution that satisfies
!  the invariance equation \varphi(xi+rho) - P( \varphi(xi) ) = 0, then for any xi0 lies in R
!  \varphi(xi + xi0) also does. 

!  We choose to compute xi0 such that we can make c1(or S1)_J = 0 

!  The Foureier coefficients CK and SK  are transformed to 

!  [ck_hat, sk_hat] ^ T = |  cos(k*xi0)   sin(k*xi0) |  * [ck, sk] ^ T
!                         | -sin(k*xi0)   cos(k*xi0) |

!  For k = 1, we have 
!  c1_hat =  cos(xi0) * c1 + sin(xi0) * s1  =  cos(xi0 - theta)
!  s1_hat = -sin(xi0) * c1 + cos(xi0) * s1  = -sin(xi0 - theta)

! which is actually a clockwise rotation through an angle xi0  around the origin 

! So we have the conclusion that: 
!   if is_ck = 0, we need to make s1_hat = 0,  so xi0 = theta 
!   where theta is the phase anlge of [c1, s1]  in polar coordinate 

!   else if  is_ck = 1, we need to make C1_hat = 0,  so xi0 = pi/2 + theta  
 
 
! Now the problem becomes the computation of xi0 by c1 and s1, be careful with the 
! quadrant 

!  where r1 = sqrt(c1^2 + s1^2), (c1, s1) is the first element in ck and sk respectively

! 
!       Input  Variables 
!  nf        dimension of ck and sk 
 
!       Input and Output Variables 
!  ck, sk    Foureier coefficients CK and SK and the updated ones such that sk(1) = 0          

!  Routine Used:
!     

!  Finally Revised by Yu -- 20160903
!***********************************************************************
subroutine rot_fcs( ck, sk, ind_cs1, is_ck, n, nf)

use dp_mod
use pi_mod
implicit none
 
! Input  and Output Declaration   
integer, intent(in)      :: n,  nf, ind_cs1, is_ck  
real(kind=dp), intent(inout)   ::  ck(nm, nf), sk(nm, nf)  
 
! Local Variable
integer :: i, k 
real(kind=dp)  :: ckcopy(nm, nf), skcopy(nm, nf), c1, s1, r1, theta, xi0, cn(nf), sn(nf)
real(kind=dp)  :: datan_2pi

  ! make a copy 
  ckcopy = ck
  skcopy = sk
    
  ! extract cos(alpha) and sin(alpha) from ck(1) and sk(1)
  ! assume we fix s1_1, from the first components 
  c1 = ck(ind_cs1, 1)
  s1 = sk(ind_cs1, 1)
  r1 = dsqrt(c1**2 + s1**2)
  
  theta = datan_2pi( c1, s1, 1) ! [0, 2pi]
  
  !  compute the phase angle of (c1, s1), and be careful with the quadrant
  print*, 'theta:', theta;  print*;  !read* 
    
  ! compute xi0, according to the  value of is_ck
  if(is_ck == 1) then 
    if(dabs(c1) < 1.d-14) theta = 0.d0
    xi0 = theta + pi/2.d0
  else 
    if(dabs(s1) < 1.d-14) theta = 0.d0  
    xi0 = theta 
  endif  
  
  print*, 'Angle of (c1, s1) = ', theta, c1, s1
  print*, 'clockwise rotation through an angle: ', xi0 

  xi0_rot = xi0 
    
  call trigrec_c(xi0, nf, cn)  
  call trigrec_s(xi0, nf, sn)  
  
  ! test with element-by-element computation  --- keep this for the moment!
  ! ck_hat =  cos( k*xi0 ) * ck + sin(k*xi0 ) * sk 
  ! sk_hat = -sin( k*xi0 ) * ck + cos(k*xi0 ) * sk 
  
  ! we only need to rotate A_0 (C0_0, CK_0, SK_0), the first segment 
  do i = 1, n 
    do k= 1, nf, 1
      ck(i, k) =  cn(k) * ckcopy(i, k) + sn(k) * skcopy(i, k)
      sk(i, k) = -sn(k) * ckcopy(i, k) + cn(k) * skcopy(i, k)
    enddo    
  end do
  
  print*, ind_cs1, '-th coefficient = 0, is_ck=', is_ck; print*
  
  if(is_ck == 1) then 
    print*, 'is_ck = 1 , set c(', ind_cs1, ',1) = 0' 
    print*, ckcopy(ind_cs1, 1), ck(ind_cs1, 1); print*
  else 
    print*, 'is_ck = 0 , set s(', ind_cs1, ',1) = 0' 
    print*, skcopy(ind_cs1, 1), sk(ind_cs1, 1); print* 
    
  endif 
 
  return  
end subroutine rot_fcs

!***********************************************************************
!     ****   assign_ind   ****
! The assignment of the two most important index arrays: ind_dg and ind_dg_pre
!  by specify the value of irow1, irow2 and irow3  and the approach by opt 

!  Finally Revised by Yu -- 20170615
!***********************************************************************
subroutine assign_ind
use dp_mod
implicit none
 
! Local Variable
integer :: i, j, k, irow_c0, irow_c1_shift, irow_c1_back, maxl(n/2), ind_c1_fix
real(kind=dp) ::  c1(n/2)
                  
    
  call fix_index(ck, sk, n, nf, ind_cs1, is_ck)
  
  print*, 'To avoid phase shift indetermination, fix', ind_cs1, '-th Four coef, is_ck = ', is_ck 
  print*, ck(ind_cs1, 1), sk(ind_cs1, 1)
!  print*; read* 
  
  ! select the two columns to remove from dg_all 
  irow1 = 0; irow2 = 0; irow3 = 0
  
  
  ! invariant curve indetermination, fix one component of c0
  irow_c0 =  (ind_c0 -1)*nf2p1 + 1
  
  ! phase shift  -- keep as zero 1 coordinate of c1,s1
  irow_c1_shift = (ind_cs1-1)*nf2p1 +(1-is_ck)**nf + 2 
   
  
  ! falling back to p.o. -- keep fixed one non-zero component of c1,s1 
  c1 = ck(1:n/2, 1)
  maxl = maxloc( dabs(c1) )
  ind_c1_fix = maxl(1)
  
 
   !   opt == 10  11   -- succeed
  ! fix the energy h, and put rho and tp as unknowns, the dimension should be the same, but the DG matrix are not 
  
    
!   opt == 20  21  
  ! fix rho to some irrational number, and put h and tp as unknowns, the dimension should be the same, but the DG matrix are not 
  
      
  irow_c1_back = (ind_c1_fix-1)*nf2p1 + 2 
  
  if(opt == 10 .or. opt == 20) then 
    irow1  = irow_c0
    irow2  = irow_c1_shift
    irow3  = irow_c1_back
    
    print*, 'Avoid falling back to p.o., fix ', ind_c1_fix, '-th c1', c1;
!    print*; read*
    
  elseif(opt == 11  .or. opt == 21) then   
  
    irow1  = irow_c0
    irow2  = irow_c1_shift
    
  elseif(opt == 100) then  
    ! -- Alex's suggestion  --- discard, even not converge at all 
    ! invariant curve indetermination, fix one component of c0
    irow1  = irow_c0
    
    irow2  = irow_c1_shift 
  endif  
    
     
  print*, 'irow1, irow2, irow3 ', irow1, irow2, irow3
!  print*; read*;  
  
  
  if(irow2 /= 0)   call rot_fcs(ck, sk, ind_cs1, is_ck,  n, nf)
  
  ! ************* the key part for Jacobi matrix *********************** 
  ! assignment of the index of unknonws and Jacobi matrix (ind_dg) to be used 
  ! for the prediction, only c0_k in Four Coefs is fixed to avoid curve indetermination
  j = 1; k = 1 
  do i = 1, na_all, 1
    if(i == irow1) cycle 
    ind_dg_pre(j) = i
    j = j + 1
     
    if(i == irow2 .or. i == irow3)  cycle 
    ind_dg(k) = i 
    k = k + 1 
  enddo 
  
!  print*, 'check the dimension of ind_dg_pre and ind_dg = = ?', len_ind_dg_pre, len_ind_dg
!  print*, j-1, k-1 
!  print*; read*
  
  return  
end subroutine assign_ind



!***********************************************************************
!     ****   fix_index   ****
!  Specify the index  of one coefficient associated with the first coordinate (k=1) to be fixed 
!  to avoid indetermination in phase shift 

!  this should be done only once for each curve, so call this routine before the refinement loop

!  Assume \varphi(xi) is an invariant curve we computed, then \varphi(xi+alpha)
!  for any value of alpha will also be a solution, this will cause problem in continuation 
 
!  As suggested by Josep-Maria, we could fix C1 (or S1)  with index ind1 
!    to this purpose, more for consideration of easier programming, we remove one unknown instead of 
!    adding one more equation. 

! -- as a consequence, we have   n*(2*nf+1) equations and  n*(2*nf+1) unknowns ...
    
!** Comment by Josep-Maria: Advanced course on long time integrations, P59, Title: Computation of a torus
!   This is not a problem, as long as we use the general rutine we have described, specifying the kernel dimension to be zero.

!-- PhD thesis, chapter 6, Methodolody, Section 6.2.5,   P100 

!-- Strategy to eliminate unknowns to avoid indeterminations

!   - 2. To avoid the indetermination in phase shift, mostly for continuation cosideration
!       2.1. Eliminate C1_J or C1_J, where J is chosen in order to have || (C1_J, S1_J) ||_2  =  max_{i=1,..., n} || (C1_i, S1_i) ||_2
!            The reason for this choice is heuristic and tries to prevent asking conditions on harmonics of low amplitudes.

!       2.2  Between C1_J and S1_J, we select the coordinate for which min(|C1_J|, |S1_J| ) is achieved. This avoids the situation:
!            assume before the correction, C1_J=0 and |S1_J| = alpha. If we eliminate S1_J, we are forcing  
!            || (C1_J, S1_J) ||_2 > alpha and therefore we will not be able to continue if || (C1_J, S1_J) ||_2 > decreases
!            along the family. If |C1_J| is different from zero but small, we will be able to continue but the continuation step
!            will be artificially reduced.  

!** Another strategy is: Set C(k)_J = 0 as long as S(k)_J remains different from zero, that is  ( C(k)_J, S(k)_J) \= 0 
!   Usually , we take k=1 (the x component) and J has been chosen in order to maximize || ( C(k)_J, S(k)_J) ||_2. 

!   Comment: the first strategy is better, because it doesn't depend on a Fourier coefficients being different from zero 
!   and C1_1 = C1_x = 0, we need to introduce 2 variables as the indice of unknowns to fix the value 
    

!       Input Variables 
!  ck,sk    see gdg_curve 

!       Output Variables 
!  ind_cs1      the index of component corresponding to k=1 to be fixed as zero.  
!  is_ck        flag, 1: fix C1_{ind_cs1}, 0:  S1_{ind_cs1}    
      

!       Module-based Varaibles
!  ndim, nf   for purpose of size declaration

!  routine Used:
!    none  

!  Finally Revised by Yu -- 20160906
!***********************************************************************
subroutine fix_index(ck, sk, n, nf, ind_cs1, is_ck)

implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration 
integer, intent(in)            :: n,  nf
integer, intent(out)           ::  ind_cs1, is_ck   
real(kind=dp), intent(in)      ::  ck(nm, nf), sk(nm, nf)  
 
! Local Variable
integer :: i, ind_loc(1)
real(kind=dp)  :: cs_norm(n)
  
  ! evaluate the coefficients of the x componet || (C1_J, S1_J) ||_2 for J = 1, ..., nf 
  ! but better to take the one in position coordinates  than velocity 
  do i = 1, n, 1
    cs_norm(i) = dsqrt( ck(i, 1)**2 + sk(i, 1)**2 )
  end do  
  
!  print*, 'n=',n, 'nf=', nf, 'sqrt|ck^2+sk^2| =', cs_norm
!  print*, cs_norm; read*
  
  ! take J that achieves  max { || (C1_J, S1_J) ||_2 } 
  ind_loc = maxloc(cs_norm)
  
  ! to extract the value of ind_loc 
  ind_cs1 = ind_loc(1)
  
  print*, 'ind_cs1 =',  ind_cs1
!  read*    
  
  ! between C1_J and S1_J, select the one that achieves min ( |C1_J|, |S1_J|)
  if( dabs( ck(ind_cs1, 1) ) .le. dabs( sk(ind_cs1, 1) ) ) then 
    is_ck = 1
  else 
    is_ck = 0
  endif     
  
  return  
end subroutine fix_index

!***********************************************************************
!     ****   dvphi_da_debug   ****
! Debug dvphi_da,  dp_da  

!  ****** Output Variables ******
!dvphi_da_diff, dp_da_diff, dg_da_diff
!     
! TODO: something is wrong, debug !! 

!  Finally Revised by Yu -- 20170612
!***********************************************************************
subroutine dvphi_da_debug(xi, im, step, gr_map, dvphi_da_diff, dp_da_diff)
implicit none
 
! Input  and Output Declaration  
integer, intent(in) :: im  
real(kind=dp), intent(in)      ::  xi, step  
real(kind=dp), intent(out)     ::  dvphi_da_diff(n,na1), dp_da_diff(n,na1) 
 
external :: gr_map 

! Local Variable
integer :: i, j, k,  icol,   nrow 

! debug: (x-h)  -- (x+h)
real(kind=dp)  :: c0mh(nm), ckmh(nm,nf), skmh(nm,nf),    &  ! check dvphi_da_dif
                  c0ph(nm), ckph(nm,nf), skph(nm,nf),    &   
  
                  ! ******************  debug **************************
                  ! check with center difference f(x+h)-f(x-h)  / 2h
                  ! so the error is of order h^2, but the coefficient in front of h^2 might be big
                  ! to avoid its effect, use Richardson extrapolation explained by Alex 
                  ! ( diff(h)*k^2 - diff(k*h) ) / (k^2-1), is much more accurate 
                  
                  funmh(n), pvfmh(n),  funph(n), pvfph(n), fun1(n),  & 
                  dp_dxmh(n, n), dp_dxph(n, n)  
                  
    
      ! ----------- check d \varphi / dA -----------------------------
      !  ( \varphi(A_i + h) -  \varphi(A_i + h) ) / h 
      dvphi_da_diff = 0.d0
      dp_da_diff    = 0.d0 
   
   ! the change on A^im, affect \varphi_im (xi) and   P ( \varphi_im (xi) )
   ! but \varphi_im (xi) is only for m > 1
   
   
      do i = 1, n
      
         nrow = n*im + i 
         
        do k = 1, nf2p1, 1 
        
          icol = (i-1)*nf2p1 + k  
          
          c0ph = c0; c0mh = c0
          ckph = ck; ckmh = ck
          skph = sk; skmh = sk
           
          if(m > 1 .and. im < m-1)  call varphi( xi, n, im+1, nf, c0, ck, sk, fun1)
           
          ! c0, ck, sk 
          if(k == 1 ) then 
            c0ph(nrow) = c0ph(nrow) + step
            c0mh(nrow) = c0mh(nrow) - step
        
          elseif( k <= nf+1) then 
            j = k -1
            ckph(nrow, j) = ckph(nrow, j) + step
            ckmh(nrow, j) = ckmh(nrow, j) - step
          else 
            j = k - nf - 1 
            skph(nrow, j) = skph(nrow, j) + step
            skmh(nrow, j) = skmh(nrow, j) - step
          endif 
    
          call varphi( xi, n, im, nf, c0ph, ckph, skph, funph) 
          call varphi( xi, n, im, nf, c0mh, ckmh, skmh, funmh) 

          
          dvphi_da_diff(:, icol) = (funph - funmh) / step / 2
      
          ! -------- check d P / d A ----------------
          ! -- x+h  --
          call gr_map( funph,  pvfph, 1, dp_dxph)
          call gr_map( funmh,  pvfmh, 1, dp_dxmh)               
        
          ! --  check d P / d A by center difference----------------
          dp_da_diff(:,  icol) =  ( pvfph  - pvfmh  ) / step  / 2
      
        enddo   
      end do
  
  return  
end subroutine dvphi_da_debug



!***********************************************************************
!     ****   dg_da_debug   ****

!  Finally Revised by Yu -- 20170612
!***********************************************************************
subroutine dg_da_debug(xi, step, gr_map, dg_da_diff)
implicit none
 
! Input  and Output Declaration   
real(kind=dp), intent(in)      ::  xi, step  
real(kind=dp), intent(out)     ::  dg_da_diff(nm,na_all)  
 
external :: gr_map 

! Local Variable
integer :: i, j, k,  icol, row_st, nrow, im 
! debug: (x-h)  -- (x+h)
real(kind=dp)  :: c0mh(nm), ckmh(nm,nf), skmh(nm,nf),    &  ! check dvphi_da_dif
                  c0ph(nm), ckph(nm,nf), skph(nm,nf),    &   
  
                  ! ******************  debug **************************
                  ! check with center difference f(x+h)-f(x-h)  / 2h
                  ! so the error is of order h^2, but the coefficient in front of h^2 might be big
                  ! to avoid its effect, use Richardson extrapolation explained by Alex 
                  ! ( diff(h)*k^2 - diff(k*h) ) / (k^2-1), is much more accurate 
                  
                  funmh(n), pvfmh(n),  funph(n), pvfph(n),   & 
                  dp_dxmh(n, n), dp_dxph(n, n),  funph_rho(n), funmh_rho(n) !, dfun_rho(n) 
                  
    
      ! ----------- check d \varphi / dA -----------------------------
      !  ( \varphi(A_i + h) -  \varphi(A_i + h) ) / h 
      dg_da_diff    = 0.d0 
   
   ! the change on A^im, affect \varphi_im (xi) and   P ( \varphi_im (xi) )
   ! but \varphi_im (xi) is only for m > 1
   
          
   do im = 0, m-1 
   
      row_st = n*im 
        
      do i = 1, n
      
        nrow = row_st + i
         
        do k = 1, nf2p1, 1
        
          icol = k +  (i-1)*nf2p1 + im*na1
          
          c0ph = c0; c0mh = c0
          ckph = ck; ckmh = ck
          skph = sk; skmh = sk
           
!          if(im > 0 ) call varphi( xi, n, im+1, nf, c0, ck, sk, fun1)
           
          ! c0, ck, sk 
          if(k == 1 ) then 
            c0ph(nrow) = c0ph(nrow) + step
            c0mh(nrow) = c0mh(nrow) - step
        
          elseif( k <= nf+1) then 
            j = k -1
            ckph(nrow, j) = ckph(nrow, j) + step
            ckmh(nrow, j) = ckmh(nrow, j) - step
          else 
            j = k - nf - 1 
            skph(nrow, j) = skph(nrow, j) + step
            skmh(nrow, j) = skmh(nrow, j) - step
          endif 
    
          call varphi( xi, n, im, nf, c0ph, ckph, skph, funph) 
          call varphi( xi, n, im, nf, c0mh, ckmh, skmh, funmh) 
        
          if(im == 0) then 
            ! to compute d G / d A by finite difference, we compute \varphi(xi+rho) as well 
            call varphi( xi+rho, n, 0, nf, c0ph, ckph, skph, funph_rho)
            call varphi( xi+rho, n, 0, nf, c0mh, ckmh, skmh, funmh_rho)
!            dfun_rho = funph_rho - funmh_rho  
!          else 
!            dfun_rho = 0.d0  
          endif 
          
          ! -------- check d P / d A ----------------
          ! -- x+h  --
          call gr_map( funph,  pvfph, 0, dp_dxph)
          call gr_map( funmh,  pvfmh, 0, dp_dxmh)               
        
          ! -- check d G / d A by center difference ----
          ! this part is dependent on the value of m 
          
               
          if(m > 1 ) then 
            ! \varphi_0(xi+rho) / d A^0 computed with im==0
            if(im == 0)  dg_da_diff(nm-n+1:nm, icol) = (funph_rho - funmh_rho ) / step / 2 
              
            if(im <= m-1) then 
              ! if m > 1, the last n-rows with \varphi_im (\xi)  - P (\varphi_im-1 (\xi) )
              ! that is dependent on A^im is missing 
              if(im > 0) dg_da_diff(row_st-n+1:row_st, icol ) = (funph-funmh ) / step / 2
              dg_da_diff(row_st+1:row_st+n, icol )            = (-pvfph  + pvfmh ) / step / 2
            endif   
            
          else  ! m==1  
              
            dg_da_diff(row_st+1:row_st+n, icol) = (funph_rho - funmh_rho - pvfph  + pvfmh ) / step / 2
                
          endif 
 
        enddo   
      end do
  enddo 
  
  return  
end subroutine dg_da_debug



