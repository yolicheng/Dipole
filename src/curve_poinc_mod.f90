!  **** My curve Module*******
!  Compute and refine the invariant curve obtained by Poincare map method,  with the energy and rotation number are both fixed
!  Take PRTBP for example, choosing the Poincare section to be x = x_L4 

! ** Note ** 
! for the initialization, we need follow: init_curve --> alloc_arr (init_nf included)

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

module curve_mod

! TODO: what if we put nf as a variable instead of a constant 
!       does it cause problem in decaring the arrays? 

! TODO: what if dp and pi, pi2 conflict with the ones in the main routine?
use dp_mod ! dp and pi, pi2 
use pi_mod 
use time_mod
use poinc_mod, only : isv   ! to avoid the many pass of parameter isv 

implicit none
save 

! subroutines from libraries, declare external attribute, 
! this is only necessary for module to call external subroutines or functions 
external      :: deltx,  trigrec_c, trigrec_s   ! -- libgr.a 
!real(kind=dp), private :: dlange, dnrm2 ! donesn't work here 

real(kind=dp), private ::  tol, tol_err 

! make this public 
integer, dimension(:), allocatable ::  ind_dg  ! the used one and zeros ones 
real(kind=dp), dimension(:),    allocatable ::  c0
real(kind=dp), dimension(:, :), allocatable ::  ck, sk 
real(kind=dp) :: h0_curv  ! energy of the curve, can be used as an unknown 


! --- For declaration of array ---- 
integer, private   ::  n, nf, nf2,  nf2p1,  na_all,       & ! basic dimension 
                       nc_dg, nr_dg, nc_dgall, nr_dgall,  & ! dimension
                       nitmax,  ind_c0, & ! Newton method  
                       debug, opt ! for debug 
                       
!integer, private   ::   nc_dgxi1,  ! possible use 

! avoid phase shift indetermination, fix C1_J(or S1_J)
integer :: ind_cs1, is_ck 


! --- For save -- 
! to write the output of refinment process
integer, private :: fwork, fcurve_approx  


contains
!  include  'refine_curve_poinc.f90'
  
  
  
!***********************************************************************
! Assignment of the energy for curve_mod, might be used as an unknonw 
! For the Time T2 method, adding one more equation for energy. 

subroutine init_curve_h0(h00) 

implicit none
real(kind=dp), intent(in)    :: h00

h0_curv = h00 
end subroutine init_curve_h0


!***********************************************************************
! Fix one component of C0 to avoid curve indetermination  

subroutine init_ind_c0(ind0) 

implicit none
integer, intent(in)  :: ind0

ind_c0 = ind0 
end subroutine  init_ind_c0



!***********************************************************************
! The initialization for control  parameters of Newton method, and which approach to take 
! in fact, we will always choose opt = 3, where we remove one unknowns to get rid of the phase shift unknowns  
 
subroutine init_curve(n0, nitmax0, tol0, tol_err0, opt0) 

implicit none
integer, intent(in)          :: n0, nitmax0, opt0  
real(kind=dp), intent(in)    :: tol0, tol_err0 
logical :: isopen  

  debug = 0
!  print*, 'Debug (1) or not for the module curve_mod?'
!  read(*,*)  debug
  
  ! the dimension of the map, for Poincare map method n = ndim-2,
  ! for Time T map method, n = ndim  
  n     = n0 
  opt   = opt0 
  
  ! termination control for Newton method 
  nitmax  = nitmax0 ! maximum time of iterations allowed
  tol     = tol0
  tol_err = tol_err0
    
!  for save - TODO: for the moment, not used
   fwork = 6
!   fwork = 222 
!   open(fwork, file='curve_poinc.work')
 
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
! -- opt, the strategy for Newton Method 
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
  na_all   = n*nf2p1
  nr_dgall = na_all  
  nc_dgall = na_all  
  
  ! if we want to add more equations or remove certain unknowns, change the value 
  ! of nr_dg and nc_dg. 
  
  ! **NOTE** we always take opt=3 
  
  ! opt = 1, 2, 3 are for Poincare map method 
  ! opt = 11, 12 are for Time Map method 
  if( opt == 1) then ! 
 !  TODO :keep constant 2 unknown,  --- perhaps we can remove this one ...
 ! that is non zero c1,or s1, to make it simple, we take the one with ind_fxd, two much  
   
    nr_dg = nr_dgall   
    nc_dg = nc_dgall - 2 

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
  
    nr_dg = nr_dgall  
    nc_dg = nc_dgall - 1
  
  elseif( opt == 3) then  ! -- only use this one  
    ! -- approach 3 : remove one unknows to make  c1_1  = 0 to avoid phase shift indetermination 
    ! so we have one less unknonws, and the same number of equtaion
       
    nr_dg = nr_dgall
    nc_dg = nc_dgall - 1
     
  elseif (opt == 11 .or. opt == 12 .or. opt == 13) then 
  ! add one more equation as the last one  H(\varphi(0)) - h = 0 
  !    unknown: rho, C0, CK,SK     n*nf2p1 + 1  - 2   = n*nf2p1 - 1 
  !             remove C0_k for curve indetermination and (C1 or S1) for  phase shift indetermination.  

  ! opt = 11 : fix the energy and T2, keep rho as an unknonw 
  !            we add one more column as the differential w.r.t. rho  
  
  ! opt = 12 : fix rho and T2, keep h as an unknown, 
  !            only the right-bottom entry of DG is non-zero and equals -1 
  
  ! opt = 13 : fix rho and h,  keep T2 as an unknown 
  !            we add one more column as the differential w.r.t. T2   

    nc_dgall = na_all + 1 ! one more column for d../d rho 
    nr_dg    = nr_dgall + 1
    nc_dg    = nc_dgall - 2 
    
  endif 
    
  print*, 'size(ind_dg), nc_dg', size(ind_dg), nc_dg; print*; read*
  
  if (allocated(ind_dg) .eqv. .false.) then
    allocate(ind_dg(nc_dg)) 
  else 
    if( size(ind_dg) /=  nc_dg) then 
     deallocate(ind_dg);  allocate(ind_dg(nc_dg))
    endif   
  endif
   
  
  print*, 'number of rows and columns in DG: ' 
  print*, 'nr_dg, nc_dg, nr_dgall, nc_dgall:', nr_dg, nc_dg, nr_dgall, nc_dgall
  read*

  return  
end subroutine init_nf

! This file includes the routines that are related to the Poincare map method, 
! and with a suffix '_poinc'   

! This file will be included in  the module  'curve_mod'

! The ones related to the Time-T map will be named with suffix '_time', 
! and will be saved in a standalone file 'refine_curve_time', also included in 'curve_mod '

!  --1. refine_curve_poinc: the one update all the unknonws
!  --2. gdg_curve_poinc:    compute the differential of the system w.r.t. the unknonws 
!  --3. gdg_inv_xi1_poinc:  the differential of one value of xi 

! **********************************************************************
! This routine is to refine the invariant curve using the  the Poincare map method,
! with the initial guess of the coefficients computed by gr_four   

!**NOTE: Alex's suggestion
!    do not deal with the indetermination of phase shift at this moment  
! but there r problems here!!!!
! during the iterations in Newton Method, the change of Fourier coefficients could 
! lead to vy**2 < 0, not possible to continue...

!     1. eliminate two unknonws h and rho, in order to fix a particular curve,
!       by looking at the Poincare map representation of the torus, we can see we have bunch of 
!       invariant curve, which can be specified by the rotation number rho 
!       so, we eliminate h and rho       

!  So we end up with a  nr (= n*nf2p1) - by -  nc (= n*nf2p1) system of linear equations 
!  we have the same number of equations and unknowns 

!   we have to 
!   The linear equations to be solved:
!                  dg *  dX = - ferr
! 
!  1.  the whole diffrential dg: D F/ d A (C_0, C's, S's) 
!  2.  the error: ferr,  to be refined by Newton method 

!  -- the Linear equation F to slove for Newton Method : 
!     F =  \varphi(xi_i + rho) - P (xi_i) = 0 ],  i = 1,..., Nf    ! the invariance condition 

! TODO: since everytime we are going to update the value of nf according to the value of the Fourier coefficients 
!       so, we will do the initialization of (nf, nf2p1, nf2...., here..) 

!       Input Variables 
!  h        perscribed energy level 
!  rho      the prescribed rotation number           

!       Input-Output Variables 
! c0, ck, sk:  All the Fourier coefficients to be refined 
      
!       Output Variables 
!  isref    flag to show if the refinement is successful or not  
!  iter     number of iterates, used for automatic step control 
         

! Finally Revised by Yu -- 20160826
!----------------------------------------
subroutine refine_curve_poinc(n, nf, rho, c0, ck, sk, isref, iter, gr_poinc, deriv, gr_cj, gr_dcj_dx)

implicit none

! Input  and Output Declaration   
integer, intent(in)     ::  n, nf
real(kind=dp), intent(inout)    ::   rho  ! the prescribed energy, might be updated when opt = 11 
real(kind=dp), intent(inout)    ::   c0(n), ck(n,nf), sk(n,nf) 
integer, intent(out)     ::  isref, iter  
  
! Local Variable
integer        :: info, num_incr, nrhs, ncol, k, irow, i, irow2   ! iter, 
 
                      
real(kind=dp)  :: dg(nr_dg, nc_dg), ferr(nr_dg), dx(nc_dg), & ! Netwon method 
                  dc0(n), dck(n,nf), dsk(n, nf), dx_all(nc_dgall), & ! correction  
                  dxm_pre, dxm, dferrm, nsqrt, &     ! termination control
                  dxi, xi, pt_approx(n) !check the new curve
                  
real(kind=dp) ::  dcj_dx(n) ! the diffrential of CJ w.r.t. PV_0, for Time T2 method 
                  
!integer :: i, j, kcol, kcol_end 
!real(kind=dp)  :: dgaux(nr_dg, nc_dg), wr(nr_dg), wi(nr_dg), vr(nr_dg, nr_dg)  !check the eigenvalue                   
                    
real(kind=dp), external :: dnrm2 
external :: gr_poinc,  deriv, eigrg,  gr_cj, gr_dcj_dx


   debug = 0
   isv  =  1
   call init_nf(nf)
   
!   print*, 'Check refine_curve_poinc:'
!   print*, 'n,nf, nr_dg, nc_dg, nc_dgall, rho: ', n, nf, nr_dg, nc_dg, nc_dgall, rho
!   print*;  read*
   
!  avoid the indetermination in phase shift, opt == 

!  subroutine fix_index(ck, sk, ind_cs1, is_ck)
  call fix_index(ck, sk, n, nf, ind_cs1, is_ck)
  print*, 'To avoid phase shift indetermination, fix', ind_cs1, '-th Four coef, is_ck = ', is_ck 
!  print*; read* 
   
  irow = 0; irow2 = 0 
  
  ! do the index array for dg to keep from dg_all 
  if( opt == 11 .or. opt == 12) then 
  
    irow  = (ind_c0 -1)*nf2p1 + 1
    
    ! phase shift 
    irow2 = (ind_cs1-1)*nf2p1 +(1-is_ck)**nf + 2  ! keep zero 1 coordinate of c1,s1
    call rot_fcs(ck, sk, ind_cs1, is_ck, n, nf)
    
!    irow2 = (ind_cs1-1)*nf2p1 + is_ck*nf + 2  ! keep constant 1 coordinate of c1,s1
    
  elseif(opt == 13) then 
    irow  = (ind_c0 -1)*nf2p1 + 1
    irow2 = (ind_cs1-1)*nf2p1 +(1-is_ck)**nf + 2  ! keep zero 1 coordinate of c1,s1
    call rot_fcs(ck, sk, ind_cs1, is_ck, n, nf)
    
  elseif(opt == 2 ) then  
     
     irow2 = (ind_cs1-1)*nf2p1 + is_ck*nf + 2
      
  elseif(opt == 3) then    
     irow = (ind_cs1-1)*nf2p1 + (1-is_ck)*nf + 2
     call rot_fcs(ck, sk, ind_cs1, is_ck, n, nf)
     
  elseif(opt == 1) then     ! -- discard this one!
  
    irow  = (ind_cs1-1)*nf2p1 + (1-is_ck)*nf + 2  ! phase shift
    call rot_fcs(ck, sk, ind_cs1, is_ck, n, nf)
    
    irow2 = (ind_cs1-1)*nf2p1 + is_ck*nf + 2
    
  endif    
     
  print*, 'irow, irow2, ', irow, irow2
  print*, 'nc_dgall, nc_dg , size(ind_dg)', nc_dgall, nc_dg , size(ind_dg)
  print*; read*; print*; read*
  
  k = 1 
  do i = 1, nc_dgall,1
    if(i == irow .or. i == irow2) cycle 
    ind_dg(k) = i 
    k = k +1 
  enddo 
  
    
  nrhs = 1  ! only one column in right-hand side 
  
  isref = 1 ! by default, we succeed the refinement  
   
  ! counter for iteration and increase of the modulus of correction 
  iter = 0          ! index of iteration 
  num_incr = 0      ! number of times when the modulus correction increases
  
 
  do 
    ! initialize all the correcitons to be zeros  
    dc0 = 0.d0
    dck = 0.d0
    dsk = 0.d0
    
    iter = iter + 1 
    
    if(iter > nitmax) then 
      write(fwork,*) 'Too many iterations! ', iter
      isref = 0 
      return
    endif 
    
    write(fwork,*) '***********', iter, 't-th iteration ***************'
    
    ! ---- differential + error ---- 
!    subroutine  gdg_curve_poinc( rho, c0, ck, sk, dg, ferr, gr_poinc )
    print*, 'befor gdg_curve_poinc, rho, n, nf,  size(c0, ck, sk):'
    print*, rho, n, nf, size(c0), size(ck), size(sk); print*;read*
     
    call  gdg_curve_poinc( rho,  c0, ck, sk, dg, ferr, gr_poinc, deriv, gr_cj, gr_dcj_dx)
    
    
    if(isv == 0) then 
      print*, 'Enter forbidden region: v**2 < 0!'
      isref = 0;       return
    endif 
!   
!    ! -- check ferr and dg after gdg_curve_poinc 
!     print*, 'check dg and ferr after gdg_curve_poinc: '
!     print*, 'nr_dg, nc_dg, n, nf', nr_dg, nc_dg, n, nf; print*; read*

!     if( debug == 1) then 

!      print*, 'dg = '
!      do i = 1, nr_dg, 1
!        write(*,'(10e20.10)') dg(i, :)
!      end do
!      print*; read*  
!      
!      print*, 'ferr = '
!      write(*,'(10e20.10)') ferr
!    endif 

    ! use dxm_pre to detect the modulus of the error, in principle, it should decrease 
    ! we allow it increases for two times, if it exceeds, the refinement fails 
    
    if( iter > 1)  dxm_pre = dxm

!      print*; print*, 'Check dg and ferr as input before calling deltx: '
!      print*, 'dg = '
!      do i = 1, nr_dg, 1
!        write(*,'(10e20.10)') dg(i, 1:10)
!      end do
!      print*; read*  
!      
!      print*, 'ferr = '
!      write(*,'(10e20.10)') ferr(1:10)
!      print*; read*
      
    call deltx( nr_dg, nc_dg, nrhs, dg, ferr, dx, info)
    
!     TODO debug! dg are zeros!
!     print*; print*, 'Check dg and ferr as output( should not change) after calling deltx: '
!      print*, 'dg = '
!      do i = 1, 3, 1
!        write(*,'(9e20.10)') dg(i, :)
!      end do
!      print*; read*  
!      
!      print*, 'ferr = '
!      write(*,'(9e20.10)') ferr(:)
!      print*; read*
!      
!!    if( debug == 1) then 
!      print*, 'the correction : dx'   
!      write(*, '(9e24.14)') dx_all(:) 
!      print*; read* 
!!    endif      
    
    if( info /= 0 ) then 
      write(fwork,*)  'SSL Solver Failed to converge!';  !read*  
      isref = 0
      return
      
    else   
 
    ! Terminate the iteration, if the change trend of the modulus of the dx fluctuates, stop
    ! in principle, the modulus of dx will decrease, we allow it increases twice as suggested by Gerard 
    
    ! Instead of the L2 norm, we look at the maximum norm ?? 
      print*, 'Maximu |dx|:', maxval(dabs(dx))
      
      dxm    = dnrm2(nc_dg, dx,   1) / dsqrt( dble(nc_dg) )
      dferrm = dnrm2(nr_dg, ferr, 1) / dsqrt( dble(nr_dg) )
      
      print*, 'dxm=',dxm, 'derr=', dferrm, 'tol=', tol;    read*  ; !!read*    !ck
      
      ! error less that tolerance
      if( dferrm < tol_err ) then
        print*, 'Error less than tol_err', dferrm, tol_err
        return  
      end if
      
      ! correction is small enough, stop the iteration and return also.....      
      if(dxm < tol) then 
        write(fwork,*) 'Modulus of correction less than tol! Succeed!', dxm 
        return
      endif
      
      
      ! allow the correction to increase at most two times 
      if(iter > 1 .and. dxm > dxm_pre) then 
        num_incr = num_incr + 1 
        
        if(num_incr >= 4 ) then ! allow the correction  dx to increase no more than twice... 
          isref = 0
          print*, 'Terminate refinement! |dx| has increased twice!'; !read*  
          return 
        endif 
        
      endif  
        
      ! -- print the old values of the unknowns to check
      if(debug == 1) then 
        write(fwork,*) 'Old: X0 -- C0  \\  ck, | sk  for \\ (k = 1, ', nf , ')'
      
        ! write the coefficients C's and S's into file fcs.dat 
        write(fwork, *)  c0 
        print*; !read*  
       
        do k = 1, n, 1
!          write(fwork, *)  k, '-th Fourier mode'
          write(fwork, '(8e24.14)')  ck(k, :)   
          write(fwork, '(8e24.14)')  sk(k, :) 
        enddo
      
        ! -- print the correction to check
        write(fwork,*) 'Correction: dX --  dc0 \\  dck, | dsk \\ (k = 1, ', nf , ')'
      endif 
      ! ------------- debug = 1 -----------------
      
      ! -- extract d C's and d S's----
      dx_all          = 0.d0
      dx_all(ind_dg)  = dx 
      
      
      ! ---- update the correction dx for the corresponding control variables 
      do k = 1, n, 1
      
        ! the start index of k-th block (component), for C0
        ncol = 1 + (k - 1) * nf2p1 
 
        if(debug == 1)  print*, 'start column is ', ncol, 'for', k, '-th component'
        
        ! -- d C0  
        dc0(k)    = dx_all(ncol)
        
        ! -- d ck, k=1_nf
        dck(k, :) = dx_all(ncol+1 : ncol+nf)
         
        ! -- d sk, k=1_nf   
        dsk(k, :) = dx_all(ncol+nf+1 : ncol+2*nf)
          
      end do

      !      ! *********** Time Map approach *************
!      ! add one more equtaion H(\varphi(0)) - h = 0 in the last row  
!      
!      ! opt = 11: fix energy h, rho as one more unknown   
!      ! opt = 12: fix rho,    h as one more unknown 

      if(opt == 11 ) then 
      
        rho = rho + dx(nc_dg)
        print*, 'rho, drho', rho, dx(nc_dg)
        print*; read*
        
      elseif(opt == 12) then 
      
        h0_curv = h0_curv + dx(nc_dg)  
        call init_curve_h0(h0_curv)
        
        print*, 'h0, dh', h0_curv,  dx(nc_dg); print*; read*
        
      elseif (opt == 13) then 
        tp = tp + dx(nc_dg)
        call init_time(tp, h0_curv)
      endif  
      
      
      ! --- write the correcitons in the working file to check ?? TODO 
               ! write the coefficients C's and S's into file fcs.dat 
      if(debug == 1) then 
        write(fwork, '(10e24.14)')   dc0 
        do k = 1, n, 1
          write(fwork, '(10e24.14)')  dck(k, :)
          write(fwork, '(10e24.14)')  dsk(k, :) 
        enddo
        print*;  read*  
      endif 
      
      ! --- update the unknowns:  c0 + ck + sk
      c0  = c0 + dc0
      ck  = ck + dck
      sk  = sk + dsk 
      
       
      if(iter > 1) open(fcurve_approx, file='curve_approx.dat', access='append', status = 'old')
    
      ! --- save the indermediate corrected curves 
      dxi = pi2 / 100 
      
      do k = 1, 100, 1
        xi = (k-1)*dxi 
        call varphi(xi, n, nf, c0, ck, sk,  pt_approx)
        write(fcurve_approx, *) xi, pt_approx
      end do
      write(fcurve_approx, *); write(fcurve_approx, *)
      
      close(fcurve_approx)
    endif   ! info == 1
   
  enddo  ! iteration   
  
end  subroutine refine_curve_poinc


!*************************  gdg_curve_poinc **************************** 
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
! dg_all      the differential of F w.r.t. the unkowns, 
!             dimension: n*(2*nf+1) - by -   n*(2*nf+1) -1       
! 
!        Module-Based Varaibles
! ndim, nf, nf2, nf2p1,  pi2 
! nr_dg = n * ( 2*nf+1 )  , nc_dg = n * ( 2*nf+1 ) - 1


!  Finally Revised by Yu -- 20160731
! ---- TODO : check 20160826
!------------------------------------------------------------------------
subroutine  gdg_curve_poinc( rho, c0, ck, sk, dg, ferr, gr_poinc, deriv, gr_cj, gr_dcj_dx)
implicit none

! Input  and Output Declaration   
real(kind=dp), intent(in)       ::  rho,  c0(n)  
real(kind=dp), intent(out)      ::  dg(nr_dg, nc_dg),  ferr(nr_dg)  
real(kind=dp), intent(inout)    ::  ck(n, nf), sk(n, nf)  
! Local Variable
integer        :: i, irow, irow2, k, ncol, nrmn, nlmn 
real(kind=dp)  :: dxi, xi, dg_xi1(n, nc_dgall),  ferr_xi1(n), dg_all(nr_dgall, nc_dgall) 
real(kind=dp)  :: xi0, pv_xi0(n), h_xi0, dcj_dx(n), dcj_da(na_all), dvphi_da0(n, na_all), &
                  dcj_da2(na_all), & ! for debug dcj_da  
                  temp, pvf2(n), dp_dx2(n,n), dh_dt, dh_dt2, dpv_xi0(n),h2 ! debug dcj_dt 


external ::   gr_poinc, deriv, gr_cj, gr_dcj_dx   
  
  ! Initialize as zeros, then we only need to deal with nonzero components 
  dg = 0.d0 
  nrmn = min0(nr_dg, nr_dgall)
  nlmn = min0(nc_dg, nc_dgall)
  
!  ----- deal with argument xi_i one by one  ---------   

! time step size for discretisize the parameter space, dxi = 2*pi/(1+2*nf), 
! and we need  to evaluate xi_i = i * dxi, where i = 0, ..., 2*nf 

  ! step size for xi 
  dxi   = pi2 / nf2 
  
  print*, 'before gdg_inv_xi1_poinc:'; print*; read* 
  
  do i = 1, nf2p1, 1
  
    ! the current xi to be evaluated
    xi = (i-1) * dxi
    
    ! the starting row for each block of xi_i, i = 0, ..., 2*nf 
    irow = 1 + (i-1) * n
    
!     print*, 'i=', i, 'xi=', xi, 'starting row:', irow; print*; read*
        
    !  subroutine gdg_inv_xi1_poinc( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, gr_poinc, deriv)
    call gdg_inv_xi1_poinc( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, gr_poinc, deriv)
    
    ! TODO: the call of cj2v within gr_poinc may return isv=0, 
    ! that is vy^2 < 0, how to deal with this case?
    
!    if(debug == 1)  print*, 'after gdg_inv_xi1:  ferr_xi1 = ', -ferr_xi1;  !read* 
     
    ! the sub-block in dg associated with one value of xi, ndim rows starting from irow
    dg_all(irow : irow+n-1,  :) =  dg_xi1
    
    ! error in F (the invariance equations)
    ferr(irow : irow+n-1 )      =  -ferr_xi1
    
!    print*, 'dg_xi1,  dg_all:'
!    do k =  1, nc_dgall, 1
!      write(fwork, '(12f18.8)') dg_xi1(:,k), dg_all(irow:irow+n-1, k)
!    end do  
!    print*; read* 
  end do 
  

! *********** Time Map approach *************
!      ! add one more equtaion H(\varphi(0)) - h = 0 in the last row  
!      
!      ! opt = 11: fix  h  and t2,  rho as one more unknown   
!      ! opt = 12: fix rho and t2,  h as one more unknown 
!      ! opt = 13: fix h   and rho, t2 as one more unknonw 
   
  if( opt == 11 .or. opt == 12 .or. opt == 13) then
    ! the last row is the energy : H(\varphi(0)) - h  = 0 
    ! add one more row for the extra energy equation, to correct ferr_h =  h - H..

    xi0 = 0.d0 
    call varphi(xi0, n, nf, c0, ck, sk, pv_xi0) 
    call gr_cj(pv_xi0, h_xi0)
    ferr(nr_dg) = h0_curv - h_xi0 
    
    print*, 'h0, H(\varphi(0)),dh',  h0_curv,h_xi0, ferr(nr_dg)
    print*;  read*
    
    ! --- last row for dg_all, d H/ D C0, C1  = d H / d (x,y,z,vx,vy,vz)  * d (x,y,z,vx,vy,vz) / D C0, C1
    ! subroutine deriv_cjlf(pv, dcj)
    
    call gr_dcj_dx(pv_xi0,     dcj_dx)
    call dvphi_da_cal(xi0,  dvphi_da0)
    dcj_da = matmul(dcj_dx, dvphi_da0)

!    subroutine dcj_da_debug(xi, n, nf, c0, ck, sk, gr_cj, dcj_da)
    call dcj_da_debug(xi0, n, nf, c0, ck, sk, gr_cj, dcj_da2)
  
!   ! ************** check dcj_da ********************    
    if(nf == 2) then     
      print*, 'check dcj_da by matmul:'
      write(*,'(9e20.10)') dcj_da
      print*; read*  
    
      print*, 'check dcj_da by centre difference:'
      write(*,'(9e20.10)') dcj_da2
      print*; read*
    endif    
    
!   ! ************** check dcj_da ********************
    ! opt: 11,  d H/ d rho = 0, no need to do assignment since it's initialized to be zero 
    ! opt: 12,  d (-h)/ dh = -1
    if(opt == 12)  dg(nr_dg, nc_dg) = -1.d0 
    
    ! opt: 13,  dH / dT2 = DH/  * DX_0 / DT2
    if(opt == 13) then 
      call deriv(0.d0, pv_xi0, n, dpv_xi0)
      dh_dt = dot_product(dcj_dx, dpv_xi0)
      print*, 'by routine, dh_dt:', dh_dt
      
      
      ! --- debug dh_dt ----- 
      temp = 1.d-4 
      call gr_poinc( pv_xi0, pvf2, 0, dp_dx2)
      call gr_cj(pvf2, h2)
      dh_dt2 = (h2 - h0_curv) / temp 
      print*, 'dh_dt by centre difference:', dh_dt2 
      print*; read*
      
      ! --- debug dh_dt -----
    endif 
    
    
    ! TODO -- 2017-05-27 18:49:08 
    !  dcj_da 
    print*, 'size(ind_dg) == nc_dg-1, na_all', size(ind_dg), nc_dg, na_all
    print*; read*
    
    dg(nr_dg, 1:nc_dg-1)  =  dcj_da(ind_dg(1:nc_dg-1))   
    
    
  endif 
  
  ! only keep the column ind_dg  from  dg_all
  
  dg(1:nrmn, 1:nc_dg-1) = dg_all(1:nrmn, ind_dg)

   
  if(debug == 1) then   
    print*, 'finish dg computation!'
    print*, 'DG, dimension: ', shape(dg)
  
  
    do i = 1, 2, 1
      write(fwork, '(10e20.10)') dg(i,:) ! dg(i,1:20)
    end do
  
    print*, 'ferr, dimension: ', size(ferr)
    write(fwork, '(10e20.10)')  ferr     ! ferr(1:30) 
  
  endif 
  
  return 
end subroutine gdg_curve_poinc



!***********************************************************************
!     ****   gdg_inv_xi1_poinc  ****
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

! TODO: 
!    put a general name for the routine to compute Poincare map and the differential 
!   gr_poinc, and gr_diffpoinc
                    
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
subroutine gdg_inv_xi1_poinc( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, gr_poinc, deriv )
implicit none

! Input  and Output Declaration   
!integer, intent(in)            ::  nf 
real(kind=dp), intent(in)      ::  xi, rho, c0(n), ck(n, nf), sk(n, nf) 
                                   
real(kind=dp), intent(out)     ::  dg_xi1(n, nc_dgall), ferr_xi1(n) 
 
! Local Variable
integer        :: i, j, k     

!real(kind=dp)  :: cn(nf2), sn(nf2), cn_rho(nf2), sn_rho(nf2) ! d \varphi(xi+rho) / d rho


real(kind=dp)  :: fun0(n),  pv0_rho(n),  &  ! the initial state 
                  pvf(n), dp_dx(n, n),   &  ! poinc
                  dp_da(n, na_all),      &  ! dp_da
                  dvphi_da(n, na_all), dvphi_rho_da(n, na_all),  &  ! for dg_xi1
                  dgdrho(n),             &  ! for dvphi_dro_xi1
                  dfun0(n), & ! for dvphi_dt
                                    
                  ! check with center difference f(x+h)-f(x-h)  / 2h
                  ! so the error is of order h^2, but the coefficient in front of h^2 might be big
                  ! to avoid its effect, use Richardson extrapolation explained by Alex 
                  ! ( diff(h1)*n^2 - diff(h2) ) / (n^2-1) , assume h2 = n*h1, is much more accurate 
                  
                  step, funmh(n), pvfmh(n),  funph(n),  pvfph(n),  & 
                  dp_dxmh(n, n), dp_dxph(n, n),  & 

                  dp_da_diff(n, na_all), &   ! check dp_da 
                  funph_rho(n), funmh_rho(n), dg_da_diff(n, na_all), & ! check dg_da the final check!!!!!
                  
                  ! (x-h)  -- (x+h)
                  c0mh(n), ckmh(n,nf), skmh(n,nf), det, &  ! check dvphi_da_dif
                  c0ph(n), ckph(n,nf), skph(n,nf), dvphi_da_diff(n, na_all) 
                  
                  
                  

!  External subroutines                  
external ::  gr_poinc  

!TODO:  check the dimension, there seems to be problems here 
!print*,'dimension of dg: ', nr_dg, nc_dg 
!print*, 'dimension of dg_xi1:', shape(dg_xi1); print*; read* 
! 
 
!  print*,'Debug: gdg_inv_xi1_poinc'
!  print*,'xi, rho,c0: ', xi, rho, c0; print*; read* 
!  
!  
  debug = 0      ! -- ckd 
  step  = 1.d-4  ! to check with center difference
  
  dg_xi1 = 0.d0 ! initialization
  
  ! ---- The Jacobi Matrix of the general function F (the ndim coordinates) w.r.t. A---
  
  ! Evaluate the invariant curve \varphi(xi), with the initial guess of coefficients (c0, ck, sk)
  call varphi( xi, n, nf, c0, ck, sk, fun0) 
  
  ! -- check \varphi, compare fun0 with the output of subroutine four_seri   -- 
  
  ! use general routines to compute poincare map and its differential 
               
!  subroutine std_map( ptin,  ptf, isdp,  dp)  isdp = 1, compute the differential             
!  call std_map( fun0,  pvf, 1, dp_dx)  

!    subroutine PoincMap_plf( pvin, pvf, isdp, dp_dx)
  call gr_poinc( fun0,  pvf, 1, dp_dx)    

    ! just to check 
    if(debug == 1) then  !--ckd, put condition to 2
      print*, 'Finish gr_poinc!'
      print*, 'dp_dx, dimension:', n, ' X ', n
      do i = 1, n, 1
        write(fwork, '(10f22.14)') dp_dx(i, 1:30)
      end do
    
    endif 
           
  ! \varphi( xi+rho ) = pv0_rho
  call varphi( xi+rho, n, nf, c0, ck, sk, pv0_rho) ! the truncated Fourier series 
  
  ! the error to be corrected by Newton Method:  \varphi(xi+rho) - P( \varphi(xi) ) = 0
  ! for the components (y, vx, vy)
  ! first deal with the columns that are related to A_T for the phase components 

  ferr_xi1 = pv0_rho - pvf
 
 !  --- compute d \varphi / dA  ---------
  !  d \varphi(xi) / d A 
  call dvphi_da_cal( xi, dvphi_da)
  
  
 !---------------------------- debug ----------------------------
    if(debug == 1) then ! -- ckd -- put condition to 2
      ! print dvphi_da on screen to compare with the one computed by center difference 
      print*, 'pv0 =', fun0
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
    
          call varphi( xi, n, nf, c0ph, ckph, skph, funph) 
          call varphi( xi, n, nf, c0mh, ckmh, skmh, funmh) 
  
      ! error here, dimension of dvphi_da_diff,  TODO: debug     
          dvphi_da_diff(:, (i-1)*nf2p1 + k) = (funph - funmh) / step / 2
      
          ! to compute d G / d A by finite difference, we compute \varphi(xi+rho) as well 
          call varphi( xi+rho, n, nf, c0ph, ckph, skph, funph_rho)
          call varphi( xi+rho, n, nf, c0mh, ckmh, skmh, funmh_rho)
      
          ! -------- check d P / d A ----------------
      
          ! -- x+h  --
          call gr_poinc( funph,  pvfph, 1, dp_dxph)
    
          call gr_poinc( funmh,  pvfmh, 1, dp_dxmh)               
       
          ! --  check d P / d A by center difference----------------
          dp_da_diff(:,  (i-1)*nf2p1 + k ) =  ( pvfph  - pvfmh  ) / step  / 2
      
          ! -- check d G / d A by center difference ----
          dg_da_diff(:, (i-1)*nf2p1 + k ) = (funph_rho-funmh_rho  - pvfph  + pvfmh ) / step / 2
 
        enddo   
      end do
  
      ! -- d vphi / d A by routine  -- ckd, put into comment
      print*, 'check dvphi_da, dimension', size(dvphi_da), 'should be', n, n*nf2p1
  
      do i = 1, n !ck
        write(fwork, '(10f16.8)') dvphi_da(i, :)
      enddo   
      read*  
  
      ! --- check d vphi / d A  ---
      print*, 'check d vaphi / d A by difference:  '  ! -- ckd 
      do i =  1, n, 1
        write(fwork, '(10f16.8)'), dvphi_da_diff(i, :)
      end do
      print*; read*
  
      print*, 'check ddvphi_da_diff - dvphi_da:  '  ! -- ckd 
      do i =  1, n, 1
        write(fwork, '(10f16.8)'), dvphi_da_diff(i, :)-dvphi_da(i, :)
      end do
      print*; read*
  
      ! --- check d P  / d A ----
      print*, 'check d P  / d A by difference:  ' 
      do i =  1, n, 1
        write(fwork, '(12f18.8)')  dp_da_diff(i, :)
      end do
      print*; read* 

    endif ! if debug == 1   


 
  ! d P / dA  =  d pv_pc{y,vx,vy} / d \varphi(xi) * d \varphi(xi) / d A
  !           = Differential of poincare map  * dvphi_da
  
  ! dimension:  n-by-n  X  n-by-n*nf2p1
  dp_da = matmul(dp_dx, dvphi_da) 
      
  !  d \varphi( xi+rho ) / d A 
  call dvphi_da_cal( xi+rho, dvphi_rho_da)
  
  
  ! Here is the solution: 
  ! If rho is unknown, opt = 11, dg_xi1 has one more column(the last one)
  ! that is the differential w.r.t. rho, calculated by dvphi_drho_xi1 
  
  dg_xi1(:, 1:na_all) = dvphi_rho_da - dp_da  ! dg_da 
  
  

  
    ! --- check d G  / d A ----
    ! TODO: problem in this part for opt=3
    if(debug == 1) then 
 
      print*, 'dp_da = dp_dx X  dvphi_da' 
      do i = 1, n, 1
        write(fwork, '(12f18.8)')  dp_da(i, :)
      enddo 
      print*;  read* 
  
      print*, 'dp_da_diff - dp_da' 
      do i = 1, n, 1
        write(fwork, '(12f18.8)')  dp_da_diff(i,:)-dp_da(i, :)
      enddo 
      print*;  read*
  
      print*, 'check dvphi_rho_da, dimension: ', size(dvphi_rho_da) !ck
      do i = 1, n !ck
        write(fwork, '(12f18.8)') dvphi_da(i, :)
      enddo   
      print*; read*  
 
 
      print*, 'check d G  / d A by difference:  ' 
      do i =  1, n, 1
        write(fwork, '(12f18.12)')  dg_da_diff(i, :)
      end do
      print*; read* 
 
      print*, 'd G / d A = dg_dx1 = dvphi_rho_da - dp_da ' 
      do i = 1, n, 1
        write(fwork, '(12f18.12)')  dg_xi1(i, :)
      enddo 
      print*; read* 
  
    endif  

!  print*, 'd G / d A = dg_dx1 = dvphi_rho_da - dp_da ' 
!  do i = 1, n, 1
!    write(fwork, '(12f18.12)')  dg_xi1(i, 1:24)
!  enddo 
!  print*; read* 
  
  return  
end subroutine gdg_inv_xi1_poinc  

!***********************************************************************
!     ****   upd_four   ****
!  Decide the approriate value of nf, and declare the Fourier coefficients CK, SK, 
!  and update the array size (if necessary) for Poincare map method. 

! The criteria is: 
!  the maximum norm of the last half of the Fourier Coef should be of one order 
!  of magnitude less than the tolrence, < tol / 1.d1, otherwise, we double nf=2nf 
!  and set the new ones as  zeros  
! 
!  ****** Input  Variables ******
!  isupd            0: Initialize c0, ck, sk, declare zero arrays of dimension n-by-nf.
!                   1: update ck,sk 

!  ****** Input and Output Variables ******
!  nf0              the old value of number of harmonics, to be update  is isupd=1
!                   ** NOTE **  Remerber to update nf0! since it is not a global value 

!  c00, ck0, sk0    the old fourier coefficients, dimension n X nf0


!  Module-base Varaibles:
!  private: n, tol 
!  public: c0, ck, sk 

!  Finally Revised by Yu -- 20170520
!***********************************************************************
subroutine upd_four(isupd, nf0) !, c00, ck0, sk0) !, c0, ck, sk)

implicit none
! Input  and Output Declaration  
integer, intent(in)          :: isupd 
integer, intent(inout)       :: nf0
!real(kind=dp), allocatable, intent(inout) :: c00(:), ck0(:,:), sk0(:,:)   
 
! Local Variable
integer        :: infb2, inf, nfmn, i 
real(kind=dp)  :: ckm, skm, tolb10, ck_copy(n,nf0), sk_copy(n,nf0), c0_copy(n)
real(kind=dp)  :: dlange 
 
  nf = nf0
  
  if(isupd == 0) then 
    print*, 'The initialization for C0, CK, SK.'; 
    print*, 'upd_four:, n, nf', n, nf; print*; read*
   
    call init_nf(nf)
    
    call alloc_arr_1d(n, c0)
    call alloc_arr_2d(n, nf0, ck); 
    call alloc_arr_2d(n, nf0, sk); 
    
    c0 = 0.d0; ck = 0.d0; sk = 0.d0
    return
  endif
  
!  write(*, *) '# Old c0, ck, sk : n, nf = ', n, nf 
!  do i = 1, n, 1
!    write(*, '(6e24.14)')   c0(i);    ! write(*,*) 
!    write(*, '(10e24.14)')  ck(i, 1:10); ! write(*,*) 
!    write(*, '(10e24.14)')  sk(i, 1:10) 
!   write(*,*);  ! write(*,*) ! two blank line2 to seperate components 
!  end do
!  print*; read*
  
  ck_copy = ck; sk_copy = sk; c0_copy = c0 
  
! -- check the Maximum norm of csf and sif, and choose an appropriate value for nf 
  ! FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
  inf = nf0 / 2
  ckm = dlange('m', n, inf, ck(:, inf+1 : nf0), n, 0.d0)
  skm = dlange('m', n, inf, sk(:, inf+1 : nf0), n, 0.d0)
  
  print*, 'maximum norm of CK and SK for the last ', nf0/2, 'Fourier harmonics'
  print*,  ckm, skm 
  
  print*; read*

  tolb10 = tol / 1.d1
  if(dmax1(ckm, skm)  < tolb10 ) then
    ! check the maximum norm of the last half of the remainning ck0(1:nf0/2), sk0(1:nf0/2),  if it's still within the tolerance, decrease nf by a half.
    infb2 = inf/2 
    ckm = dlange('m', n, inf, ck(:, infb2+1 : inf), n, 0.d0)
    skm = dlange('m', n, inf, sk(:, infb2+1 : inf), n, 0.d0)
    
    if(dmax1(ckm, skm)  < tolb10 .and. nf0>16) nf = nf0 / 2
    print*, 'Last ', infb2, 'Coefs:', ckm, skm
  else 
    ! double the number of harmonics
    nf = 2 * nf0
    
  endif 
  
  if(nf /= nf0) then ! if we need to declare arrays with new dimension
  ! NOTE ** if the input ck0, sk0 are ck, sk such that we want to update the ck,sk 
  ! the deallocate will make ck0, sk0 = 0  
    ! -- declare memory for  ck, sk 
!    print*, 'Number of fourier modes used, nf=  ', nf, 'Original nf0 = ', nf0 
    call alloc_arr_2d(n, nf, ck); call alloc_arr_2d(n, nf, sk)
    ck = 0.d0; sk = 0.d0 
  endif 
         
  nfmn = min0(nf, nf0)
  
  ! assignment for public module variables
  print*, 'nf=', nf, '  nf0=', nf0, '  nfmn=', nfmn ; print*; read*
  c0            = c0_copy                   ! the constant term 
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
!     ****   alloc_arr_1d   **** 
!  allocate array csf(1:nf), sif(1:nf) for the FFT to obtain Fourier coefficients
!***********************************************************************
subroutine alloc_arr_1d( nf0, csf)  

implicit none
! Input  and Output Declaration  
integer, intent(in)       :: nf0 
real(kind=dp), allocatable, intent(inout) :: csf(:) !, sif(:)   

  if (allocated(csf) .eqv. .false.) then
    allocate(csf(nf0)); ! allocate(sif(nf0))
  else 
    deallocate(csf);    ! deallocate(sif)
    allocate(csf(nf0)); ! allocate(sif(nf0)) 
  endif
    
  return
end subroutine alloc_arr_1d

!***********************************************************************
!     ****   allococate 2D array  of dimention n-by-nf0
!   use 0 to avoid conflicts  with the variables from the module 
!***********************************************************************
subroutine alloc_arr_2d(n0, nf0, ck0)!, sk0)  

implicit none
! Input  and Output Declaration  
integer, intent(in)       ::n0,  nf0 
real(kind=dp), allocatable, intent(inout) :: ck0(:,:) !, sk0(:,:)   
  
  if (allocated(ck0) .eqv. .false.) then
    allocate(ck0(n0, nf0));!  allocate(sk0(n, nf0)) 
  else 
    deallocate(ck0);       ! deallocate(sk0) 
    allocate(ck0(n0, nf0)); ! allocate(sk0(n, nf0)) 
  endif 
    
  return
end subroutine alloc_arr_2d


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

!  routine Used:
!     trigrec_c, trigrec_s

! TODO: probably we will not need to compute pv using a standalone subroutine 
!       because at the point when pv is needed, we have already computed  cn, sn 

!  Finally Revised by Yu -- 20160730
!***********************************************************************
subroutine varphi( xi, n, nf, c0, ck, sk, pv)
implicit none

! Input  and Output Declaration 
integer, intent(in)     ::  n, nf
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
!     ****   dvphi_da_cal   ****
!   for a certain value of xi, the Jacobi Matrix of the invariant curve \varphi 
!   w.r.t. the  Fourier coefficients A.  
!   
!   so dvphi_da is a matrix with ndim same block in diagonal line 
!   we store by Fourier modes, k = 0, 1, ..., nf 

!   ** NOTE ** 
!   if rho is an unknown, dgda has one more column, which is the differential w.r.t. rho 

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
!  n, nf, nf2, nf2p1

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
integer        :: i,  jcol
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
      write(fwork, '(12f18.8)') cn
     
      print*, 'sn='
      write(fwork, '(12f18.8)') sn 
    endif   
     
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
subroutine dvphi_drho_xi1( xi, rho, dgdrho)
implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration   
real(kind=dp), intent(in)      ::  xi, rho 
real(kind=dp), intent(out)     ::  dgdrho(n) 
 
! Local Variable
integer        :: i,   k
real(kind=dp)  :: cn(nf), sn(nf), xi_rho, temp 
  

  dgdrho = 0.d0 ! initialization
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
real(kind=dp), intent(inout)   ::  ck(n, nf), sk(n, nf)  
 
! Local Variable
integer :: i, k 
real(kind=dp)  :: ckcopy(n, nf), skcopy(n, nf), c1, s1, r1, theta, xi0, cn(nf), sn(nf)
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
  print*, 'theta:', theta; print*; read* 
    
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
!  read*
    
  call trigrec_c(xi0, nf, cn)  
  call trigrec_s(xi0, nf, sn)  
  
  ! test with element-by-element computation  --- keep this for the moment!
  ! ck_hat =  cos( k*xi0 ) * ck + sin(k*xi0 ) * sk 
  ! sk_hat = -sin( k*xi0 ) * ck + cos(k*xi0 ) * sk 
  
  do i = 1, n 
    do k= 1, nf, 1
      ck(i, k) =  cn(k) * ckcopy(i, k) + sn(k) * skcopy(i, k)
      sk(i, k) = -sn(k) * ckcopy(i, k) + cn(k) * skcopy(i, k)
    enddo    
  end do
  
  print*, 'check carefully if ', ind_cs1, '-th coefficient = 0? ';
!  print*, 'if is_ck = 0, S(', ind_cs1, ',1) = 0' 
!  print*, 'if is_ck = 1, C(', ind_cs1, ',1) = 0'
  print*;   print*, 'is_ck = ', is_ck; ! read* 
  
  if(is_ck == 1) then 
    print*, 'is_ck = 1 , set c(', ind_cs1, ',1) = 0' 
    print*, ckcopy(ind_cs1, 1), ck(ind_cs1, 1); print*
  else 
    print*, 'is_ck = 0 , set s(', ind_cs1, ',1) = 0' 
    print*, skcopy(ind_cs1, 1), sk(ind_cs1, 1); print* 
    
  endif 
      
  print*; read*
  
!  print*, ' ck \\ sk'
!  do i = 1, n, 1
!    write(*, '(10e20.10)') ck(i, :); print* 
!    write(*, '(10e20.10)') sk(i, :)
!    print*; print*; read*
!  end do
 
  return  
end subroutine rot_fcs




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
real(kind=dp), intent(in)      ::  ck(n, nf), sk(n, nf)  
 
! Local Variable
integer :: i, ind_loc(1)
real(kind=dp)  :: nm(n)
  
  ! evaluate the coefficients of the x componet || (C1_J, S1_J) ||_2 for J = 1, ..., nf 
  do i = 1, n, 1
    nm(i) = dsqrt( ck(i, 1)**2 + sk(i, 1)**2 )
  end do  
  
!  print*, 'n=',n, 'nf=', nf, 'sqrt|ck^2+sk^2| =', nm
!  print*, nm; read*
  
  ! take J that achieves  max { || (C1_J, S1_J) ||_2 } 
  ind_loc = maxloc(nm)
  
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



end module curve_mod
