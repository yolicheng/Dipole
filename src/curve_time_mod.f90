module curve_time_mod

!  **** My curve Module based on Time T map method *******
!  Compute and refine the invariant curve obtained by stroboscopic map (Time T map) 
!  Among h, T2 and rho, we need to fix two to specify one torus 

! opt=10, we fall back to the p.o., so we have to apply JM's trick, to fix one component of c1(s1)

!  Figure out the dimension of the problem, the number of unknonws and equations 
!  which one should be fixed, or leave free for the Newton Method

!  ** Note ** 
!  for the initialization and the update of nf and c0,ck,sk, we use upd_four, 
!  set isupd = 0 and = 1 respectively.

!  the ideas are different between the Time t map and Poincare map method.  
!  so do a module seperately, although they have some similar part, but do not mix the two methods

!  Around a periodic orbit with central part, we have plenty 2D tori  
 
!    1. take a periodic orbit we have already computed, and compute the MM, associated eigenvalues and eigenvectors
!    2. take the linear flow as the initial guess for invariant curve \phi: x0 + epsilon * (vr*cos(xi) - vi*sin(xi))
!       -- epsilon is the distance from the periodic orbit, vr+ivi is the associated eigenvector, 
!    5. do a simple Fourier analysis to compute the Fourier coefficients as the initial guess \bm X_0
!    6. Impose the invariance equations, and specify the energy level, to compute the invariant curve 
!    7. Globalize the curve to an invariant torus, it is easy, we only need to integrate for one period for each point  


! ** NOTE **
!    3.1  there are two indeterminations: the curve indetermination and  the phase shift indetermination, 
!         it can be avoided by fixing C1_1 = 0, this is done by remove C1_1 from the unknows,
!         although in the system equations

!         it is expressed by adding another equation  C1_1 = 0. 
!         since if we keep C1_1 as one unknown, then the Jacobi Matrix will be singular, which makes it difficult for Newton Method 

!    3.2  by fixing the energy H0 and the rotation number \rho, we fix a particular invariant curve 

!    3.3  the unknowns are only the Fourier coefficients: C0, CK, SK 

! First try: fix rho and T2, let free of h, that means no additional equation for the energy. 


!    so, we have the number of unknowns : n * (2*nf+1) - 1 
!    where n is the dimension of the phase space to be evaluated,  n = ndim -2
!    since we fix h0 and x0, decrease dimension by 2
!    nf is the number of Fourier modes. 

!    the number of equations 
!           n * (2*nf+1)  

! ** Note ** 
!    1. we have 2 more equations than unknows, but Josep Maria claims the kernel dimension is zero... 
!    2. use the general name for vector field and energy (Jacobi costant)
!    3. subroutines from libraries, need to declare external attribute **
!         external  ::  deltx, trigrec_c, trigrec_s      ! ----  libgr.a 


use dp_mod   
use pi_mod 

implicit none
save 

! subroutines from libraries, declare external attribute, 
! this is only necessary for module to call external subroutines or functions 
external      :: deltx,  trigrec_c, trigrec_s   ! -- libgr.a 

real(kind=dp), private ::  tol, tol_err, xi0_rot
integer       ::  irow1, irow2, irow3, m, nf2p1, na_all, nc_dg, nr_dg ! index of unknonws removed 

integer, private       :: nf_min, nf_max 
real(kind=dp), private :: ds_min, ds_max 

! make this public 

integer, dimension(:), allocatable ::  ind_dg, ind_dg_pre  ! the used one and zeros ones 
real(kind=dp), dimension(:),    allocatable ::  c0, dx_pre
real(kind=dp), dimension(:, :), allocatable ::  ck, sk 
real(kind=dp)  :: tp, rho,  h0  ! make the three parameters of the torus public 


! --- For declaration of array ---- 
integer, private   ::  n,   nm, na1,  nf, nf2,  len_ind_dg, len_ind_dg_pre, &
!                  nf2p1,  na_all, nc_dg, nr_dg, & ! basic dimension 
                       nr_dgall, nc_dgall, len_dx_pre, & ! for dg_all, including all the unknonws and equations
                       nitmax,   ind_c0, &  ! Newton method  
                       debug, opt ! for debug 
                       

! avoid phase shift indetermination, fix C1_J(or S1_J)
integer :: ind_cs1, is_ck 


! --- For save -- 
! to write the output of refinment process
integer, private :: fwork, fcurve_approx  


contains
  
  ! common routine between time map method and Poinc map method 
  include  'curve_mod.f90' 
  include  'curv_angle.f90'
  
!********************* init_time_h0 ************************************
! Assignment of the energy for curve_mod, might be used as an unknonw 
! For the Time T2 method, adding one more equation for energy. 

subroutine init_time_h0(h00) 

implicit none
real(kind=dp), intent(in)    :: h00

  h0 = h00 
end subroutine init_time_h0


!************************* init_time_tp *********************************
subroutine init_time_tp(tp0) 
implicit none

real(kind=dp), intent(in)   :: tp0

  tp = tp0
  return  
end subroutine init_time_tp

!************************* init_ind_c0**********************************
! Fix one component of C0 to avoid curve indetermination  

subroutine init_ind_c0(ind0) 

implicit none
integer, intent(in)  :: ind0

  ind_c0 = ind0 
  return
end subroutine  init_ind_c0


!***********************************************************************
!     ****   tori_cont_time   **** ! TODO  
! The continuation of tori family, using predictor-correcter scheme  for Time T map method
! all the variables are declared global, all modifications of the values are updated

! opt = 10 for the first torus and opt = 11 for the following ones. 

!  *** step control ****

!  We check the angle formed by the three last tori 
!  (to compute this angle we take into account h, δ, ρ and the Fourier coefficients of order zero). 
!  If it is greater than a given tolerance (we have typically used 15◦), 
!  we divide by 2 the continuation step and restart
!  In some situations, it is necessary to restart the process from the first of the three last tori.


!  ****** Input Variables ******
!  The old values of unknonws,      rho, tp, h0, c0, ck, sk 
!  ntori    number of tori to be computed 
!  ds       the step size for continuation 

!  Finally Revised by Yu -- 20170715
!***********************************************************************
subroutine tori_cont_time( ds, dir, ntori, opt0, opt1, itori, fcs, fcurve,  & 
       fpara_curve, ftori, gr_map, deriv, gr_cj, gr_dcj_dx)

use dp_mod
implicit none
 
! Input  and Output Declaration  
integer, intent(in)          ::  ntori, opt0, opt1, fcs, fcurve, fpara_curve, ftori 
integer, intent(out)         ::  itori  ! number of tori computed 
real(kind=dp), intent(inout) ::  ds, dir   
external :: gr_map, deriv, gr_cj, gr_dcj_dx

real(kind=dp) :: dlange  

 
! Local Variable, it is not even possible to store the previous three curves
integer       :: nf_pre(3), nf0, nf_old, isincr
real(kind=dp) :: c0_pre(3,nm), rho_pre(3), tp_pre(3), h0_pre(3),    &  
                 dxi, xi,  pv_curv(n), pv_curv_rho(n), cj, pvf(n), &
                 ckm, skm, ds0, csangle, &
                 rho_old, tp_old, h0_old, c0_old(nm), dir0 
                
! the previous three curves  
real(kind=dp), allocatable, dimension(:,:)  ::  ckpre1, skpre1, ckpre2, skpre2,  &
                           ckpre3, skpre3, ck_old, sk_old

   
integer :: isrestart, isref, i,  niter, debug, np, incr_time 

real(kind=dp) :: dnrm2 


debug = 0

isrestart = 0
itori     = 1 
ds0  = ds 
dir0 = dir 

nf_pre(1) = nf
rho_pre = 0.d0; tp_pre = 0.d0; h0_pre = 0.d0
 csangle = 0.d0 


do 
  
  if(itori > ntori) return 
  
  ! The approach (opt) is different for the first and the following tori, opt0 and opt1, respectively
  if(itori == 1) then 
    call init_opt(opt0) 
  else 
    call init_opt(opt1) 
  endif   
  
  
  ! ******** restart from the previous curve *************        
    
  if(isrestart == 1) then  
    print*, 'ds == ?', ds 
    ! -- keep nf, change ds, it doesn't have to be ds/2.d0 
    if(dabs(ds) < ds_min) then 
      ds = ds * 2.d0 
      if(nf == nf_max) then 
        print*, 'ds needed is smaller than the  minimum! Stop'
        stop
        
      else 
        nf = nf*2;  isrestart = 2
        cycle 
      endif
      
    elseif(dabs(ds) > ds_max) then
      print*, 'Cannot increase ds, so increase nf!' 
        ds = ds_max*dsign(1.d0, ds)
        isrestart = 2; 
        cycle
    endif
     
    
    call upd_four(0, nf);
    rho = rho_pre(1); tp = tp_pre(1); h0 = h0_pre(1) 
    c0  = c0_pre(1, :);   ck(:, 1:nf_pre(1))  = ckpre1;  sk(:, 1:nf_pre(1))  = skpre1
    
    isrestart = 0 
    
    if(nf == nf_pre(1) ) then 
       goto  201
    else 
       goto  200 ! recompute all the dx_pre, and do init_nf and assign_ind 
    endif     
    
    
  elseif(isrestart == 2) then  
    ! -- only update the value of nf, using the old values before refinement 
      
    ! -- only one case, we do this, that is niter > nitmax 
    print*; print*, '************ Increase Nf ************';print*;
    nf = 2*nf
    
    print*, 'nf, nf_old: ', nf, nf_old; print*; ! !read* 
    
    if(nf > nf_max) then 
      print*, 'Nf needed is too big! nf > nf_max', nf, nf_max 
      stop
    endif 
     
           
    call init_nf(nf);   call upd_four(0, nf) 
    
    rho = rho_old;  tp = tp_old;   h0 = h0_old;  
    c0  = c0_old;
    ck(:, 1:nf_old) = ck_old;  sk(:, 1:nf_old) = sk_old;
    
    
  endif 
  
  ! **********************************************************
  
  print*; print*, '*********', itori, '-th torus:  nf0, nf, ds, csangle, dir =', nf_pre(1), nf, ds, csangle, dir ; print* 
!  if(itori > 55) !read*
  
  ! --- save the variable before refinement 
  nf_old  = nf;        c0_old  = c0
  rho_old = rho; tp_old = tp;  h0_old = h0;  
  call alloc_arr_2d(nm, nf, ck_old); ck_old = ck 
  call alloc_arr_2d(nm, nf, sk_old); sk_old = sk
  
  call refine_curve_time(n, nf, rho, c0, ck, sk, isref, niter, gr_map, deriv, gr_cj, gr_dcj_dx)
  
!  
  if(isref == 2) then 
    ! |d F| > tol with the Newton method converge, we need to increase NF 
    isrestart  = 2 
    print*, 'Converge with |d F| > tol! Increase nf and restart again!'; !read*
    cycle

  elseif (isref == 0) then 
    print*, 'Newton method fails to converge! Restart with ds/2, current ds=', ds; !read*
    
    ! first we tried to increase ds, but failed to converge, so we try increase nf 
    if(isincr == 1) then 
!      nf = nf * 2  
      ds = ds / 4.d0;   isrestart = 1;   cycle 
    endif 
      
      if(itori == 1) then    
        print*, 'Fail for the first torus! Increase nf'
        isrestart = 2;  
        cycle 
      endif

     ds = ds / 2.d0; 
     isrestart = 1 
     cycle
         
  else 
  
  ! -------------- succeed with a new curve ----------------------- 
    !  check if NF is accurate enough, if not, restart again with 2nf 
    !  check if we fall back to a p.o., ck_k, sk_k with k>=1 has zero norm 
    
    print*, 'succeed with the refinement, check if we get a new curve!'; print* 
    
    ckm = dlange('m', nm, nf, ck, n, 0.d0)
    skm = dlange('m', nm, nf, sk, n, 0.d0)
    print*,  ckm, skm 
    if(dsqrt(ckm*ckm+skm*skm) < tol/1.d1) then 
      print*, 'We fall back to p.o.!, ckm, skm =', ckm, skm; print*; !read*
      stop 
    endif 
    
    nf0 = nf;  call check_tail(isincr)
    
    if(isincr == 1) then 
      ! in principle we need to double it 
      if(nf >= 2*nf_pre(1) ) then 
         print*, 'Last 1/4 beyond the tolerance! Decrease ds!';  !!read*
         ds = ds / 2.d0; 
         isrestart = 1
         cycle
      else 
         
        ! if we have already the maximum number of Fouerier Modes, we decrease the step size and restart. 
        if(nf == nf_max) then
           ds = ds / 2.d0 
           isrestart = 1
           cycle
         endif   
         
        ! we only allow update three times, if we are not able to jump over it, increase nf. 
        incr_time = incr_time + 1 
        if(incr_time > 3)    then 
           if(nf <= nf_max/2 ) then 
              nf = 2*nf 
            else 
              ds = ds / 2.d0
              isrestart = 1; 
              cycle 
            endif 
        endif       
         
        print*, 'Try increase ds and jump over the gap! incr_time = ', incr_time  !!read*, 
        
        ! TODO:  we can keep the old nf and use 2*ds to jump over the gap 
         
        ds        = ds* 1.3d0;  ! for 2, it jump too fast,  ! using 1.6 we are able to jump 1/8 resonance for tori=22
        isrestart = 1 ;  
        cycle
      endif 
    
    elseif(isincr == -1 .and. itori > 1) then 
        
      call upd_four(-1, nf)  ! we need less Fourier modes 
    
    endif 
    
    
    ! ------- succeed with a new curve --------
    print*, 'Succed with a new curve! ntori = ', itori; print*;  !!read*
    isrestart = 0
    
    call init_nf(nf);   
  endif   
  
  
  if(debug == 1) then 
    print*, 'If isrestart = 1 or isref = 0, error!', isrestart, isref; 
    print*; !read*
  endif 
  
  
  if(ntori == 1)  goto  100
  
  
!  **********  step control ************  
!  We check the angle formed by the three last tori 
!  (to compute this angle we take into account h, δ, ρ and the Fourier coefficients of order zero). 
!  If it is greater than a given tolerance (we have typically used 15◦), 
!  we divide by 2 the continuation step and restart
!  In some situations, it is necessary to restart the process from the first of the three last tori.
  
  ! the previous curve is refined successfully, save the refined values 
  ! copy the old c0, ck, sk, nf, if for the upd_four, it may happen that
  ! we need to double nf and restart again 
  
!  print*, 'rho, tp, h0', rho, tp, h0; print* 
  
  if(itori == 1) then 
    nf_pre(1)  = nf;   c0_pre(1,:) = c0 
    rho_pre(1) = rho;  tp_pre(1)   = tp;  h0_pre(1) = h0
    call alloc_arr_2d(nm, nf, ckpre1); ckpre1 = ck 
    call alloc_arr_2d(nm, nf, skpre1); skpre1 = sk 
    
      
  elseif(itori == 2 ) then 
    if(debug == 1) then 
      print*, 'Save the second curve after refinement!'
      print*; !read*
    endif 
    
    nf_pre(2) = nf_pre(1);    c0_pre(2,:) = c0_pre(1,:)
    rho_pre(2) = rho_pre(1);  tp_pre(2)   = tp_pre(1);   h0_pre(2) = h0_pre(1)
    call alloc_arr_2d(nm, nf_pre(1), ckpre2); ckpre2 = ckpre1 
    call alloc_arr_2d(nm, nf_pre(1), skpre2); skpre2 = skpre1    
    
    nf_pre(1) = nf;    c0_pre(1,:) = c0 
    rho_pre(1) = rho;  tp_pre(1) = tp;   h0_pre(1) = h0
    call alloc_arr_2d(nm, nf, ckpre1); ckpre1 = ck
    call alloc_arr_2d(nm, nf, skpre1); skpre1 = sk    
      
  else 
  
!   print*, 'More than three curves, check the angle !'
    ! itori > 3, save itori as the last one, ckpre1...also update the previous two
    ! check if the angle is small (less than 15 degree), if so, we keep the result 
    ! if not, we start from the previous curve with smaller ds.
    
    ! subroutine check_angle(c0_pre, rho_pre, tp_pre, h0_pre, isrestart, dir)
    c0_pre(3, :) = c0;  rho_pre(3) = rho; tp_pre(3) = tp; h0_pre(3) = h0  
    call check_angle(c0_pre, rho_pre, tp_pre, h0_pre, isrestart, csangle)
    
    if (isrestart == 1)  then 
      print*, 'ds = ', ds 
!        decr_time = decr_time + 1 
!        if(decr_time > 3 .or. dabs(ds) < ds_min * 2.d0) then 
          ! we restart from the second previous curve pre(2), instead of pre(1) 
!          ds = ds / 2.d0 
!          cycle 
!        endif 
      ds = ds / 2.d0 
!      isrestart = 1
      cycle 
    endif   
     
    
    if(debug == 1) then 
      print*, 'check the angle between the previous tori, itori =', itori, 'angle=', csangle, 'dir=', dir
      do i = 1, 3, 1
       print*, rho_pre(i), tp_pre(i), h0_pre(i)
      enddo  
      print*; !read*
      print*, '*********************************************************'
    endif 
    
    if(csangle < 0) then 
      print*, 'The new refined curve reverses the sense! Restart!' 
      dir = -dir 
      ds = ds / 2.d0 
      isrestart = 1
      cycle 
    endif 
    
    nf_pre(2:3)  = nf_pre(1:2);  c0_pre(2:3,:) = c0_pre(1:2,:) ! fixed dimension 
    rho_pre(2:3) = rho_pre(1:2); tp_pre(2:3)   = tp_pre(1:2); h0_pre(2:3) = h0_pre(1:2) 
    ! update pre3 to be the old pre2 
    call alloc_arr_2d(nm, nf_pre(2), ckpre3); ckpre3 = ckpre2
    call alloc_arr_2d(nm, nf_pre(2), skpre3); skpre3 = skpre2
     
    ! update pre2 to be the old pre1
    call alloc_arr_2d(nm, nf_pre(1), ckpre2); ckpre2 = ckpre1
    call alloc_arr_2d(nm, nf_pre(1), skpre2); skpre2 = skpre1
      
    ! update the latest tori 
    nf_pre(1)  = nf;     c0_pre(1,:) = c0 
    rho_pre(1) = rho;  tp_pre(1) = tp; h0_pre(1) = h0
    call alloc_arr_2d(nm, nf, ckpre1); ckpre1 = ck 
    call alloc_arr_2d(nm, nf, skpre1); skpre1 = sk
  endif 
  
  ! --- save the Fourier Coefs in fcs.dat and the Fourier approximated curve in fcurve_refn_time.dat
  !     and the associated parameters in para_curve_cont, the tori in  ---- 
    
100 print*, 'Got a new curve! save fcs and the curve!' ! only the A^0 is needed 
  
  do i = 1, n, 1
    write(fcs, '(6e24.14)')   c0(i);     write(fcs,*) 
    write(fcs, '(10e24.14)')  ck(i, :);  write(fcs,*) 
    write(fcs, '(10e24.14)')  sk(i, :) 
    write(fcs,*);    write(fcs,*) ! two blank lines to seperate the components in phase vector
  end do
  
  ! Evaluate the truncated Fourier series for \hat T(xi)
  np  = nf2p1
  dxi = pi2/(np-1)
  do i = 1,  np, 1  
    xi = (i-1) * dxi
     
    call varphi( xi,     n, 0, nf, c0, ck, sk, pv_curv)
    call varphi( xi+rho, n, 0, nf, c0, ck, sk, pv_curv_rho)
    call gr_cj(pv_curv, cj)
    write(fcurve,  '(14e24.14)')   xi, pv_curv, pv_curv_rho, cj 
    
    ! if(i > 1) open(ftori, file=fntori_time_refn, access='append', status='old')
    if(i < np) call plob(pv_curv, 0.d0, tp, n, n, dsign(1.d0, tp), 1, ftori, deriv, gr_cj, pvf)
  enddo 
  
  write(ftori,*)
  write(fcurve,*);   write(fcurve,*);  

  write(fpara_curve, '(I5, 11e24.14, 1f5.0)', advance='no')  nf, rho, tp, h0, pi2/tp, rho/tp,  c0(1:n), dir
  do i = 1, n, 1
    write(fpara_curve, '(12e24.14)', advance='no')  ck(:, i),  sk(:, i)
  end do 
  write(fpara_curve,*);
    
  print*, 'Finish para_curve!';  print*; ! !read*
  
  if(ntori == 1) return 
  
  ! update ds is we have  a refined curve 
   if(niter > 8 ) then 
  
    if(dabs(ds) >= ds_min*2.d0)  ds = ds / 2.d0
    
  elseif(niter <= 3 .and. itori > 3) then 
  
    if(dabs(ds) <= ds_max /2.d0)  ds = 2.d0 * ds
    
  endif  
  
  
  ! update the index of torus and the step size 
  incr_time = 0
!  decr_time = 0
  itori   = itori + 1 
  
    call init_opt(opt1);    
    ! subroutine  gdg_pre_time( rho, c0, ck, sk, dg_pre, gr_map, deriv, gr_cj, gr_dcj_dx)
200  call init_nf(nf); call assign_ind  

     call  gdg_pre_time( rho, c0, ck, sk, dx_pre, gr_map, deriv, gr_cj, gr_dcj_dx)  
      
!    print*, 'dir, ds, |dx_pre|', dir, ds, dnrm2(len_dx_pre+1, dx_pre, 1); print*
!    write(*,'(10f20.10)') dx_pre; print*   

201 if( itori >= 3) then 
    ! we only need new c0, rho, tp, and h0 
    ! check the angle and the direction! 
    
    c0_pre(3, :) = c0  + dir*ds*dx_pre(1:na_all:nf2p1);
    if( opt == 11) then  
      rho_pre(3)   = rho + dir*ds*dx_pre(na_all+1); 
      tp_pre(3)    = tp  + dir*ds*dx_pre(na_all+2)
      h0_pre(3)    = h0 
    endif  
  
    call check_angle(c0_pre, rho_pre, tp_pre, h0_pre, isrestart, csangle)
    if(isrestart == 1)  then 
      ds = ds / 2.d0
      isrestart = 0 
    endif    
  endif 
 
  ! subroutine curve_pre_time( ds, dir)
  call curve_pre_time( ds, dir)
  
  print* 
  
enddo   


  return  
end subroutine tori_cont_time

!***********************************************************************
!     ****  check_angle   ****
!  check the angle of the last two curves, and the predicted one 
!  if it is small (less than 15 degree), we keep the result 
!  if not, we start from the previous curve with smaller ds
!  if the angle is greater than 90 degree, we reverse the direction of the guess dx 

!  ****** Input Variables ******
!  the parameters used to compute the angle between the last two curves and 
!  the current predicted or refined one 

!  c0_pre(3,n), rho_pre(3), tp_pre(3), h0_pre(3)

!  ****** Output Variables ******
!  isrestart     flag to restart or not      
!  dir           the continuation direction  


!  Finally Revised by Yu -- 20170616
!***********************************************************************
subroutine check_angle(c0_pre, rho_pre, tp_pre, h0_pre, isrestart, csangle)

use dp_mod
implicit none
 
! Input  and Output Declaration   
real(kind=dp), intent(in)    :: c0_pre(3, nm), rho_pre(3), tp_pre(3), h0_pre(3)
integer, intent(out)         :: isrestart     
real(kind=dp), intent(out)   :: csangle  
 
! Local Variable
integer :: i, k 
real(kind=dp)  :: curv3(3, 1+n)  
  
  isrestart = 0
!  debug = 0 
  
  ! the third rows for curv3 are the current refined curve should be put as the first row
  ! so we need to resort the 
  
  if(opt == 11) then 
    curv3(2:3, 1) = rho_pre(1:2) 
    curv3(2:3, 2) = tp_pre(1:2) 
    curv3(1, 1:2)   = (/rho_pre(3), tp_pre(3)/)
!    curv3(1:2, 3) = h0_pre(1:2)
  endif 
  
  
  k = 1
  do i = 1, n
    if(i == ind_c0) cycle
    curv3(2:3, 2+k) = c0_pre(1:2,i)
    curv3(1,   2+k) = c0_pre(3,i)
    k = k+1
  enddo

!  if(debug == 1) then  
    print*, 'check curv3 for angle computation befor curv_angle:'
    do i = 1, 3
      print*, curv3(i, :)
    enddo 
    print*; ! !read*
!  endif      
    
!  cos(14.984 ^o) = 0.966
  call curv_angle(curv3, 1+n, csangle) ! TODO: to debug csangle 
!  print*, 'csangle: ', csangle 
 
    
  if(dabs(csangle) < 0.96593) then 
    print*, 'The angle between the last three curves are two large > 15, restart!'; print*
    isrestart = 1
    
  elseif(csangle < 0.d0) then 
    dx_pre = -dx_pre   
    print*, 'reverse the direction of the kernel!'; print* 
  endif 

  return  
end subroutine check_angle


!***********************************************************************
!     ****   curve_pre_time   ****
!  Update the unknowns using the prediction  dx_pre by step size ds 
!  can be used for restart with a new ds 

!  ****** Input Variables ******
! ds      the step size, and the prediction vector is global module variable 

!  Finally Revised by Yu -- 20170616
!***********************************************************************
subroutine curve_pre_time( ds, dir ) !,  gr_map, deriv, gr_cj, gr_dcj_dx)
implicit none
 
! Input  and Output Declaration 
!integer, intent(in) :: opt1  
real(kind=dp), intent(in)      ::  ds, dir   

! Local Variable
integer    ::  i, nrow, k, ncol
real(kind=dp)  ::  dc0(nm), dx_pre0(len_dx_pre+1)  
real(kind=dp), allocatable, dimension(:,:)  ::  dck, dsk   

!  print*, 'check if |dx_pre| == 1', dnrm2(na_all+2, dx_pre,1); print*; !read*
  
  dx_pre0 = dir*ds*dx_pre 
  
!  call init_opt(opt1)
!  print*, 'opt == ', opt; print*; !read*
   
  call alloc_arr_2d(nm, nf, dck);  call alloc_arr_2d(nm, nf, dsk); 
  
  ! update c0, ck, sk and rho, tp according to the value of opt 
  do i = 0, m-1  
    do k = 1, n, 1
      nrow = i*n + k 
      ncol = 1 + (k - 1) * nf2p1 + i*na1 
      dc0(nrow)    = dx_pre0(ncol)                   ! -- d C0 
      dck(nrow, :) = dx_pre0(ncol+1    : ncol+nf)    ! -- d ck, k=1_nf
      dsk(nrow, :) = dx_pre0(ncol+nf+1 : ncol+2*nf)  ! -- d sk, k=1_nf 
    end do
  enddo   
      
!  print*, 'dx_pre(irow1) == 0?', dx_pre(irow1); print*
  
  print*, 'nf=', nf, 'shape, c0, dc0,ck,dck', shape(c0), shape(dc0), shape(ck), shape(dck);print*
  
  c0 = c0 + dc0; ck = ck + dck; sk = sk + dsk 

  if(opt == 11) then     
!    print*, 'Old rho, tp:', rho, tp 
    
    rho = rho + dx_pre0(na_all+1)
    tp  = tp  + dx_pre0(na_all+2)
      
!    print*, 'rho, tp, drho, dtp', rho, tp,  dx_pre(na_all+1:na_all+2)
!    print*; !read*
  endif

  print*, ' -----  New prediction ----------- '; print*
  print*, rho, tp, h0, c0(1:n);   print* 
  
  return  
end subroutine curve_pre_time




! ************************ refine_curve_time ***************************
! This routine is to refine the invariant curve using the  Time T map method,
! with the initial guess of the coefficients computed by gr_four

!  ** NOTE ** The refinement is successful either  |dx| < tol or |dF|  < tol 
!  for the case isref=1 but |dF|  > tol, we need to increase Nf, otherwise 

! opt=10  --  2017-06-04 10:47:00 
! No additional energy equation

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
subroutine refine_curve_time(n, nf, rho, c0, ck, sk, isref, iter, gr_map, deriv, gr_cj, gr_dcj_dx)

implicit none

! Input  and Output Declaration   
integer, intent(in)     ::  n, nf
real(kind=dp), intent(inout)    ::   rho  ! the prescribed energy, might be updated when opt = 11 
real(kind=dp), intent(inout)    ::   c0(nm), ck(nm,nf), sk(nm,nf) 
integer, intent(out)            ::  isref, iter  
  
! Local Variable
integer        :: info, num_incr, nrhs, ncol,nrow,  k,  i ! , j  
 
                      
real(kind=dp)  :: dg(nr_dg, nc_dg), ferr(nr_dg), dx(nc_dg), & ! Netwon method 
                  dc0(nm), dck(nm,nf), dsk(nm, nf), dx_all(na_all), & ! correction  
                  dxm_pre, dxm, dferrm,  &     ! termination control
                  dxi, xi, pt_approx(n)        !check the new curve
                  
                    
real(kind=dp), external :: dnrm2 
external :: gr_map,  deriv, eigrg,  gr_cj, gr_dcj_dx

  call init_nf(nf); call assign_ind
  
   
  debug = 0
  
!  print*, 'finish the assignment of ind_dg and ind_dg_pre!'; print*; !read*
   
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
    if(debug == 1) then 
      print*, 'befor gdg_refn_time, rho, n, m, nf,  size(c0, ck, sk):'
      print*, rho, n, m, nf, size(c0), size(ck), size(sk); print*;!read*
    endif 
    
     
    call  gdg_refn_time( rho, c0, ck, sk, dg, ferr, gr_map, deriv, gr_cj, gr_dcj_dx)
!    print*, 'maxnorm(ferr):', maxval(dabs(ferr))

    ! use dxm_pre to detect the modulus of the error, in principle, it should decrease 
    ! we allow it increases for two times, if it exceeds, the refinement fails 
    
    if( iter > 1)  dxm_pre = dxm

    if(opt == 10) then 
      ! the first torus 
      call deltx( nr_dg, nc_dg, nrhs, dg, ferr, dx, info)
      call deltx_qr(nr_dg, nc_dg, dg, ferr,  dx, info)


    elseif(opt == 11) then
     
      call deltx_qr(nr_dg, nc_dg, dg, ferr,  dx, info)
      
    endif 
    
    if( info /= 0 ) then 
      write(fwork,*)  'SSL Solver Failed to converge!';  !!read*  
      isref = 0
      return
      
    else   
 
    ! Terminate the iteration, if the change trend of the modulus of the dx fluctuates, stop
    ! in principle, the modulus of dx will decrease, we allow it increases twice as suggested by Gerard 
    
    ! Instead of the L2 norm, we look at the maximum norm ?? 
!      print*, 'Maxnorm|dx|:', maxval(dabs(dx))
      
      dxm    = dnrm2(nc_dg, dx,   1) / dsqrt( dble(nc_dg) )
      dferrm = dnrm2(nr_dg, ferr, 1) / dsqrt( dble(nr_dg) )
      
      print*, 'dxm=',dxm, 'derr=', dferrm, 'tol=', tol;   ! !read*  ; !!!read*    !ck
!      if(iter == 1) !read*

      ! error less that tolerance
      if( dferrm < tol_err) then
        print*, 'Error less than tol_err', dferrm, tol_err
        return  
      end if
      
      ! correction is small enough, stop the iteration and return also.....      
      if(dxm < 1.d-14) then 
        if(dferrm > tol_err*1.d1)  isref = 2 ! deal with this case by increasing nf and restart again.... 
        write(fwork,*) 'Modulus of correction less than tol! Succeed!', dxm 
        return
      endif
      
      
      ! allow the correction to increase at most two times 
      if(iter > 1 .and. dxm > dxm_pre) then 
        num_incr = num_incr + 1 
        
        if(num_incr >= 3) then ! allow the correction  dx to increase no more than twice... 
          isref = 0
          print*, 'Terminate refinement! |dx| has increased three times!'; !!read*  
          return 
        endif 
        
      endif  
        
      ! -- print the old values of the unknowns to check
      if(debug == 1) then 
        write(fwork,*) 'Old: X0 -- C0  \\  ck, | sk  for \\ (k = 1, ', nf , ')'
      
        ! write the coefficients C's and S's into file fcs.dat 
        write(fwork, *)  c0 
        print*; !!read*  
       
        do k = 1, nm, 1
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
      
      ! TODO: for the other cases... 
      dx_all(ind_dg)  = dx(1:len_ind_dg) 
     
      
      ! ---- update the correction dx for the corresponding control variables 
    do i = 0, m-1 
      
      
      do k = 1, n, 1
      
        ! the start index of k-th block (component), for C0
        nrow = n*i + k
        ncol = 1 + (k - 1) * nf2p1 + i*na1
 
        if(debug == 1)  print*, 'start column is ', ncol, 'for', k, '-th component'
        
        
        ! be careful with this part.... 
        ! -- d C0  
        dc0(nrow)    = dx_all(ncol)
        
        ! -- d ck, k=1_nf
        dck(nrow, :) = dx_all(ncol+1 : ncol+nf)
         
        ! -- d sk, k=1_nf   
        dsk(nrow, :) = dx_all(ncol+nf+1 : ncol+2*nf)
          
      end do
    enddo 
  
      !      ! *********** Time Map approach *************
!      ! add one more equtfaion H(\varphi(0)) - h = 0 in the last row  
!      
!      ! opt = 11: fix energy h, rho as one more unknown   
!      ! opt = 12: fix rho,    h as one more unknown 

      if(opt == 10 .or. opt == 11) then 
!        print*, 'Old rho, tp:', rho, tp 

        rho = rho + dx(nc_dg-1)
        tp  = tp  + dx(nc_dg)
        
!        print*, 'rho, tp, drho, dtp', rho, tp,  dx(nc_dg-1:nc_dg)
!        print*;  !!read*
      endif  
      
      
      ! --- write the correcitons in the working file to check ?? TODO 
               ! write the coefficients C's and S's into file fcs.dat 
      if(debug == 1) then 
        write(fwork, '(10e24.14)')   dc0 
        do k = 1, n, 1
          write(fwork, '(10e24.14)')  dck(k, :)
          write(fwork, '(10e24.14)')  dsk(k, :) 
        enddo
        print*;  !!read*  
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
        call varphi(xi, n, 0, nf, c0, ck, sk, pt_approx)
        write(fcurve_approx, *) xi, pt_approx
      end do
      write(fcurve_approx, *); write(fcurve_approx, *)
      
      close(fcurve_approx)
    endif   ! info == 1
   
  enddo  ! iteration   
  
end subroutine refine_curve_time


!*************************  gdg_all_time **************************** 
!  This routine computes the whole Jocobix matrix and the error matrix for systems of equations 
!  to be solved  by Newton Method, as input of deltx 
!                  dg *  dX = - ferr
!                 
!  For the nf values of xi_i, i = 0, ..., 2*nf 
!  
!  1.  the whole diffrential dg: D F/ d   A(C_0, C's, S's) + D F / rho, tp, h   
!  2.  the error:   ferr = - F   
!                   F(x) =  \varphi(xi_i + rho) -  \varphi(xi_i)    
!                            +   H(xi_0) - h0 = 0 


!  -- the Linear equation F to slove for Newton Method : 
!     F:  \varphi(xi_i + rho) - \varphi(xi_i) = 0 ],  i = 1,..., Nf   ! the invariance condition  
!          Note that the invariance conditon is for all the n components  
  
!  the unknowns are written in column vector, of dimension n *(2*Nf+1):
!  ---  X = [A_x, A_y, A_z, A_vx, A_vy, A_vz]  
  
!   for each xi, dg = d F / d A  is compute by subroutine gdg_inv_xi1 
!        
  
! ** NOTE ** 
!    TO make everything clear, we compute the full matrix to aviod confusion
!    Do the substraction of corresponding submatrix in the routine refine_curve_poinc

!       Input Variables 
!  h            the presribed vaule of energy 
! c0, ck, sk   the initial unkowns to be refined, in total is of dimension  n * (2*nf+1)


!       Output Variables 
! ferr_all    the error in F to be refined, dimension:   n*(2*nf+1) 
! dg_all      the differential of F w.r.t. the unkowns, 
!             dimension: n*(2*nf+1) - by -   n*(2*nf+1) -1       
! 
!        Module-Based Varaibles
! ndim, nf, nf2, nf2p1,  pi2, nr_dg, nc_dg, na_all, nr_dgall, nc_dgall 


!  Finally Revised by Yu -- 20170609


!------------------------------------------------------------------------
subroutine  gdg_all_time( rho, c0, ck, sk, dg_all, ferr_all, gr_map, deriv, gr_cj, gr_dcj_dx)
implicit none

! Input  and Output Declaration   
real(kind=dp), intent(in)       ::  rho,  c0(nm), ck(nm, nf), sk(nm, nf) 
real(kind=dp), intent(out)      ::  dg_all(nr_dgall, nc_dgall),  ferr_all(nr_dgall)  
external ::   gr_map, deriv, gr_cj, gr_dcj_dx   

! Local Variable
integer        :: i, im, row_st, debug 

real(kind=dp)  :: dxi, xi, dg_xi1(nm, na_all),  ferr_xi1(nm), &                 !dvphi_da(n, na1)
                  pvf(n), dpvf(n), dcj_dx(n), dcj_da(na1), dvphi_da(n, na1), & ! drho, d H/ d A 
                  dgdrho(n), cj,  pv(n), dp_dx(n,n)
                  
real(kind=dp)  :: step, dcj_da_h(na1), dcj_da_h2(na1), dcj_da_diff(na1), & ! debug dcj_da only dependent on A^0 !-- ckd
                  dg_dtp_h(nm), dg_dtp_h2(nm), dg_dtp_diff(nm),&  ! debug dcj_dt  -- ckd
                  dg_drho_h(n), dg_drho_h2(n), dg_drho_diff(n)   ! debug dcj_drho  -- ck

  debug = 0
  
  if(debug == 1) then 
    print*, 'check gdg_all_time: nr_dgall, nc_dgall: ', nr_dgall, nc_dgall
    print*, 'size dg_all, ferr_all', shape(dg_all), size(ferr_all) ;    print*; !read*
  endif 

  
  ! Initialize as zeros, then we only need to deal with nonzero components
  dg_all   = 0.d0 
  ferr_all = 0.d0 
  
!  ----- deal with argument xi_i one by one  ---------   

! From JM's paper: 
! discretisize the parameter space,  dxi = 2*pi/(2*nf+1), 
! and we need  to evaluate xi_i = i * dxi, where i = 0, ..., 2*nf 

  ! step size for xi 
  dxi   = pi2 / nf2p1 

  ! --- compute the full Jacobi Matrix and error ------ 
  
  do i = 1, nf2p1, 1
  
    ! the current xi to be evaluated
    xi = (i-1) * dxi
    
    ! the starting row for each block of xi_i, i = 0, ..., 2*nf 
    row_st =  (i-1) * nm
    
    
    !  subroutine gdg_da_xi1( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, gr_map) 
    !-- ckd !! 2017-06-12 20:48:40 
    call gdg_da_xi1( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, gr_map)
    
    ! error in F, for all the segment with the last on as invariance equations  
    ferr_all(row_st+1 : row_st+nm)  =  -ferr_xi1
    

    ! -- Jacobi Matrix for Newton Method --  
    ! the sub-block in dg associated with one value of xi, ndim rows starting from row_st
!    dg_xi1(nm, na_all), ferr_xi1(nm) 
    dg_all(row_st+1 : row_st+nm,  1:na_all) =  dg_xi1
    
    
    ! ************** deal with extra unknowns and equations ************ 
    
    ! the third last column  is  the differential w.r.t. rho, only the last column 
    ! dimension: n  for the last segment  
    call dvphi_drho_xi1(xi, rho, dgdrho) !-- ckd  
    dg_all(row_st+nm-n+1 : row_st+nm, nc_dgall-2) = dgdrho
     
     
    if(debug == 1) then  
    !--- debug dg_rho for the last segment --  ckd 2017-06-13 11:51:14 
      step = 1.d-3
      call dg_drho_debug(xi, step/2.d0, dg_drho_h)
      call dg_drho_debug(xi, step, dg_drho_h2)
      dg_drho_diff = (4.d0*dg_drho_h - dg_drho_h2)  / 3.d0 
      print*, 'dg_drho, centre difference, routine:'
      print*, dg_drho_diff;       print*, dgdrho 
      print*, 'diff in dg_drho, maxnorm:', maxval(dabs(dgdrho-dg_drho_diff)), &
               maxval(dabs(dgdrho-dg_drho_h)), maxval(dabs(dgdrho-dg_drho_h2))
      print*; !read*
    endif 
    !-------------------------------------------------------------
    
     
    ! the second last column is  the differential w.r.t. tp, associated to the map P( \varphi_m-1(\xi) ) 
    ! row_st for one xi 
    
    ! dg_dtp ! TODO: debug, something is wrong here. 
    do im = 0, m-1
      call varphi(xi, n, im, nf, c0, ck, sk, pv)
      call gr_map(pv,  pvf, 0, dp_dx) 
      call deriv(0.d0, pvf, n, dpvf)  
      ! since the derivative is w.r.t. TP/m, so we need to divde dpvf by m to obtain dg_dtp 
      dg_all(row_st+im*n+1 : row_st+im*n+n, nc_dgall-1) =  -dpvf/m 
    enddo    
    
     
     ! TODO : debug, something wrong here !
     if(debug == 1) then        
      ! -- debug dg_dtp ---   ckd! 2017-06-13 11:37:46 -------------------
      !    1.d-4, 4.1791015092940142E-012
      !    1.d-3, 3.5746960946880790E-012
      step = 1.d-4
      call dg_dtp_debug(xi, step/2.d0, gr_map, dg_dtp_h) 
      call dg_dtp_debug(xi, step,      gr_map, dg_dtp_h2) 
      dg_dtp_diff = (4.d0*dg_dtp_h - dg_dtp_h2) / 3.d0 

      print*, 'dg_dtp, centre difference, routine:'
      print*, dg_dtp_diff;     print*, dg_all(row_st+1 : row_st+nm, nc_dgall-1) 
      print*, 'diff in dg_dtp, maxnorm:', maxval(dabs(dg_dtp_diff-dg_all(row_st+1 : row_st+nm, nc_dgall-1))), &
               maxval(dabs(dg_dtp_h-dg_all(row_st+1  : row_st+nm, nc_dgall-1))), & 
               maxval(dabs(dg_dtp_h2-dg_all(row_st+1 : row_st+nm, nc_dgall-1)))
      print*; 
      
      ! TODO: something wrong there, debug!          
    endif 
    !-------------------------------------------------------------------
    
    
    ! if i==1, xi = 0, evaluate the energy equation here and the differential. 
    ! --- last row .. d H/ D C0, C1  = d H / d (x,y,z,vx,vy,vz)  * d (x,y,z,vx,vy,vz) / D C0, C1
    !     d H / drho = 0,  d H / d TP = 0 
    if(i == 1 ) then 
!      print*, 'check xi==0?, i', xi, i; print*; !read*  ! -- ck 
      
      ! energy equation 
      call varphi(xi, n, 0, nf, c0, ck, sk,  pv)
      call gr_cj(pv, cj)
      ferr_all(nr_dgall) = h0 - cj 
      
      print*, 'dcj, cj0, cj', h0-cj, h0, cj 
      print*; !!read*
       
      call gr_dcj_dx(pv,     dcj_dx)
      call dvphi_da_cal(xi,  dvphi_da)
      
      dcj_da = matmul(dcj_dx, dvphi_da)
      
      if(debug == 1) then 
        ! check dcj_da  !! --- ckd 2017-06-13 07:40:03 
        step = 1.d-4
        call dcj_da_debug(xi, step/2.d0, n, m, nf, c0, ck, sk, gr_cj, dcj_da_h) 
        call dcj_da_debug(xi, step,      n, m, nf, c0, ck, sk, gr_cj, dcj_da_h2) 
        dcj_da_diff = (4.d0*dcj_da_h - dcj_da_h2 ) / 3.d0 
        print*,'maxnorm in diff of dcj_da:', maxval(dabs(dcj_da_diff-dcj_da))
        print*;!read*
        
      endif 
      
      
      ! the last row associated to c0, ck, sk -- and eliminating three unknowns
      dg_all(nr_dgall, 1:na1) = dcj_da  ! the C0,  ck, sk the initial curve 
        
      ! the last two entries: dg /d rho = 0, dg / dtp = 0
        
    endif ! xi=0, for equation and differential w.r.t. h0 
     
    ! the last column is  the differential w.r.t. h0, only the last entry is nonzero as -1
    dg_all(nr_dgall, nc_dgall) = -1.d0
      
  enddo    
  
  return 
end subroutine gdg_all_time


!*************************  gdg_refn_time **************************** 
!  This routine extract from dg_all, ferr_all to obtain the coefficient matrix 
!  and error for the linear equations to be solved  by Newton Method, as input of deltx 
!        dg *  dX = - ferr  
 
! **Note** 
!  -ferr is returned as ferr_all from routine gdg_all_time 


!  According to the approach we select, which is specified by the value of opt:
!--  opt = 100
!    Alex's suggestions:  Remove energy equation, keep h free and T, rho fixed
!                         unknonws: c0,ck,sk  minus two for indeterminations
                          
!--  opt = 10:   ---- works --- 
!     JM's approach for first torus:  
!        include energy equation, eliminate c0_i, c1(s1)_j, c1(s1)_k
!        unknowns: tp, rho, h0(fixed)


!       Input Variables 
!  dg_all, ferr_all   the full matrix for dg and ferr  


!       Output Variables 
! dg, ferr   the Jacobi and error matrix to be refined by Newton's method 
!            dimension specified according to the value of opt  
! 
!        Module-Based Varaibles
! ndim, nf, nf2, nf2p1,  pi2, nr_dg, nc_dg

!  Finally Revised by Yu -- 20170611
! ---- TODO : check  
!------------------------------------------------------------------------
!subroutine  gdg_all_time( rho, c0, ck, sk, dg_all, ferr_all, gr_map, deriv, gr_cj, gr_dcj_dx)

subroutine  gdg_refn_time( rho, c0, ck, sk, dg, ferr, gr_map, deriv, gr_cj, gr_dcj_dx)
implicit none

! Input  and Output Declaration   
real(kind=dp), intent(in)       ::  rho,  c0(nm), ck(nm, nf), sk(nm, nf)  
real(kind=dp), intent(out)      ::  dg(nr_dg, nc_dg),  ferr(nr_dg)  
external ::   gr_map, deriv, gr_cj, gr_dcj_dx   

! Local Variable
integer :: debug
real(kind=dp)  ::   dg_all(nr_dgall, nc_dgall), ferr_all(nr_dgall) 
  
  ! Initialize as zeros, then we only need to deal with nonzero components 
  dg   = 0.d0 
  ferr = 0.d0 
  debug = 0
  
!  subroutine  gdg_all_time( rho, c0, ck, sk, dg_all, ferr_all, gr_map, deriv, gr_cj, gr_dcj_dx)
  call gdg_all_time( rho, c0, ck, sk, dg_all, ferr_all, gr_map, deriv, gr_cj, gr_dcj_dx)
  
  ! extract the necessary columns and rows from dg_all and ferr_all according to opt 
  
  ! columns related to c0,ck,sk  
!  print*, 'check the dimension of dg and dg_all:', shape(dg), shape(dg_all)
!  print*; !read*  !-- ckd 
  
  ! 2017-06-13 07:51:02  -- ckd! dg is full rank, rank(dg) = nc_dg 
  dg(1:nr_dg, 1:len_ind_dg)  =  dg_all(1:nr_dg, ind_dg)
  ferr(1:nr_dg)              =  ferr_all(1:nr_dg)
  
  ! columns related to rho, tp and h0     
  if(opt == 10 .or. opt == 11) then 
    ! unknonws rho and tp, take the third and second last columns in dg_all
!    print*, 'nc_dg == len_ind_dg +2', nc_dg, len_ind_dg; print*; !read*
     
    dg(1:nr_dg, nc_dg-1:nc_dg) = dg_all(1:nr_dg, nc_dgall-2:nc_dgall-1)
     
  endif 
   
  return 
end subroutine gdg_refn_time


!*************************  gdg_pre_time **************************** 
! The prediction of DG for the new curve to do continuation 

!  According to the approach we select, which is specified by the value of opt:
!--  opt = 100
!    Alex's suggestions:  Remove energy equation, keep h free and T, rho fixed
!                         unknonws: c0,ck,sk  minus two for indeterminations
  
 ! For the prediction, since we already know the rank is 2, do a little trick to fix the rank
 ! instead of using TOL as the precision  to decide the rank 
                          
!--  opt = 10:   ---- works --- 
!     JM's approach for the continuation:  
!      -- include energy equation  
!      -- eliminate c0_i to specify the curve 
!      -- fix energy h0 
       
!    the kernel of the new system is DF_p is two dimensional associated to  c1(s1)_j, c1(s1)_k
!    we select the one orthogonal to the direction corresponding to phase shift indetermination

!       Input Variables 
!  dg_all, ferr_all   the full matrix for dg and ferr  

!       Output Variables 
! dg_pre    the prediction of the particial system 
! 

!  Finally Revised by Yu -- 20170611
!------------------------------------------------------------------------

subroutine  gdg_pre_time( rho, c0, ck, sk, dx_pre0, gr_map, deriv, gr_cj, gr_dcj_dx)
implicit none

! Input  and Output Declaration   
real(kind=dp), intent(in)       ::  rho,  c0(nm), ck(nm, nf), sk(nm, nf)  
real(kind=dp), intent(out)      ::  dx_pre0(na_all+2)
external ::   gr_map, deriv, gr_cj, gr_dcj_dx   

! Local Variable
integer :: i,  dg_rank, jcol, k, dim_ker, ind_null, debug 

real(kind=dp)  ::  dg_all(nr_dgall, nc_dgall), ferr_all(nr_dgall), dg_pre(nr_dg, len_dx_pre), &
                   dg_base(len_ind_dg_pre+2, len_ind_dg_pre+2) ! pay attention to the dimension of dg_pre  
 
real(kind=dp)  ::  dcs_tlt(na_all), dck_tlt(nf), dsk_tlt(nf), cn(nf), sn(nf), pdct, pdct_mn
                  
real(kind=dp), allocatable :: dg_null(:,:)
real(kind=dp)  :: dnrm2

  debug = 0 
  
!  subroutine  gdg_all_time( rho, c0, ck, sk, dg_all, ferr_all, gr_map, deriv, gr_cj, gr_dcj_dx)
  call gdg_all_time( rho, c0, ck, sk, dg_all, ferr_all, gr_map, deriv, gr_cj, gr_dcj_dx)
  
  ! extract the needed columns related to c0, ck,sk, for the prediction of DF_P
  if(opt == 10 .or. opt == 11) then 
   
     dg_pre(1:nr_dg, 1:len_ind_dg_pre)  =  dg_all(1:nr_dg, ind_dg_pre)
    
     ! rho, tp as the last two columns
     dg_pre(1:nr_dg, len_ind_dg_pre+1:len_ind_dg_pre+2)  =  dg_all(1:nr_dg, nc_dgall-2:nc_dgall-1)
  endif 
  
  if(debug == 1) then 
    print*,'Check, the last 3 columns dg_pre, dg_all:'
    do i = 1, nr_dg,  1
      print*,  dg_pre(i,len_ind_dg_pre+1:len_ind_dg_pre+2), dg_all(i, nc_dgall-2:nc_dgall-1)
    end do
    print*;!read*
  endif   
      
  
  ! here dg_pre includes the columns associated to rho and tp 
  ! since we already know the rank is 2, do a fixed value for dg_rank here. 
  
  call null_cal(nr_dg, len_ind_dg_pre+2, dg_pre, dg_rank, dg_base)
  dim_ker = 2  ! fix it, we do not need to modify tol for different tori 
  
!  dim_ker = len_ind_dg_pre+2 - dg_rank
  
!  print*, 'dimension of the kernel:', dim_ker, 'dg_rank=', dg_rank
!  print*, 'if dim_ker ==2? if not, tune the value of tol in null_cal'; print*; !!read*
    
  call alloc_arr_2d(len_ind_dg_pre+2, dim_ker, dg_null)  
  dg_null = dg_base(:, len_ind_dg_pre+3-dim_ker:len_ind_dg_pre+2)
  
  
  ! take the one orthogonal to coefficient modulation, which can be obtained 
  ! by differentiating (1) w.r.t. \xi_0
  ! \tilta CK  = CK*cos(k\xi_0) - SK*sin(k\xi_0)
  ! \tilta SK  = CK*sin(k\xi_0) + SK*cos(k\xi_0)
  
  ! d \tilta CK / dxi_0 = -CK*sin(k\xi_0)*k - SK*cos(k\xi_0)*k
  ! d \tilta SK / dxi_0 =  CK*cos(k\xi_0)*k - SK*sin(k\xi_0)*k
  
  !  \xi_0 is the angle used rotated by rot_fcs to avoid phase shift indetermination for last torus 
  
  print*, 'for the last curve, xi0_rot=', xi0_rot  
  call trigrec_c(xi0_rot, nf, cn)  
  call trigrec_s(xi0_rot, nf, sn)  
  
  ! the same rule of memery saving forr dcs_tlt and d 
  do i = 1, n, 1
    ! the starting  index of i-th block  for i-th component
    jcol = (i-1)*nf2p1 + 1 
    dcs_tlt(jcol) = 0.d0  
    do k = 1, nf 
      dck_tlt(k) = - ck(i, k)* sn(k) * k  - sk(i, k)* cn(k) * k  
      dsk_tlt(k) =   ck(i, k)* cn(k) * k  - sk(i, k)* sn(k) * k
    enddo
    
    dcs_tlt(jcol+1    : jcol+nf)   =  dck_tlt
    dcs_tlt(jcol+nf+1 : jcol+2*nf) =  dsk_tlt
  end do
  
!  print*, 'check the dimension of differential of phase shift and dg_null:'
!  print*, 'na_all == len_ind_dg_pre?', na_all, len_ind_dg_pre; print*; !read*
  
  print*, 'check which kernel is orthogonal to the differential of phase shift'
  
  ind_null = 1; pdct_mn = 1.d3
  do i = 1, dim_ker, 1
!    print*, i, '-th kernel'
    
    pdct  = dot_product(dcs_tlt(ind_dg_pre), dg_null(1:len_ind_dg_pre,i))
    if(dabs(pdct) < dabs(pdct_mn)) then 
      pdct_mn = pdct; ind_null = i  
    endif  
     
    print*, pdct 
  enddo 
  
!  print*,'finish the pdct! choose ', ind_null, '-th column'; print*
!  !read*;  
  
  ! for better array operation, we return dx_pre0 as a full matrix, with the removed c0_k as zero 
  ! and leave the indeterminations avoidance in refine_curve_time 
  dx_pre0 = 0.d0 
  
  ! TODO, only deal with opt=10 and 11 at the moment  --  2017-06-15 09:35:30 
  if(opt == 10 .or. opt == 11) then 
    dx_pre0(ind_dg_pre) = dg_null(1:len_ind_dg_pre, ind_null)
    
    ! the third last and second last: drho, dtp
!    print*,'size(dx_pre0) = na_all+2?', size(dx_pre0), na_all+2 !-- ckd 
!    print*; !read*
    
    dx_pre0(na_all+1:na_all+2) = dg_null(len_ind_dg_pre+1:len_ind_dg_pre+2, ind_null)
    
!    print*, 'the update on rho and tp:'
!    print*, dg_null(len_ind_dg_pre+1:len_ind_dg_pre+2, ind_null)

  endif   
  
  dx_pre = dx_pre /dnrm2(na_all+2, dx_pre,1)
  
  
end subroutine gdg_pre_time


!***********************************************************************
!     ****   dg_dtp_debug   ****
!  check dg_dtp 

!  Finally Revised by Yu -- 20170613
!***********************************************************************
subroutine dg_dtp_debug(xi, step, gr_map, dg_dtp)
implicit none
 
! Input  and Output Declaration  
!integer,  intent(in)      :: im 
real(kind=dp), intent(in)      ::   xi, step 
real(kind=dp), intent(out)     ::   dg_dtp(nm)  
 
external :: gr_map 
 
! Local Variable
integer :: i, row_st
real(kind=dp) ::  tp0,tpp, tpm, dp_dx(n,n), pvfp(n),pvfm(n), pv0(n)  

  tp0 = tp 
  
  do i = 0, m-1  
  
    row_st = i*n 
    
    call varphi(xi, n, i, nf, c0, ck, sk, pv0)
    
    ! tpp, tpm 
    tpp = tp0 + step 
    tpm = tp0 - step 
  
    call init_time_tp(tpp)
    call gr_map(pv0,  pvfp, 0, dp_dx)  
  
    call init_time_tp(tpm)
    call gr_map(pv0,  pvfm, 0, dp_dx)
   
    dg_dtp(row_st+1 : row_st+n) = (- pvfp + pvfm ) / step / 2.d0
    
  enddo 
  
   
  call init_time_tp(tp0)    
  return  
end subroutine dg_dtp_debug


!***********************************************************************
!     ****   dg_drho_debug   ****
!***********************************************************************
subroutine dg_drho_debug(xi, step, dg_drho)
implicit none
 
! Input  and Output Declaration  
real(kind=dp), intent(in)      ::   xi, step 
real(kind=dp), intent(out)     ::   dg_drho(n)  
 
 
! Local Variable
real(kind=dp) ::  pv_rhop(n), pv_rhom(n) 

  call varphi(xi+rho+step, n, 0, nf, c0, ck, sk, pv_rhop)
  call varphi(xi+rho-step, n, 0, nf, c0, ck, sk, pv_rhom)
      
  dg_drho = (pv_rhop - pv_rhom) / step / 2.d0
   
  return  
end subroutine dg_drho_debug




end module curve_time_mod
