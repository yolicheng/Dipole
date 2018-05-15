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
subroutine refine_curve_poinc(n, nf, rho, c0, ck, sk, isref, iter, gr_poinc, gr_cj, gr_dcj_dx)

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
external :: gr_poinc,  eigrg,  gr_cj, gr_dcj_dx

   debug = 1
   isv  =  1
   call init_nf(nf)
   
!   print*, 'Check refine_curve_poinc:'
!   print*, 'n,nf, nr_dg, nc_dg, nc_dgall, rho: ', n, nf, nr_dg, nc_dg, nc_dgall, rho
!   print*;  read*
   
   
!   write(*, *) '# Before refinemetn c0, ck, sk : n, nf = ', n, nf 
!    do i = 1, n, 1
!      write(*, '(6e24.14)')   c0(i);       ! write(*,*) 
!      write(*, '(10e24.14)')  ck(i, 1:10); ! write(*,*) 
!      write(*, '(10e24.14)')  sk(i, 1:10) 
!      write(*,*);  
!    end do
!    print*; read*
    
    
!  The first thing to do is to choose which parameter to fix in order to 
!  avoid the indetermination in phase shift, opt == 

!  subroutine fix_index(ck, sk, ind_cs1, is_ck)
  call fix_index(ck, sk, n, nf, ind_cs1, is_ck)
  print*, 'To avoid phase shift indetermination, fix', ind_cs1, '-th Four coef, is_ck = ', is_ck 
!  print*; read* 
  
  ! Do phase shift for opt=3 poincare map approach to avoid phase shift. 
  if(opt == 3) then 
    call rot_fcs(ck, sk, ind_cs1, is_ck, n, nf)
  endif 
    
  nrhs = 1  ! only one column in right-hand side 
  
  isref = 1 ! by default, we succeed the refinement  
   
  ! counter for iteration and increase of the modulus of correction 
  iter = 0          ! index of iteration 
  num_incr = 0      ! number of times when the modulus correction increases
  
  nsqrt =  dsqrt( dble(nc_dg) )
  
 
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
     
    call  gdg_curve_poinc( rho,  c0, ck, sk, dg, ferr, gr_poinc, gr_cj, gr_dcj_dx)
    
    
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
    
!     print*; print*, 'Check dg and ferr as output( should not change) after calling deltx: '
!      print*, 'dg = '
!      do i = 1, nr_dg, 1
!        write(*,'(10e20.10)') dg(i, :)
!      end do
!      print*; read*  
!      
!      print*, 'ferr = '
!      write(*,'(10e20.10)') ferr
!      print*; read*
      
    if( debug == 1) then 
      print*, 'the correction : dx'   
      write(*, '(8e24.14)') dx 
      print*; read* 
    endif      
    
    if( info /= 0 ) then 
      write(fwork,*)  'SSL Solver Failed to converge!';  !read*  
      isref = 0
      return
      
    else   
 
    ! Terminate the iteration, if the change trend of the modulus of the dx fluctuates, stop
    ! in principle, the modulus of dx will decrease, we allow it increases twice as suggested by Gerard 
    
    ! Instead of the L2 norm, we look at the maximum norm ?? 
    
      dxm    = dnrm2(nc_dg, dx,   1) / nsqrt
      dferrm = dnrm2(nr_dg, ferr, 1) / nsqrt
      
      print*, 'dxm=',dxm, 'derr=', dferrm, 'tol=', tol;  ! read*  ; !!read*    !ck
      
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
      ! 
      ! we have to recover to this same format to apply the following assignment
      if(opt == 3) then  
      ! **** For poincare map approach ******* 
       ! the index of C1(s1)_{ind_cs1} in the vector of all the unknowns 
        irow   = (ind_cs1-1)*nf2p1 + (1-is_ck)*nf + 2 
        dx_all(irow ) = 0.d0 
        
        dx_all(1: irow-1)        = dx(1: irow-1)
        dx_all(irow+1: nc_dgall) = dx(irow : nc_dg)
        
      elseif(opt == 1) then 
        ! remove  of C1(s1)_{ind_cs1}, one to avoid phase shift, one the avoid falling back to p.o. 
        irow   = (ind_cs1-1)*nf2p1 + (1-is_ck)*nf + 2 
        irow2  = (ind_cs1-1)*nf2p1 + is_ck*nf + 2
        
        dx_all(irow ) = 0.d0 
        dx_all(irow2) = 0.d0
        
        if(irow > irow2) then 
          dx_all(1: irow2-1)        = dx(1: irow2-1)
          dx_all(irow2+1: irow-1)   = dx(irow2: irow-2)
          dx_all(irow+1: nc_dgall)  = dx(irow-1: nc_dg)
        else 
          dx_all(1: irow-1)         = dx(1: irow-1)
          dx_all(irow+1: irow2-1)   = dx(irow: irow2-2)
          dx_all(irow2+1: nc_dgall) = dx(irow2-1: nc_dg)
        endif 
       
      elseif (opt == 2) then 
      
        irow2  = (ind_cs1-1)*nf2p1 + is_ck*nf + 2
        dx_all(irow2) = 0.d0
        
        dx_all(1: irow2-1)         = dx(1: irow2-1)
        dx_all(irow2+1: nc_dgall)  = dx(irow2: nc_dg)
      endif 
      
      
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
      
        rho = rho + dx(nc_dgall)
        print*, 'rho, drho', rho, dx(nc_dgall)
        print*; read*
        
      elseif(opt == 12) then 
        h0_curv = h0_curv + dx(nc_dgall)  
        call init_curve_h0(h0_curv)
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


!  routine Used: -- TODO 
!     None

!  Finally Revised by Yu -- 20160731
! ---- TODO : check 20160826
!------------------------------------------------------------------------
subroutine  gdg_curve_poinc( rho, c0, ck, sk, dg, ferr, gr_poinc, gr_cj, gr_dcj_dx)
implicit none

! Input  and Output Declaration   
real(kind=dp), intent(in)       ::  rho,  c0(n)  
real(kind=dp), intent(out)      ::  dg(nr_dg, nc_dg),  ferr(nr_dg)  
real(kind=dp), intent(inout)    ::  ck(n, nf), sk(n, nf)  
! Local Variable
integer        :: i, irow, irow2, k, ncol 
real(kind=dp)  :: dxi, xi, dg_xi1(n, nc_dgall),  ferr_xi1(n), dg_all(nr_dgall, nc_dgall) 
real(kind=dp)  :: xi0, pv_xi0(n), h_xi0, dcj_dx(n) 


external ::   gr_poinc, gr_cj, gr_dcj_dx 
  
  ! Initialize as zeros, then we only need to deal with nonzero components 
  dg = 0.d0 
  
  
!  ----- deal with argument xi_i one by one  ---------   

! time step size for discretisize the parameter space, dxi = 2*pi/(1+2*nf), 
! and we need  to evaluate xi_i = i * dxi, where i = 0, ..., 2*nf 

  ! step size for xi 
  dxi   = pi2 / nf2p1
  
!  print*, 'before gdg_inv_xi1_poinc:'; print*; read* 
  
  do i = 1, nf2p1, 1
  
    ! the current xi to be evaluated
    xi = (i-1) * dxi
    
    ! the starting row for each block of xi_i, i = 0, ..., 2*nf 
    irow = 1 + (i-1) * n
    
     print*, 'i=', i, 'xi=', xi, 'starting row:', irow; print*; read*
        
    !  subroutine gdg_inv_xi1_poinc( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, gr_poinc)
    call gdg_inv_xi1_poinc( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, gr_poinc)
    
    ! TODO: the call of cj2v within gr_poinc may return isv=0, 
    ! that is vy^2 < 0, how to deal with this case?
    
    if(debug == 1)  print*, 'after gdg_inv_xi1:  ferr_xi1 = ', -ferr_xi1;  !read* 
!    
!    print*, 'dg_xi1 = '
!    do k =  1, n, 1
!      write(fwork, '(12f18.8)') dg_xi1(k,:)
!    end do  
!    print*; read*    
    
    ! the sub-block in dg associated with one value of xi, ndim rows starting from irow
    dg_all(irow : irow+n-1,  :) =  dg_xi1
    
    ! error in F (the invariance equations)
    ferr(irow : irow+n-1 )      =  -ferr_xi1
  end do 


! *********** Time Map approach *************
!      ! add one more equtaion H(\varphi(0)) - h = 0 in the last row  
!      
!      ! opt = 11: fix energy h, rho as one more unknown   
!      ! opt = 12: fix rho,    h as one more unknown 

  if( opt == 11 .or. opt == 12 ) then
    ! the last row is the energy : H(\varphi(0)) - h  = 0 
    ! add one more row for the extra energy equation, to correct ferr_h =  h - H..

    xi0 = 0.d0 
    call varphi( xi0, n, nf, c0, ck, sk, pv_xi0) 
    call gr_cj(pv_xi0, h_xi0)
    ferr(nr_dgall) = h0_curv - h_xi0 
    
    ! --- last row for dg_all, d H/ D C0, C1  = d H / d (x,y,z,vx,vy,vz)  * d (x,y,z,vx,vy,vz) / D C0, C1
!    subroutine deriv_cjlf(pv, dcj)
    call gr_dcj_dx(pv_xi0, dcj_dx)
    
    print*, 'dcj_dx', dcj_dx 
    do k = 1, n, 1
      ncol = nf2p1*(k-1)+1 
      dg_all(nr_dgall,ncol:ncol+1) =  dcj_dx(k)
    enddo 
    
    ! opt: 11,  d H/ d rho = 0, no need to do assignment since it's initialized to be zero 
    ! opt: 12,  d (-h)/ dh = -1
    
    if(opt == 12)   dg_all(nr_dgall, nc_dgall) = -1.d0 
    
  end if 
  
  if(opt == 2) then 
    ! fix a constant value  
    irow2 = (ind_cs1-1)*nf2p1 + is_ck*nf + 2 
    
    dg(:, 1 : irow2-1)     = dg_all(:, 1 : irow2-1)
    dg(:, irow2 : nc_dg)   = dg_all(:,  irow2+1 : nc_dgall)
    
  else if( opt == 3) then 
  ! delete the irow-th unknowns from dx, and irow-th column from dg  
      
    ! the index of C1(s1)_{ind_cs1} in the vector of all the unknowns 
    irow = (ind_cs1-1)*nf2p1 + (1-is_ck)*nf + 2 
    
    dg(:, 1 : irow-1)    = dg_all(:, 1 : irow-1)
    dg(:, irow : nc_dg)  = dg_all(:,  irow+1 : nc_dgall)
  
  else if(opt == 1) then 
   ! remove another unknown  
    irow  = (ind_cs1-1)*nf2p1 + (1-is_ck)*nf + 2 
    irow2 = (ind_cs1-1)*nf2p1 + is_ck*nf + 2 
    
    print*, 'irow, irow2:', irow, irow2; print*; read*
    
    if(irow > irow2) then 
      dg(:, 1 : irow2-1)     = dg_all(:, 1 : irow2-1)
      dg(:, irow2 : irow-2)  = dg_all(:,  irow2+1 : irow-1)
      dg(:, irow-1 : nc_dg)  = dg_all(:,  irow+1 : nc_dgall)
    
    else 
      dg(:, 1 : irow-1)      = dg_all(:, 1 : irow-1)
      dg(:, irow : irow2-2)  = dg_all(:,  irow+1 : irow2-1)
      dg(:, irow2-1 : nc_dg) = dg_all(:,  irow2+1 : nc_dgall)
    
    endif   
  endif 
  
  
  if(debug == 1) then   
  print*, 'finish dg computation!'
  print*, 'DG, dimension: ', shape(dg)
  
  
  do i = 1, nr_dg, 1
    write(fwork, '(10e20.10)') dg(i,:)
  end do
!  print*; read*
  
  print*, 'ferr, dimension: ', size(ferr)
  write(fwork, '(10e20.10)') ferr 
!  print*; read*
  
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
subroutine gdg_inv_xi1_poinc( xi, rho, c0, ck, sk, dg_xi1, ferr_xi1, gr_poinc )
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
  debug = 1      ! -- ckd 
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
        write(fwork, '(10f22.14)') dp_dx(i, :)
      end do
    
!      det = dp_dx(1, 1)*dp_dx(2, 2) - dp_dx(1, 2)*dp_dx(2, 1)
!      print*, 'det (D P / D X)  = ', det 
!      print*; read*  
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
  
  dg_xi1(:, 1:na_all) = dvphi_rho_da - dp_da 
  
  if(opt == 11) then 
    call dvphi_drho_xi1(xi, rho, dgdrho)
    dg_xi1(:, nc_dgall) = dgdrho
  endif 


  
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
!    write(fwork, '(12f18.12)')  dg_xi1(i, :)
!  enddo 
!  print*; read* 
!  
  return  
end subroutine gdg_inv_xi1_poinc  
