!  **** My tori Module*******
!  Compute and refine the invariant curve obtained by time T2 map, following several steps : 

!  -- 1. take a stroboscopic map on the invariant tori with a time interval t2 = 2*pi/omega2 
!  
!  -- 2. do Fourier analysis to obtain 2 dominant frequecies 

!  -- 3. using the general Fourier.f routine to compute the Fourier coefficients, which provide the initial seed for refinement

!  -- 4. Use Newton method to refine the coefficients obtained from the previous steps 
 
! ** Note ** 
!    1. use the general name for vector field and energy (Jacobi costant)
!    2. subroutines from libraries, need to declare external attribute **
!         external  ::    trigrec_c, trigrec_s   ! -- libgr.a 


module tori_mod

use dp_mod
use pi_mod 
implicit none
save 

! subroutines from libraries, declare external attribute, 
! this is only necessary for module to call external subroutines or functions 

external  ::   trigrec_c, trigrec_s   ! -- libgr.a 

! public ones - 
! shared by both methods : nf, ndim, nitmax, tol0, tol_err 
integer       ::  nf, nidm, nitmax, nf2, nf2p1 
real(kind=dp) ::  tol, tol_err 

! --- For declaration of array ----  
integer, private   ::  nc_dg, nr_dg, nc_dgall, nr_dgall, nc_dgxi1  !dimension
                        
! --- For save 
! to write the output of refinment process
integer, private :: fwork, fout 

! ---Module private variables--- 
!    ** must be initialized before any call of this module
! nf        number of Fourier modes
! nitmax    maximum of iteration for Newton method 
! tol       tolerance of correction used to terminate Netwon method 
! tol_err   tolerance of error usde to terminate Newton method

contains

!***********************************************************************
!     ****   init_fft  ****
! Initialize the commonly shared variables
! to use assignment, add 0 as the suffix to distinguish with the module private ones 
!***********************************************************************
subroutine init_tori( nf0, nitmax0, ndim0, tol0, tol_err0)

implicit none

! Input  and Output Declaration  
integer, intent(in)        ::  nf0 , nitmax0, ndim0   
real(kind=dp), intent(in)  ::  tol0, tol_err0   
 


! -- public  ---   
  ! commonly used value for dimension declaration 
  nf    = nf0
  nf2   = 2*nf  
  nf2p1 = nf2 + 1
  
  ! termination control for Newton method 
  tol     = tol0
  tol_err = tol_err0
  nitmax  = nitmax0 

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
! -- Local Variable -- 
  ! number of columns of rows in dg,  the lineare quation: dg * dx = -ferr
  nc_dg = 6*nf2p1
  nr_dg = nc_dg+1
  
  ! dimension of dgall, the number of column is  nc_dg + 3 (have 3 more unknowns)
  nc_dgall = nc_dg + 3
  nr_dgall = nr_dg 
  
  ! dimension of dg_xi1, the Jacobi matrix for one value of xi
  ! without h, we only have nc_dg + 2  columns, (two more unknowns: rou, t2)
  nc_dgxi1 = nc_dg + 2 
!  nr_dgxi1 = 6  ! discard, this is easy to understand, so it's not necessary to introduce another variable
  

!  for save - TODO: for the moment, not used
   fwork = 6
!  fwork = 222 
  open(fwork, file='fur.work')
  
  
  return  
end subroutine init_tori

  
! **********************************************************************
subroutine refine_curve( h, rou, t2, c0, ck, sk, isref )


! This routine is to refine the invariant curve, with the initial guess of the coefficients computed by   
! the subroutine main_refncurve. 

!**NOTE: Advanced course on long time integrations, P59, Title: Computation of a torus
!  The tori we are looking for are embedded in 2–parametric families, which can be paramtrized by 2 parameters among h, rou, T2. 
!  Therefore, in order to compute a torus, we  need to 
!   
!    1. eliminate one coordinate of C0, in order to fix a curve on the torus. Here we take C0_x 
!    2. set a coordinate of C1 (S1) equal to zero and eliminate it, in order to get rid of the phase shift indetermination. 
!       Here we take C1_ind(S1_ind), ind is specified by routine  'fix_index', according to the strategy in Josep-Maria's thesis
!    3. eliminate two unknonws among h, T2 and rou, in order to fix a particular torus.
!       here, we eliminate h and T2     
!       Comments: there are a lot of invariant tori around the torus we get, we  
!  
!  So we end up with a   nr ( = 1+6*nf2p1 ) - by -  nc (= 6*nf2p1 - 1) system of linear equations, 
!  with unique solution but which has more equations than uknowns. 
!  This is not a problem, as long as we use the general rutine we have described, specifying the kernel dimension to be zero

! TODO: why the kernel dimension is zero? 

!   we have to 
!   The linear equations to be solved:
!                  dg *  dX = - ferr
! 
!  1.  the whole diffrential dg: D F/ d [h, rou, t2, A(C_0, C's, S's) ]
!  2.  the error: ferr,  to be refined by Newton method 

!  -- the Linear equation F to slove for Newton Method : 
!     F =    [ H( \varphi(xi_0)) - h = 0;                                   ! energy 
!         F1: \varphi(xi_i + rou) - \phi_t2 (xi_i) = 0 ],  i = 1,..., Nf    ! the invariance condition 

!       Input Variables 
!  h       energy level to fix           

!       Input-Output Variables 
!  rou, t2, c0, ck, sk :  All the Fourier coefficients and rotation number / period 
      
!       Output Variables 
!  isref    flag to show if the refinement is successful or not           

! TODO-NO1: put the mostly often used parameters as the private parameters to share within all the subroutines inclued

! Finally Revised by Yu -- 20160801
!----------------------------------------

  implicit none

! Input  and Output Declaration   
  real(kind=dp), intent(in)       ::  h  ! the prescribed energy 
  real(kind=dp), intent(inout)    ::  rou, t2, c0(6),  ck(6,nf), sk(6,nf) 
  integer, intent(out)     ::  isref  
  
! Local Variable
  integer        :: ind, is_ck, & ! the index of unknowns C1 to eliminate 
                    iter, i, j, info, num_incr, nrhs,  &
                    ncol, k, kcol, kcol_end  
                       
  real(kind=dp)  :: dg(nr_dg, nc_dg), ferr(nr_dg), dx(nc_dg), & ! Netwon method 
                    drou, dt2, dc0(6), dck(6,nf), dsk(6,nf),  & ! correction  
                    dxm_pre, dxm, dferrm, &   ! termination control
                    dg_all(nr_dgall, nc_dgall)
                    
  real(kind=dp), external :: dnrm2 
  
  nrhs = 1  ! only one column in right-hand side 
   
! refine iteratively 
  iter = 0          ! index of iteration 
  num_incr = 0      ! number of times when the modulus correction increases

  !  subroutine plob(y0,t0,tf, n, tdir, ftag, ispl, deriv, gr_cj,  y) 
  !  TODO without plot, only compute the variational matrix 
  !      but we could save the orbit as the first try and plot the orbit to check the segments 
  !  TODO: note that to call gr_lf and gr__cjlf, we need to call lf_mod, and initialize 
  
  ! TODO: there are two possible way to get \phi_t: the time T map or the Poincare map 
  open(100, file = './dat/tori/curvall.dat' )  !ck 
  call  plob(pv0, 0.d0, t2, 42, 1, 100, 1, gr_lf, gr_cjlf,  pv_var) 
  
  pvf  = pv_var(1:6) ! \phi_t2( \varphi(xi)) = \phi_t2( pv_vphi)
  stm  = reshape(pv_var(7:42), (/6,6/)) 
  
  
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
    
    ! -- differential + error 
    ! subroutine  gdg_curve( nf, h, rou, t2, c0, ck, sk, ferr,  dg)  
    call  gdg_curve( h, rou, t2, c0, ck, sk, ferr, dg_all)  
    
    ! -- extract the corresponding sub-matrix from ferr, and dg_all 
    !    we fix the coordinates of c0_1( c0_x ) and c1_ind (is_ck=1) or s1_ind (is_ck=0)
    !    the actual Jacobi Matrix should be dg_all eliminate : 
    !     1.  the first column -- D F/ dh  -- set to a presribed energy
    
    !     2.  the 4-th column c0_1( c0_x ) -- set to zero 
    
    !     3.  if is_ck==1,  the 5 + (ind-1)*(2*nf+1) -th column  -- set to zero
    !            is_ck==0,  the 6 + (ind-1)*(2*nf+1) -th column  
    
    !         ind is the component of pv, the block (c0,ck,sk) starts from  4+(ind-1)*(2*nf+1),
    !         and c1 (s1) is the second (third) component in each block, index is 4+(ind-1)*(2*nf+1) + 1 (2)
    ! 
    ! subroutine fix_index( ck, sk, ind, is_ck)
    call fix_index(ck, sk, ind, is_ck) 
    
    ! we fix the energy, and eliminate h from the unknowns.  We also fix 
    dg(:, 1:2) = dg_all(:, 2:3)
    
    ! if is_ck = 1, keep s1_ind; elseif is_ck = 0,  keep c1_ind 
    ! TODO: be careful with this assignment, it is not always dg(:, 3), unless ind = 1
    
    if( ind == 1 ) then
      dg(:, 3) = dg_all(:,  5+is_ck)
    else
    ! we need to copy the all the columns from  5 to  4 + (ind-1)*(2*nf+1)
      ! number of columns to copy, ind-1 blocks + c0_ind 
      ncol = (ind-1)*nf2p1 
      
      dg(:, 3: 3+ncol)  = dg_all(:,  5:5+ncol)
      
      ! for the next block, deside if we keep s1_ind or c1_ind 
      dg(:, 4+ncol) = dg_all(:,  6+is_ck+ncol)
      
      ! copy the rest columns, till the last one
      dg(:, 5+ncol : nc_dg) = dg_all(:,  8+ncol : nc_dgall)      
      
    end if
    
    ! use dxm_pre to detect the modulus of the error, in principle, it should decrease 
    ! we allow it increases for two times, if it exceeds, the refinement fails 
    
    if( iter > 1)  dxm_pre = dxm
    
!    subroutine deltx( nr, nc, nrhs, a, b, x, info)
    call deltx( nr_dg, nc_dg, nrhs, dg, ferr, dx, info)

    if( info /= 0 ) then 
      write(fwork,*)  'SSL Solver Failed to converge!';  read*
      isref = 0
      return
      
    else   
 
    !  Terminate the iteration, if the change trend of the modulus of the dx fluctuates, stop
    !  in principle, the modulus of dx is decreasing, we allow it increases twice as suggested by Gerard 

      dxm    = dnrm2(nc_dg, dx,   1)
      dferrm = dnrm2(nc_dg, ferr, 1)
      
      print*, 'dxm=',dxm, 'derr=', dferrm, 'tol=', tol;  read* !ck
      
      !TODO: if ferr is small enough, stop the iteration and return also.....
      if( dferrm < tol_err ) then
        print*, 'Error less than tol_err', dferrm, tol_err
        isref = 1
        return  
      end if
      
      if(iter > 1 .and. dxm > dxm_pre) then 
      
        num_incr = num_incr + 1 
        
        if(num_incr >= 3 ) then ! allow the correction  dx to increase no more than twice... 
          isref = 0
          print*, 'Terminate refinement! |dx| starts increasing!'; read*
          return 
        endif 
          
      endif  

      if(dxm < tol) then 
        write(fwork,*) 'Modulus of correction less than tol! Succeed!', dxm 
        return
      endif  
    
      ! -- print the old values of the unknowns to check
      !    note that c0_x is fixed.... 
      write(fwork,*) 'Old: X0 --  rou, t2, c0_{y,z} \\  ck, sk_{x,y,z} \\ (k = 1, ', nf , ')'
      
      ! write the coefficients C's and S's into file fcs.dat 
      write(fwork, *)  rou, t2, c0(1), c0(2),  c0(3)
      do k = 1, nf, 1
        write(fwork, *)  ck(1,k), sk(1, k), ck(2,k), sk(2, k), ck(3,k), sk(3, k) 
      enddo
      
      ! -- print the correction to check
      write(fwork,*) 'Correction: dX --  drou, dt2, dc0_{x,y,z} \\  dck, dsk_{x,y,z} \\ (k = 1, ', nf , ')'
      
      drou = dx(1)
      dt2  = dx(2)
      
      ! -- C's and S's
      ! ** Note: one coordinate of C1 is fixed (of ind ) -- the same rule as dg 
      
      ! if is_ck = 1, keep s1_ind; elseif is_ck = 0,  keep c1_ind 
      ! TODO: be careful with this assignment, it is not always dg(:, 3), unless ind = 1
    
      if( ind == 1 ) then
         ! -- the x- component 
         if(is_ck == 1) then 
           dck(1, 1) = dx(3)
         else 
           dsk(1, 1) = dx(3)
         endif 
         
         ! if ind = 1,  for the first block (x component) in dX, c0_x and c1_x (or s1_x) are fixed 
         ! this block is of size (nf-1)*2, and c1_2 start at column 4 of dX  ( 4+(nf-1)*2 -1 = nf2 + 1)  
         dck(1, 2: nf) = dx(4: nf2p1 : 2)  
         dsk(1, 2: nf) = dx(5: nf2p1 : 2)
         
         ! -- the rest component 
         ! the rest blocks are of dimension nf2p1 (c0, ck, sk), starts at column nf2p1+1 + (k-2)*nf2p1 of dX 
         do k = 2, nf, 1
           kcol = nf2p1 + 1 + (k-2) * nf2p1
           dc0(k)    = dx(kcol)
           dck(k, :) = dx(kcol+1: nf2p1 : 2)
           dsk(k, :) = dx(kcol+2: nf2p1 : 2)
         end do
         
      else
        ! -- the x- component 
        !  c0_1 is already fixed, if ind \= 1, this block has 2*nf components, 
        !  from c1_1 to c1_nf, starts at column 3 of dX  ( 3+ nf*2 -1 = nf2+2) 
         
        dck(1, :) = dx(3: nf2+2 : 2)  
        dsk(1, :) = dx(4: nf2+2 : 2)
        
        ! -- the rest component before ind-th block, deal with as a whole block, starts at column nf2+3 of dX   
        do k = 2, ind-1, 1
           kcol     = nf2 + 3 + (k-2) * nf2p1
           kcol_end = kcol + nf2p1 - 1  
           dc0(k)    = dx(kcol)
           dck(k, :) = dx(kcol+1: kcol_end : 2)
           dsk(k, :) = dx(kcol+2: kcol_end : 2)
        end do 
        
        ! -- the  ind-th block, of size nf2, starts at column nf2 + 3 + (k-1) * nf2p1 of dX 
        kcol     = nf2 + 3 + (ind-2) * nf2p1
        kcol_end = kcol + nf2 - 1
        
        dc0(ind) = dx(kcol)
        
        if(is_ck == 1) then 
          dck(ind, 1) = dx(kcol+1)
        else 
          dsk(ind, 1) = dx(kcol+1)
        endif 
        
        dck(ind, 2:nf) = dx(kcol+2: kcol_end: 2)
        dsk(ind, 2:nf) = dx(kcol+3: kcol_end: 2)        
        
        !  --- the rest blocks: ind+1 : nf-th block 
        ! for each block, we have nf2p1 components to assign(c0+ck+sk) 
        ! starts at column   nf2 + 3 + (ind-2) * nf2p1 + nf2  = nf2p1*ind + 1
        
        ! --- TODO: do a subroutine for the standard assignment of a whole block, 
        do k = ind+1, nf, 1
           kcol     = nf2p1*ind + 1 + (k-ind-1) * nf2p1
           kcol_end = kcol + nf2p1 - 1  
           dc0(k)    = dx(kcol)
           dck(k, :) = dx(kcol+1: kcol_end : 2)
           dsk(k, :) = dx(kcol+2: kcol_end : 2)
        end do
        
        ! --- update the unknowns: rou + t2 + c0 + ck + sk
        rou = rou + drou 
        t2  = t2 + dt2
        c0  = c0 + dc0
        ck  = ck + dck
        sk  = sk + dsk 
        
        ! --- write the correcitons in the working file to check ?? TODO 
               ! write the coefficients C's and S's into file fcs.dat 
        write(fwork, *)  drou, dt2, dc0(1), dc0(2),  dc0(3)
        do k = 1, nf, 1
          write(fwork, *)  dck(1,k), dsk(1, k), dck(2,k), dsk(2, k), dck(3,k), dsk(3, k) 
        enddo
        
        end if
    
    endif   
   
  enddo    

  return  
end subroutine refine_time

!***********************************************************************
!     ****   fix_index   ****
!  Specify the index of unknowns to be fixed to avoid indeterminations 
  
!  There are two indeterminations: 1--invariant curve and 2--phase shift,
!  As suggested by Josep-Maria, we could fix one coordinate of C0 and C1, 
!  with index as ind0 and ind1, respectively 
    
! -- as a consequence, we have 3 + 6*(2*nf+1) - 3 =  6*(2*nf+1) unknowns ...
!    but we still have   1 + 6*(2*nf+1) equations.... 
    
!** Comment by Josep-Maria: Advanced course on long time integrations, P59, Title: Computation of a torus
!   This is not a problem, as long as we use the general rutine we have described, specifying the kernel dimension to be zero.

!-- PhD thesis, chapter 6, Methodolody, Section 6.2.5,   P100 

!-- Strategy to eliminate unknowns to avoid indeterminations
!   - 1. To avoid the indetermination in invariant curve
!        keep constant (to eliminate) one coordinate of C0 (exp, C0_x=0). For continuation, take care if this condition is valid. 

!   - 2. To avoid the indetermination in invariant curve, mostly for continuation cosideration
!       2.1. Eliminate C1_J or C1_J, where J is chosen in order to have || (C1_J, S1_J) ||_2  =  max_{i=1,...,n} || (C1_i, S1_i) ||_2
!            The reason for this choice is heuristic and tries to prevent asking conditions on harmonics of low amplitudes.

!       2.2  Between C1_J and B1_J, we select the coordinate for which min(|C1_J|, |S1_J| ) is achieved. This avoids the situation:
!            assume before the correction, C1_J=0 and |S1_J| = alpha. If we eliminate S1_J, we are forcing  
!            || (C1_J, S1_J) ||_2 > alpha and therefore we will not be able to continue if || (C1_J, S1_J) ||_2 > decreases
!            along the family. If |C1_J| is different from zero but small, we will be able to continue but the continuation step
!            will be artificially reduced.  

!** Another strategy is: Set C(k)_J = 0 as long as S(k)_J remains different from zero, that is  ( C(k)_J, S(k)_J) \= 0 
!   Usually , we take k=1 (the x component) and J has been chosen in order to maximize || ( C(k)_J, S(k)_J) ||_2. 

!   Comment: the first strategy is better, because it doesn't depend on a Fourier coefficients being different from zero 
!   and C1_1 = C1_x = 0, we need to introduce 2 variables as the indice of unknowns to fix the value 
    

!       Input Variables 
!  nf       see gdg_curve            
!  ck,sk    see gdg_curve 

!       Output Variables 
!  ind      the index of component in C1 to be fixed as zero.             

!  Routine Used:
!    none  

!  Finally Revised by Yu -- 20160731
!***********************************************************************
subroutine fix_index(ck, sk, ind, is_ck)

implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration 
!integer, intent(in)            ::  nf
integer, intent(out)           ::  ind, is_ck   
real(kind=dp), intent(in)      ::  ck(6, nf), sk(6, nf)  
 
! Local Variable
integer :: i, ind_loc(1)
real(kind=dp)  :: nm(nf)
  
  ! evaluate the coefficients of the x componet || (C1_J, S1_J) ||_2 for J = 1, ..., nf 
  do i = 1, nf, 1
    nm(i) = dsqrt( ck(1,i)**2 + sk(1,i)**2 )
  end do  
  
  ! take J that achieves  max { || (C1_J, S1_J) ||_2 } 
  ind_loc = maxloc(nm)
  
  ind = ind_loc(1)
  
  ! between C1_J and S1_J, select the one that achieves min ( |C1_J|, |S1_J|)
  if( dabs( ck(1,ind) ) .le. dabs( sk(1,ind) ) ) then 
    is_ck = 1
  else 
    is_ck = 0
  endif     
  
  return  
end subroutine fix_index

!*********************************************************************** 
!  This routine computes the coefficients and error for the linear equations to solve
!  by Newton Method, as input of deltx_lnx 
!                  dg *  dX = - ferr
!                 
!  For the nf values of xi_i, i = 0, ..., 2*nf 
! 
!  1.  the whole diffrential dg: D F/ d [h, rou, t2, A(C_0, C's, S's) ]
!  2.  the error: ferr,  to be refined by Newton method 

!  -- the Linear equation F to slove for Newton Method : 
!     F = [ H( \varphi(xi_0)) - h = 0;                                   ! energy 
!       F1: \varphi(xi_i + rou) - \phi_t2 (xi_i) = 0 ],  i = 1,..., Nf   ! the invariance condition  
!          Note that the invariance conditon is for all the  6 components (x,y,z,  vx,vy,vz)
  
!  the unknowns are written in column vector, of dimension 3 + 6 *(2*Nf+1) :
!  ---  X = [  h, rou, t2, A_x,  A_y, ... , A_vz ] 
  
! the component of dg =     
!     -- the first row,  d H /  d (h, rou, t2,  C_0 ) 
!        where, the nonzero components are:  d H/ dh = -1  and d H / d C_0
  
!     -- the rest row,  d F1 / d [h, rou, t2, A(C_0, C's, S's) ]
!        where d F1 / d h = 0, and the rest are computed by subroutine gdg_f1inv

!           [  d \varphi(xi_i+rou) / d rou,  -d \phi_t2 ( \varphi(xi_i) ) / dt2 ]
  
! -- 2. TO make everything clear, we will compute the full matrix to aviod confusion
!       Do the substraction of corresponding submatrix outside this subroutine. 

!       Input Variables 
!  nf           number of Fourier modes, take as initial guess 5 suggested by Gerard 
!  h            the presribed vaule of energy 
! rou           rotation number: 2pi / omega2 * omega1
!  t2           period along the frequecy omega2:   2pi / omega2 
!  c0, ck, sk   the initial unkowns to be refined, in total is of dimension 6 * (2*nf+1)


!       Output Variables 
!  ferr     the error in F to be refined, dimension:   1 + 6*(2*nf+1) 
! dgall     the differential of F w.r.t. the unkowns, 
!           dimension:  [ 1 + 6*(2*nf+1) ] - by -  [ 3 + 6*(2*nf+1) ]        
! 
!        Module-Based Varaibles
! pi2, nf, nf2p1, 
! nr_dg = 6 * ( 2*nf+1 ) + 1, nc_dg = 6 * ( 2*nf+1 ) 


!  Routine Used: -- TODO 
!     None

!  Finally Revised by Yu -- 20160731 ---- TODO : check
!----------------------------------------
subroutine  gdg_curve( h, rou, t2, c0, ck, sk, ferr_all, dg_all)

implicit none

! Input  and Output Declaration   
real(kind=dp), intent(in)    ::  h, rou, t2, c0(6), ck(6, nf), sk(6, nf) 
real(kind=dp), intent(out)   ::  ferr_all(nr_dgall),  dg_all(nr_dgall, nc_dgall) 
 
! Local Variable
integer        :: i, j, k, m, nf2p1, ncol, irow  
real(kind=dp)  :: dxi, xi, dg_xi1(6, nc_dgxi1 ),  ferr_xi1(6), &
                  h0, dh(6) 

real(kind=dp), parameter :: pi2 = 8.d0 * datan(1.d0)

  
  ! Initialize as zeros, then we only need to deal with nonzero components 
  dg_all = 0.d0 
  
  ! ---  The Jacobi Matrix of the general function F ----
 
  ! -- the first row of F is  H( \varphi(xi_0)) - h = 0 
  !    where xi_0   = C_0 ([x,y,z,vx,vy,vz])
  !    if we fix the presribed energy, just remove the first column 
 
  ! -- d H / d h
  dg_all(1, 1) = -1.d0  
  
  ! -- d H / d C_0 ^(x,y,z, vx,vy,vz)
  !    the index of C_0^(x,y,z, vx,vy,vz) :    4 + i* (2*nf+1), where i = 0, 5 
  
  !    C_0 =  H( \varphi(0)) ) is obtained  by substituding xi=0 into the truncated Fourier Series \varphi
  !    which is the actual value of (x,y,z, vx, vy,vz)_0 to compute the energy h
  
  call deriv_cjlf(c0, dh)
  dg_all(1, 4 : 4 + 5*nf2p1: nf2p1)  =  dh  
  
  !---- H0 - h 
  call gr_cjlf(c0, h0)
  ferr_all(1) = h0 - h
  
!  ----- deal with argument xi_i one by one  ---------   

! time step size for discretisize the parameter space, dxi = 2*pi/(1+2*nf), 
! and we need  to evaluate xi_i = i * dxi, where i = 0, ..., 2*nf 

  ! step size for xi 
  dxi   = pi2 / nf2p1
  
  ! the number of columns in dg, also the size od ferr 
  ncol  = 3 + 6*nf2p1
  
  do i = 1, nf2p1, 1
    ! the current xi to be evaluated
    xi = (i-1) * dxi
    
    ! the starting row for each block of xi_i, i = 0, ..., 2*nf 
    irow = 2 + (i-1)*6
    
!   subroutine gdg_inv( nf, xi, rou,  c0,  ck, sk, dg, ferr)  
    call gdg_inv_xi1( xi, rou, c0, ck, sk,  dg_xi1, ferr_xi1)
    
    ! the sub-block in dg associated with one value of xi
    dg_all(irow:irow+5, 2:ncol) = dg_xi1
    
    ! error in F (the invariance equations)
    ferr_all(2+(i-1)*6 : i*6+1 ) =  ferr_xi1 
    
  end do 
  
  return 
end subroutine gdg_curve


!***********************************************************************
!     ****   dvphi_da   ****
!   for a certain value of xi, the Jacobi Matrix of the invariant curve \varphi 
!   w.r.t. the  Fourier coefficients A.  
!   
!   we denote A = [ C0, C's, S's] as the Fourier coefficients, then A is of dimenson 6 * (2*Nf+1), write as components A = [A_x, A_y, A_z, A_vx, A_vy, A_vz] 
  
!  where \varphi as a truncated Fourier series.
!        \varphi (xi) = C_0 + Sigma_ k from 1 to Nf of { ( C_k cos(k*xi) + S_k cos(k*xi)) }
  
!       Input Variables 
!  xi       the argument for the invariant curve \varphi 
!  nf       number of Fourier modes
                    
!       Output Variables 
!  dgda      the differential (Jacobi Matrix) of \varphi w.r.t. the Fourier coefficients A            
 
!  Routine Used:
!   trigrec_c, trigrec_s 

!  Finally Revised by Yu -- 20160729
!***********************************************************************
subroutine dvphi_da_cal( xi, dgda)
implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration   
real(kind=dp), intent(in)      ::  xi  
real(kind=dp), intent(out)     ::  dgda(6, nf2p1) 
 
! Local Variable
integer        :: nf2p1, i, j, k, jcol
real(kind=dp)  :: cn(nf2), sn(nf2)

  print*;  print*, '-------- check dvphi/da ----------'
  
  dgda = 0.d0 ! initialization
  
 ! compute the trigometric function  xi 
  call trigrec_c( xi, nf2, cn)  
  call trigrec_s( xi, nf2, sn)
   
  
  ! The Jacobi Matrix of the general function F (the 6 coordinates) w.r.t. the A
  ! --- for 6 components, each is of dimension 6 * ( 2*Nf+1 ), the same block in the diagonal line, the rest components are zeros
  
  do i = 1, 6, 1
  
     ! the starting  index of i-th block  for i-th component
     jcol = (i-1)*nf2p1 + 1  
     
     dgda( i,  jcol ) = 1.d0       ! C_0
     
     ! the 2*nf components are cos(k*xi), sin(k*xi), k =1, ...,  nf
     do k = 1,  nf,  1 
       j = jcol + 1  + (k-1)*2
      
       print*, 'i,nf, jcol, k, j', i,nf, jcol, k, j  !ck
       read*                                   !ck 
       
       dgda( i,  j )   =  cn(k)  ! cos( k*xi ), k = j/2
       dgda( i,  j+1 ) =  sn(k)  ! sin( k*xi )
       
     enddo 
     
  end do
 
  return  
end subroutine dvphi_da_cal


!***********************************************************************
!     ****   gdg_inv_xi1  ****
!  for a certain value of xi, the Jacobi Matrix of the invariance condition (for all 6 components):  
!       F1:  \varphi(xi + rou) - \phi_t2 ( \varphi(xi) ) = 0 ],  i = 1,..., Nf
!  w.r.t. the  rotation number rou + period t2 along freq omega2 + the Fourier coefficients A                 

!  the unknowns are written in column vector, of dimension 2 + 6 *(2*Nf+1) :
!  ---  X = [ rou, t2, A_x,  A_y, ... , A_vz ] 
  
!  the component of dg =
!  -- the first 2 columns:   d F1 /  d ( rou, t2)  = 
!     [ d \varphi(xi+rou) / d rou,  -d \phi_t2 ( \varphi(xi_i) ) / dt2 ]
!       where the second component is actutally the vector field at  point \varphi(xi_i)  
!       (since the flow is the solution of the vector field, and the vector field is the derivative of the flow w.r.t. time )
  
!   --  the 3h column,  d F1 / d A_x ---  ( as a block of dimension 2*Nf + 1) 
!       -- we need to apply chain rule to compute  d \phi / d A 
!
!  d \varphi(xi+rou) / d A_x  - [ d \phi_t2 ( \varphi(xi) ) / d \varphi(xi) ] * d \varphi(xi) / d A_x
!            

! ** Note: we compute  dvphi(xi+rou) /drou and dvphi(xi+rou) /da together
!          then we only need to compute [ cos(k*xi+rou), sin(k*xi+rou) ] once  
!   
!   we denote A = [ C0, C's, S's] as the Fourier coefficients, then A is of dimenson 6 * (2*Nf+1), write as components A = [A_x, A_y, A_z, A_vx, A_vy, A_vz] 
  
!  where \varphi as a truncated Fourier series.
!        \varphi (xi+rou) = C_0 + Sigma_ k from 1 to Nf of { ( C_k cos(k*xi+rou) + S_k cos(k*xi+rou)) }
  
!       Input Variables  (initial guess)
!!  nf      number of Fourier modes -- module-based 
!  xi       the argument for the invariant curve \varphi 
!  rou      rotation number:  2pi / omega2 * omega1
!  t2       the period along the frequecy omega2 :   2pi / omega2
!  c0, ck, sk   dimension 6-by-nf,  the unknowns C_k and S_k (k=1, ..., nf), for each component
!           (x,y,z,vx,vy,vz), the values are different. 

                    
!       Output Variables 
!  dg      dimension: ndim(6), ndim * ( 2*nf+1 ) + 2. 
!          the differential (Jacobi Matrix) of invariance equations w.r.t.  rou + t2 + A
 
! ferr  dimesion 6, the error to be refined by Netwon method 
!          F1 = \varphi(xi_i + rou) - \phi_t2 (xi_i)  

!     Module based Varaibles
! ndim      dimension of the phase space  

!  Routine Used:
!   trigrec_c, trigrec_s 

! In order to make it more general, instead of 6, we put ndim
!  Finally Revised by Yu -- 20160730
!***********************************************************************
subroutine gdg_inv_xi1( xi, rou, c0, ck, sk, pvf, stm, dg_xi1, ferr_xi1)
implicit none
!integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration   
!integer, intent(in)            ::  nf 
real(kind=dp), intent(in)      ::  xi, rou, t2, c0(ndim), ck(ndim, nf), sk(ndim, nf), &
                                   pvf(ndim), stm(ndim, ndim)  ! the final state of phi_t(xi) 
                                   
real(kind=dp), intent(out)     ::  dg_xi1(ndim, nc_dgxi1 ), ferr_xi1(ndim) 
 
! Local Variable
integer        :: i, k, jcol
real(kind=dp)  :: cn(nf2), sn(nf2), cn_rou(nf2), sn_rou(nf2), &
                  pv0(ndim), vf0(ndim), pv0_rou(ndim),  &  
!                  pvf(ndim), stm(ndim,ndim), pv_var(42),  ! the integration of flow to compute the curve 
                  dvphi_da(ndim, nf2p1), dvphi_rou_da(ndim, nf2p1), dphi_da(ndim, nf2p1)

!  External subroutines                  
!external :: gr_cjlf, gr_lf  
  
  dg_xi1 = 0.d0 ! initialization
  
 ! the trigometric function at k*xi+rou 
 
 !  ---- discard ----- because my misunderstanding
 !  cos(k*xi + rou) = cos(k*xi) * cos(rou) - sin(k*xi) * sin(rou) ! discard
 !  sin(k*xi + rou) = cos(k*xi) * sin(rou) + sin(k*xi) * cos(rou) ! discard 
!  crou = dcos(rou)
!  srou = dsin(rou)
!  cn_rou = cn*crou - sn*srou 
!  sn_rou = cn*srou + sn*crou

  call trigrec_c( xi+rou, nf2, cn_rou)  
  call trigrec_s( xi+rou, nf2, sn_rou)
 
  ! ---- The Jacobi Matrix of the general function F (the 6 coordinates) w.r.t. the [rou, t2, A] ---
  
  ! -- 2nd column: d F / d t2 = d \phi ( \varphi(xi) ) / dt
  !    which is the vector field at \varphi(xi):  pv_vphi(6) -- position + velocity
  
  ! Evaluate the invariant curve \varphi(xi), with the initial guess of coefficients (c0, ck, sk)
  ! TODO: \varphi(xi) is given value or not... I think it's not, we need to recompute for each pv0 
  call varphi( xi, c0, ck, sk, pv0) 
  
  ! \varphi( xi+rou ) = pv0_rou
  call varphi( xi+rou, c0, ck, sk, pv0_rou) ! the truncated Fourier series 

         
  !  - d \phi_t2 (  \varphi(xi) ) / dt2 = - vector field at pv0
  call gr_lf(0.d0, pv0, 6, vf0)
  
  dg_xi1(:,2)  = -vf0 
 
  !  --- d F1 / d A =  d \varphi( xi + rou) / d A -  d \phi / d A 
  !  where we apply the chain rule to get: 
  !   d \phi / d A =  stm *  d \varphi(xi) / d A 
  !      where  stm  =  [ d \phi( xi ) / d \varphi(xi) ],  is the state transition matrix  
    
  ! the error to be corrected by Newton Method:  \phi_t2( pv_vphi) - \varphi(xi+rou)  
  ferr_xi1 = pv0_rou - pvf  
  
  !  d \varphi( xi+rou ) / d A 
  call dvphi_da_cal( xi+rou, dvphi_rou_da)
  
  !  --- compute d \phi / dA  ---------
  !  d \varphi(xi) / d A 
  call dvphi_da_cal( xi, dvphi_da)
  
  ! d \phi / dA =  *  d \varphi(xi) / d A
  dphi_da = matmul(stm, dvphi_da) 
  
  ! --- for 1st column + 3rd ~ rest columns, we need to handle component-by-component
  
  do i = 1, 6, 1  ! 6 components (x,y,z, vx, vy, vz)
  
     ! 1st column: d F1 / d rou  = d \varphi ( xi + rou) / d rou  = 
     !                           = Sigma from k=1 to nf of { -C_k * k*sn_rou(k) + S_k * k*cn_rou(k) }
     do k = 1, nf, 1
       dg_xi1(i, 1) = dg_xi1(i, 1) - ck(i, k) * k*sn_rou(k) + sk(i, k) * k*cn_rou(k) 
     end do
     
     ! d F1 / d A(x,y,z,vx,vy,vz):  start from the 3rd column 
     ! each is of dimension 6 * ( 2*Nf+1 ), the same block in the diagonal line, the rest are zeros
     
     ! the starting  index of i-th block, the dimension of each block is nf2p1  
     jcol = (i-1)*(nf2+1) + 3  
     
     print*, 'i, jcol', i, jcol   ! ck 
     read*  !ck
     
     ! ----- d \varphi( xi + rou) / d A - d \phi / d A -----
     dg_xi1( i,  jcol: jcol+2*nf ) =   dvphi_rou_da(i, :) - dphi_da(i, :)    
     
     ! -- ck --
     print*, 'dvphi_drou_da', dvphi_rou_da(i, :)
     print*, 'dphi_da', dphi_da(i, :) 
     print*; read* 
       
  end do
 
  return  
end subroutine gdg_inv_xi1 


!***********************************************************************
!     ****   varphi   ****
!  Compute the truncated Fourier series of the invariant curve \varphi(xi) 
!  given the coefficients and initial phases

!       Input Variables 
!  nf           number of Fourier modes
!  xi           the argument for the invariant curve \varphi 
!  c0, ck, sk   dimension 6-by-nf,  the unknowns C_k and S_k (k=1, ..., nf), for each component
!               (x,y,z,vx,vy,vz), the values are different. 

!       Output Variables 
!  pv           the state vector: position + velocity            

!  Routine Used:
!     trigrec_c, trigrec_s

! TODO: probably we will not need to compute pv using a standalone subroutine 
!       because at the point when pv is needed, we have already computed  cn, sn 

!  Finally Revised by Yu -- 20160730
!***********************************************************************
subroutine varphi( xi, c0, ck, sk, pv)
implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration 
real(kind=dp), intent(in)    :: xi, c0(6), ck(6, nf), sk(6, nf)   
real(kind=dp), intent(out)   :: pv(6) 
 
! Local Variable
integer :: i, k 
real(kind=dp)  :: cn(nf), sn(nf)
  
 ! compute the trigometric function at k*xi 
  call trigrec_c( xi, nf, cn)  
  call trigrec_s( xi, nf, sn)
  
  do i = 1, 6, 1
    pv(i) = c0(i)
    do k = 1, nf, 1
       pv(i) = pv(i) + ck(i,k) * cn(k) + sk(i,k) * cn(k)
    end do
  end do  

  return  
end subroutine varphi


end module tori_mod
