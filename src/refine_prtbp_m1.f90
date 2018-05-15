! **********************************************************************
!  This two routines  refine_curve_poinc_m1 + fix_index, are used to avoid 
!  the indetermination in phase shift, fix c1_ind or s1_ind based on the strategy 
!  in Josep-Maria's thesis

!  since in this way we will have 1 less unknows than the equations, so we discard this approach at the moment.

! This routine is to refine the invariant curve, 
! with the initial guess of the coefficients computed by gr_four   
 

!**NOTE: Alex's note about the Poincare map method 
!  The tori we are looking for are embedded in 2–parametric families, which can be paramtrized by 2 parameters among h, rho   
!  The curve on the torus is fixed by thePoincare section, so we do not have curve indetermination

!  Therefore, in order to compute a torus, we  need to 
!    1. set a coordinate of C1 (S1) equal to zero and eliminate it, in order to get rid of the phase shift indetermination. 
!       Here we take C1_ind(S1_ind), ind is specified by routine  'fix_index', according to the strategy in Josep-Maria's thesis

!    2. eliminate two unknonws h and rho, in order to fix a particular curve,
!       by looking at the Poincare map representation of the torus, we can see we have bunch of 
!       invariant curve, which can be specified by the rotation number rho 
!       so, we eliminate h and rho       

!  So we end up with a   nr ( n*nf2p1 ) - by -  nc (= n*nf2p1 - 1) system of linear equations, 

!  with unique solution but which has more equations than uknowns. 
!  This is not a problem, as long as we use the general rutine we have described, specifying the kernel dimension to be zero

! TODO: why the kernel dimension is zero? 

!   we have to 
!   The linear equations to be solved:
!                  dg *  dX = - ferr
! 
!  1.  the whole diffrential dg: D F/ d A (C_0, C's, S's) 
!  2.  the error: ferr,  to be refined by Newton method 

!  -- the Linear equation F to slove for Newton Method : 
!     F =  \varphi(xi_i + rho) - P (xi_i) = 0 ],  i = 1,..., Nf    ! the invariance condition 

!       Input Variables 
!  h        perscribed energy level 
!  rho      the prescribed rotation number           

!       Input-Output Variables 
! c0, ck, sk:  All the Fourier coefficients to be refined 
      
!       Output Variables 
!  isref    flag to show if the refinement is successful or not           

! Finally Revised by Yu -- 20160826
!----------------------------------------
subroutine refine_curve_poinc_m1( h, rho, c0, ck, sk, isref, deriv, gr_cj, cj2v, dvind_dx)

implicit none

! Input  and Output Declaration   
real(kind=dp), intent(in)       ::  h  ! the prescribed energy 
real(kind=dp), intent(inout)    ::  rho, c0(n),  ck(n,nf), sk(n,nf) 
integer, intent(out)     ::  isref  
  
! Local Variable
integer        :: ind_cs1,  is_ck, & ! the index of unknowns C1 to eliminate 
                  iter, i, j, info, num_incr, nrhs,  &
                  ncol, k, kcol, kcol_end  
                       
real(kind=dp)  :: dg(nr_dg, nc_dgm1), ferr(nr_dg), dx(nc_dgm1), & ! Netwon method 
                  dc0(n), dck(n,nf), dsk(n, nf),  & ! correction  
                  dxm_pre, dxm, dferrm, &   ! termination control
                  dg_all(nr_dg, nc_dg+1)
                    
real(kind=dp), external :: dnrm2 
external :: deriv, gr_cj, cj2v, dvind_dx, poinc, deltx   

  nrhs = 1  ! only one column in right-hand side 
   
! refine iteratively 
  iter = 0          ! index of iteration 
  num_incr = 0      ! number of times when the modulus correction increases

  !  TODO: note that to call gr_lf and gr_cjlf, we need to call lf_mod, and initialize in the main routine 
  
  ! we do not need to compute ind_cs1 for every iteration, so compute it before the do loop   
  ! subroutine fix_index( ck, sk, ind_cs1, is_ck)
  call fix_index(ck, sk, ind_cs1, is_ck) 

  ! fix the corresponding coefficient to 0, and do not update during the Newton method 
  if(is_ck == 1) then 
    ck(ind_cs1, 1) = 0.d0
  else 
    sk(ind_cs1, 1) = 0.d0
  endif 
  
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
    
    ! -- differential + error - 
    ! subroutine  gdg_curve_poinc( h, rho, c0, ck, sk, ferr_all, dg_all, deriv, gr_cj)
    call  gdg_curve_poinc( h, rho,  c0, ck, sk,  ferr, dg_all, deriv, gr_cj, cj2v, dvind_dx)
    
    ! -- extract the corresponding sub-matrix from ferr_all, and dg_all 
    
    !    we fix the coordinates of  c1_ind (is_ck=1) or s1_ind (is_ck=0)
    !    the actual Jacobi Matrix should be dg_all eliminate : 
    
    !     1.  if is_ck==1,  the 2      + (ind-1)*nf2p1    -th column  -- set to zero
    !            is_ck==0,  the 2 + nf + (ind-1)*nf2p1    -th column  
    
    !     where ind is the component of pv for poincare section pv(ind) = p0, 
    !     the block (c0,ck,sk) starts from  1+(ind-1)*nf2p1, and has nf2p1 elements
    !     and c1 (s1) is the second (third) component in each block, index is 1+(ind-1)*nf2p1 + 1 (1+nf)
    ! 
    !  -- eliminate the column related to c1_{ind_cs1} if is_ck = 1 from the unknowns 
    !     and fix the value to be 0 
    
    ! when we update the Fourier coefficients, we need to do the opposite for the corresponding column
    ! if is_ck = 1, keep s1_{ind_cs1};  elseif is_ck = 0,  keep c1_{ind_cs1} 
    
    ! we need to copy the all the columns from  1 to  1+(ind_cs1-1)*nf2p1 
    ! if is_ck = 1, we eliminate c1_{ind_cs1},  copy the columns before and after this one, 
    ! just to skip this very column 
    ! 
  
    ! index of the first column for ind-th block - 1
    ncol = (ind_cs1 - 1) * nf2p1 + 1
    
    ! copy all the column before c1_{ind_cs1}
    if(is_ck == 1) then 
      dg(:, 1 : ncol)        = dg_all(:, 1 : ncol)
      dg(:, 1+ncol : nc_dg)  = dg_all(:, 2+ncol : nc_dg+1)
         
    else ! (is_ck == 0), eliminate the column numbered ncol+nf+1 
      dg(:, 1 : ncol+nf)       = dg_all(:, 1 : ncol+nf)
      dg(:, 1+ncol+nf : nc_dg) = dg_all(:, 2+ncol+nf : nc_dg+1)
    endif 
    
    ! use dxm_pre to detect the modulus of the error, in principle, it should decrease 
    ! we allow it increases for two times, if it exceeds, the refinement fails 
    if( iter > 1)  dxm_pre = dxm
    print*, 'nr_dg, nc_dg, ndim, nf', nr_dg, nc_dg, ndim, nf
    !read*  
!    subroutine deltx( nr, nc, nrhs, a, b, x, info)
     print*, 'dg '
    do i = 1, nr_dg, 1
      write(*,'(12f18.8)') dg(i, :)
    end do
    print*; !read*  
    
    call deltx( nr_dg, nc_dg, nrhs, dg, -ferr, dx, info)

    write(*, '(12f18.10)') dx 
    !!read*    
    
    if( info /= 0 ) then 
      write(fwork,*)  'SSL Solver Failed to converge!';  !read*  
      isref = 0
      return
      
    else   
 
    !  Terminate the iteration, if the change trend of the modulus of the dx fluctuates, stop
    !  in principle, the modulus of dx is decreasing, we allow it increases twice as suggested by Gerard 

      dxm    = dnrm2(nc_dg, dx,   1) / dsqrt( dble(nc_dg) )
      dferrm = dnrm2(nr_dg, ferr, 1) / dsqrt( dble(nr_dg) )
      
      print*, 'dxm=',dxm, 'derr=', dferrm, 'tol=', tol;   read*  ; !!read*    !ck
      
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
          print*, 'Terminate refinement! |dx| starts increasing!'; !read*  
          return 
        endif 
          
      endif  

      if(dxm < tol) then 
        write(fwork,*) 'Modulus of correction less than tol! Succeed!', dxm 
        return
      endif  
    
      ! -- print the old values of the unknowns to check
      !    note that c0_x is fixed.... 
      write(fwork,*) 'Old: X0 -- C0  \\  ck, | sk  for \\ (k = 1, ', nf , ')'
      
      ! write the coefficients C's and S's into file fcs.dat 
      write(fwork, *)  c0 
      print*; !read*  
       
      do k = 1, nf, 1
        write(fwork, *)  ck(:, k), sk(:, k) 
      enddo
      
      ! -- print the correction to check
      write(fwork,*) 'Correction: dX --  dc0 \\  dck, | dsk \\ (k = 1, ', nf , ')'
      
      
      ! -- Update C's and S's
      ! ** Note: one coordinate of C1 is fixed (of ind ) -- the same rule as dg 
      !          so do not update the corresponding coefficient, but should we put it to 0? when? 
      ! if is_ck = 1, keep s1_ind; elseif is_ck = 0,  keep c1_ind 
    
      print*, 'ind_cs1 = ', ind_cs1; !read*  
      
      do k = 1, n, 1
        ! the start index of k-th block (component), for C0
        ncol = 1 + (k - 1) * nf2p1 
          
        print*; print*, 'start column is ', ncol, 'for', k, '-th component'
        
        if(k == ind_cs1) then  ! we skip the c1 or s1
          
          if(is_ck == 1) then 
            dck(k, 1) = 0.d0
            if(nf >= 2) dck(k, 2:nf) = dx(ncol+1 : ncol+nf-1)
            
            dsk(k, :) = dx(ncol+nf : ncol+2*nf-1)
            
            print*, 'is_ck == 1' !ck

            
          else  
            dck(k, :) = dx(ncol+1: ncol+nf)
            
            dsk(k, 1) = 0.d0
            if(nf >= 2) dsk(k, 2:nf) = dx(ncol+nf+1 : ncol+2*nf-1)
            
            print*, 'is_ck == 0'  !ck
            
          endif  
          
           print*, 'dx(',ncol,':', ncol+2*nf-1,') = ', dx(ncol: ncol+2*nf-1)!ck 
          
          
               
        else
          ! for the components after ind_cs1-th one, copy the whole block
          if(k > ind_cs1)  ncol = ncol - 1  
          
          ! this part is the same for the blocks not equal to ind_cs1
          dck(k, :) = dx(ncol+1:  ncol+nf)
          dsk(k, :) = dx(ncol+1+nf:  ncol+2*nf)
          
          print*, 'dx(',ncol,':', ncol+2*nf,') = ', dx(ncol: ncol+2*nf)
        
        endif 
        !!read*    
        ! update C0  
        dc0(k) =  dx(ncol)
        
        print*, 'k-th component: c0, ck,sk: '
        print*, dc0(k), dck(k,:), dsk(k, :)
        !!read*    
        
      end do

      ! --- write the correcitons in the working file to check ?? TODO 
               ! write the coefficients C's and S's into file fcs.dat 
      write(fwork, *)   dc0 
      do k = 1, n, 1
        write(fwork, *)  dck(k, :), dsk(k, :) 
      enddo
      print*; !read*  
      
      ! --- update the unknowns:  c0 + ck + sk
      c0  = c0 + dc0
      ck  = ck + dck
      sk  = sk + dsk 
        
    endif   ! info == 1
   
  enddo  ! iteration   

  return  
end subroutine refine_curve_poinc

!***********************************************************************
!     ****   fix_index   ****
!  Specify the index of unknowns to be fixed to avoid indeterminations 
  
!  There are two indeterminations:    2--phase shift,
!  As suggested by Josep-Maria, we could fix one coordinate of C0 and C1, 
!  with index as ind0 and ind1, respectively 
    
! -- as a consequence, we have  ndim*(2*nf+1) - 1 =  ndim*nf2p1- 1 unknowns ...
!    but we still have   1 + ndim*(2*nf+1) equations.... 
    
!** Comment by Josep-Maria: Advanced course on long time integrations, P59, Title: Computation of a torus
!   This is not a problem, as long as we use the general rutine we have described, specifying the kernel dimension to be zero.

!-- PhD thesis, chapter 6, Methodolody, Section 6.2.5,   P100 

!-- Strategy to eliminate unknowns to avoid indeterminations

!   - 2. To avoid the indetermination in phase shift, mostly for continuation cosideration
!       2.1. Eliminate C1_J or C1_J, where J is chosen in order to have || (C1_J, S1_J) ||_2  =  max_{i=1,..., ndim} || (C1_i, S1_i) ||_2
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
!  ind      the index of component in C1 to be fixed as zero.             

!       Module-based Varaibles
!  ndim, nf   for purpose of size declaration

!  routine Used:
!    none  

!  Finally Revised by Yu -- 20160731
!***********************************************************************
subroutine fix_index(ck, sk, ind_cs1, is_ck)

implicit none
integer, parameter :: dp = kind(1.d0)   

! Input  and Output Declaration 
!integer, intent(in)            ::  nf
integer, intent(out)           ::  ind_cs1, is_ck   
real(kind=dp), intent(in)      ::  ck(n, nf), sk(n, nf)  
 
! Local Variable
integer :: i, ind_loc(1)
real(kind=dp)  :: nm(n)
  
  ! evaluate the coefficients of the x componet || (C1_J, S1_J) ||_2 for J = 1, ..., nf 
  do i = 1, n, 1
    nm(i) = dsqrt( ck(i, 1)**2 + sk(i, 1)**2 )
  end do  
  
  write(*,*) nm
  
  ! take J that achieves  max { || (C1_J, S1_J) ||_2 } 
  ind_loc = maxloc(nm)
  
  ! to extract the value of ind_loc 
  ind_cs1 = ind_loc(1)
  
  
  print*, 'ind_cs1 =',  ind_cs1
  !!read*    
  
  ! between C1_J and S1_J, select the one that achieves min ( |C1_J|, |S1_J|)
  if( dabs( ck(ind_cs1, 1) ) .le. dabs( sk(ind_cs1, 1) ) ) then 
    is_ck = 1
  else 
    is_ck = 0
  endif     
  
  return  
end subroutine fix_index

