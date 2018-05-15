!***********************************************************************
!     ****   trace   ****
! Compute the trace of the monodromy matrix  M, given the six eigenvalues  
! Follow Alex's notes: stability for trace 
!  We know the basic property of Momodromy matrix, 

! -- 1. symplectic: det(M) = 1, M^T * \Omega * M = 1, where \Omega = [0, In;  -In, 0]
! -- 2.  1 / \lambda, \bar \lambda is also the eigenvalues  (the inverse and the conjugate)
! -- 3. define trace (\lambda) = \lambda + 1 / \lambda

! by the value of trace, we can get the stability behaviour of the P.O. 
!   if we have 
! -1. two couples of real eigenvalues,   |tr1|, |tr2| >= 2, (if lambda = 1, tr = 2)  
! -2. 1 pair of real, and 1 pair of complex(modulus 1) a+bi,
!       for real,  |tr1| > 2;  
!       for complex, tr2 = 2*a  ==> |tr2| = 2|a| <= 2 (cross x-axis, b=0 ==> |tr2| = 2) 
! -3. 2 couples of complex (of modulus 1), |tr1|, |tr2| <= 2
! -4. complex quadruple, tr1 = \bar tr2, what is the modulus of trace in this case???  
!  > 2, I suppose... but with the same real part  TODO: figure this out.... 
!     
! So, we do  
! -- 1. pre-process to move the pair of eigenvalue as 1 of multiplicity  2 
!   --- forget about this one, at the moment we try to compute all 6 traces, and have a look... 
! -- 2. compute the trace of the rest four eigenvalues and save them, forget the order at the moment 


!  ****** Input Variables ******
!  wr1, wi1      dimension n X 1, the real and imagnary part of the current point to sort 
!  n             dimension 

!  ****** Output Variables ******
!  tr            the trace of the 

 
!  Routine Used: 
!     sort_1d 

!  Finally Revised by Yu -- 20161226, the first day after Christmas
!***********************************************************************
subroutine trace( wr, wi, n,  tr_re, tr_im)

use dp_mod
!use sort_mod
implicit none
 
! Input  and Output Declaration   
integer, intent(in)            :: n 
real(kind=dp), intent(in)      :: wr(n), wi(n)    
real(kind=dp), intent(out)     :: tr_re(n), tr_im(n)
 
! Local Variable
integer        :: i, debug  
real(kind=dp)  :: a, b, temp, prsc

  ! set the debug flag 
  debug = 0  
  
  if(debug == 1)   then 
     print*, 'n=', n
     print*, 'Input of eigvalues: real // imag'
     print*,  wr
     print*,  wi
  endif 
    
  prsc    = 1.d-8
  
  do i = 1, n, 1
!    if(dabs( wi(i) ) > prsc) then  
      ! complex   (lambda = a+bi) 
      !  tr = lambda  + 1 / lambda = a+bi + 1 /(a-bi) 
      !     = a(1+ ( 1 / (a^2+b^2) ) ) + b(1 - ( 1 / (a^2+b^2) ) ) i
      
      a    = wr(i); b = wi(i)
      temp = 1.d0 / (a*a + b*b)
      tr_re(i) = a * ( 1.d0 + temp )
      tr_im(i) = b * ( 1.d0 - temp )
      
      if( dabs(tr_im(i)) < prsc )  tr_im(i) = 0.d0  ! if too small, treat as real 
!    else 
!      ! real 
!      ! tr = 2*a 
!      tr_re(i) =  wr(i) + 1.d0/wr(i)
!      tr_im(i) = 0.d0
!    endif 
    
  end do
  
  return
  
end subroutine trace 
  
!  ! --- romove lambda = 1 (1), deal with the rest n-2 eigenvalues
!  ! compute the distance of each eigenvalues from (1, 0i), take the 
!  ! two smallest as the two repeat lambda = 1 
!  do i = 1, n 
!    dst1(i) = dsqrt( (wr(i) - 1.d0 )**2 + wi(i)**2 )
!  enddo   
!  call sort_1d( dst1, n, 1, ind_dst1)
!  
!  ! remove the first 2, take them as the last 2 eigenvalues in wr/wi_st and ind_st 
!  
!  
!  ! the rest eigenvalues to be sorted
!  ind_stm2 = ind_dst1(3:n)
!  wrm2     = wr( ind_stm2 )
!  wim2     = wi( ind_stm2 )  
!  
!  ! -- sort by real(modulus) > complex (modulus), insertion sort --
!  ! ** NOTE ** CKD -- 2016-12-27 18:59:17 
!  if(is_one == 1) then 
!    
!    if(debug == 1) then 
!       print*, 'is_one == 1'; print*; read*
!    endif 
!    
!    wr_stm2 = wrm2;  wi_stm2 = wim2 ! make a copy 
!    ind_work = (/(i, i = 1, nm2)/)
!    
!    i = 0
!   
!    do  
!    
!      i = i + 1  
!      if( i > nm2) exit
!      
!      ind_cur = ind_work(i) ! the appropriate location for the current component
!      
!      j = i 
!      
!      if( i>1 .and. ind_cur > i)  j = ind_cur 
!      
!      do  
!        j = j + 1
!       
!        if(j > nm2) exit 
!        
!        if(debug == 1) then 
!          print*, 'Pre-sorting wr_st, wi_st: ' 
!          write(*, '(4f24.14)') wr_stm2(1:i) 
!          write(*, '(4f24.14)') wi_stm2(1:i) 
!       
!          print*; print*, '*** i, j, ind_cur***  ', i, j, ind_cur; print*
!          print*, 'i-/j- th component of wr, wi:' 
!          write(*, '(4f24.14)') wr_stm2(i), wr_stm2(j) 
!          write(*, '(4f24.14)') wi_stm2(i), wi_stm2(j)  
!          print*; read*
!        endif 
!        
!        ! *** i-th component is real ***
!        if( dabs( wi_stm2(i) )   <= prsc  ) then    
!          if ( dabs( wi_stm2(j) ) > prsc)  cycle           ! real > complex
!          
!          !-- both real, compare the modulus
!          if( dabs( dabs( wr_stm2(i) ) - dabs( wr_stm2(j) ) ) < prsc  ) then 
!            ! --if modulus the same, consider also the sign 
!            if(debug == 1) print*, 'Both real, same modulus!'
!            
!            if( wr_stm2(i) > wr_stm2(j) )   cycle  
!          
!          else 
!            if(debug == 1) print*, 'Both real, different modulus!'
!            if ( dabs( wr_stm2(i) ) >= dabs( wr_stm2(j) ) ) cycle
!            
!          endif  
!          
!        else 
!        
!          ! *** i-th component is complex ***   
!          if ( dabs( wi_stm2(j) )  > prsc )  then    !  complex < real 
!           
!            ! compare the modulus 
!            nmi = dsqrt( wr_stm2(i)**2 + wi_stm2(i)**2 ) 
!            nmj = dsqrt( wr_stm2(j)**2 + wi_stm2(j)**2 ) 
!            
!            if(debug == 1) print*, 'nmi, nmj, ds', nmi, nmj, nmi-nmj
!            
!            if (dabs(nmi- nmj) < prsc )  then 
!            ! modulus the same, compare the imaginary part,  take the smaller modulus as the bigger one, 
!            ! since it's closer to real axis.
!           
!              if(debug == 1)  print*, 'both complex, Same modulus' 
!              
!              ! check the real part, bigger modulus wins... 
!              if( dabs( wr_stm2(i) ) >  dabs( wr_stm2(j) ) )  cycle 
!              
!              if( dabs( dabs(wr_stm2(i))- dabs(wr_stm2(j)) ) < prsc ) then 
!                ! if modulus of real is the same, take positive real > negative
!                if(debug == 1) print*, 'Same modulus of real, check the sign'
!                if( wr_stm2(i) > wr_stm2(j)  )  cycle 
!                
!                if ( wr_stm2(i) == wr_stm2(j) ) then 
!                  ! if real is also the same, compare the imaginary part 
!                  if(debug == 1)  print*, 'We need to compare imaginary part' 
!              
!                  ! same modulus and same real part,  conjugate: positive real imaginary > negative imaginary 
!                  ! compare the imaginary part 
!                  if( dabs(wi_stm2(i)) < dabs( wi_stm2(j)) )  then 
!                    print*, 'Both complex, |wi(i)| <= |wj(i)|' 
!                
!                  elseif( dabs( dabs(wi_stm2(i)) - dabs( wi_stm2(j) ) ) < prsc) then 
!                    if(debug == 1) print*, 'imaginary part of same modulus' 
!                
!                    if( wi_stm2(i) > wi_stm2(j) )  then 
!                     if(debug == 1) print*, 'positve > negative imaginary part ' 
!                      cycle 
!                   endif
!                  endif  ! compare image 
!                endif 
!              endif  
!                 
!            else  
!               ! modulus different, bigger one wins  
!               if(nmi > nmj ) cycle
!            endif ! compare modulus
!              
!          endif 
!        endif 
!       
!        ! For the rest cases in real or complex, switch the components of 
!        ! All the previous elements from i to j-1, since they are in order  

!        ind_temp        = ind_work(j)   ! the appropriate location for the current component
!        ind_work(i+1:j) = ind_work(i:j-1)
!        ind_work(i)     = ind_temp
!        
!        if(debug == 1) print*, 'ind_work', ind_work
!          
!        wr_stm2 = wrm2(ind_work)
!        wi_stm2 = wim2(ind_work)
!        
!        if(debug == 1) then 
!          print*, wr_stm2
!          print*, wi_stm2 
!        endif 
!         
!        ind_cur = j 
!      end do
!    end do
!    
!    goto 110
!       
!  endif   
!    
!  ! ----------------------------------------------------------------  
!  ! -- sort by distance with the previous point   
!  ! the default order: 1, 2, ..., n-2, will be updated later 
!  ind_work = 0
!  
!  ! compute the distance matrix 
!  do i = 1, nm2, 1
!    ! the distance to the i-th component
!    do j = 1, nm2, 1
!      ds_1c(j) = dsqrt(  ( wrm2(j)-wr0(i) )**2 + ( wim2(j)-wi0(i) )**2  )
!    end do
!    
!    !subroutine sort_1d( a, na, order, ind_sorted) ! increasing order 
!    call sort_1d( ds_1c, nm2, 1, ind_1c)
!    
!    
!    ! save the sorted ds and in in array for future use 
!    ds_arr(:,  i)  = ds_1c
!    ind_arr(:, i)  = ind_1c 
!    
!  enddo 
!  
!  ! Still use insertion method to locate the index of i-th component in ind_st
!  ! In principle, take the closet eigenvalue in ind_1c(1), but in case of bifuraction, 
!  ! we can have ds_1c(1) = ds_1c(2), in this case, we take the next value ind_1c(1+1) 
!  
!  ! ** NOTE **  This can only be done after we compute the complete array ds_arr and ind_arr 
!  do i = 1, nm2, 1
!    ! the i-th column for the closet element to the i-th one in (wr0, wi0)
!    ind_1c = ind_arr(:, i)
!    ds_1c  = ds_arr(:, i)
!    
!    if(debug == 1) then   ! to check  --- no use with output this part 
!      print*;  print*, 'Sort ds and ind:'
!      print*, 'ind_1c', ind_1c 
!      print*, 'ds_1c',  ds_1c
!    endif    
!    
!    ! TODO, debug, iter = 25, the previous one is 1,2,3,4,5,6    
!    ! but this one we have:    1           5           6           5           2           3
!!   and 4-th is missing.
!! original eig: real // imag
!!   3.9284000350000001      -0.99583357480000001      -0.99583357480000001       0.25455656020000000        1.0000000000000000        1.0000000000000000     
!!   0.0000000000000000        9.1189315349999994E-002  -9.1189315349999994E-002   0.0000000000000000        2.6777886680000001E-007  -2.6777886680000001E-007



!   ! debug, iter = 49, 6 appears twice.... TODO:debug 
!   ! only one loop of check 
!!     5           6           3           4           2           1
!!      1           2           6           5           3           6     
!    ! Take the index of first closest element that is different from the sorted ones 
!    ! If ind_1c(1) is alreay taken by element in [1:i-1], we take the second one 
!    ! If the element in [i+1:n] is equal to ind_1c(1), we need to compare ds to 
!    ! determine if we need to take the successive one instead
!    
!    ind_temp = 1
!    isdup    = 0  ! the i-th component in ind_st is duplicated with the previous one[1:i-1] ?
!    
!    do 
!      isdup = 0 
!      
!      do j = 1, nm2, 1
!        
!        if(j == i) cycle 
!        
!        if(debug == 1) print*, 'i, j, ind_temp, ind_work(j), ind_1c(ind_temp)', i, j, ind_temp, ind_work(j), ind_1c(ind_temp)
!      
!        if( ind_1c(ind_temp) == ind_work(j) ) then 
!        
!          isdup = 1
!          if(j < i ) then 
!        
!            ind_temp = ind_temp + 1 
!          
!          else  
!            ! j > i  TODO: this part, I do not understand 
!            if( ds_1c(ind_temp) > ds_arr(ind_temp,j) )    ind_temp = ind_temp + 1
!            if(debug == 1) then 
!              print*, 'ds_1c', ds_1c
!              print*, 'ds_arr(ind_temp,*)',  ds_arr(ind_temp,:)
!            endif
!         
!          ! TODO: what if the two distance are equivalent?  --- there is not too much point in this part 
!          ! 2016-12-29 16:22:43   Discard this part, just keep the the distance comparision
!!          if( dabs( ds_1c(1)- ds_arr(1,j) ) < prsc )  then 
!!            ! if they have the same real part, we compare the imaginary part, by modulus? 
!!            ! otherwisze, we compare the real part 
!!            if( dabs( wr_st(ind_st(ind)) - wr_st(ind_st(j)) ) < prsc  ) then 
!!              if ( dabs( wi_st(ind_st(ind)) ) <  dabs( wi_st(ind_st(j))) ) ind = ind +1 
!!            else 
!!              if ( dabs( wr_st(ind_st(ind)) ) <  dabs( wr_st(ind_st(j))) ) ind = ind +1
!!            endif   
!!          
!!          else 
!!            if( ds_1c(1) > ds_arr(1,j) )    ind = ind + 1
!!          endif 
!          
!          endif  
!         
!        endif   
!      end do
!      
!      if (isdup == 0) exit 
!    enddo 
!    
!    ind_work(i) = ind_1c(ind_temp)
!    
!    if(debug == 1) then  ! TO check --- this one helps... I don't know why 
!      print*; print*, 'after sort:'
!      print*, 'ind_temp =', ind_temp, 'ind_work', ind_work 
!      print*;  read*
!    endif 
!    
!  end do
!  
!  
!   ! update the eigenvalues by the sorted index 
!110   ind_stm2      = ind_stm2(ind_work)  ! TODO
!      ind_st(1:nm2) = ind_stm2

!      wr_st(1:nm2) = wr(ind_stm2)                            
!      wi_st(1:nm2) = wi(ind_stm2)
!  
!  if(debug == 1) then  ! to check 
!    print*; print*, 'after sort:'
!    print*, 'ind_temp =', ind_temp, 'ind_st', ind_st 
!    print*;  read*

!    print*, 'Original wr, wi'
!    print*, wr
!    print*, wi
!    print*; read* 
!   
!    print*, 'sorted wr, wi'
!    print*, wr_st 
!    print*, wi_st 
!   endif 
!   
!!  print*,  'ind_st', ind_st   
!  return  
!end subroutine trace


!!***********************************************************************
!!     ****   pair   ****
!! check two eigenvalues to see if one is the inverse of the other 
!! if  

!!  ****** Input Variables ******


!!  ****** Output Variables ******


!!  Routine Used:
!!     

!!  Finally Revised by Yu -- 20160
!!***********************************************************************
!subroutine pair( in, out  optional)

!use dp_mod
!implicit none
! 
!! Input  and Output Declaration   
!real(kind=dp), intent(in)      ::    
!real(kind=dp), intent(out)     ::    
! 
!! Local Variable
!real(kind=dp)  ::
!  
!    
!  source

!  return  
!end subroutine pair











