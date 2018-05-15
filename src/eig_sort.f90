!***********************************************************************
!     ****   eig_sort   ****
! Sort the eigenvalues (real+imag in complex form) of dimension n-by-n  
! although we deal with the Hamiltonian system, where the eigenvalues come in pairs 
! we do not apply this property here 

! ** NOTE ** 
! Initialize the ind_st array, to avoid random assignment with value 1,...,6
! 

! if np = 1, we sort by the standard that    --- ind_st = 1,..,6
!   pure real first (by modulus ), complex(imag/=0) (amplitude )

! if np >= 2, sort by the nearest distance with the eigenvalues of the previous point 
!   use a n-by-n distance matrix, and an index vector 
!  --- ind_st = 0

!  ****** Input Variables ******
!  wr0, wi0      dimension n X 1, the real and imagnary part of the previous point
!  wr1, wi1      dimension n X 1, the real and imagnary part of the current point to sort 
!  is_one        flag for the number of points to sort, 
!                1: only 1 point (np=1 approach),  others: more than 2 points  (np>=2 option) 

!  ****** Output Variables ******
!  wr1_st, wi1_st    dimension n X 1, the real and imagnary part of the sorted current point 
!  ind_st            the ordered index, output for further check 
!  Routine Used: 
!     sort_1d 

!  Finally Revised by Yu -- 20161226, the first day after Christmas
!***********************************************************************
subroutine eig_sort( wr0, wi0, wr, wi, n, is_one, wr_st, wi_st, ind_st)

use dp_mod
use sort_mod
implicit none
 
! Input  and Output Declaration   
integer, intent(in)     :: n, is_one 
integer, intent(out)     :: ind_st(n)
real(kind=dp), intent(in)      ::  wr0(n), wi0(n), wr(n), wi(n)    
real(kind=dp), intent(out)     ::  wr_st(n), wi_st(n)  
 
! Local Variable
integer  :: i, j, ind_arr(n,n), ind_1c(n), ind_temp, ind_cur, debug, isdup 
real(kind=dp)  :: ds_1c(n), nmi, nmj, ds_arr(n,n), prsn 


  ! set the debug flag 
  debug = 0  
  
  if(debug == 1)   then 
     print*, 'is_one = ', is_one , 'n=', n
     print*, 'Previous one:'
     print*,  wr0 
     print*,  wi0
  endif 
  
    
  ind_cur = 1 
  prsn    = 1.d-8
  
  ! -- sort by real(modulus) > complex (modulus), insertion sort --
  ! ** NOTE ** CKD -- 2016-12-27 18:59:17 
  if(is_one == 1) then 
    
    if(debug == 1) then 
       print*, 'is_one == 1'; print*; read*
    endif 
    
    wr_st = wr;  wi_st = wi ! make a copy 
    ind_st = (/(i, i=1,n)/)
    
    i = 0
   
    do  
    
      i = i + 1  
      if( i > n) exit
      
      ind_cur = ind_st(i) ! the appropriate location for the current component
      
      j = i 
      
      if( i>1 .and. ind_cur > i)  j = ind_cur 
      
      do  
        j = j + 1
       
        if(j > n) exit 
        
        if(debug == 1) then 
          print*, 'Pre-sorting wr_st, wi_st: ' 
          write(*, '(6f24.14)') wr_st!(1:i) 
          write(*, '(6f24.14)') wi_st!(1:i) 
       
          print*; print*, '*** i, j, ind_cur***  ', i, j, ind_cur; print*
          print*, 'i-/j- th component of wr, wi:' 
          write(*, '(6f24.14)') wr_st(i), wr_st(j) 
          write(*, '(6f24.14)') wi_st(i), wi_st(j)  
          print*; read*
        endif 
        
        ! *** i-th component is real ***
        if( dabs( wi_st(i) ) <= prsn  ) then    
          if ( dabs( wi_st(j) ) > prsn)  cycle           ! real > complex
          
          !-- both real, compare the modulus
          if( dabs( dabs( wr_st(i) ) - dabs( wr_st(j) ) ) < prsn  ) then 
            ! --if modulus the same, consider also the sign 
            if(debug == 1) print*, 'Both real, same modulus!'
            
            if( wr_st(i) > wr_st(j) )   cycle  
          
          else 
            if(debug == 1) print*, 'Both real, different modulus!'
            if ( dabs( wr_st(i) ) >= dabs( wr_st(j) ) ) cycle
            
          endif  
          
        else 
        
          ! *** i-th component is complex ***   
          if ( dabs( wi_st(j) )  > prsn )  then    !  complex < real 
           
            ! compare the modulus 
            nmi = dsqrt( wr_st(i)**2 + wi_st(i)**2 ) 
            nmj = dsqrt( wr_st(j)**2 + wi_st(j)**2 ) 
            
            if(debug == 1) print*, 'nmi, nmj, ds', nmi, nmj, nmi-nmj
            
            if (dabs(nmi- nmj) < prsn )  then 
            ! modulus the same, compare the imaginary part,  take the smaller modulus as the bigger one, 
            ! since it's closer to real axis.
           
              if(debug == 1)  print*, 'both complex, Same modulus' 
              
              ! check the real part, bigger modulus wins... 
              if( dabs( wr_st(i) ) >  dabs( wr_st(j) ) )  cycle 
              
              if( dabs( dabs(wr_st(i))- dabs(wr_st(j)) ) < prsn ) then 
                ! if modulus of real is the same, take positive real > negative
                if(debug == 1) print*, 'Same modulus of real, check the sign'
                if( wr_st(i) > wr_st(j)  )  cycle 
                
                if ( wr_st(i) == wr_st(j) ) then 
                  ! if real is also the same, compare the imaginary part 
                  if(debug == 1)  print*, 'We need to compare imaginary part' 
              
                  ! same modulus and same real part,  conjugate: positive real imaginary > negative imaginary 
                  ! compare the imaginary part 
                  if( dabs(wi_st(i)) < dabs( wi_st(j)) )  then 
                    print*, 'Both complex, |wi(i)| <= |wj(i)|' 
                
                  elseif( dabs( dabs(wi_st(i)) - dabs( wi_st(j) ) ) < prsn) then 
                    if(debug == 1) print*, 'imaginary part of same modulus' 
                
                    if( wi_st(i) > wi_st(j) )  then 
                     if(debug == 1) print*, 'positve > negative imaginary part ' 
                      cycle 
                   endif
                  endif  ! compare image 
                endif 
              endif  
                 
            else  
               ! modulus different, bigger one wins  
               if(nmi > nmj ) cycle
            endif ! compare modulus
              
          endif 
        endif 
       
        ! For the rest cases in real or complex, switch the components of 
        ! All the previous elements from i to j-1, since they are in order  

        if(debug == 1)  print*, 'ind_st', ind_st
        
        ind_temp = ind_st(j)   ! the appropriate location for the current component
        ind_st(i+1:j) = ind_st(i:j-1)
        ind_st(i) = ind_temp
        
        if(debug == 1) print*, 'ind_st', ind_st
          
        wr_st = wr(ind_st)
        wi_st = wi(ind_st)
        
        if(debug == 1) then 
          print*, wr_st
          print*, wi_st 
        endif 
         
        ind_cur = j 
      end do
    end do
    
    return
       
  endif   
    
  ! ----------------------------------------------------------------  
  ! -- sort by distance with the previous point   
  ! the default order: 1, 2, ..., n, will be updated later 
  ind_st = 0
  
  ! compute the distance matrix 
  do i = 1, n, 1
    ! the distance to the i-th component
    do j = 1, n, 1
      ds_1c(j) = dsqrt(  ( wr(j)-wr0(i) )**2 + ( wi(j)-wi0(i) )**2  )
    end do
    
    !subroutine sort_1d( a, na, order, ind_sorted) ! increasing order 
    call sort_1d( ds_1c, n, 1, ind_1c)
    
    ! save the sorted ds and in in array for future use 
    ds_arr(:,  i)  = ds_1c
    ind_arr(:, i)  = ind_1c 
    
  enddo 
  
  ! Still use insertion method to locate the index of i-th component in ind_st
  ! In principle, take the closet eigenvalue in ind_1c(1), but in case of bifuraction, 
  ! we can have ds_1c(1) = ds_1c(2), in this case, we take the next value ind_1c(1+1) 
  
  ! ** NOTE **  This can only be done after we compute the complete array ds_arr and ind_arr 
  do i = 1, n, 1
    ! the i-th column for the closet element to the i-th one in (wr0, wi0)
    ind_1c = ind_arr(:, i)
    ds_1c  = ds_arr(:, i)
    
    if(debug == 1) then   ! to check  --- no use with output this part 
      print*;  print*, 'Sort ds and ind:'
      print*, 'ind_1c', ind_1c 
      print*, 'ds_1c',  ds_1c
    endif    
    
    ! TODO, debug, iter = 25, the previous one is 1,2,3,4,5,6    
    ! but this one we have:    1           5           6           5           2           3
!   and 4-th is missing.
! original eig: real // imag
!   3.9284000350000001      -0.99583357480000001      -0.99583357480000001       0.25455656020000000        1.0000000000000000        1.0000000000000000     
!   0.0000000000000000        9.1189315349999994E-002  -9.1189315349999994E-002   0.0000000000000000        2.6777886680000001E-007  -2.6777886680000001E-007



   ! debug, iter = 49, 6 appears twice.... TODO:debug 
   ! only one loop of check 
!     5           6           3           4           2           1
!      1           2           6           5           3           6     
    ! Take the index of first closest element that is different from the sorted ones 
    ! If ind_1c(1) is alreay taken by element in [1:i-1], we take the second one 
    ! If the element in [i+1:n] is equal to ind_1c(1), we need to compare ds to 
    ! determine if we need to take the successive one instead
    
    ind_temp = 1
    isdup    = 0  ! the i-th component in ind_st is duplicated with the previous one[1:i-1] ?
    
    do 
      isdup = 0 
      
      do j = 1, n, 1
        
        if(j == i) cycle 
        
        if(debug == 1) print*, 'i, j, ind_temp, ind_st(j), ind_1c(ind_temp)', i, j, ind_temp, ind_st(j), ind_1c(ind_temp)
      
        if( ind_1c(ind_temp) == ind_st(j) ) then 
        
          isdup = 1
          if(j < i ) then 
        
            ind_temp = ind_temp + 1 
          
          else  
            ! j > i  TODO: this part, I do not understand 
            if( ds_1c(ind_temp) > ds_arr(ind_temp,j) )    ind_temp = ind_temp + 1
            if(debug == 1) then 
              print*, 'ds_1c', ds_1c
              print*, 'ds_arr(ind_temp,*)',  ds_arr(ind_temp,:)
            endif
         
          ! TODO: what if the two distance are equivalent?  --- there is not too much point in this part 
          ! 2016-12-29 16:22:43   Discard this part, just keep the the distance comparision
!          if( dabs( ds_1c(1)- ds_arr(1,j) ) < prsn )  then 
!            ! if they have the same real part, we compare the imaginary part, by modulus? 
!            ! otherwisze, we compare the real part 
!            if( dabs( wr_st(ind_st(ind)) - wr_st(ind_st(j)) ) < prsn  ) then 
!              if ( dabs( wi_st(ind_st(ind)) ) <  dabs( wi_st(ind_st(j))) ) ind = ind +1 
!            else 
!              if ( dabs( wr_st(ind_st(ind)) ) <  dabs( wr_st(ind_st(j))) ) ind = ind +1
!            endif   
!          
!          else 
!            if( ds_1c(1) > ds_arr(1,j) )    ind = ind + 1
!          endif 
          
          endif  
         
        endif   
      end do
      
      if (isdup == 0) exit 
    enddo 
    
    ind_st(i) = ind_1c(ind_temp)
    
    if(debug == 1) then  ! TO check --- this one helps... I don't know why 
      print*; print*, 'after sort:'
      print*, 'ind_temp =', ind_temp, 'ind_st', ind_st 
      print*;  read*
    endif 
    
    ! TODO: if there any exceptions? extreme case to handle? 
  end do
  
  ! update the eigenvalues by the sorted index 
  wr_st = wr(ind_st)                            
  wi_st = wi(ind_st)
  
  if(debug == 1) then  ! to check 
    print*; print*, 'after sort:'
    print*, 'ind_temp =', ind_temp, 'ind_st', ind_st 
    print*;  read*

    print*, 'Original wr, wi'
    print*, wr
    print*, wi
    print*; read* 
   
    print*, 'sorted wr, wi'
    print*, wr_st 
    print*, wi_st 
   endif 
   
!  print*,  'ind_st', ind_st   
  return  
end subroutine eig_sort













