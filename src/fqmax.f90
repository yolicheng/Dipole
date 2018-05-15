subroutine fqmax(xin, nfq, xout)
!  this subroutine is to compare the current frequency and magnitude with the maximum array
!  and return the updated maximum array, which always have the first nfq maximum amplitude

!  but the constant term, which equals to fq = 0, doesn't have any sense, so we will discard this term 

! 	Input Parameter
!  xin 		:	the current frequency + amp + real + imag
!  nfq		:	number of  frequencies to keep, we only consider the number less than 6(included) 
!                       actually we take nfq = 4

! Input & Output 
!  xout(nfq,4)	: 	the current and updated frequency(the first column), and amplitude(the second column) 

implicit none 
 
integer, parameter:: dp = kind(1.d0) 

integer, intent(in)::  nfq  
real(kind=dp), intent(in):: xin(4) ! freq + amp + real + imag
real(kind=dp), intent(inout):: xout(nfq, 4) ! As the final state
 
! local parameter 
real(kind=dp) :: x(nfq, 4),  amp  
integer :: j
!print*, 'nfq', nfq
 
x = xout

amp = xin(2) ! the current amplitude 


if ( amp > x(1,2) ) then  
 ! replace the maximum one, the first array 
  xout(1,:)    = xin
  if (nfq == 1) return 
 
  xout(2:nfq,:) = x(1:nfq-1,:)

else
  if (nfq == 1) return 
  if ( amp > x(2,2) ) then  
 ! replace the second maximum one, the second array 
  xout(2,:)    = xin 
  if (nfq == 2) return   
  
  xout(3:nfq,:) = x(2:nfq-1,:)
  
  else
    if (nfq == 2) return
    if ( amp > x(3,2) ) then  
    ! replace the third maximum one, the third array 
    xout(3,:)    =  xin 
    if (nfq == 3) return   
  
    xout(4:nfq,:) = x(3:nfq-1,:)


    else
      if (nfq == 3) return 
      if ( amp > x(4,2) ) then  
 ! replace the 4-th maximum one, the 4-th array 
      xout(4,:)    =  xin  
      if (nfq == 4) return   
  
      xout(5:nfq,:) = x(4:nfq-1,:)

else
  if (nfq == 4) return
  if ( amp > x(5,2) ) then  
 ! replace the 5-th maximum one, the 5-th array 
  xout(5,:)    =  xin 
  
  if (nfq == 5) return   
  
  xout(6:nfq,:) = x(5:nfq-1,:)
    
else
  if (nfq == 5) return
  if ( amp > x(6,2) ) then  
 ! replace the 6-th maximum one, the 6-th array 
  xout(6,:)    =  xin  
  if (nfq == 6) return   
  
  xout(7:nfq,:) = x(6:nfq-1,:)
  
else
  if (nfq == 6) return
  if ( amp > x(7,2) ) then  
 ! replace the 6-th maximum one, the 6-th array 
  xout(7,:)    =  xin 
  if (nfq == 7) return   
  
  xout(8:nfq,:) = x(7:nfq-1,:)

else
  if (nfq == 7) return 
  if ( amp > x(8,2) ) then  
 ! replace the 6-th maximum one, the 6-th array 
  xout(8,:)    =  xin 
  if (nfq == 8) return   
  
  xout(9:nfq,:) = x(8:nfq-1,:)

else
  if (nfq == 8) return
  if ( amp > x(9,2) ) then  
 ! replace the 6-th maximum one, the 6-th array 
  xout(9,:)    =  xin  
  if (nfq == 9) return   
  
  xout(10:nfq,:) = x(9:nfq-1,:)

else
  if (nfq == 9) return  
  if ( amp > x(10,2) ) then  
 ! replace the 6-th maximum one, the 6-th array 
  xout(10,:)    =  xin 
  if (nfq == 10) return   
  
  xout(11:nfq,:) = x(10:nfq-1,:)
               
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif




return
end subroutine fqmax

  
   
