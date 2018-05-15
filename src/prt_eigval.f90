!***********************************************************************
!     ****   prt_eigval   ****
! Print the eigenvalue by (real +  imaginary) to file, or on the screen(ftag == 6)

!  ****** Input Variables ******
!  n       dimension of the  phase space 
! ftag     file tag 
! wr       real part of the eigenvalues
! wi       imaginary part of the eigenvalues
!
!  Routine Used:
!    None  

!  Finally Revised by Yu -- 20161014
!***********************************************************************
subroutine prt_eigval( n, ftag, wr, wi)

use dp_mod
implicit none
 
! Input  and Output Declaration   
integer, intent(in)         ::  n, ftag 
real(kind=dp), intent(in)   ::  wr(n), wi(n)
 
! Local Variable
integer ::  i
  
!  write(ftag, *)  wr , wi 
!  if(ftag /= 6) close(ftag) 

  do i = 1, n, 1
    write(ftag, '(2e20.10)', advance='no')  wr( i ), wi( i ) 
  end do  
  
  write(ftag, *)
  
!  print*, 'check the file to save eigenvalues'; print*; read*


  return  
end subroutine prt_eigval

