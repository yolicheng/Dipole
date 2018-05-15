subroutine swap2gr(x0, cs, xg)

! ***********  swap the state to general form ************

! rule of     x_i, x_i+1, x_i+2
! x -> y  -> z  -> x -> y   
! case 1: N =  [1 0 0]  x  -> y  -> z 
! case 2: N =  [0 1 0]  y  -> z  -> x 
! case 3: N =  [0 0 1]  z  -> x  -> y 
! 

!    Input
!  x0	the state of point of which to compute
!  sgn	the sign of q/m 

! 	output 
!  xg   the general form of the state


implicit none
integer, parameter :: dp = kind(1.d0) , n = 6

integer, intent(in) ::   cs
real(kind=dp), intent(in) ::  x0(n)
real(kind=dp), intent(out) :: xg(n)


if(cs == 1) then 
! copy the initial one
  xg = x0 
!  write(*,'(6f8.4)') xg
  
elseif(cs == 2) then 
!(C1, C2, C3, C4, C5, C6) → (C2, C3, C1, C5, C6, C4).
!  print*, '(C1, C2, C3, C4, C5, C6) → (C2, C3, C1, C5, C6, C4)'
  xg = x0((/2,3,1,5,6,4/) )

elseif( cs== 3) then 
! (C1, C2, C3, C4, C5, C6) → (C3, C1, C2, C6, C4, C5).
!  print*, '(C1, C2, C3, C4, C5, C6) → (C3, C1, C2, C6, C4, C5)'
  xg = x0((/3,1,2,6,4,5/) )

endif


end subroutine swap2gr
