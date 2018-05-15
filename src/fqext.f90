subroutine fqext(feqm, nr, lt, isft, f_fft, nfq,  fqmx)
!  this subroutine is extract the frequency and amplitude from the data returned by twofft, which is complex numbers
!  and save them to file ftag if isft == 1, and return the first nfq dominant frequencies, to be written to file later

! the structure of fqmx 
!  1st row:             0, a0, a0, 0 
!  the other rows:      f, amp, real, imag

! 	Input
! nr 	:	number of  rows in feqm
! feqm	:	 feq - mod - real - imag 
! lt	: 	length of signal, used to compute frequency
! isft  :       flag to indict to save all the results from dtwofft or not ? 
! f_fft :	file tag to save the frequency and amplitude for fft1.dat
! nfq	:	number of  frequencies to keep, we only consider the number less than 10(included) 
!                       normally we take nfq = 4

!      Onput
! fqmx(nfq+1, 4):  the selected frequency(the first column), and amplitude(the second column) - real -imag


implicit none 
integer, parameter:: dp = kind(1.d0) 
real(kind=dp), parameter :: pi = 4.d0*datan(1.d0)

integer, intent(in) ::  nr, isft, f_fft, nfq    

real(kind=dp), intent(in)   ::  lt ,   feqm(nr, 4)
real(kind=dp), intent(out)  ::  fqmx(nfq+1, 4) 

! local variables
integer :: i, j
real(kind=dp) :: f, amp, iamp, ramp, x_curr(4), fqmxm1(nfq, 4) 


 CHARACTER(LEN=*), PARAMETER  :: fmt  = "(4f22.12)" ! the format for magnitude output


! the initial value should be zero
fqmxm1 = 0.d0  

fqmx(1, :) =  feqm(1,:)

do i = 1, nr

  if(isft == 1) write(f_fft, fmt)  feqm(1, :)
  
! the first row is the constant term  
  if( i > 1) then 
    x_curr = feqm(i, :)
    call fqmax(x_curr, nfq, fqmxm1)
  endif 
  
enddo 

! the second and latter rows are the frequencies with largest amplitude  
fqmx(2:nfq+1, :) = fqmxm1

print*, nfq, 'dominant frequencies';  print*; 
do j = 1, nfq+1 
  write(*,'(4f20.10)') fqmx(j, :)
enddo 
print*     
  
! add a blank line to seperate orbit as block
if(isft == 1)  write(f_fft, *) 
 
return
end    
    
    
