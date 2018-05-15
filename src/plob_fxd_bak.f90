subroutine plob_fxd(y0, t0, h, nt, ftag, deriv, yall) 
!********************************************************************
!  this subroutine is to plot a segment of trajectory, the intial state y0(6), 
!  time span for integration is [0, tf], with fixed stepsize as h, in this case we can't 
!  obtain the end time to be exact tf

! stop conditon: end time tf? or steps n? 


!  And save all the date in file ftag, with 2 blank lines between 2 different orbits
   
!  better to do the write loop outside the subroutine--- 
!  in case to change the name of the data file, and write during the do loop

!  Input variables  
!    y0(6) :  inital state of trajectory 
!    t0    :  start time for the integration
!    tf    :  end time for the integration  
!    ftag  :  data file to save the result of integration 
   
   
! Output Variables
!    yall(nt,6)  :  all the state obtain from the integration 
  
!  ROUTINE USED: GR_RK78     GR_RTBP(A,B,N,F)
!********************************************************* ***********

implicit none 

integer, parameter:: dp = kind(1.d0) 

real(kind=dp) :: r(13,6),b(6),f(6), hmin, hmax, e1, t, y(6)  ! rk78
		 

integer, intent(in):: nt, ftag   
real(kind=dp), intent(in):: y0(6),t0, h
real(kind=dp), intent(out):: yall(nt,6) ! As the final state
 
integer :: ni  
   
external deriv

! use hmin=hmax=h, to tell the integrator is fixed-stepsize
hmin = h  
hmax = h  
e1   = 1.d-13
 
y = y0      
 
! Integrate the periodic orbits by gr_rk78, write the x, y in txt file with ftag
t =  t0 !0.d0

ni = 1
print*, 'ni, nt', ni,nt
read(*,*)  

!  call gr_rk78(t,y,6,h,hmin,hmax,e1,r,b,f, deriv)
  
do while( ni .lt. nt) 
  write(ftag,'(7e20.10)') t, y
  write(*,*) ni, t, y
!  read(*,*)
  call gr_rk78(t,y,6,h,hmin,hmax,e1,r,b,f, deriv)
  ni = ni+1
enddo


!do while( ni < nt ) 
!  yall(ni,:) = y
!  write(ftag,'(7e20.10)') t, y
!  write(*,'(1I5, 8e20.10)') ni, t, h, y
!  
!  print*, 'ck rk78, t, y,h, hmin, hmax'
!  print*, t, y,h, hmin, hmax 
!  read(*,*)
!  
!  call gr_rk78(t,y,6,h,hmin,hmax,e1,r,b,f, deriv)
!  ni = ni+1
!  
!enddo

! save the last point
yall(ni,:) = y
write(ftag,'(7e20.10)') t, y     


!write(*,*) 'after plob_fxd, t0, t, tf, yf', t0, t, tf, y !ck
read(*,*)  !ck

write(ftag,*)  ! one blank line use as block seperator
!write(ftag,*) 

return
end subroutine plob_fxd


