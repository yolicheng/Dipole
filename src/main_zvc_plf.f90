! Produce the data for zero velocity curve plot of PLF problem  

!  --- data saved in ./dat/zvc_plf/   
!  1. zvc_plf_3d.dat       -- x-y-h 

!  -- to process data by gnuplot, finally we get table 

!   use script zvc_plf.pl in ./dat/ directory 
! --2.1  produce zvc_plf_3d.tab 
! --2.2  fix energy level, plot  
! --2.3  plot directly with zvc_plf_2d.tab
!         p 'zvc_plf_2d.tab' w l 
         
! 
use dp_mod
use plf_mod 
implicit none

integer, parameter  ::  ndim   = 4  ! for lf
! for zvc 
real(kind=dp) ::  xlim(3), ylim(3)


call init_plf     

! ----- Produce the data for zero velocity curve to plot 
xlim = (/-2d0, 5.d-3, 2d0/)
ylim = (/-2d0, 5.d-3, 2d0/)
call zvc_plf(xlim, ylim)

stop 

end
