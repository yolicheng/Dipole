program main_eq
! 2017-02-28 11:37:10 
! Put the equilibria of all the three cases in one file to plot, separated by block
! finished -already

use dp_mod
use lf_mod  
implicit none

! Local Variables

open(100, file='./dat/eq_all.dat', access='append',status='replace')
call  eq_dat( 100 )

stop
end program main_eq



















  
