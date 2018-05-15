program fft_prtbp

!Do a Fourier analysis of the refined torus computed, to check if the fourier analysis is reliable
!>> rho=1.1140282949439548;         t2 = 20.929626993648455
!>> w2 = 2*pi/t2;  w1 = rho / t2

!--- Original torus
! Unit: radian/unit time
!w1 =  0.0532273363152640
!w2 =  0.300205317041071

! Unit : cycle / unit time, to keep consistent with fftsc 
!w1 = 0.00847139368218902
!w2 = 0.0477791601495561

! --- by fft_ob(furian), we have the first two freq with the maximal amplitudes
!  unit:  cycle/unit time
! w1 = 4.7779160148975076E-002
! w2 = 0.15180887412977237

! but after we do some trick that sets 8.4713933651788618E-003 as basic frequecy
! ** NOTE ** the basic frequecies are not uniquely defined, they do not have to have some physical meanings
!            very nice, the precision is quite high!  We can trust fft_ob!!!
! w1 =   8.4713933651788618E-003
! w2 =   4.7779160148919413E-002



! This is a point on the torus
!  -0.48419417582513136       0.87534078610935684        9.7355202644228162E-003   4.7682020962989688E-003   3.0000047683715820

use dp_mod
use pi_mod
use emconst_mod 
use rtbpconst_mod
implicit none


integer, parameter  :: ndim0   = 4  ! for PRTBP
                         


! Local Variables
 
! -- for fft_ob 
integer       :: nfmax, nsm, nitmax 
real(kind=dp) :: pvi(ndim0),  xmax, tf, cj0  

external :: gr_cjprtbp, gr_prtbp ! Jacobi constant + vector field 


! -- assign the value of mass ratio for EM system and put the values into the common module --- 
call init_emconst
print*, 'emrat = ', emrat 
call init_rtbpconst(emrat, ndim0)

print*, 'check the mass ration, mu =', mu, 'ndim=', ndim 
read*


!pvi = (/xl4 + 1.d-3, yl4 + 2.d-3, 1.d-4, -1.d-4/) ! add a small deviation on L4 to obtain an initial point
! cj0 = 3.0000046704015602  

!  -0.48419417582513136       0.87534078610935684        9.7355202644228162E-003   4.7682020962989688E-003   3.0000047683715820

pvi = (/-0.48419417582513136d0,  0.87534078610935684d0, &
       9.7355202644228162d-3,   4.7682020962989688d-3/)
 cj0 =   3.0000047683715820

tf     = 1.d4 
nitmax = 20  
nsm    = 2**16
nfmax  = 20 
xmax   = 4.d0


open(60,file = './dat/tori_prtbp/fft_ob.dat', access ='append', status = 'replace' )
write(60 ,*)  ' # Refined torus to do fft : t (x y vx vy)   cj'

open(61  ,file = './dat/tori_prtbp/fft_fcs.dat', access ='append',status = 'replace')
open(62  ,file = './dat/tori_prtbp/fft_fcsbas.dat', access ='append',status = 'replace')

! -- do fourier analysis on this orbit ---
print*, 'DO fourier analysis!'; read*
!subroutine fft_ob(x0, n, lt, np, nitmax, nfmax, xmax, fob, ffc, ffbas, deriv, gr_cj) 
call fft_ob(pvi, ndim, tf, nsm,  nitmax, nfmax, xmax, 60, 61, 62, gr_prtbp, gr_cjprtbp)  

read* 


stop

end 


 



