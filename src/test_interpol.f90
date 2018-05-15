program test_interpol

!  Study the planar Restricted Three Body Problem, and choose L4, which is of type 
!  center X center, so around L4 we have a plenty 2D tori. 
!    1. start from a point that is close to L4, and integrate the orbit 
!    2. compute the Poincare maps, npoinc = 10000, with Poincare section taken as x = p0 (or y= p0)
!    3. compute the rotation number of numerically computed invariant curve 
!    4. use the linear interpolation to obtain equally spaced points w.r.t. the angle between the relative position vector from L4 and the x-axis!
!    5. do a simple Fourier analysis to compute the Fourier coefficients as the initial guess to refine the curve 
!    6. Impose the invariance equations, and specify the energy level, to compute the invariant curve 
!    7. Globalize the curve to an invariant torus 


!  --- data saved in ./dat/curve_prtbp/  subfolder 
!  1. torus.dat         -- the original orbit of the 2D torus 
!  2. curve.dat         -- the invariant curve obtained by Poincare map 
!  3. ob_interpol. dat  -- linear interpolated points for Fourier analysis
!  4. fun.dat           -- approximated orbits by the Fourier representation 
!  
use dp_mod
use emconst_mod 
use rtbpconst_mod
use sort_mod

implicit none
integer, parameter  ::  npoinc =  100 , & !20 !1000  ! number of Poincare maps, to start test, 11 is enough?
                        ndim0   = 1, & !   !!!! for PRTBP
                        npvar   = 4, & ! 20 for variational matrix, 4 for state vector 
                        np_interpol = 201, &  ! in order to use fourier, we need odd number of points 
                        nf   =  10 ! how many Fourier modes ???
                         
real(kind=dp)       ::  pi2 = 8.d0*datan(1.d0),  day = 24.d0 * 60.d0 * 60.d0  


! Local Variables
real(kind=dp) :: p0, tmax, xmax, pvi(npvar), tf, pvf(npvar), tic, toc, cj0,   &  ! torus
                 xl4, yl4, & ! L4 equilibirum point
                 tpc(npoinc), xpc(npoinc, ndim0), cj, t, pv(ndim0), & !curve 
                 pt(npoinc, 2),  rho, phi(npoinc), & ! rotnum
                 ptnew(np_interpol, ndim0), pv_in(np_interpol),  darg,   funn(ndim0), &  ! interpol + fourier
                 arg1, ptnew1(ndim0)   ! -- test f(x) = 1 + 0.25sin(x) + 0.5cos(x)
  
real(kind=dp), dimension(0:nf, ndim0) ::  csfn, sifn ! fourier coefs of  all the compoments
real(kind=dp), dimension(0:nf)        ::  CSF, SIF   ! fourier coefs of 1 component
 
                  
integer ::  i, j, k, l,  ncrs, ind_sorted(npoinc), &
            ftorus, ispl, fcurve, fcurve_interpol, ffc, ffun, ffcbas, & 
            ind, dir, imax, pt_col, isvy, np, dircurve, & ! poinc 
            iscurve, isinterpol, isfourier !  control which subroutine to call
            
            
character(len=70) :: fntorus, fncurve, fncurve_interpol, fnfc, fnfun,  fnfcbas   


external :: gr_cjprtbp, gr_prtbp ! Jacobi constant + vector field 

real(kind=dp) ::  four_seri  !funcs

iscurve    = 0
isinterpol = 1
isfourier  = 1

! ** NOTE ** we have to take a look at how the invariant curve moves, and set the value of dircurve
!            1: counterclockwise;  -1: clockwise
! we have -1 in this case... 
dircurve =  -1

! -- assign the value of mass ratio for EM system and put the values into the common module --- 
call init_emconst
print*, 'emrat = ', emrat 
call init_rtbpconst(emrat, ndim0)

print*, 'check the mass ration, mu =', mu, 'ndim=', ndim 
read*

! take a point that is close the L4 point,  with the big primary at (mu, 0, 0), 
! the triangular libration point L4 form a equilateral triangle with the two primaries. 
! so the location of L4 is ( -1/2+mu, sqrt(3)/2  ) and  L5 located at ( -1/2+mu, -sqrt(3)/2  )
xl4 = -.5d0 + mu 
yl4 = dsqrt(3.d0) / 2.d0

! the initial point --- instead of take random points, we fix the energy and poincare section 
!pvi = (/xl4 + 1.d-3, yl4 + 2.d-3, 1.d-4, -1.d-4/) ! add a small deviation on L4 to obtain an initial point
! cj0 = 3.0000046704015602  

!this is a test, to fix the energy level, and x0 = xl4,  and change the value of y, to have a look of a lot tori. 

! files to save the original torus + invariant curve + interpolated one + approximated one 
!ftorus = 21;  fntorus    =  './dat/curve_prtbp/torus_test.dat'

fcurve   = 25;  fncurve  =  './dat/curve_prtbp/curve_test.dat'
open(fcurve, file = fncurve, status = 'replace')
! test a simple funtion f= 1 + + .5cos(x) + .25sin(x), with period 2pi,
! take 100 sample points, and do interpolation for 200 points, and check the result    
darg = pi2 / npoinc
do i = 1, npoinc, 1
  tpc(i) = (i-1)*darg
  xpc(i, 1) = 1.d0 + .5d0*dcos(tpc(i)) + .25d0*dsin(tpc(i))
  write(fcurve, * ) tpc(i), xpc(i,1)
end do

! -- interpolated curve (equally spaced in angle within [0, 2pi]) --
fcurve_interpol = 26;   fncurve_interpol  =  './dat/curve_prtbp/curve_interpol_test.dat' 
if( isinterpol == 1) then 
  open(fcurve_interpol  ,file = fncurve_interpol)
  write(fcurve_interpol, *)  ' # Interpolated curve: argument    f'  ! the explanation line 
  call  interpol( xpc, npoinc, 1, tpc, np_interpol, fcurve_interpol, ptnew)
  close(fcurve_interpol)
endif
  

! -- Fourier analysis : Fourier coefficients + approximated curve --    
! the coefficients computed by general Fourier analysis:  fourier.f + fun.f 
ffc = 111;  fnfc = './dat/curve_prtbp/fcs.dat'
open(ffc, file = fnfc) 
write(ffc,*) '# For each column: (x,y, vx,vy),  for each line: the Four Coef ck(i), sk(i), i=0,...,', nf 

! the approximated Fourier function, to check if they match the original data 
ffun = 222;  fnfun = './dat/curve_prtbp/fun_test.dat'
open(ffun,  file = fnfun) 
write(ffun, *) ' # argument    (x,y, vx,vy)' ! TODO -- check the energy? cj'

! -- general Fourier analysis for nf Fourier modes, with np_interpol points  
! --  read from fcurve_interpol for the array ptnew --- 
if( isinterpol == 0) then 
  open(fcurve_interpol, file = fncurve_interpol) 
  read(fcurve_interpol,*) ! skip the comment line 
  do i = 1, np_interpol, 1
    read(fcurve_interpol,*) arg1, ptnew1
    ptnew(i,:) = ptnew1
  end do
endif 

! check ptnew 
!call system( 'gnuplot /home/yu/Dropbox/dipole/dat/curve_prtbp/curve_interpol.pl' )

!  --- Fourier analysis:  deal with the six components one by one 
do k = 1, ndim, 1
  
  pv_in =  ptnew(:, k) ! one component, the k-th column  
 
!  write(*,'(10f16.10)')  pv_in ! ck!!
!  read* 
  
  call gr_foun( pv_in, np_interpol, nf, csf, sif)
  
  ! check gr_foun and gr_four(my routine)
!  print*, 'csf and sif by gr_foun'
!  write(*,'(11f16.10)') csf 
!  write(*,'(11f16.10)') sif 
!  read* 
  
  csf = 0.d0;  sif = 0.d0
!  subroutine gr_four(f, n, m, csf, sif)
  call gr_four(pv_in, np_interpol, nf, csf, sif)
!  print*, 'csf and sif by my own routine gr_four'
!  write(*,'(11f16.10)') csf 
!  write(*,'(11f16.10)') sif 
!  read*   
  
  csfn(:, k) = csf
  sifn(:, k) = sif
  
end do
  
  
! write the coefficients C's and S's into file fcs.dat (x,y,z,vx,vy,vz)  cf // sf !x  // cfy-sfy -- ......
do k = 0, nf, 1
   write(ffc, *)  csfn(k, :) 
   write(ffc, *)  sifn(k, :) 
   write(ffc,*)
enddo


! Period = 2*pi in radian w.r.t. the angle of the points 

darg = pi2/np_interpol ! take this value to keep coherent with the original points  

! Evaluate the truncated Fourier series with the above coefficients
do l = 1,  np_interpol, 1  
 
    arg1 = (l-1) * darg
    
    ! the 6 components
    do k = 1, ndim, 1
!      function four_seri( phi, flag, m, cof, sif)  result(fun)
      funn(k) = four_seri( arg1, 1,  nf, csfn(:, k), sifn(:, k) ) 
    end do
    
    write(ffun, *) arg1, funn  !, cj 
    
!    print*, 'Points (arg - f ):  original  -interpolated:'
!    write(*, *)  arg1,  ptnew(l,:), funn
!    
!    write(*, *)  arg1, funn 
!    read*
  enddo 
  
  write(ffc, *)
  write(ffun, *) 

stop

end 


 



