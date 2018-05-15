program  main_rotnum
!  Compute the rotation number of the invariant curve saved in ./dat/curve_poinc/ob.dat 
!  Use the method from Carles' class note P63

implicit none
integer, parameter  :: dp = kind(1.d0)
real(kind=dp), parameter :: pi2 = 8.d0*datan(1.d0)
integer, parameter  :: npoinc = 1000, &   ! number of Poincare maps, to start test, 11 is enough?
                       np_interpol = 1001, &  ! in order to use fourier, we need odd number of points 
                       nf   =  100 ! how many Fourier modes ???
                       
! Local Variables
real(kind=dp) :: tpc(npoinc), xpc(npoinc, 6), cj, t, pv(6), & ! curve
                 pt(npoinc, 2), rho, arg(npoinc), & ! rotnum
                 ptnew(np_interpol, 6), pv_in(np_interpol), dt, fun6(6)  ! interpol + fourier
                 
integer ::  i,  k, l,  ncrs, ind_sorted(npoinc),  & 
            fob,  fob_interpol, ffcs, ffun  ! file tag  
            
character(len=70) ::  fnob, fnob_interpol, fnfcs, fnfun 

real(kind=dp), dimension(0:nf, 6) ::  CSF6, SIF6 ! fourier coefs of  all six compoments
real(kind=dp), dimension(0:nf)    ::  CSF, SIF   ! fourier coefs of 1 component
 
! claim the self-defined function 
real(kind=dp) ::  funcs 
            
! --- deal with the file open to save or read ---- 
fob   = 25;   fnob  =  './dat/curve_poinc/ob.dat'
open(fob  ,file = fnob)
read(fob, *) ! skip  the first comment line 


! -- interpol --
fob_interpol = 26;   fnob_interpol  =  './dat/curve_poinc/ob_interpol.dat'
open(fob_interpol  ,file = fnob_interpol)
write(fob_interpol, *)  ' # argument      (x,y,z,vx,vy,vz)'  ! the comment line 

! -- Fourier analysis --
! the coefficients computed by general Fourier analysis:  fourier.f + fun.f 
ffcs = 111;  fnfcs = './dat/curve_poinc/fcs.dat'
open(ffcs, file = fnfcs) 
write(ffcs,*) '# For each column: (x,y,z,vx,vy,vz),  for each line: the Four Coef ck(i), sk(i), i=0,...,', np_interpol 

! the approximated Fourier function, to check if they match the original data 
ffun = 222;  fnfun = './dat/curve_poinc/fun.dat'
open(ffun,  file = fnfun) 
write(ffun, *) ' # argument    (x,y,z, vx,vy,vz)' ! TODO -- check the energy? cj'

print*, 'Files names for data read and write'
print*, fnob, fnob_interpol, fnfcs, fnfun;  print*   !ck


! read the Poincare maps from the data file 
do i = 1, npoinc, 1  
!  write(fpc, '(8e24.14, 1I5)')  tpci, xpc, cj, ncrs !-- from poinc_n.f90
  read(fob, '(8e24.14, 1i5)') t, pv, cj, ncrs
  tpc(i) = t
  xpc(i, :) = pv
enddo   
  
  
! -- compute the rotation number -- 
pt = xpc(1:npoinc,  1:2) ! only use x-y

!subroutine rotnum( pt, np, rho)
call  rotnum( pt, npoinc, rho, arg, ind_sorted)

! update the points by an increasing order of arg 
xpc = xpc(ind_sorted+1, :) 

! ----- linear interpolation to get equally spaced points in angle -----------
! do the linear interpolation by the value of arg, to obtain npoinc points equally spaced in arg within the interval [0, 2pi]

call  interpol( xpc, npoinc, 6, arg, np_interpol, fob_interpol, ptnew)
 
! -- save all the points to file ob_interpol.dat, better to do this in the routine interpol  


! -- general Fourier analysis for nf Fourier modes, with np_interpol points  

! deal with the six components one by one 
do k = 1, 6, 1
  pv_in =  ptnew(:, k) ! one component, the k-th column  
 
  call gr_foun( pv_in, np_interpol, nf, csf, sif)
  
  csf6(:, k) = csf
  sif6(:, k) = sif
  
end do
  
  
! write the coefficients C's and S's into file fcs.dat (x,y,z,vx,vy,vz)  cf // sf !x  // cfy-sfy -- ......
do k = 0, nf, 1
  
   write(ffcs, *)  csf6(k, :) 
   write(ffcs, *)  sif6(k, :) 
   write(ffcs,*)
!     sif6(1,k), csf6(2,k), sif6(2,k), csf6(3,k), sif6(3,k), &
!                     csf6(4, k), sif6(4,k), csf6(5,k), sif6(5,k), csf6(6,k), sif6(6,k)
enddo


! Period = 2*pi/rho, the relation between the period and the rotation number 
dt = pi2/rho/np_interpol ! take this value to keep coherent with the original points  

! Evaluate the truncated Fourier series with the above coefficients
do l = 1, np_interpol*2 
 
    t = (l-1) * dt
    
    ! the 6 components
    do k = 1, 6, 1
      fun6(k) = funcs( t, rho, nf, csf6(:, k), sif6(:, k) )
    end do
    
!    call gr_cjlf(fun6, cj)
    write(ffun, *) t, fun6  !, cj 
    
  enddo 
  
  write(ffcs, *)
  write(ffun, *) 

end 


 



