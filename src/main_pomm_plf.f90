program pomm_plf

! Given the I.C. of the periodic orbits, compute the MM and eigenvalues for stability analysis
!  --- data saved in ./dat/map_plf/h/  subfolder 

! Using the p.o. of Type III,   ./po3_initall.dat     ./po3_oball.dat 
 
  
use dp_mod
use pi_mod 
use poinc_mod
use plf_mod 
use po_plf_mod
implicit none

integer, parameter  ::  ndim   = 4,  &  ! for lf
                        npvar  = 20, &  ! 20 for variational matrix, 4 for state vector 
                        npp    = 100  
! Local Variables
integer       ::  debug, i, k
real(kind=dp) ::  beta0, cj_eq

! eigenvalues and eigenvectors around the equilibria                 
real(kind=dp) ::  dlf(ndim, ndim), wr(ndim), wi(ndim), vr(ndim, ndim), phi(ndim, ndim)


! Po 
integer       ::  ispl,  ipo, fob, finit, feig, npo
real(kind=dp) ::  cj0, y0(ndim), TP, yf(ndim), dg(1,1), cj, kern_dg , x0,  &
                  pvi(npvar), pvf(npvar), tf
                  
character(len=70) :: fnob, fninit, fneig


! take eq1  in planar normal case as an example, type center X saddle  equilibria
! with beta = 1.d0, there is resonance when beta is around 1, avoid this for the moment 
debug = 0

    
! to assign the value of  eq     
call init_plf     
 
beta0 = 2.d0   ! use 2 to check 
call init_beta(beta0)   
     

call gr_cjplf(eq, cj_eq)
print*, 'check energy! ', cj_eq, eq

! -- the energy level we r taking, smaller than the one at the equilibria, so we have a 
! closed UFO 

 cj0 = cj_eq + 0.05d0
 print*, 'cj0 = ', cj0, '= 4.3767487109222252 ' ; read*


fob = 100; finit = 200; feig = 88
! better to use the general name, and copy to the specified folder later
! npo = 11
!write (fninit, fmt='(a)')  './dat/map_plf/h3/po3_initall.dat' ! I.C. for Type III p.o.s 
!write (fninit, fmt='(a)')  './dat/map_plf/po_init1.dat' ! I.C. for Type III p.o.s 

! 632 orbits....
npo = 632
write (fninit, fmt='(a)')  './dat/map_plf/po/pox_init_all.dat' ! I.C. for Type I p.o.s 

open(finit, file = fninit, status='old')
read(finit, *)  ! blank comment line 


write(fnob,    fmt='(a)')  './dat/map_plf/po/pox_ob.dat'   ! general name for   p.o. orbits
write (fneig,  fmt='(a)')  './dat/map_plf/po/pox_eig.dat'  ! eigenvalues values

open(fob,   file = fnob,   status='replace', access='append')
open(feig,  file = fneig,  status='replace', access='append')

! read the I.C. of the p.o., and Integrate for one period to compute the
! Monodromy Matrix and the eigenvalues
do ipo = 1, npo, 1
  read(finit, *) TP, y0, cj
  write(*,*) TP, y0, cj
!  print*; read*

  
  ! --- plot the orbit for one revolution and compute the Eigenvalues of the Monodromy Matrix to check the stability
  pvi(1:4)    = y0
  pvi(5:20)   = 0.d0
  pvi(5:20:5) = 1.d0
  
!  subroutine plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 
  call plob(pvi, 0.d0, TP, ndim, npvar, 1, 1, fob, gr_plf, gr_cjplf,  pvf) 
    
  phi = reshape(pvf(5:20), (/4,4/))
   
!  print*; print*, 'MM'  
!  do i = 1, ndim, 1
!    print*, phi(i,:)
!    end do
!  print*; read*
  
    
  call eigrg(phi, ndim,1, wr, wi, vr)
 
!  print*, 'Finish one p.o.! Check the stability!'; read* 
    
  !  Save the eigenvalues to files, avoid 444 and 6(by default the screen)
  ! Save the eigenvalues to files + x0 as the first parameter  
  write(feig, '(2e24.14)', advance='no')   pvi(1), cj  ! x
  do i = 1, ndim, 1
    write(feig, '(2e24.14)', advance='no')  wr( i ), wi( i ) 
  end do  
  write(feig,*)
  
!  call prt_eigval( ndim, feig, wr, wi )
    
enddo 

 close(feig)
stop 
 
end 
 


