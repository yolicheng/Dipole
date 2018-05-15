program pomfd_plf

! This is to compute the periodic orbits, as well the associated manifold if it is unstable,
! in the planar LF problem of the normal case 

! I just take 100 initial points equally-spaced on the  manifold, integrate and check if the stable and unstable always overlap??

! TODO:
! --1. refine the periodic orbit
! --2. compute the manifold

!  --- data saved in ./dat/pomfd_plf/  subfolder 

!  1. po.dat         -- the refined periodic orbits
!  2. poit.dat       -- the intermediate orbits during refinement 
! 
use dp_mod
use pi_mod 
use plf_mod 
implicit none

integer, parameter  ::  ndim   = 4,  &  ! for lf
                        npvar  = 20, &  ! PLF: 20 for variational matrix, 4 for state vector 
                        nmf    = 100    ! number of orbits on the manifold
                        
! to compute the manifold, ta   

                      
! Local Variables
integer       ::  debug, i 
real(kind=dp) ::  beta0 

! eigenvalues and eigenvectors around the equilibria                 
real(kind=dp) ::  wr(ndim), wi(ndim), vr(ndim, ndim), phi(ndim, ndim)

! P.O. 
integer       ::  fob, fpo_init, finit, fmfd, col
real(kind=dp) ::  y0(ndim), TP, cj0, yf(ndim) 

! mfd
real(kind=dp) ::  pvi(npvar), pvf(npvar), ymfd(nmf, 4), yi(ndim),  tf, epsl, vep(ndim)
                  
character(len=70) :: fnob, fnpo_init, fninit


! take eq1  in planar normal case as an example, type center X saddle  equilibria
! with beta = 1.d0, there is resonance when beta is around 1, avoid this for the moment 
debug = 0
    
! to assign the value of  eq     
call init_plf     
 
beta0 = 2.d0   ! use 2 to check 
call init_beta(beta0)   
     
!call gr_cjplf(eq, cj)
!print*, 'check energy! ', cj, eq

! ** NOTE **
! here we only deal with h3, if we change the energy level, be careful with the name of the folder 
! no to overwrite the available data 

! Files open and comment
fob = 100; fpo_init = 200;  finit = 300; fmfd = 888
! better to use the general name, and copy to the specified folder later
!write(fnob, fmt='(a)')      './dat/pomfd_plf/po1_ob.dat' !general name for 1 p.o.
!write (fninit,  fmt='(a)')  './dat/pomfd_plf/po1_init.dat'
!write (fnmfd,  fmt='(a)')  './dat/pomfd_plf/mfd_init.dat'

write(fnpo_init,  fmt='(a)')  './dat/map_plf/h3/poy_init.dat'

! use new name  not to overwrite the old files
write(fnob,       fmt='(a)')  './dat/map_plf/h3/poy_ob1.dat' 
write(fninit ,    fmt='(a)')  './dat/map_plf/h/mfd_init1.dat'


open(fob, file = fnob, status='replace')

! read from the refined periodic orbits.
open(fpo_init, file = fnpo_init, status='old')
read(fpo_init, *) tp, y0, cj0 
write(*, *) tp, y0, cj0  ; read* ! to check

open(finit, file = fninit, status='replace', access='append')
write(finit,*) '# I.C. for manifolds:   (x, y, vx, vy)   cj'

! --- compute the manifold, using the eigenvector of Monodromy Matrix as linear approximation
  ! --- plot the orbit for one revolution and compute the Eigenvalues of the Monodromy Matrix to check the stability
pvi(1:ndim) = y0  
pvi(5:20)   = 0.d0
pvi(5:20:5) = 1.d0
! do not plot, only use this for Monodromy matrix computation
call plob(pvi, 0.d0, TP, npvar, 1, fob, 0, gr_plf, gr_cjplf,  pvf) 
 close(fob) 
 
 
! check the final state 
print*, 'I.C. :', pvi(1:4)
print*, 'F.C. :', pvf(1:4) 
print*; read*

! Monodromy Matrix  
phi = reshape(pvf(5:20), (/4,4/))
do i = 1, ndim, 1
   print*, phi(i,:)
end do
print*; read*
    
call eigrg(phi, ndim, 1, wr,wi, vr)
print*, 'Pick the dominant unstable eigenvector ( lambda_max > 1)! Input the column'
read(*, *) col 
vep = vr(:, col)
print*, 'eigenvector :', vep; print*; read*

!print*, 'Check if the p.o. is unstable? ' 
!print*; read* 

! Compute the manifold 
!subroutine gr_mfdinst(ypo, n, tp, mmat, nmf, epsl, ftag, ymfd, deriv, gr_cj)
! the integration time for the manifolds
tf = 3.d3 

! --- right half of unstable manifold ---- 
epsl = 1.d-6
call gr_mfdinst(y0,  4, tp, vep, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 


! -- ck
!open(345, file = './dat/map_plf/h3/ck_po.dat', status='replace', access='append')
!call plob(y0, 0.d0, tp, ndim, 1, 345, 1, gr_plf, gr_cjplf,  yf) 
! close(345)

call ck_mfdinst(y0,  4, tp, vep, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 
print*, 'Check mfd finished!'; read*

! --- plot the manifold 
open(fmfd, file = './dat/map_plf/h/mfd_obup.dat', status='replace', access='append')
write(fmfd,*) '# manifold:  t   (x, y, vx, vy)   cj'

do i = 1, nmf, 1
  yi = ymfd(i, 1:ndim)
  call plob(yi, 0.d0, tf, ndim, 1, fmfd, 1, gr_plf, gr_cjplf,  yf) 
end do
 close(fmfd)

! --- left half of unstable manifold ---- 
!tf = 8.d2
open(fmfd, file = './dat/map_plf/h/mfd_obum.dat', status='replace', access='append')
write(fmfd,*) '# manifold:  t   (x, y, vx, vy)   cj'
epsl = -1.d-6
call gr_mfdinst(y0,  4, tp, vep, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 

do i = 1, nmf, 1
  yi = ymfd(i, 1:ndim)
  call plob(yi, 0.d0, tf, ndim, 1, fmfd, 1, gr_plf, gr_cjplf,  yf) 
end do
 close(fmfd)


! --- left half of  stable manifold ---- 
print*, 'Pick the stable eigenvector ( 1 / lambda_max < 1)! Input the column'
read(*, *)  col 
vep = vr(:, col)
print*, 'eigenvector :', vep; print*; read*

open(fmfd, file = './dat/map_plf/h/mfd_obsp.dat', status='replace', access='append')
write(fmfd,*) '# manifold:  t   (x, y, vx, vy)   cj'
epsl = 1.d-6
call gr_mfdinst(y0,  4, tp, vep, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 

do i = 1, nmf, 1
  yi = ymfd(i, 1:ndim)
  call plob(yi, 0.d0, tf, ndim, -1, fmfd, 1, gr_plf, gr_cjplf,  yf) 
end do
 close(fmfd)
 

open(fmfd, file = './dat/map_plf/h/mfd_obsm.dat', status='replace', access='append')
write(fmfd,*) '# manifold:  t   (x, y, vx, vy)   cj'
epsl = -1.d-6
call gr_mfdinst(y0,  4, tp, vep, nmf, epsl, finit, ymfd, gr_plf, gr_cjplf) 

do i = 1, nmf, 1
  yi = ymfd(i, 1:ndim)
  call plob(yi, 0.d0, tf, ndim, -1, fmfd, 1, gr_plf, gr_cjplf,  yf) 
end do
 close(fmfd)
 
 
end 
 
 
