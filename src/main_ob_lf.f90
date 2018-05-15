program main_ob_lf
! 2018-04-22 17:46:23 
!# the integration of curve for torus.
! given the initial state of one point, do the integration. 

!       
use dp_mod
use lf_mod ! the system related parameters of lorentz force problem
use pi_mod

implicit none

integer, parameter ::  ndim = 6       ! the dimension of LF problem
!integer, parameter ::  np   = 5       ! 4 points on the curve    ---curve15_4pt.dat
!integer, parameter ::  np   = 16       ! 4 points on the curve    ---curve15_4pt.dat

integer, parameter ::  np  = 65      ! the points on the curve --- curve15.dat 
!integer, parameter ::  nepoch  = 256     ! the points on the curve 

integer, parameter ::  nepoch  =64     ! the points on the curve 

! Global  Variables, the assignment in the main routine can use an appendix 0 to avoid confliction. 
!lf_mod:      beta, cs, sgn, eq    
   

real(kind=dp) ::  beta0, tp2 
 

! Local Variables
 
real(kind=dp) :: y0(np*ndim), yf(np*ndim), pt(ndim), xi, cj, &
                 t0, tf, pt0(ndim), ptf(ndim), tdir
 

! 4 initial points from invariant curve and the integration as torus  
integer :: fcurve, ftori, fcurve_map, fcurve_map_nospace,ipt, i, j,  ispl, debug     
 character(len=70) :: fncurve, fntori, fncurve_map,fncurve_map_nospace

 debug = 0

! Intialize the state vector of the equilibirum points 
call init_lf  

print*, 'Please input beta (q/m):'
read(*,*) beta0 
call init_beta(beta0)   
    
print*, 'check, cs, ieq, beta', cs, ieq, beta 
read* 
    
!call gr_cjlf(eq, cj)
!print*,'check energy!, cj, ieq,', cj, ieq, eq 
!read*

fncurve     =  './dat/curve_lf/time/torus21_end/curve15.dat'
!fncurve     =  './dat/curve_lf/time/torus21_end/curve15_4pt.dat'
fntori      =  './dat/curve_lf/time/torus21_end/torus15.dat' ! remember to modify the name 
fncurve_map =  './dat/curve_lf/time/torus21_end/curve15_map.dat'  
fncurve_map_nospace =  './dat/curve_lf/time/torus21_end/curve15_map_nospace.dat' 

!fncurve     =  './dat/curve_lf/time/torus21_end/curve15.dat'
!fntori      =  './dat/curve_lf/time/torus21_end/torus15_all.dat' ! remember to modify the name 
!fncurve_map =  './dat/curve_lf/time/torus21_end/curve15_map_all.dat' 
!fncurve_map_nospace =  './dat/curve_lf/time/torus21_end/curve15_map_all_nospace.dat' 

fcurve = 20; ftori = 21; fcurve_map = 22; fcurve_map_nospace=23

open(fcurve, file = fncurve,  status='old')  !read from old file 
read(fcurve,*) ! one comment line 

open(ftori, file = fntori, access  ='append', status='replace') 
write(ftori, *) '# Tori: t (x,y,z,vx,vy,vz) cj  ', np, 'points at one epoch' 

open(fcurve_map, file = fncurve_map, access  ='append', status='replace') 
write(fcurve_map, *) '#: t (x,y,z,vx,vy,vz) cj -- The Time-T map of the curves within T2 ', np, 'points at ', nepoch, 'epoches'  

open(fcurve_map_nospace, file = fncurve_map_nospace, access  ='append', status='replace') 
write(fcurve_map_nospace, *) '#:The Time-T map of the curves for data processing', np, 'points at ', nepoch, 'epoches'  


ipt = 0
y0 = 0.d0

open(333, file = './dat/ob_guess.dat', access ='append',status='replace')

 
!subroutine plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 
  
tdir = 1.d0; ispl = 1  
t0 = 0.d0;  tf = 0.d0; 
tp2 = 2.15052849370480d0 ! T2 = 2.15052849370480 from para_curve31.dat 
 
!plot orbit for 4 initial points
! plob is okay! but plob_n seems wrong.


!subroutine plob_n(y0, ndim, np, t0,tf, tdir, ispl, fob, yf, deriv, gr_cj) 


! save the initial points 
do j = 1, np 
    
     
   read(fcurve, *) xi, pt, cj ! read from
   
   
   write(fcurve_map,*) tf, pt, cj
   if(j<np) write(fcurve_map_nospace,*) tf, pt, cj
   
   do i = 1, 64/np-1 
     read(fcurve, *)
   enddo 
   
   if(debug == 1) print*, tf, pt,  cj
   call gr_cjlf(pt, cj) 
   
   
   ! check the orbit of the first point 
!   subroutine plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 

!   call plob(pt, t0, tp2, ndim, ndim, tdir, ispl, 333, gr_lf, gr_cjlf,  ptf)
!   print*, 'Finish the orbit of first point on curve, pt0, ptf'
!   print*, pt,cj
!   call gr_cjlf(ptf, cj) 
!   print*, ptf, cj
!   write(333,*);   write(333,*)
!----------------------------------------------

   
   y0( (j-1)*ndim+1 : j*ndim)  = pt  
   if(debug == 1) print*, tf, pt,  cj
enddo 


if(debug == 1) print*; 
 
!   close(333)
!   stop
  
write(fcurve_map,*); write(fcurve_map,*)
print*  
  
  
do i = 1, nepoch

  tf = tp2/nepoch*i
  
  ! plot the four points one-by-one 
  ! save the points 
  do j = 1, np 
    
    pt = y0( (j-1)*ndim+1 : j*ndim) 
    call plob(pt, t0, tf, ndim, ndim, tdir, ispl, ftori, gr_lf, gr_cjlf,  ptf)
    
    yf( (j-1)*ndim+1 : j*ndim)  = ptf 
    
    call gr_cjlf(ptf, cj)
    write(fcurve_map,*) tf, ptf,  cj 
   
    if(j<np) write(fcurve_map_nospace,*) tf, pt, cj
    if(debug == 1) print*, tf, ptf,  cj
    
  enddo  
   
  y0 = yf;  t0 = tf
  
  write(fcurve_map,*);   write(fcurve_map,*)
  if(debug == 1) print*
enddo 


 close(fcurve_map);  close(fcurve_map_nospace) 


 
end program main_ob_lf



















  
