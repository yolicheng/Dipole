program main_po_asymt

! 2017-02-21 18:43:56 
! For asymmetric periodic orbit, use the general routine
! Integrate for one period, and return to the 
        
use dp_mod
use lf_mod ! the system related parameters of lorentz force problem
use pi_mod
use gr_po_mod

implicit none

integer, parameter ::  ndim = 6       ! the dimension of LF problem 
integer, parameter ::  neq  = 42,  &  ! compute also the variational matrix, 6+36=42
                       npo  = 10       !  

! Global  Variables, the assignment in the main routine can use an appendix 0 to avoid confliction. 
!lf_mod:      beta, cs, sgn, eq    
!poinc_mod:   p0, tmax, xmax, h0 
!             n,  ind_vel, ind, dir, imax  
  
  
!  fix one component: x(ind) = p0

real(kind=dp) ::  beta0, p0  
real(kind=dp) ::  tol0, tol_err0      ! error control 

integer      ::   ind,  dir_curve,  idcor    

! Local Variables
integer :: debug, i,  ivr,  ifam, ispo,  &
           tdir, ipo, isgn, idsgn, niter

real(kind=dp) ::  wr(6), wi(6), vr(6, 6),   &   ! eigenspace of variational matrix  
                  poinst(6), vrchs(6), epsl_po, cj, tp0, pof(6), & !pofam
                  dhmax, dhmin, dpara, cj0,  para1(2),  para0, &   ! continuation
                  dlf(6, 6), wr_mm(6), wi_mm(6), vr_mm(6, 6), &  ! MM 
                  pvi(42), pvf(42), phi(6,6), yf(42), tp,  y0(6)
                   
real(kind=dp), allocatable  :: x0_pre(:, :), x0(:), dg(:,:)  ! for the two previous points
                    
                  
integer :: fpo, fpoinit, fmmegv, fmmegp         
character(len=70) :: fnpo, fnpoinit, fnmmegv, fnmmegp              

real(kind=dp), external :: dnrm2 ! from  library 

! debug or not 
debug = 0

idcor = 1

! already tried the '-' branch 
dir_curve = 1  ! plus family for continuation 


! Intialize the state vector of the equilibirum points 
call init_lf  

print*, 'Please input beta (q/m):'
read(*,*) beta0 
call init_beta(beta0)   
    
print*, 'check, cs, ieq, beta', cs, ieq, beta 
read* 
    
! At the moment --- skip this real unit --- TODO  
! Also compute the unit of distance, time and velocity, which are all function of beta
! but at this moment, we will not use the real quantity of the units. 
! call lfunit(beta, runit, tunit, vunit)
    
print*, 'Do you want to change the sign of equilibirum point (Yes = 1) ?'
read(*,*)  isgn
    
! If the equilibirum points with different signs are not symmetric,  
! the periodic orbits will also be different, we have to compute seperately 
if(isgn == 1) then 
  idsgn     = mod(ind, 3) + 1 
  eq(idsgn) = -eq(idsgn)
endif     
 
! Jacobi matrix of the lorentz force with respective to the state 
call dflrtz(eq, dlf)

! check the energy 
call gr_cjlf(eq, cj)
print*,'check energy!, cj, ieq,', cj, ieq, eq 
read*

! limits of step size for continuation
dhmax = 1.d0  ! for y
dhmin = 1.d-10 

! the variation in energy, we need to decrease to go from the closed region to open one
dpara = dir_curve*1.d-3  

! Compute the eigenvalues and eigenvectors of dlf
call eigrg(dlf, ndim, 1, wr,wi, vr)

print*, 'Which column to use for initial guess of p.o.?'
read(*,*) ivr 

vrchs = vr(:, ivr)
print*, vrchs 
vrchs = vrchs/dnrm2(6,vrchs,1)
print*, dnrm2(6,vrchs, 1), vrchs
read* 

! commomly used  varaibles
epsl_po = 1.d-3     
tol0    = 1.d-10 ;  tol_err0 = 1.d-10
call  init_errctr(debug, tol0, tol_err0)


! the initial guess for po, move along the corresponding eigenvector for a small distance epsl_po 
poinst =  eq  +  epsl_po * vrchs !  
print*, 'poinst=', poinst 
call gr_cjlf(poinst, cj0)

!  ***************** compute the 3 families of PO ***************************
! save the data of po and  its initial state 
fpo = 20;  fpoinit = 21;  fmmegv = 23;  fmmegp = 24 

! use new name, to avoid overwritten the good ones, and the type of approach selected is used as suffix
write(fnpo,    fmt='(a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/apo.dat'  
write(fnpoinit,fmt='(a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/apoinst.dat'
write(fnmmegv, fmt='(a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/apommegv.dat'
write(fnmmegp, fmt='(a,i0,a,i0,a)')   './dat/po/case', cs, '/eq', ieq, '/apommegp.dat' 

print*, fnpo, fnpoinit, fnmmegv, fnmmegp

open(fpo,file=fnpo, access ='append',status='replace')
write(fpo,*)  '# t      (x y z vx vy vz)        cj'

open(fpoinit,file=fnpoinit, access ='append',status='replace')
write(fpoinit,*) '# TP  I.C.(x y z vx vy vz)    CJ    dpara   ds'

open(fmmegv, file=fnmmegv, access  ='append',status='replace') 
open(fmmegp, file=fnmmegp, access ='append',status='replace')

! Initialization of parameters to do numercial continuation
tdir  = 1 ! integrate forward

print*, 'check wi to determine which family to study'
print*, 'wi=',  wi

! Use the period as initial guess, TP = 2*pi/wi, wi is the pure imaginary eigenvalue
if( mod(ivr,2) == 1) then 
  ifam = ivr/2+1
else 
  ifam = ivr/2
endif  

tp0 =  2*pi / dabs(wi(2*ifam))  
print*, 'tp0=', tp0, 'wi=', wi(2*ifam)
read*
 
print*, 'Input ind,  add constraint:  p0 = poinst(ind)'
read(*,*) ind
p0 = eq(ind)
    
call init_gr_po(idcor, ind, p0) 
allocate(x0_pre(2, nctr)) 
allocate(x0( nctr)) 
allocate(dg(ntar, nctr)) 


! plot the initial guessed orbit  to see if it is a periodic orbit
open(333, file = './dat/po_guess.dat', access ='append',status='replace')
call plob(poinst, 0.d0,  tp0,  6, 6, 1, 1, 333 , gr_lf, gr_cjlf, pof)  
print*, 'check the initial guessed orbit!'; close(333); read*

!  one parameter continuation, along peroid or energy, here, we take energy.... 
! the second p.o. is obtained by change a little bit the energy, 
! and the succesive ones are  obtained by extrapolation using the previous two. 

ipo = 0

x0(1:nc)   = poinst(ctr) 
para0      = tp0  

if(tfree == 1) x0(nctr) = tp0 
x0_pre(1,:) = x0 ;  para1(1) = tp0 

! x0 inclue all the control variables, while y0 only refer to the state vector 

do  
!  if(ipo > 500) exit
  if(ipo > npo) exit !
  
  ! the new guess for periodic orbit 
  y0(ctr) = x0(1:nc); tp = para0 
  if(idcor == 7 .or. idcor == 1)  y0(ind) = p0 
  if(tfree == 1)  tp = x0(nctr)
  
  call gr_cjlf(y0, cj)
  print*, 'Check the engery: ', cj
  print*, 'tp0, y0';  print*, tp, y0;  read*
  
  ! refinement using Newton method 
  call refn_po(y0, tp, yf, dg, niter, ispo, gr_lf, gr_cjlf) 
  print*, 'The refined period', tp 
  print*, 'IC   \\ FC'
  print*, y0
  print*, yf(1:6)
  print*; read*
  
  if(ispo == 0)  then 
    dpara = dpara / 2.d0
    print*, 'Too big step in continuation, fail! Decrease by half', dpara
    
    if (dabs(dpara) < 1.d-8) then 
      print*, 'Too small step of H for continuation! <1.d-8, stop the continuation', dpara;
      read* 
      exit
    endif  
    
    ! use the two previous good p.o.s for predicition with smaller stepsize
    call pred_h(x0_pre, para1, dpara, ipo, x0, para0)
    cycle 
    
  else 
    
    ! we obtain a new periodic orbit
    ipo      = ipo + 1
    x0(1:nc) = y0(ctr)
    para0    = tp 
    
    if(tfree == 1)   x0(nctr) = tp 
    
    print*, 'before update.  for ', ipo, '-th p.o.'
    print*, 'tp, x0', tp, x0 
    print*, 'y0', y0 
    print*, 'para0=', para0
    print*, 'para1=', para1   
    print*, 'x0', x0 
    print*, 'x0_pre'
    do i = 1, 2 
      print*, x0_pre(i,:)
    enddo 
    print*; read*
    
    ! -- update x0_pre and para1 
    if(ipo > 2) then 
      x0_pre(1, :)  = x0_pre(2,:)  
      para1(1)    = para1(2)
      
      x0_pre(2, :)  = x0
      para1(2)      = para0 
   
    else 
      x0_pre(ipo, :) = x0
      para1(ipo)     = para0
    endif 
    
    print*, 'after update'  ! 0--- the peroid is not updated
    print*, para1   
    print*, 'x0_pre'
    do i = 1, 2 
      print*, x0_pre(i,:)
    enddo 
    print*; read*
    
    ! -- modify the stepsize based on the iterations of Newton method 
    if(niter < 3)  dpara = dmin1( dabs(2*dpara), dhmax) * dpara / dabs(dpara)
    if(niter > 5)  dpara = dpara / 2.d0 
    
    ! --- plot po for 1 period -----
    open(fpoinit, file = fnpoinit, status='old', access='append')
    write(fpoinit, *) tp,  x0, cj0;  close(fpoinit)
    write(*, *) tp,  x0, cj0 ; print*; read*  ! print to screen 
    
    open(fpo, file = fnpo, status= 'old', access='append')
    
    ! --- plot the orbit for one revolution and compute the Eigenvalues of the Monodromy Matrix to check the stability
    pvi(1:6)    = y0
    pvi(7:42)   = 0.d0
    pvi(7:42:7) = 1.d0
    
!    subroutine plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 
    call plob(pvi, 0.d0, tp, 6, 42,  1, 1,  fpo, gr_lf, gr_cjlf,  pvf) 
    
    print*; print*, 'Check pvi, pvf:'
    print*, pvi(1:6)
    print*, pvf(1:6)
    print*; read*
    
    phi = reshape(pvf(7:42), (/6,6/))
    
    print*;  print*, 'Check MM:'
    do i = 1, 6, 1
      print*, phi(i,:)
    end do
    print*;! read*
    
    call eigrg(phi, 6, 1, wr_mm, wi_mm, vr_mm)

    print*, 'eigenvalues, real part'
    print*,  wr_mm 
  
    print*, 'eigenvalues, imaginary part'
    print*,  wi_mm 
   
    call prt_eigval( 6, fmmegv, wr_mm, wi_mm )
    call prt_eigvec( 6, fmmegp, wi_mm, vr_mm )
  
    write(fmmegp,*) ! add a blank line to seperate eigenvector matrix
   
    print*, 'Check the type of P.O.:' 
    print*, 'Finish one p.o.!';  print*; read* 
    
    close(fpo) 
    
    print*, 'Call pred_h to compute the new guess!'
    
    ! ---- continuation ------------
    ! since we have only 1 free parameter to continue along: x 
    ! we take a step forward and update the energy level  
  
    ! ---- prediction for the guess of a new p.o -------
    !  using linear extrapolation based on the previous 2 p.o.s
!    subroutine pred_h( x0_pre, h_pre, dh, ipo, x0, h)
    call pred_h(x0_pre, para1, dpara, ipo, x0, para0)
    
!    if( debug == 1 ) then
      print*, 'Initial guess for', ipo, '-th new p.o.:'
      print*,  x0, para0 ; read*
!    end if
    
  endif  ! ispo == 1
 
enddo   


stop
end program main_po_asymt



















  
