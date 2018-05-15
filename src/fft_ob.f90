subroutine fft_ob(x0, n, lt, np,  nitmax, nfmax, &
                  xmax, fob, ffc, ffbas, nf_bas, fa_bas, isrpt, deriv, gr_cj) 
! freq unit:  cycle / unit  time  


! 2016-08-03 12:30:07  
!    use do loop to deal with all the n components, this is especially for the Fourier analysis of an orbit .
!    the orbit is saved in file with tag fob, 
!    and detected frequecies and coefficients are saved in ffbas 

                 
! for a given p.o.,   with X0 + TP0, use furian from fft_mod to detect domimant freqs
 
! number of sample better to be a power of 2 
! 2^12 is enough (only need to satisfy Nyquist Criteria)

! fs = (NP-1) / LT > R_N = 2 * B      ===>      NP > 2 * LT * B + 1
! for the tested function, B = 27, take time length LT = 10 
!    NP > 2 * 10 * 27 + 1 =  541,  with 2**10 = 1052,  2**9 = 512,  2**8 = 256

! with the refinement algorithm we proposed, only Nyquist Criteria is needed to be satisfied....
 
!       Input Variables 
! x0           initial state
! n            number of dimension of the phase space(x0)
! lt           the time length for sampling: for po, it is usually 1 period  
! np           nuber of sampling points 

!***********  Output Variables  **************
! nf_bas       number of basic frequecies detected 
! fa_bas         basic frequecy + amplitude
!              use fixed dimension 5*2n for the convenience of passing parameters

!  --- For Fourier frequencies extraction  
! tol          Tolerance in correction for refinement 
! tolres       Maximum amplitude of the residual 
! nitmax       maximum iteration for Newton method to do the refinement 

! nfmax        maximal number of frequencies to keep      
! xmax         detect escape, non-sense to do analysis for the escape segment

!  ---- the results are saved in these files ----
! fob          file tag for orbit 
! ffc          file tag for the frequencies + coefs 
! ffbas        file tag for the basic frequencies and the error in the linear combination  


!  Routine Used:
!    furian   from module fft_mod 

! Finally Revised by Yu -- 20160614
!----------------------------------------
use fft_mod
use dp_mod
implicit none

  
! Input  and Output Declaration
integer, intent(in)             ::  n, np, nitmax,  fob, ffc, ffbas
integer, intent(inout )         ::  nfmax ! probably we will update the value of this 

real(kind=dp), intent(in)       ::  x0(n), lt, xmax !,  tol, tolres
 
integer, intent(out)            ::  nf_bas, isrpt 
real(kind=dp), intent(out)      ::  fa_bas(5, n) ! freq+amp


! Local variables
integer::  i, j, sep, fob_ref, nf,  isbas,  ordmax, nfbas, isescp 
           

!fx_fft, fy_fft, fz_fft, fvx_fft, fvy_fft, fvz_fft !  --fft

real(kind=dp) ::  xf(n), xl(np,n), x(np), h,   & !fft
                  tst, tend, &  ! evaluate time  
                  fcout(100,5), resmax, fc(100,5), fbas(20,5), delmax, tolf ! frebas -- to keep coherent with the size of array frebas
                  
external :: deriv, gr_cj 


! Initialize for module fft_mod
call init_fft( np, lt, nfmax, nitmax)

! For basic frequecies detection
nf_bas = 0 
fa_bas = 0.d0 

fob_ref = 37
open(fob_ref, file = './dat/ob_ref.dat', access='append', status='replace')

!ffbas   = 38
!open( ffbas, file = './dat/fbas.dat',    access='append', status='replace')


! the fixed stepsize for the integration of the orbit
h = lt / np ! time interval  for integration 
 
print*, 'h, lt, np = ', h, lt, np
read* 

! ****************************** orbit sampling ******************************
! use the automatic step contol as a reference 
print*, 'Before plob, x0 = '
print*, x0
!subroutine plob(y0, t0, tf, ndim, nvar, tdir, ispl, ftag, deriv, gr_cj,  y) 
call plob(x0, 0.d0, lt, n, n, 1,  1, fob_ref, deriv, gr_cj,  xf) 

 close(fob_ref)
print*, 'lt, xf', lt, xf; print*   
read*
  
xl   = 0.d0 ! to avoid mistake, for every orbit, reinitialize the state vector
!xmax = 1.d0 ! the maximum value to detect escape, 1 is enough i think...
      
call cpu_time(tst)    
 
print*, 'Before plob_fxd, x0='
print*,  x0 ; print*

!subroutine plob_fxd(y0, n, t0, h0, np, xmax, ftag, yall,  isescp, deriv, gr_cj) 
call plob_fxd(x0, n, 0.d0, h, np, xmax, fob, xl, isescp, deriv, gr_cj) 

if(isescp == 1) return

call cpu_time(tend)
  
print*, 'Elapsed time for orbit sampling is ', tend-tst; print*;

! check the final state and stop time? 
print*, 'lt, yf', np*h, xl(np,:)
print*
read*
    
!! use the components of the state vector from the integration with fixed stepsize    
!x = xl(:, 1) ! x
!y = xl(:, 2) ! y 
!z = xl(:, 3) ! z 
 
! ************* Furier Analysis  *************************

! ---- Detect the basic frequencies ----------- 
print*; print*, '****Detect the basic frequencies(1 = Yes)****'
read(*,*) isbas 

! -- component by component  -- 
! -- Initialize of ordmax and tolf
if(isbas == 1) then
!-- general 
!  ordmax = 10
!  tolf   = 1.d-6

! for detection of frequecies
  ordmax = 40
  tolf   = 1.d-3
endif 

! we only do for position, forget about the velocity at this moment
do i = 1, n/2, 1
  nfmax = 5
  x = xl(:, i ) 
  
  write(ffbas,*); write(ffc,*) '# ****************', i, '-th component   ****************' 
  call furian( x, ffc, nf, fcout, sep, resmax)     
  
      
  if(isbas == 1) then
  
    fc(1:nf,:) = fcout(2:nf+1, :) 
  
    print*, nf, 'Frequencies to be detected :'
    do j = 1, nf 
      write(*, '(I5, 4e24.14, f10.0)')  j, fc(j, :)
    enddo   
    write(*, *)  ! add a blank line 
   
   print*, 'Not satisfied with', nfmax, 'frequecies, use more frequecies ? (1 = Yes)' 
   read*, isrpt
   if(isrpt == 1) then 
     print*, 'Input the new value of nfmax'
     read*, nfmax 
     return
   endif   
      
    !  subroutine frebas( fc, nf,  nfbas, fbas, ffbas)
    call frebas( fc, nf, ffc, ordmax, tolf,  nfbas, fbas, delmax)
    nf_bas = max0(nfbas, nf_bas)
    
    write(ffbas,*); write(ffbas,*) '****************', i, '-th component   ****************' 
    write(ffbas,*) '# Total number of basic freqs:', nfbas
    write(ffbas,*) '# Max error of the linear combinations= ', delmax
    
    do j = 1, nfbas
       write(ffbas,'(4f26.14, f10.0)') fbas(j, :)
       fa_bas(j, 2*i-1:2*i) = fbas(j, (/1,4/)) ! freq+amp 
    enddo
      
  endif
  
  read*
end do

print*, 'check fa_bas! nf_bas = ', nf_bas
do i = 1, nf_bas, 1
  print*, fa_bas(i, :)
end do
print*; read* 

return  
end subroutine fft_ob

















  
