subroutine fft_ob(x0, lt, np, tol, tolres, nitmax, maxnf, &
                  xmax, fob, ffc, deriv, gr_cj) 
                  
! for a given p.o.,   with X0 + TP0, use furian from fft_mod to detect domimant freqs
 
! number of sample better to be a power of 2 
! 2^12 is enough (only need to satisfy Nyquist Criteria)

! fs = (NP-1) / LT > R_N = 2 * B      ===>      NP > 2 * LT * B + 1
! for the tested function, B = 27, take time length LT = 10 
!    NP > 2 * 10 * 27 + 1 =  541,  with 2**10 = 1052,  2**9 = 512,  2**8 = 256

! with the refinement algorithm we proposed, only Nyquist Criteria is needed to be satisfied....
 
!       Input Variables 
!  x0           initial state
!  lt           the time length for sampling: for po, it is usually 1 period  
!  np           nuber of sampling points 
! tol, tolres, nitmax      control of error for refinement


!  maxnf        maximal number of frequencies to keep      
!  xmax         detect escape, non-sense to do analysis for the escape segment
! fob           file tag for orbit 
! ffc           file tag for the frequencies + coefs 


!       Output Variables 
!  a            

!  Routine Used:
!    furian  

! Finally Revised by Yu -- 20160614
!----------------------------------------
use fft_mod

implicit none

integer, parameter :: dp = kind(1.d0) 
  
! Input  and Output Declaration
integer, intent(in)             ::  np, nitmax, maxnf, fob, ffc 
real(kind=dp), intent(in)       ::  x0(6), lt,  tol, tolres, xmax 
 

! local variables
integer::  i, j, isft, fx_fft, fy_fft, fz_fft, fvx_fft, fvy_fft, fvz_fft, fob_ref !  --fft

real(kind=dp) ::  xf(6), xl(np,6), x(np),y(np), z(np),  & 
                  vx(np), vy(np), vz(np), h,  & ! fft  
                  tst, tend, &  ! evaluate time  
                  
                  !the domimant frequency which has the first maxnf maximum amplitude
                  fqmx_x(maxnf+1, 4), fqmx_y(maxnf+1, 4), fqmx_z(maxnf+1, 4), fqmx_vx(maxnf+1, 4), fqmx_vy(maxnf+1, 4), fqmx_vz(maxnf+1,4) , &
                  fourx, foury, fourz,  fourvx, fourvy, fourvz,   t 
      
DOUBLE COMPLEX   fft1(np),fft2(np), fft3(np),fft4(np), fft5(np),fft6(np) 
! character(len=70) ::  fnpo, fnx, fny, fnz, fnvx, fnvy, fnvz, fna, fnmx  

external :: deriv, gr_cj 

call init_fft( np, lt, tol, tolres, nitmax)
      
!fpo_fft = 20; fx_fft = 30; fy_fft = 31;  fz_fft = 32;   ! save x-y-z and  vx-vy-vx
!fvx_fft = 33; fvy_fft = 34;  fvz_fft = 35;
! 
!fqmx = 36; ! save the dominant magnitude and frequency

fob_ref = 37
open(fob_ref, file = './dat/ob_ref.dat', access='append', status='replace')


! I think there is no need to save each data.... 
!write(fnpo, fmt='(a)') './dat/poformfly/fft/po_fft.dat'
!write(fnx,  fmt='(a)') './dat/poformfly/fft/x_fft.dat' 
!write(fny,  fmt='(a)') './dat/poformfly/fft/y_fft.dat' 
!write(fnz,  fmt='(a)') './dat/poformfly/fft/z_fft.dat' 
!write(fnvx,  fmt='(a)') './dat/poformfly/fft/vx_fft.dat' 
!write(fnvy,  fmt='(a)') './dat/poformfly/fft/vy_fft.dat' 
!write(fnvz,  fmt='(a)') './dat/poformfly/fft/vz_fft.dat' 

!! open an empty file using the fname

!open(fpo_fft, file= fnpo, access ='append',status='replace')  
!open(fx_fft, file= fnx, access ='append',status='replace')  
!open(fy_fft, file= fny, access ='append',status='replace')
!open(fz_fft, file= fnz, access ='append',status='replace')

!open(fvx_fft, file= fnvx, access ='append',status='replace')  
!open(fvy_fft, file= fnvy, access ='append',status='replace')
!open(fvz_fft, file= fnvz, access ='append',status='replace')


!write(fnmx, fmt='(a)') './dat/poformfly/fft/fqmx.dat' 
!open(fqmx,  file= fnmx, access ='append',status='replace')  
!write(fqmx,*) '# fx	ax 	fy 	ay 	fz 	az'

! the fixed stepsize for the integration of the orbit
h = lt / np ! time interval  for integration 
 
print*, 'h, lt, np=', h, lt, np
read* 

! ****************************** orbit sampling ******************************
! use the automatic step contol as a reference 
print*, 'Before plob, x0 = '
print*, x0
call plob(x0, 0.d0, lt, 6, 1,  fob_ref, 0, deriv, gr_cj,  xf) 
 close(fob_ref)
print*, 'lt, xf', lt, xf; print*   
read*
  
xl   = 0.d0 ! to avoid mistake, for every orbit, reinitialize the state vector
!xmax = 1.d0 ! the maximum value to detect escape, 1 is enough i think...
      
call cpu_time(tst)    
 
print*, 'Before plob_fxd, x0='
print*,  x0 ; print*

!  subroutine plob_fxd(y0,t0, h0, np, xmax, ftag, yall, deriv, gr_cj)
call plob_fxd(x0, 0.d0, h, np, xmax, fob, xl, deriv, gr_cj) 

call cpu_time(tend)
  
print*, 'Elapsed time for orbit sampling is ', tend-tst; print*;

! check the final state and stop time? 
print*, 'lt, yf', np*h, xl(np,:)
print*
read*
    
    
!! use the components of the state vector from the integration with fixed stepsize    
x = xl(:, 1) ! x
y = xl(:, 2) ! y 
z = xl(:, 3) ! z 
  
  
!vx = xl(:, 4) ! vx
!vy = xl(:, 5) ! vy 
!vz = xl(:, 6) ! vz 
 
! ************* Furier Analysis  *************************
! -- for x ---
call furian( x,  nf, fcout, sep)     
write(ffc, *) ' #  For x:     feq    cc    ss    amp'
write(ffc, *) ' #  detect ', nf, ' frequecies! '
write(ffc, *) ' #  Maximal modulus of residual = ', resmax
write(ffc, *)
 
do i = 1, nf+1 
  write(ffc, *) fcout(i, 1:4)
enddo   
write(ffc, *)  ! add a blank line  

! -- for y --- 
call furian( y,  nf, fcout, sep, resmax)   
  
write(ffc, *) ' #   For y:     feq    cc    ss    amp'
write(ffc, *) ' #  detect ', nf, ' frequecies! '
write(ffc, *) ' #  Maximal modulus of residual = ', resmax
write(ffc, *)
 
do i = 1, nf+1 
  write(ffc, *) fcout(i, 1:4)
enddo   
write(ffc, *)  ! add a blank line  

! -- for z --- 
call furian( z,  nf, fcout, sep, resmax)   
  
write(ffc, *) ' #   For z:     feq    cc    ss    amp'
write(ffc, *) ' #  detect ', nf, ' frequecies! '
write(ffc, *) ' #  Maximal modulus of residual = ', resmax
write(ffc, *)
 
do i = 1, nf+1 
  write(ffc, *) fcout(i, 1:4)
enddo   
write(ffc, *)  ! add a blank line  

  
! ---- discard dtwofft due to low accuracy ------    
call cpu_time(tst)    
call dtwofft(x, y,fft1,fft2, np)   ! for x and y
call dtwofft(z, vz,fft3,fft6, np) ! for z and vz
call dtwofft(vx, vy,fft4,fft5, np) ! for vx and  vy

call cpu_time(tend)  
print*, 'Elapsed time for  3 * dtwofft is ', tend-tst
  

! Extract real magnitude and frequency from fft1 returned by twofft, save them in seperate files   
print*,'#********* for x **********'  
write(*, *) '# ********* x *********'
 
isft = 0 ! do not save all the frequencies and corresponding amplitude

!subroutine fqext(fft, np, lt, f_fft, maxnf, fqmx)
call fqext(fft1, np, lt, isft, fx_fft, maxnf, fqmx_x) 
  
!write(*, *) '# ********* y *********'
call fqext(fft2, np, lt, isft,  fy_fft, maxnf, fqmx_y) 
  
  
!  write(*, *) '# ********* z *********'
call fqext(fft3, np, lt, isft, fz_fft, maxnf, fqmx_z) 
  
  
!!  write(*, *) '# ********* vx *********'
!call fqext(fft4, np, lt, isft, fvx_fft, maxnf, fqmx_vx)  
!  
!!  write(*, *) '# ********* vy *********'
!call fqext(fft5, np, lt, isft, fvy_fft, maxnf, fqmx_vy) 
!  
!!write(*, *) '# ********* vz *********'

!call fqext(fft6, np, lt, isft, fvz_fft, maxnf, fqmx_vz) 

! use the first maxnf frequencies that has the largest magnitude  to recover the Fourier series  
do i =  1, np, 1
  t = (i-1) * h 
  call fourseri(t,  maxnf,  fqmx_x, fourx )   !x
  call fourseri(t,  maxnf,  fqmx_y, foury )   !y
  call fourseri(t,  maxnf,  fqmx_z, fourz )   !z
  
! save the difference between the real one and the  Fourier Series      
!  print*, 't, x, y, z',  t, x(i), y(i), z(i), vx(i), vy(i), vz(i)
!  print*, 'Fourier: ',   t, fourx, foury, fourz, fourvx, fourvy, fourvz
!  print*, 'dx=fourx-x',  t, fourx-x(i),  foury-y(i),  fourz-z(i), fourvx-vx(i),  fourvy-vy(i),  fourvz-vz(i) 
!  
  write(fdrmx, '(12f14.8)'),  t, x(i), fourx,  dabs(fourx-x(i)), y(i), foury,  dabs(foury-y(i)), z(i), fourz,  dabs(fourz-z(i))
!  read* 
end do
  

!write(fqmx, *) '#  fx  - ax  - rex  - imx  -- fy  - ay  - rey  - imy  -- fz  - az  - rez  - imz '
!write(fqmx, *) '#  fvx - avx - revx - imvx -- fvy - avy - revy - imvy -- fvz - avz - revz - imvz '

do j = 1, maxnf+1
   write(fqmx,'(6f14.8)')   fqmx_x(j,1), fqmx_y(j,1),fqmx_z(j,1), fqmx_x(j,2), fqmx_y(j,2), fqmx_z(j,2)
   
!  write(fqmx,'(12f14.8)')   fqmx_x(j,1), fqmx_y(j,1),fqmx_z(j,1),fqmx_x(j,2:4), fqmx_y(j,2:4),  fqmx_z(j,2:4)
!  write(fqmx,'(12f14.8)')   fqmx_vx(j, :), fqmx_vy(j, :),  fqmx_vz(j, :) 
  
!  write(*,'(12f14.8)')   fqmx_x(j, :),  fqmx_y(j, :),  fqmx_z(j, :)
!  write(*,'(12f14.8)')   fqmx_vx(j, :), fqmx_vy(j, :),  fqmx_vz(j, :) 
!!  write(fqmx, *);    print* 
enddo   
  
  
! Pick the first maxnf=4 dominant frequencies and then save them in the form 
!! fx- ax - fy - ay - fz - az, keep 4 digits after the decimer, for convenient use in latex tabular

write(fqmx, *);    print*   

return  
end subroutine fft_ob

















  
