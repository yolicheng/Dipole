program test
! this is for test of  simple fortran functions 
use dp_mod
use pi_mod
implicit none

! Variables
real(kind=dp) ::  a, alog2 ! logrithm base 2 

! -- for multidimensional array ----
integer       :: i, j, iseed(4), idist, &
                 nf, ind , &
                 nr, nc, nrhs, info
                 
real(kind=dp) ::  b(2,2,2), b2(8), bb(2), bvec(2,2,2), b1d(8), &
                  a1(3), a2(3), a12(3), &
                  arg, cs, si, & ! test trigometric functions   
                  rho, x, dx, f, f2, & ! test the phase shift of a Fourier curve   
                  tann(20000), cotn(20000), & !
                  la(2,2), lb(2), lx(2), invla(2,2) ! test deltx + lapack 
 
integer :: rank 
real(kind=dp) ::  arr0(4,6), arr0_base(6,6), arr1(3,4), arr1_base(4,4), arr1T(4,3),arr1T_base(3,3)                 
!  do zero velocity curves                  
real(kind=dp) ::  h0, xf                   
 
              
! we don't even need to declare external for the routines from the library '
!external :: dlarnv!, cross_product ! from LAPACK to generate normal distribution random data 
! check cotangent(n*x) in recurrent method  ! 

! ------ check the null space ------
!A=\left[{\begin{array}{cccccc}1&0&-3&0&2&-8\\0&1&5&0&-1&4\\0&0&0&1&7&-9\\0&0&0&0&0&0\end{array}}\,\right].
arr0(1,:) = (/1, 0,-3,0,2,-8/)
arr0(2,:) = (/0, 1,5, 0,-1,4/)
arr0(3,:) = (/0, 0,0, 1,7, -9/)
arr0(4,:) = 0.d0 


call null_cal(4, 6, arr0, rank, arr0_base)
do i = 1, 6
 print*, arr0_base(i,:)
enddo 
print*; read*


arr1(1,:) = (/1,1,1,1/); arr1(2,:) = (/1,2,3,4/); arr1(3,:) = (/4,3,2,1/); 
call null_cal(3, 4, arr1, rank, arr1_base)

print*, 'the null space of A!'
do i = 1, 4
 print*, arr1_base(i,rank+1:4)
enddo 
print*; read*


arr1T= transpose(arr1)
call null_cal(4, 3, arr1T, rank, arr1T_base)

print*, 'the null space of A!'
do i = 1, 3
 print*, arr1T_base(i,rank+1:3)
enddo 
print*; read*


! do zero velocity curves
h0 = 4.3767487109222252 
xf = 2.d0
!call zvc_plf( h0, xf)

stop 


!subroutine deltx( nr, nc, nrhs, a, b, x, info) !check this routine 
la(1, :) = (/2.d0, 1.d0/); la(2, :) = (/3.d0, -2.d0/)
lb = (/4.d0, -1.d0/)
nr = 2;  nc = 2; nrhs = 1
!subroutine inv(A, nrow, ncol, Ainv)

call inv(la, nr, nc, invla)
lx = matmul(invla, lb) 
print*, 'x = ', lx ; print*;read* 

print* , 'by Lapack'
la(1, :) = (/2.d0, 1.d0/); la(2, :) = (/3.d0, -2.d0/)
call deltx(nr, nc, nrhs, la, lb, lx, info) 
print*, 'x = ', lx ; print*
stop 

x = 0.167743d0 
print*, 'cot(x) = 1/tan(x) = ',  1.d0 / dtan(x) 
read*
open(101, file = './dat/cotn.dat')
call trigrec_cot(x, 20000, cotn)
do i = 1, 20000, 1
  write(101, '(3e24.14)')  x*i, cotn(i),  1.d0 / dtan(x*i)
end do

! check tangent(n*x) in recurrent method  !ckd 
open(101, file = './dat/tann.dat')
x = 0.167743d0 
call trigrec_tan(x, 20000, tann)
do i = 1, 20000, 1
  write(101, '(3e24.14)')  x*i, tann(i), dtan(x*i)
end do
stop

 
! -- ckd, a periodic function with different angles are only a shift in phase.  
! test f(x) = 1 + 0.5sin(x) + 0.25*cos(x) 
! and do a phase shift: f(x+rho) =  1 + 0.5sin(x+rho) + 0.25*cos(x+rho)
! where x, rho \in [0, 2pi], and check the two curves 
x = -3*pi+0.321d0 
print*, 'mod(x>-3pi, 2pi)', x, dmod(x, pi2)
read*
stop

open(100, file = './dat/phaseshift.dat')
rho = 0.177764d0*pi2
dx = pi2/100
do i = 1, 100, 1
  x = (i-1)*dx 
  f = 1.d0 + 0.5d0*dsin(x) + 0.25d0*dcos(x)
  x = x+rho
  f2 = 1.d0 + 0.5d0*dsin(x) + 0.25d0*dcos(x)
  write(100,*) x-rho, f, x+rho, f2
end do

stop
 
! --- test trigrec_c and trigrec_s 
do i = -2, 2, 1 
  arg = pi + i * pi/20  
  print*, 'arg=', arg
  print*, 'intrinsic cos and sin, ', dcos(arg), dsin(arg)
  
  call trigrec_c(arg, 1, cs)
  call trigrec_s(arg, 1, si)
  print*, 'trigrec_c and trigrec_s, ', cs, si 
  read*
enddo 
! test the function cross_product ! --- tested,
  a1 = (/3,2,5/)
  a2 = (/1,7,4/)  
  
  call cross_product(a1, a2, a12)
  print*, 'cross(a1, a2) = ', a12 

! discard--- failed... because function cross does not have an explicit interface.... 
! if we declare as real(kind=dp), dimension(3) :: cross, the compiler sees it as an array. 

!  a12 = cross(a1, a2)
!  print*, 'By function: cross(a1, a2) = ', a12 
! ----  test for loop with a start index less than the end index 
! --tested! if i_start < nf, the block will not be executed....
!nf = 2
!do i = 3, nf, 1
!  print*, i, '-th iteration!, nf=', nf
!  read*
!end do

!stop 


! --- multidimensional array, to extract subarray in low dimension  --- tested 

! ----------- random data in normal distribution --------------- 
idist = 1 ! normal distribution      = 1:  uniform (0,1)  = 2:  uniform (-1,1) 
iseed = (/0,0,0,1/)
call dlarnv( idist, iseed, 8,  b2)
b = reshape( b2, (/2,2,2/) )


! ---- test reshape a multi-dimension array 2-by-2-by-2 into a column vector  
! ---tested, to get 1d array, use (/8/) as the second input 
b1d = reshape(bvec, (/8/))
print*, '1d array: b1d= ', b1d

! --- test maxloc 
print*,  maxloc(b1d) 
read*


! ---- vector operation (Arithmetic) in fortran --  
bvec = b*4 - b !tested,  in fortran90 the intrinsic functions operate component-wise on arrays


print*, 'b = : '
do i = 1, 2 
  do j = 1, 2 
    print*, b(i,j,:)
  enddo 
enddo   

print*, 'a = : '
do i = 1, 2 
  do j = 1, 2 
    print*, bvec(i,j,:)
  enddo 
enddo

bb = b(1,1,:) ! extract submatrix
print*, 'bb=', bb


! logrithm base 2  --tested 
a = 2.d0**4 

alog2 = dlog(a)/ dlog(2.d0)
print*, 'a=', a, 'log_2(a)=', alog2

end program test
