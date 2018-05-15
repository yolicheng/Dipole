program main_eqeig

!  2016-06-26 08:26:55 
! compute the behavior of eigvalues as a function of beta, the problem is that we need to sort it . 

! 2016-12-26 18:48:59 
!  After half a year, finnally I succeed with the recovery of this routine 

use dp_mod
use pi_mod
use lf_mod ! the system related parameters of lorentz force problem
use sort_mod 
implicit none

integer, parameter ::  n = 6  ! the dimension of LF problem 
    
! lf_mod Module Variables
! Global :   beta, cs, sgn, eq, ieq 
real(kind=dp) ::  beta0,  dbt, dbt0


! Local Variables
integer :: i,  iter,  feig, ind, isgn, isdif, ind_st(n)

real(kind=dp) ::   dlf(n,n), wr(n), wi(n), vr(n,n), wr0(n), wi0(n), &  ! eigenvalues of variational matrix
                   wr_st(n), wi_st(n)
                 
 character(len=70) ::  fneig           


! case 1: normal;       2: radial    3: tangential  --- updated sequence 

! If the equilibirum points with different signs are not symmetric, the periodic orbits will also be different, we have to compute seperately 
call init_lf 
 
print*, 'Check the global value : cs, sgn, eq, ieq '
print*, cs, sgn, eq, ieq; print*; read*
 

feig= 28
write(fneig, fmt='(a,i0,a,i0,a)')   './dat/eq/cs', cs, '_eq', ieq, '_eig.dat'
print*,   fneig;  print*; 
open(feig, file=fneig, access  ='append',status='replace') 
write(feig, *) '# Beta --- ', n, 'eigenvalues (real / imag)'
write(feig, *) '# eq: ', eq 

print*, 'Change the sign of eq?'
read*, isgn 
isdif = 1
if(isgn == 1) then 
  print*,'Input the index to change sign!' ; read*, ind 
  eq(ind) = -eq(ind)
  isdif = -1
endif 
print*, 'isdif, New eq:', isdif, eq;  print*; read*       
  


beta0 = -10.d0   ! cs1_eq3,  eq2 

!beta0 = -15.d0  ! cs1_eq1 

!3.4520000000    0.0000000000    8.7387868653    0.0000000000   -8.7387868653   -0.0116639529    0.8285276755   -0.0116639529   -0.8285276755    0.0116639529    0.8285276755    0.0116639529   -0.8285276755
!3.4530000000    0.0000000000    8.7420201389    0.0000000000   -8.7420201389    0.0000000000    0.8174240294    0.0000000000   -0.8174240294    0.0000000000    0.8396379265    0.0000000000   -0.8396379265
!beta0 = 3.451d0  ! cs2_eq3  -- CKD! the distances are not the same:  for  0.8174240294  <  0.8396379265 , no problem in eig_sort

dbt  = 1.d-2 ; dbt0 = dbt 


iter = 0


do  
  if(beta0 > 10.d0) exit  ! cs1_eq2, eq3  ! done 
  
!  if(beta0 > 30.d0) exit   ! cs1_eq1    [-10:25] ! --done 2017-03-03 18:27:01 
  
!  print*, 'beta0=', beta0; print*; read*
  
  iter = iter+1
  
  ! cs2_eq3, smaller dbt to remove the gap between 0 lines.  3.4525 -0.9516  
  if(cs==2 .and. ieq==3) then 
    if(dabs(beta0 - isdif*3.4525d0) < 2*dbt0 .or. dabs(beta0 - isdif*(-0.9516d0)) < 2*dbt0 ) then  
      dbt = 1.d-4
    elseif( dabs(beta0 - isdif*0.0732) < 2*dbt0 .or. dabs(beta0 - isdif*(1.5326)) < 2*dbt0 ) then
      dbt = 1.d-4
    else !if !(dabs(deta0) < dbt0) then 
      dbt = dbt0
    endif 
  endif 
  
!  print*,'dbt0, isdif', dbt0, isdif; print*;read*
  
  beta0 = beta0 + dbt 
  print*, 'beta0, dbt,dbt0, isdif', beta0, dbt,dbt0, isdif; 
!   print*; read* 
  
  call init_beta(beta0)
  call dflrtz(eq, dlf)
  
!  do i = 1, 6, 1
!    print*, dlf(i, :)
!  end do
!  print*; read*
  
  ! compute the eigenvalues and eigenvectors  of dlf
  call eigrg(dlf, n, 0, wr, wi, vr)
  
  if(iter == 1) then 
    wr0 = wr 
    wi0 = wi 
    
!    subroutine eig_sort( wr0, wi0, wr, wi, n, is_one, wr_st, wi_st, ind_st)
    call eig_sort( wr0, wi0, wr, wi, n, 1,  wr_st, wi_st, ind_st)
  
  else 
  
    call eig_sort( wr0, wi0, wr, wi, n, 0,  wr_st, wi_st, ind_st) 
    
  endif   
  
  wr0 = wr_st; wi0 = wi_st
  
  write(feig, '(1f16.10)', advance='no')  beta0 
  
!  write(*, '(1f16.10)', advance='no')  beta0  ! print to screen
  
  do i = 1, n, 1
    
    write(feig, '(2f16.10)', advance='no')  wr_st(i),  wi_st(i)
    
   ! write(*, '(2f16.10)', advance='no')  wr_st(i), wi_st(i) ! print to screen
  end do
  write(feig,*)
  
!  write(*,*) ! print to screen

!  print*; read*
  
enddo 


stop
end program main_eqeig



















  
