program main_mm_eig

! sort the eigvalues of MMat, to check the stability of a periodic orbit. 

use dp_mod
use pi_mod
use sort_mod 
implicit none

integer, parameter ::  n = 6  ! the dimension of LF problem 
    
! Local Variables
integer :: i, j,  iter,  feig, feig_st, find_st, debug, iter0, ind_st(n) 

real(kind=dp) ::   wr(n), wi(n),  wr0(n), wi0(n), &  ! eigenvalues of  mm
                   wr_st(n), wi_st(n), w0(2*n)
                 
 character(len=70) ::  fneig, fneig_st, fnind_st           


feig = 18;  feig_st = 28; find_st = 38
! ----  remember to modify the name -----
!fneig =  './dat/po/case1/eq3/ifam1/bt2/pommegv.dat'          
!fneig_st =  './dat/po/case1/eq3/ifam1/bt2/mmeig_st.dat' ! remember to modify the name 
!fnind_st  =  './dat/po/case1/eq3/ifam1/bt2/mmeig_ind_st.dat'  

!fneig    =  './dat/po/case1/eq3/pommegv.dat'           
!fneig_st =  './dat/po/case1/eq3/mmeig_st.dat'  
!fnind_st  =  './dat/po/case1/eq3/mmeig_ind_st.dat'  

!fneig    =  './dat/po/case1/eq1/beta2/pommegv.dat'          
!fneig_st =  './dat/po/case1/eq1/beta2/mmeig_st.dat' ! remember to modify the name 
!fnind_st =  './dat/po/case1/eq1/beta2/mmeig_ind_st.dat'  


!fneig    =  './dat/po/case1/eq1/beta2/pommegv.dat'          
!fneig_st =  './dat/po/case1/eq1/beta2/mmeig_st.dat' ! remember to modify the name 
!fnind_st =  './dat/po/case1/eq1/beta2/mmeig_ind_st.dat'  



open(feig, file = fneig,  status='old')  !read from old file 
read(feig,*) ! one comment line 


open(feig_st, file = fneig_st, access  ='append', status='replace') 
write(feig_st, *) '# Sorted ', n, 'eigenvalues (real / imag)'  

open(find_st, file = fnind_st, access  ='append', status='replace') 
write(find_st, *) '# Sorted index of mm_eig, original: [1, 2, 3, 4, 5, 6]'

 
iter = 0
iter0 = 1 ! by default 
!iter0 =  47  !22
   ! -- 24 debug, initialize ind_st. iter = 48 bug again!!!
   
!  iter = 48, debug, 6 appears twice        1           2           6           5           3           6

debug = 0

do  

  read(feig, *, end=99)  w0 
  
  iter = iter + 1

!  if(iter <  iter0 )  cycle   
!  if(iter == iter0)  debug = 1
!  if(iter > 55 ) exit  ! to debug 
!  
   
  ! pre-process for the Monodromy matrix 
  wr = w0(1:2*n-1:2)
  wi = w0(2:2*n:2)
  
  if(debug == 1) then 
    print*; print*, '************', iter, '-th point ***************'
!    print*, w0 
    print*, wr 
    print*, wi 
    print*;  ! read* 
  endif 


  if(iter == 1 ) then 

!  if(iter == iter0) then  !! -- debug 
   
!    print*, 'iter == 1'; print*; read* 
    wr0 = wr;    wi0 = wi 
    
!  subroutine eig_sort( wr0, wi0, wr, wi, n, is_one,  wr_st, wi_st)
    call mmeig_sort( wr0, wi0, wr, wi, 6, 1,  wr_st, wi_st, ind_st)
  
  else 
  
    if(debug == 1) then 
      print*, 'Previous point:'
      print*, wr0 
      print*, wi0
      print*;
!      read*
    endif 
  
    call mmeig_sort( wr0, wi0, wr, wi, 6, 0,  wr_st, wi_st, ind_st) 
    
  endif   
  
  wr0 = wr_st; wi0 = wi_st
  
 ! save the ordered index 
  write(find_st, *)   ind_st 
  
  if(debug == 1) then 
    print*, 'After the sorting: ind_st =  ',  ind_st ; print *
    print*, wr_st 
    print*, wi_st
    print*; read*
  endif 
  
!  open(feig_st, file = fneig_st, status = 'old', access = 'append')
  do i = 1, n, 1
    write(feig_st, '(2f16.10)', advance='no')  wr_st(i),  wi_st(i)
!    if(debug == 1)  write(*, '(2f16.10)', advance='no')  wr_st(i), wi_st(i) ! print to screen
  end do
  write(feig_st, *)   
  
!  if(debug == 1)  write(*,*) ! print to screen
  
!  read*
enddo 


99 stop
end program main_mm_eig



















  
