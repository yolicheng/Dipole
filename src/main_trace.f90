program main_trace
! compute the trace of Momodromy matrix, for stability analysis 


use dp_mod
implicit none

integer, parameter ::  n = 6  ! the dimension of LF problem 
    
! Local Variables
integer :: i,  feig, ftrace, ftrace_st, find_st,  debug, iter0, iter, ind_st(n) 

real(kind=dp) ::   wr(n), wi(n), w0(2*n),  tr_re(n), tr_im(n), &
                   wr0(n), wi0(n), wr_st(n), wi_st(n), tr_str(n), tr_sti(n)  ! sort
 character(len=70) ::  fneig, fntrace, fntrace_st, fnind_st            

feig = 18;  ftrace = 28; ftrace_st = 29; find_st = 30 

! ----  remember to modify the name -----
!fneig    =  './dat/po/case2/eq3/pommegv.dat'           
!fntrace =  './dat/po/case2/eq3/trace2.dat'  
!fntrace_st =  './dat/po/case2/eq3/trace_st2.dat'  
!fnind_st =  './dat/po/case2/eq3/trace_ind_st2.dat'  

!fneig      =  './dat/po/case2/eq1/pommegv.dat'           
!fntrace    =  './dat/po/case2/eq1/trace.dat'  
!fntrace_st =  './dat/po/case2/eq1/trace_st.dat'  
!fnind_st   =  './dat/po/case2/eq1/trace_ind_st.dat'  

!fneig    =  './dat/po/case2/eq2/pommegv.dat'           
!fntrace =  './dat/po/case2/eq2/trace.dat'  
!fntrace_st =  './dat/po/case2/eq2/trace_st.dat'  
!fnind_st =  './dat/po/case2/eq2/trace_ind_st.dat'  

   
!fneig    =  './dat/po/case2/eq2/ifam2_old/pommegv5.dat'           
!fntrace =  './dat/po/case2/eq2/ifam2_old/trace.dat'  
!fntrace_st =  './dat/po/case2/eq2/ifam2_old/trace_st.dat'  
!fnind_st =  './dat/po/case2/eq2/ifam2_old/trace_ind_st.dat'  


!fneig   =  './dat/po/case1/eq3/ifam1/bt2/pommegv.dat'          
!fntrace =  './dat/po/case1/eq3/ifam1/bt2/trace.dat'   
!fntrace_st =  './dat/po/case1/eq3/ifam1/bt2/trace_st.dat'  
!fnind_st =  './dat/po/case1/eq3/ifam1/bt2/trace_ind_st.dat'  


!fneig   =  './dat/po/case1/eq3/ifam2/bt2/pommegv.dat'          
!fntrace =  './dat/po/case1/eq3/ifam2/bt2/trace.dat'   
!fntrace_st =  './dat/po/case1/eq3/ifam2/bt2/trace_st.dat'  
!fnind_st =  './dat/po/case1/eq3/ifam2/bt2/trace_ind_st.dat'  


!fneig    =  './dat/po/case1/eq3/pommegv.dat'           
!fntrace =  './dat/po/case1/eq3/trace.dat'  
!fntrace_st =  './dat/po/case1/eq3/trace_st.dat'  
!fnind_st =  './dat/po/case1/eq3/trace_ind_st.dat'  


!fneig    =  './dat/po/case1/eq1/pommegv.dat'           
!fntrace =  './dat/po/case1/eq1/trace.dat'  
!fntrace_st =  './dat/po/case1/eq1/trace_st.dat'  
!fnind_st =  './dat/po/case1/eq1/trace_ind_st.dat'  

!fneig    =  './dat/po/case1/eq1/beta2/pommegv.dat'          
!fntrace  =  './dat/po/case1/eq1/beta2/trace.dat' ! remember to modify the name 
!fntrace_st =  './dat/po/case1/eq1/beta2/trace_st.dat' ! remember to modify the name 
!fnind_st =  './dat/po/case1/eq1/beta2/trace_ind_st.dat'  

!fneig    =  './dat/po/case1/eq1/beta2m/pommegv.dat'          
!fntrace  =  './dat/po/case1/eq1/beta2m/trace.dat'     ! remember to modify the name 
!fntrace_st =  './dat/po/case1/eq1/beta2m/trace_st.dat' ! remember to modify the name 
!fnind_st =  './dat/po/case1/eq1/beta2m/trace_ind_st.dat'  

!fneig      =  './dat/po/case1/eq2/pommegv.dat'           
!fntrace    =  './dat/po/case1/eq2/trace.dat'  
!fntrace_st =  './dat/po/case1/eq2/trace_st.dat'  
!fnind_st   =  './dat/po/case1/eq2/trace_ind_st.dat'  


!fneig      =  './dat/po/case1/eq2/ifam1/bt2/pommegv.dat'           
!fntrace    =  './dat/po/case1/eq2/ifam1/bt2/trace.dat'  
!fntrace_st =  './dat/po/case1/eq2/ifam1/bt2/trace_st.dat'  
!fnind_st   =  './dat/po/case1/eq2/ifam1/bt2/trace_ind_st.dat'  

!fneig      =  './dat/po/case1/eq2/ifam2/bt2/pommegv.dat'           
!fntrace    =  './dat/po/case1/eq2/ifam2/bt2/trace.dat'  
!fntrace_st =  './dat/po/case1/eq2/ifam2/bt2/trace_st.dat'  
!fnind_st   =  './dat/po/case1/eq2/ifam2/bt2/trace_ind_st.dat'  


fneig      =  './dat/po/case1/eq2/ifam2_bif/pommegv.dat'           
fntrace    =  './dat/po/case1/eq2/ifam2_bif/trace.dat'  
fntrace_st =  './dat/po/case1/eq2/ifam2_bif/trace_st.dat'  
fnind_st   =  './dat/po/case1/eq2/ifam2_bif/trace_ind_st.dat'  



open(feig, file = fneig,  status='old')  !read from old file 
read(feig,*) ! one comment line 

open(ftrace, file = fntrace, access  ='append', status='replace') 
write(ftrace, *) '# Trace:',  n, '(real / imag)'  

open(ftrace_st, file = fntrace_st, access  ='append', status='replace') 
write(ftrace_st, *) '# Sorted ', n, 'trace (real / imag)'  

open(find_st, file = fnind_st, access  ='append', status='replace') 
write(find_st, *) '# Sorted index of trace, original: [1, 2, 3, 4, 5, 6]'


 
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
!!  
   
  ! pre-process for the Monodromy matrix 
  wr = w0(1:2*n-1:2)
  wi = w0(2:2*n:2)
  
  if(debug == 1) then 
    print*; print*, '************', iter, '-th point ***************'
    print*, wr 
    print*, wi 
    print*;  ! read* 
  endif 

  !subroutine trace( wr, wi, n,  tr_re, tr_im)
  call trace( wr, wi, n, tr_re, tr_im)
  
!  print*, tr_re;
!  print*, tr_im; print*; read*
  
  ! sort the trace
  if(iter == 1 ) then 

!  if(iter == iter0) then  !! -- debug 
!  print*, 'iter == 1'; print*; read* 
    wr0 = tr_re;    wi0 = tr_im
    call mmeig_sort( wr0, wi0,  tr_re, tr_im, 6, 1,  tr_str, tr_sti, ind_st)
  
  else 
  
    if(debug == 1) then 
      print*, 'Previous point:'
      print*, wr0;       print*, wi0;       print*;
    endif 
  
    call mmeig_sort( wr0, wi0, tr_re, tr_im, 6, 0,  tr_str, tr_sti, ind_st) 
    
  endif   
  
  wr0 = tr_str; wi0 = tr_sti
  
  ! save the ordered index 
  write(find_st, *)   ind_st
  if(debug == 1)  print*, ind_st 
  
  do i = 1, n, 1
    write(ftrace, '(2f20.10)', advance='no')  tr_re(i),  tr_im(i)
    write(ftrace_st, '(2f20.10)', advance='no')  tr_str(i),  tr_sti(i)
    if(debug == 1)  write(*, '(2f20.10)', advance='no')  tr_str(i),  tr_sti(i)! print to screen
  end do
  write(ftrace, *); write(ftrace_st, *)  
  
  if(debug == 1)  write(*,*) ! print to screen
  
enddo 

99 stop
end program main_trace



















  
