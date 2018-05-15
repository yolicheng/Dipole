! to test if there is something wrong with the computation of the variational matrix 
! try ieq = 2, beta = 1, ifam = 3    --- tested! fine 
!  because 2924 seems too big, so I have checked this by finite difference, everything is fine here.

! --------------------  from numercial continuation  -------------------- 
!         33 -th P.O. TP:    5.6436944269033864     
!  5.64369443  0.50434094  0.71676431 -0.01184629 -0.00138981 -0.00485037 -0.00004893  2.28923670
! refined initial state, tp, ynew   5.6436944269033864       0.50434093784184986       0.71676430921462631       -1.1846288791455030E-002  -1.3898124069350972E-003  -4.8503687009210589E-003  -4.8933142540624947E-005
! Eigenvalues
! 2924.8498732306    -0.8775666537      0.4794546572   -0.8775666537     -0.4794546572    0.0003418979     1.0000000000      0.0000001126    1.0000000000     -0.0000001126
! Eigenvectors
! -0.41059322  -0.17482640  0.08843288   -0.17482640 -0.08843288    0.28832085  -0.08504032  0.00000056   -0.08504032 -0.00000056  
!  0.00725561  -0.00081914 -0.22955077   -0.00081914  0.22955077    0.41580920  -0.29678602 -0.00003684   -0.29678602  0.00003684  
!  0.40754153  -0.07701194 -0.15607194   -0.07701194  0.15607194    0.29170560  -0.00299413  0.00012738   -0.00299413 -0.00012738  
! -0.57265392  -0.23995572 -0.47023360   -0.23995572  0.47023360   -0.40260257   0.00320106  0.00001548    0.00320106 -0.00001548  
!  0.01216317   0.61456877  0.00000000    0.61456877 -0.00000000   -0.58172842  -0.26546549  0.00005112   -0.26546549 -0.00005112  
!  0.58068064   0.42036473 -0.21338268    0.42036473  0.21338268   -0.39797595   0.91334278  0.00000000    0.91334278 -0.00000000 
!-------------------------------------------------------------------------------

use lf_mod

implicit none 
!  integer, parameter :: dp=kind(1.d0), n = 42
integer, parameter ::  neq = 42
    
integer ::   ftag, tdir, i, cs, ieq, idck 
real(kind=dp) :: y0(neq),t0, tf
real(kind=dp) :: y(neq) ! As the final state
  

! Local Variables
real(kind=dp) :: r(13,neq), b(neq), f(neq), t, h, hmin, hmax, e,  & ! rk78 
		 cj , po0(6), phi(6,6),  & 
		 beta, &
		 dy(6), phi1(6), yck(neq), vf(neq), &
		 dxck ! to check the finite difference of the state transition matrix

!  external :: gr_lf, gr_cjlf
 cs = 1   
ieq = 2  
beta = 1.d0; 

call init_lf_mod
call init_lf(beta, cs, ieq) 

! 0.56123418697415695        0.0000000000000000       0.56122987432829607        0.0000000000000000        9.9943356318814896E-005   0.0000000000000000 

! 0.42413555  0.56123418  0.00000000  0.56122988  0.00000000  0.00009994  0.00000000  1.88988156 ! ifam = 1
! 0.4550459126096073D+02  0.5612316936062401D+00  0.0000000000000000D+00  0.5612317343582736D+00  0.0000000000000000D+00  0.2178198390586603D-06  0.0000000000000000D+00 ! ifam = 3

tdir = 1
tf = 0.4550459126096073D+02/ 1  ! 1 period
!  tf = 1.d-0

!po0 = (/0.56123418d0,  0.d0,   0.56122988d0,  0.d0,  0.00009994d0, 0.d0/)  ! the first p.o. 
!po0 = (/0.56123169d0,  0.d0,   0.5612317343582736d0,  0.d0,  0.21781984d-6, 0.d0/)  ! the 3rd p.o. 

po0 = (/0.5612316936062401d0,  0.d0,  0.5612317343582736d0,  0.d0,  0.217819839058660d-6, 0.d0/)  ! the 3rd p.o. 
 
!po0 = (/0.56123102415468651d0, 0.d0, 0.56123102415468651d0,  0.d0,  0.d0, 0.d0/)  ! eq2
!po0 = (/0.d0, 0.d0, 1.d0,  0.d0,  0.d0, 0.d0/)  ! eq1


tf = 5.6436944269033864d0
po0 = (/ 0.50434093784184986d0, 0.71676430921462631d0, -1.1846288791455030d-2, -1.3898124069350972d-3, & 
        -4.8503687009210589d-3, -4.8933142540624947d-5/)
 
ftag = 22 
open(ftag,file='./dat/testphi.dat', access ='append',status='replace')
 
! Integrate the periodic orbits by gr_rk78, write the x, y in txt file with ftag
t =  0.d0
t0 = 0.d0
h =  1.d-3! maybe small value is better, for a p.o. with small period, the orbit could be quite corse 

! specify the error control 
hmin = 1.d-10
hmax = 1.d0
e    = 1.d-14

y = 0.d0
y(1:6) = po0 
y(7:42:7) = 1.d0  
  
 
! Integrate the orbit in interval [t0, tf]
do while( dabs(t+h-t0) .lt. dabs(tf-t0) )
  
!----  check the vector field ---------
!  call  gr_lf(t, y, neq, vf)
!  
!  print*, 'check the variatioal matrix'
!  phi =  reshape(vf(7:42), (/6,6/))
!  
!  do i = 1, 6
!    write(*,'(6f12.8)') phi(i,:) 
!  enddo
!  print*
! --------------------------------------


  call gr_cjlf(y(1:6), cj)
    
  write(ftag,'(8e20.10)')  t, y(1:6), cj
!  write(*,'(8e20.10)')  t, y(1:6), cj!ck
  call gr_rk78(t,y, neq ,h,hmin,hmax,e,r,b,f, gr_lf )
  
    
  ! check the state transition matrixf  - for every step
!  print*; print*, 'Phi- t=',  t
!  phi = reshape(y(7:42), (/6,6/)) ! better than equivalence declaration.... 
!  do i = 1, 6
!    write(*,'(6e20.10)')  phi(i,:)
!  enddo 
!    read*
  
    
enddo
     
!if we want a periodic orbit with exact 1 period, a control of tf must be made
if(dabs(t-tf) .gt. 1.d-9) then
  
  h =  tdir*tf - t ! for the stable orbit, the final time should be -tf
! to make sure the next step can be executed within allowable interval [hmin hmax]
!    hmin should be dabs(h)
  call gr_rk78(t,y, neq,h, dabs(h), hmax, e,r,b,f, gr_lf)
    
 endif 

! save the last point  
call gr_cjlf(y(1:6), cj)      
write(ftag,'(8e20.10)') t, y(1:6), cj
  
write(ftag,*)  ! better to save as a block than index
write(ftag,*)  ! better to save as a block than index  

  
! **********************  test with finite difference   *********************** 
yck = 0.d0
yck(1:6) = po0 
yck(7:42:7) = 1.d0


! choose which column to compare 
idck = 6 ! 1,2,3, 4,5,6- all checked! fine here
dxck = 1.d-6 

yck(idck) = yck(idck) + dxck
  
  
t = 0.d0
h = 1.d-3
do while( dabs(t+h-t0) .lt. dabs(tf-t0) )
   write(ftag,'(7e20.10)')  t, yck(1:6) 
!  write(*,'(7e20.10)')  t, yck(1:6) 
  call gr_rk78(t,yck, neq ,h,hmin,hmax,e,r,b,f, gr_lf )
enddo
  
!if you want the strict periodic orbit, a control of tf must be made
if(dabs(t-tf) .gt. 1.d-9) then
  h =  tdir*tf - t 
  call gr_rk78(t,yck, neq,h, dabs(h), hmax, e,r,b,f, gr_lf)
endif 
write(ftag,'(7e20.10)')  t, yck(1:6) 
  
! by integration 
phi = reshape(y(7:42), (/6,6/)) 
print*
print*, 'By integration, phi'
do i = 1, 6
  write(*,'(6f20.10)') phi(i,:)
enddo 
  
print*
  
! the difference between yck and y
dy = yck(1:6) - y(1:6)
phi1 = dy / dxck
print*, 'By finite difference, corresponding to the', idck,'-th column of phi'
print*, phi1

 close(ftag)  
read*

! after step-by-step check, there seems no problem, although the result still displays possible mistake

! possible mistake of derivative of vy_f w.r.t vx_0  .... -- to check - 20160327
! By integration, phi
!    0.1030039770E+01   -0.6828005501E-02   -0.2832988387E-01    0.2389010455E-01   -0.1274392607E-02    0.6707329287E-01
!    0.4794943215E-01    0.9998670081E+00    0.2010907635E-01    0.1274394302E-02   -0.1533595282E-02    0.3515949127E-02
!    0.8312620601E-01    0.2239813152E-02    0.9222926286E+00    0.6707335859E-01   -0.3515960995E-02    0.1840277552E+00
!    0.2811599236E+00   -0.1282478248E-02   -0.2640218413E+00    0.9943570977E+00   -0.5009005103E-02   -0.1532044038E-01
!    0.1483823487E-01   -0.1528413739E-02   -0.1394293553E-01    0.5001684431E-02    0.9993444533E+00    0.3549161455E-01
!    0.7785330840E+00   -0.3531732445E-02   -0.7282510639E+00   -0.1558511048E-01   -0.3544529919E-01    0.9574528795E+00

! By finite difference, corresponding to the           4 -th column of phi
!   2.4136843745914401E-002   1.2765367354279412E-003   6.7069962728094623E-002  0.99436804490810937        1.1945330471969200E-002  -1.5625594147282467E-002


!**  d vx_f / d vy_0 
!By integration, phi
!    0.1030039770E+01   -0.6828005501E-02   -0.2832988387E-01    0.2389010455E-01   -0.1274392607E-02    0.6707329287E-01
!    0.4794943215E-01    0.9998670081E+00    0.2010907635E-01    0.1274394302E-02   -0.1533595282E-02    0.3515949127E-02
!    0.8312620601E-01    0.2239813152E-02    0.9222926286E+00    0.6707335859E-01   -0.3515960995E-02    0.1840277552E+00
!    0.2811599236E+00   -0.1282478248E-02   -0.2640218413E+00    0.9943570977E+00   -0.5009005103E-02   -0.1532044038E-01
!    0.1483823487E-01   -0.1528413739E-02   -0.1394293553E-01    0.5001684431E-02    0.9993444533E+00    0.3549161455E-01
!    0.7785330840E+00   -0.3531732445E-02   -0.7282510639E+00   -0.1558511048E-01   -0.3544529919E-01    0.9574528795E+00

! By finite difference, corresponding to the           5 -th column of phi
!  -1.2764831414102673E-003  -1.2854470133951345E-003  -3.5244103191445220E-003  -1.1960555569909038E-002  0.99937005665714673       -3.2959943968505177E-002

  
end  
