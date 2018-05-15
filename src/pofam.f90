subroutine pofam(yi,np,imax,e2,dir,ds,id,ntp, deriv, ynew)
! compute np new p.o.s with one initial po guess
! To make sure the initial state is periodic, refine it without check

! 	Input Varialbles
!   yi(6)	initial condition(x,y=0,z, vx=0, vy, vz=0)
!   np		number of new p.o.  if np=0, just do differential correction
!   imax	number of intersections with poincare section to get p.o.
!   e2		error control for difcor: vx,vz to be zero(1.d-9)
!   dir		dirction along the vector field, +1:increase,-1:decrease
!   ds		displacement for new PO
!   id    	flag of orbit, 1: planar lyap, 2: halo, 3: vertical lyap
!   ntp 	multiplier to get a full period     

! 	Output Variables
!   ynew	the initial condition for new p.o.s, (np+1)*7
!     		save the initial p.o in the first row

!  Routine used: dfcr, fctn, champ, adams

 
! 20150408 --- the vertical lyapunov checked  
!     imax is set to 1 to get the first intersection, pay attention to
!     the computation of g for new PO

! Finally revised : 20151224  
!***********************************************************************
implicit none
  
  
!  both the state and the variational equations are desired, so n=42
integer, parameter:: n=42, dp = kind(1.d0) 
integer :: i, fout 

real(kind=dp) :: yf(n), f(2),g(2,3), cham(4,3), tp, &
   		 cj, cj2, dcj, aaa
  
integer, intent(in) 	   :: np, imax, dir,id, ntp
real(kind=dp), intent(in)  :: yi(6), e2, ds
real(kind=dp), intent(out) :: ynew(np,7)
  
external :: deriv
  
ynew = 0.

print*, 'start pofam!'  
fout = 18 
open(fout,file='./dat/PoInSt.dat',access ='append',status='replace')  

do i = 1, np
  !subroutine dfcr(y0,imax,e2,id, xfi,yf,g) 
  call dfcr(yi,imax,e2,id, tp,yf,g, deriv) ! this is the main problem
   
!  print*, 'finish dfcr!'
  
!  print*, 'g', g
    
  tp = ntp * tp ! whole peroid 
     
  ynew(i,:) = (/tp, yi/)
  
  ! write to file, PoInSt.dat
  write(fout,'(7d20.10)') tp, yi ! fout - PoInSt.dat 
!  write(*,'(7d20.10)') tp, yi ! fout - PoInSt.dat 
!  read(*,*) aaa
  
   
! using Adams predictor to get the new initial condition
  if(np.eq.0) return
   
  call champ(i,g, dir, cham) 
  
!  print*, 'champ!, dir', dir
!  print*, 'yi, cham', yi, cham(i,:)
!  print*
!  read(*,*) aaa
  
!  print*, 'before adams, i, ds, cham, yi(1:5:2)', i, ds, cham, yi(1:5:2)
!  read(*,*) aaa 
  
  call adams(i,yi(1:5:2), ds,cham) 
  
!  print*, 'numerical continuation, new state!'
!  print*, yi
!  print*
!  read(*,*) aaa
  
enddo

 close(fout)

return
end subroutine pofam
	 
	 
