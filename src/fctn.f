	subroutine fctn(x,ind,f,g,z,y,xfi,imax, deriv )
c
C++++++++++++++++++++++++++  
C
C     AUXILIARY SUBROUTINE FOR THE COMPUTATION OF THE FUNCTION F(*) 
c	AND ITS JACOBIAN MATRIX AFTER A HALF OR ONE REVOLUTION

C          INPUT PARAMETERS:
C     X(*)       INITIAL POINT (x,y=0,z,xdot=0,ydot,zdot=0)ON y=0 WITH xdot=0,zdot=0
C     Z(*)       INITIAL POINT IF IND=1
C     IND        INDEX FOR DETERMINATION OF INITIAL POINT, 
C		 if not, initialize the first 6 components to be the state x
C		 and the left 7-42, to be identity matrix
 
C     deriv 	 The subroutine to compute vector field: gr_rtbp or gr_lf
c
c 		OUTPUT PARAMETERS:
C     F(*)       The target variables(f1=0, f2=0 )
C     G(*,*)     JACOBIAN MATRIX OF F(*) (w.r.t the control variables)
C     Z(*)       FINAL POINT UNDER THE POINCARE MAP
C     Y(*)       VECTOR FIELD AT Z(*)
C     XFI        TIME TO GO FROM THE INITIAL POINT TO Z(*)
 
C
C  by Gerard, without any modification to adjust to my routine 
c    Only available for planar lyapunov and halo orbit in dfcr process.

! try to modify this subroutine to make it available for different kind of symmetry,
! for example, wrt y=0 plane,  that means  : id = 2, 
!  initial state: (x,0,z, 0,vy,0)
!  control variables: x, z, vy  	-- 1,3(1,2,3 except id ); 3+id
!  target variables: vx = 0, vz = 0 	-- 3 + 1,3(1,2,3 except id )
!  stop condition: use poincare map to get the crossing with y=0 plane

! tha State Transit Matrix:   6-6 
! if id == 2, ctr1 = 1, ctr2 = 3 
! g (the variational matrix) : 
! | phi(3+ctr1, ctr1), phi(3+ctr1, ctr2), phi(3+ctr1, 3+id) |
! | phi(3+ctr2, ctr1), phi(3+ctr2, ctr2), phi(3+ctr2, 3+id) |

c  subroutine used: poinc, deriv
C***********************************************************************
  	
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION X(6),F(2),G(2,3),Y(42),Z(42)
 
        !print*, 'start fctn, x0, ind=0!', x, ind
        !read(*,*) aaa
        
	IF (IND.EQ.0) GO TO 3
 
	Y = Z
	GO TO 5
        
3	y = 0.d0	
	y(1:6)= x  
	y(7:42:7) = 1.d0
	
  
!     imax determine if the STM is for half or full period

5	call poinc(y,imax, xfi,z, deriv) 
!	call gr_cjrtbp(z(1:6),xmu,cj2 ) 
!	dcj = cj2 - cj 
	
	f(1) = z(4)  !vx
	f(2) = z(6)  !vz
	
	call deriv(0.D0,z,42,y)
 
!        
	f1=-y(4)/y(2)  !axf/vyf
	f2=-y(6)/y(2)  !azf/vyf
 	
! G is the Jacobi matrix of f(f1=vx, f2=vz) with respect to x0,z0,vy0
! the variational eqs are stored by columns,make sure not to make mistakes
! suggestion:
c   before you use the component of the matrix, be careful of the index!
 
	g(1,1) = z(10)+f1*z(8)  ! d vx/d x = phi(4,1)- axf/vfy*phi(2,1)  
	g(1,2) = z(22)+f1*z(20)
	g(1,3) = z(34)+f1*z(32)
	
	g(2,1) = z(12)+f2*z(8)
	g(2,2) = z(24)+f2*z(20)
	g(2,3) = z(36)+f2*z(32)
	
	RETURN
	END

