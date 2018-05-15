	subroutine champ(np,g,dir,cham)
c
C++++++++++++++++++++++++++  
C
C     COMPUTATION OF THE VECTOR FIELD GIVING THE CHARACTERISTIC CURVE
C     OF THE FAMILY
C	  INPUT PARAMETERS:
C     NP	     NUMBER OF THE LAST COMPUTED P.O.
C     G(*,*)         JACOBIAN MATRIX OF F(*) (SEE SUBROUTINE TRACA)
C     dir	     SENSE ON THE CHARACTERISTIC CURVE
C	   
c 	  OUTPUT PARAMETERS:
C     cham(*,*)  VECTOR FIELD ON THE CHARACTERISTIC CURVE AT THE LAST 4
C		     POINT, ONLY THE LAST ROW OF cham(*,*) IS COMPUTED
C		     IN THIS SUBROUTINE
C
C***********************************************************************

	implicit real*8(a-h,o-z)

c Newton's method
	integer ::  dir
	dimension g(2,3),cham(4,3), a(3)
	
	!print*, 'champ, g', g
	
	
	a(1) =  g(1,2)*g(2,3) - g(1,3)*g(2,2)
	a(2) = -g(1,1)*g(2,3) + g(1,3)*g(2,1)
	a(3) =  g(1,1)*g(2,2) - g(1,2)*g(2,1)
	
	
	sm= dir*dsqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
	
	if(np.gt.4)go to 1
	
	
	do 2 i=1,3
2	cham(np,i)=a(i)/sm
 	
 !	print*, 'champ! check! np, a,dir, sm, cham'
 !	print*,  np, a, dir, sm, cham(np,:)
 !	read(*,*) aaa 
	  
	if(np.eq.1)return
	
	nnp=np
	go to 5
	
1	do 3 i=1,3
	do 3 j=1,3
3	cham(i,j)=cham(i+1,j)

	do 4 i=1,3
4	cham(4,i)=a(i)/sm
	nnp=4
	
c   this part for the reverse of the sense, didn't understand it!	
5	c=0.d0
	do 6 i=1,3
6	c=c+cham(nnp,i)*cham(nnp-1,i)


	if(c .gt. 0.d0)return
	dir = -dir
	do 7 i=1,3
7	cham(nnp,i)=-cham(nnp,i)

	write(6,*)'reversal in the sense along the characteristic curve'
	return
	end

