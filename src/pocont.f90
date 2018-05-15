! For the continuation: Adams predictor + cham update routines
!  only for continuation along the arc-length parameter

!***********************************************************************
! Adams-bashforth predictor method, Linear multistep method is to computationally 
! solve the ordinary differential equations

!  discard the idea to do the step control here, by detect the angle among the latest 3 points.
!  Instead, do the control according to the angle betwen the vector field of the lastest two points.

! 20160329  Multistep methods refer to several previous points and derivative values. 
!           Adams method uses fixed step size for all the previous points.

!    So be careful with the input cam 

!          Input parameters:
!  np         	number of computed periodic orbits
!  x(*)       	initial conditions of the last p.o.
!  hh         	integration step along arc length of solution locus
!  cam(*,*)   	the vector field of the latest four points 

!        Output parameters:
!  X(*)       new aproximated initial conditions for p.o.

!  Module-base Parameters
!  nctr 	number of control variables

!  Revised by Yu, 20160418
! --------------------------------------------------------------------------


subroutine adams(np,x, hh, cam )
implicit none 

! Input and Output declaration 
integer, intent(inout) :: np 
real(kind=dp), intent(in) ::  hh, cam(4, nctr)
real(kind=dp), intent(inout) :: x(nctr) 
  
! Local Variables 
integer :: i 

! Adams 4-steps perdictor method for each component of X  
do i = 1, nctr
   
  if (np .ge. 4) then 
    x(i) = x(i)+ hh*( 55.d0*cam(4,i) - 59.d0*cam(3,i) + 37.d0*cam(2,i) - 9.d0*cam(1,i) ) / 24.d0

  elseif(np == 1) then 
    x(i) = x(i) + hh*cam(1,i)
    
  elseif(np == 2) then 
    x(i) = x(i) + .5d0*hh*( 3.d0*cam(2,i) - cam(1,i) )
    
  elseif(np == 3) then 
    x(i) = x(i) + hh*( 23.d0*cam(3,i) - 16.d0*cam(2,i) + 5.d0*cam(1,i) ) / 12.d0
      
  endif 
enddo
 
return
end subroutine adams
 
  
!*******************************************************************************
subroutine champ(np,g, cham, dir)
!  Compute the tangent vector field along the characteristic curve by sloving the Kernel(G)

!  The kenel V is a combination of Vi, which is the determinant of submatrix Gi discarding i-th column of G
!  cham contains previous vector fields with equal step size, to be integrated by Adams predictor method

!  g is the Jacobi Matrix for the current np-th p.o. 

!  20160420 -  do not do step control here, once this subroutine is called, this p.o. is 'good '

!	  INPUT PARAMETERS:
!     NP	     	NUMBER OF THE LAST COMPUTED P.O. WITH THE SAME STEP SIZE
!     G(*,*)         	JACOBIAN MATRIX OF F(*) (SEE SUBROUTINE TRACA)
!	   
! 	  OUTPUT PARAMETERS:
!     cham(*,*)  	VECTOR FIELD ON THE CHARACTERISTIC CURVE AT THE LAST 4
!		     	POINT, ONLY THE LAST ROW OF cham(*,*) IS COMPUTED
!		     	IN THIS SUBROUTINE 

! 	Input and Output Varaibles
!     dir	     	SENSE ON THE CHARACTERISTIC CURVE


! 	Module-based Parameters
!     ntar, nctr  
! --------------------------------------------------------------------------
 
implicit none 

integer, intent(in)       :: np 
integer, intent(inout)    :: dir 
 
real(kind=dp), intent(in)    :: g(ntar, nctr) 
real(kind=dp), intent(inout) :: cham(4, nctr)  
  
! Local Variables 
integer :: i, j, nnp 
real(kind=dp) :: a(nctr),  & ! the tangent vector on the characteristic curve of the current p.o. 
 		 gi(ntar, ntar),  det  ! the submatrix of g by deleting the i-th column

  
! A is the vector field of the characteristic curve, which is its integral curve. 
! the characteristic curve is obtained by integrating the vector field by using Adams method
  
! the following is the computation of the vector field A, the components of A are the derivative of the control variable w.r.t the arc-length paramter s.
  
! d ctr(i) / d s, where s = sum( s_i**2) is the arc length
! the vector field along the characteristic curve, which is also Ker(G), such that G * v = 0 

! 1. v = | G1 G2 ... Gn |  
!    where Gi (i=1,..,n) is the determinant of the submatrix of G discarding i-th column
! 	 vi = (-1)^(i+1) * det(Gi) 

! 2. Then normalize v to obtain a unit vector field
! Check 20160407  error is 1.d-15 compared with the direct element-by-element computation, discard the latter approach    
 
! It is too hard to add the step size control -20160412 --- come on!!!
  
! The vector field for the current point, of which the Jacobi Matrix is G
  
 
! construct Gi, which is the submatrix of G discarding i-th column 
do i = 1, nctr

  if ( i== 1) then  ! g1 = g(:, 2,:nctr)
    gi = g(:, 2 : nctr) 
       
  elseif (i < nctr) then  ! discard i-th column
    gi(:, 1 : i-1)    = g(:, 1   : i-1) 
    gi(:, i : nctr-1) = g(:, i+1 : nctr) 
     
  else ! i = nctr,  g_nctr = g(:, 1,:nctr-1)
    gi = g(:, 1 : nctr-1)  
  endif 

!  compute the determinant of gi  
  call detmat(gi, ntar, det)
    
  a(i) = det * (-1)**(i+1)  ! ker ( DG ) =  (-1)**(i+1) * det (DG/i-th column)
    
enddo 
  
! normalize and take the required sense
a = dir * a / dnrm2(nctr, a, 1) 
  
if(debug ==1 )  
print*, 'In champ, check ds, np', ds, np 
read*
 
 
! update the matrix cham, which contain the vector filed of the last four points
    
if(np .gt. 4)  then 
  cham(1:3, :) = cham(2:4, :)  
  cham(4, :)   = a 
      
else 
  cham(np, :) = a  
  if (np == 1) return
      
endif
 
nnp = min0(np, 4)
  
  
! The previous part makes sure that | c = cos( <v1, v2>) | > 0.995,  this guarantees to continue along the same branch
! And the reversal of the sense is necessary when c < 0 to make sure we continue on the same sense(direction) 
       
! Have no idea, how to continue along a new branch .... not necessary at this moment, we are only interested in the starting family
! From Gerard's book: detect whether the new point obtained is a bifurcating point. 
! If the vector product of the new vector field and the previous one is < 0, change the sign of the vector field
 
 
!  return ! without reverse  --- for test to see the difference
  c = dot_product( cham(nnp,:), cham(nnp-1, :) )
if( c .gt. 0.d0 )  return  ! along the same sense
  
! Compare by reversal and non-reversal, Conclusion: -- 20160412 
!    if | cos<v1,v2> |>= 0.995, without the inversal, the continuation will jump backwards and forwards within the same family 
!    if | cos<v1,v2> |<  0.995, with or without reversal, it will jump to another family  ??? really? how?  this only happens when we do not refrain the period, using poincare map the obtain the next mirror configuration   ! --- make clean of this one... why?

! Reveral is needed to continue along the same direction
dir = -dir
 cham(nnp, :) = - cham(nnp, :)

write(*,*)  'Reversal in the Sense along the Characteristic Curve!'
print*, 'cos<v1, v2>', c 
print*, cham(nnp-1, :)
print*, cham(nnp, :)
  
if(debug == 1) read*
  
return
end subroutine champ
