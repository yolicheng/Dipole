! *****************************************************************************
! The prediction along the arc length parameter   S :   Adams predictor + cham update routines

! the good file without dealing with the bifurcation is save in pre_arclen_good.f90


! using Adams multistep method to solve differential equations of the tangent vector along th characteristic curve 
! provided the  several previous points and derivative values. 
    
!	 Input Varaible
!  g 		 Jacobi matrix, used to compute vector field 

! 	Input and Output Varaibles 
!  ipo		index of current p.o.  2
!  dir		dirction along the vector field, +1:increase,-1:decrease
!  ds		step size for the continuation (adaptive control), lower and upper bound [ds_min, ds_max] 
!  nds 		available points with fixed step size 
!  cham		the vector filed of the latest(up to 4) points, used by Adams Multistep predictor
!               to predict the vector field of next p.o.
!  cham5 	the  vector filed of the previous points (up to 5) with fixed step size, used to double the step size  
!  curv 	the points(up to 3) on the characteristic curve, also with the fixed step size, corresponding to cham  
!  yctr 	the prediction for the next (or current) p.o. ( updated, if necessary)


! 	Input Variables
!  incrds 	flag for the step control according to the iterations of Newton Method,  by half: -1, double: 1,  or remain :0. 
! 		if angle_ds /= 0, reset to 0 

!  con_stop 	flag to terminate the continuation process when the required ds < ds_min


!  	Private Module-base Varaibles
! tol		tolerance for vx,vz to be zero (ie,1.d-13 or 1.d-16 in gerard's routine)
! prsc 		presion to terminate the iteration(ie, 1.d-13)
! nctr, ntar, tar, ctr, ntp   
  
!  Local 
!  Finally revised : 20160411 

! todo: 1. chaml is only used for check,  delete this one later after a thorough debug of the routine 
!       . 

! Routine used:  dsdecr, dsincr, champ, adams 
! Finally revised by Yu -- 20160420
! --------------------------------------------------------------------------

subroutine pre_arclen(ipo, ds, nds, dir, g, cham, cham5, curv, yctr, csangle, incrds, con_stop) 

implicit none 
integer, parameter :: dp = kind(1.d0)

integer, intent(out)   :: con_stop
integer, intent(inout) :: ipo, nds, dir, incrds

real(kind=dp), intent(in)     :: g(ntar, nctr)
real(kind=dp), intent(out)    :: yctr(nctr) 
real(kind=dp), intent(inout)  :: ds, csangle, cham(4, nctr), cham5(5, nctr), curv(5, nctr)
  
! Local Variables 
integer :: i 
real(kind=dp) ::  curv2ds(3,nctr), csangle2ds
 

! Use nds_new to add more steps after each change of ds, to make sure the previous points are not too far away for the current one 
! to avoid the possible bigger angle between <v1,v2>

print*, 'Compute the prediction of new p.o. along characteristic curve!'

    
!initialize the flags for control of stepsize to be zeros
con_stop = 0 ! too small stepsize required
 
! cos<v1,v2> = c in characteristic curve, computed for 3 consecutive points 
  
!  question: check by cham or by the finite difference by curv:  for sure calling curv_angle  is better 
 
! --compare the angle computed by either the last 2 available rows in cham or the last 3 available rows in curv 
!          and compare also v1 -- curv(2,:)-curv(1,:),  v2-- curv(3,:)-curv(2,:) for nds >=3 
! --Conclusion : the angle computed by cham is not reliable!!! because curv is the actual points, sometimes there is big difference

! Remember, all the compuation is done approximately, not so accurate.... 
  

if (incrds == -1 ) then     
! the case when we need to decrease the stepsize 
   
!subroutine dsdecr( ipo, ds, nds, yctr, cham, cham5, curv, con_stop )
  call dsdecr( ipo, ds, nds,  yctr, cham, cham5, curv, con_stop)
  
  if(con_stop == 1) then 
    print*, 'Continuation terminated because required arc step for curvature control by halfing is too small!'
    read*;   return
  endif 
  
  call curv_angle(curv(1:3,:), nctr, csangle)  
  
  print*, 'check decrese of ds by half!'
  print*, 'ds=',ds, 'ipo=', ipo, ',	nds=', nds, 'csangle=', csangle ! ck
!  read*  ! ck
  
else 

  ! --------   update cham and cham5 ------------
  nds = nds + 1

  if(debug == 1) then ! ckd-- move to debug mode 
    print* ;  
    print*, 'nds=', nds, '  call champ to update cham'!ck
    print*; print*, 'check cham5, before update'
    do i = 1, 5 
      print*, cham5(i,:)
    enddo  
    read*
  endif   

  !-- how to proceed some more steps to pass the bifurcation???
  ! subroutine champ(np,g, cham, dir)   
  call champ(nds, g, cham, dir) 
  
   
  if(dnrm2(nctr, cham(min0(nds, 4),:), 1) < 1.d-13) then 
    con_stop = 1; return 
  endif   
  
  if(debug == 1) then ! ck-- TODO move to debug mode 
    print*, 'check  champ by champ'
    do i = 1, 4 
     print*, cham(i,:)
    enddo  
    read* 
  endif 
  
  if ( nds .le. 5) then 
    cham5(nds,:) = cham(min0(nds, 4), : )
  else 
    cham5(1:4,:) = cham5(2:5, : ) 
    cham5(5,:) = cham(min0(nds, 4), : )
  endif   

  if(debug == 1) then ! ckd-- move to debug mode 
    print*, 'check cham5, after update'
    do i = 1, 5 
     print*, cham5(i,:)
    enddo  
    read* 
  endif  
  

! before doubling the step, check the possible new 3 points with stepsize of 2*ds
  if ( incrds == 1) then 
!    print*, 'nds >=5 to double ds! nds=',nds ; read* !ck
    
    curv2ds(1:3, :) = curv(1:5:2,:)
    call curv_angle(curv2ds, nctr, csangle2ds)
    
!    print*, 'check angle for possible doubling of ds, csangle=', csangle2ds
!    if(debug == 1)  read* 
    
    if ( dabs(csangle2ds) < csangle_min .and. anglectr == 1)  then 
      incrds = 0 
      
    else 
      print*, 'Double ds!!'; !read*
      call dsincr( ds, nds, cham, cham5, curv )  
      csangle = csangle2ds
    endif 
  
  endif  ! incrds == 1
  
endif  ! incrds = -1 

if (nds > 2 .and. debug == 1) then 
  print*;  print*, 'Check the angle after updating cham!'
  print*, 'cos<v1,v2> =', csangle, 'only make sense for nds>2, nds=', nds
  print*; read* 
endif

    
! Integrate the characteristic curve along the arc-length  parameter  s, the differential system is 
! the last row of cam computed by subroutine champ.... 

! yctr is passed from pofam, or updated by subroutine dsdecr(angle_ds = 1)
 call adams(nds, yctr, ds, cham) 

return
end subroutine pre_arclen 

!******************************** dsincr ****************************************
! This routine deal with the case when the step size need to be doubled with niter < 3, and ds*2 < ds_max 
! if ds*2 > ds_max, keep the current ds without increasing 

! Varaibles need to be updated: cham, cham5, curv 

!  All the reassignments have been  checked carefully. -20160419     

! 	Input and Output Varaibles
!  ds 		the current and updated prediction of the stepsize for the next p.o. 
!  nds 		number of available p.o.s in cham5, with the same fixed step size 
!  cham		the vector filed of the latest(up to 4) points, used by Adams Multistep predictor
!               to predict the vector field of next p.o. 
!  cham5 	all the previous vector filed with fixed step size, used to double the step size 
!  curv 	the points(up to 3) on the characteristic curve, also with the fixed step size, corresponding to cham  

! 	Module-based Parameters
!  nctr, debug, ds_max 
!******************************** dsincr ****************************************
subroutine dsincr( ds, nds, cham, cham5, curv )

implicit none 
integer, parameter :: dp = kind(1.d0)

integer, intent(inout) :: nds 
real(kind=dp), intent(inout) :: ds, cham(4, nctr), cham5(5, nctr), curv(5, nctr)

! Local Parameters 
integer :: i

! Here, the rows of curv and cham are corresponding to each other
! if nds = 3, we keep the available vectors as the last 2 rows in curv, because only the last 2 rows will be used, if necessary (in champ)
! if nds = 1, keep only available the first row in curv
 
! check cham, cham5, curv after the dscrease  ! -- ckd 
if (debug == 1) then  
  print*; print*, 'Before dsincr, check!'
  print*, 'ds=',ds, 'nds=', nds

  print*, 'cham'
  do i = 1, 4
    print*, cham(i, :)
  enddo
  print*          

  print*, 'cham5'
  do i = 1, 5
    print*, cham5(i, :)
  enddo
  print*   

  print*, 'curv'
  do i = 1, 5
    print*, curv(i, :)
  enddo
  print*   
endif 

!   if ds*2 > ds_max, we cannot use the previous np:np-4: 2 points...., not fixed step size 

if ( ds * 2.d0 <= ds_max) then 
          
  print*, 'ds needs to be increased!'
 
  if(debug == 1) read*  !ck
        
  ! we have the previous 5 points to provide us 2 points with stepsize as 2*ds 
  ds = ds * 2.d0
  
  print*, 'nds>=5', nds 
  if(nds < 5 .and. debug == 1)  read*
  
  cham(1:3,:) = cham5(1:5:2, : ) 
  curv(1:3,:) = curv(1:5:2,  : ) 
  
  nds = 3
  cham5(1:3,:) = cham(1:3, : )

         
else 
  
  ! this part will never be called, because the condition for tha clal of dsincr is   if( nds .ge. 4 .and. 2.d0*ds .le. ds_max ) 
  print*, 'Will never reach this part,  Error!!!! ds*2.d0 > ds_max'
  read* 

! if ds == ds_max already, keep everything unchanged 
!  if(dabs(ds-ds_max) .gt. 1.d-14) then
!    cham(1,:) = cham5(nds,: ) 
!    curv(1,:) = curv(3, :)
!    nds = 1

!    cham5(1, :) = cham(1, : )
!  endif 
!  
!  ds = ds_max  
endif 
  
! check cham, cham5, curv after the dscrease  ! -- ckd  
if (debug == 1) then 
  print*; print*, 'After dsincr, check!'
  print*, 'ds=',ds, 'nds=', nds

  print*, 'cham'
  do i = 1, 4
    print*, cham(i, :)
  enddo
  print*          

  print*, 'cham5'
  do i = 1, 5
    print*, cham5(i, :)
  enddo
  print*   

  print*, 'curv'
  do i = 1, 5
    print*, curv(i, :)
  enddo
  print*
  
  print*, 'After the increase!, nds = ', nds, 'ds=', ds
  read* 
          
endif 

return      
end subroutine dsincr



!******************************** dsdecr *******************************************
! This routine is to decrease the stepsize by half when niter > 5, and update the corresponding array.
! Most importantly, compute yctr as the new prediction of initial conditon of ipo-1 -th p.o., so ipo should be decreased by 1   

! it is done before calling the subroutine champ, and because ispo=6 when niter>5, so curv has not been updated. 
! so the arrays cham, cham5 and curv all have not been updated with the current p.o. yet, all the first nds-rows are available
! and corresponding to each other  

! 20160419 -- need to check this routine 

! 	Output Varaibles 
!  yctr 	updated new predictions of the control variables for the ipo-th p.o. 
!  con_stop 	flag to terminate the continuation process when the required ds < ds_min


! 	Input and Output Varaibles
!  ipo  	the index of current p.o. 
!  ds 		the current and updated prediction of the stepsize for the next p.o. 
!  nds 		number of available p.o.s in cham5, with the same fixed step size 
!  cham		the vector filed of the latest(up to 4) points, used by Adams Multistep predictor
!               to predict the vector field of next p.o. 
!  cham5 	all the previous vector filed with fixed step size, used to double the step size 
!  curv 	the points(up to 3) on the characteristic curve, also with the fixed step size, corresponding to cham  

! 	Module-based Parameters
!  nctr, debug, ds_max 

! 	Module-based Parameters
!  nctr, debug, ds_max , issymt, ctr
!******************************** dsdecr *******************************************

subroutine dsdecr( ipo, ds, nds,  yctr, cham, cham5, curv, con_stop)

implicit none 
integer, parameter :: dp = kind(1.d0)

integer, intent(out)   :: con_stop
integer, intent(inout) :: ipo, nds 

real(kind=dp), intent(out)    :: yctr(nctr) 
real(kind=dp), intent(inout)  :: ds, cham(4, nctr), cham5(5, nctr), curv(5, nctr)

! Local Varaibles
integer :: i, ismin
real(kind=dp) :: v1(nctr), v2(nctr), ds2 ! c, c2,
 
con_stop = 0


! check cham, cham5, curv before dsdecr
if(debug ==  1) then 
  print*, 'Befor dsdecr, check!'
  print*, 'ds=',ds, 'nds=', nds
  read*

  print*, 'cham'
  do i = 1, 4
    print*, cham(i, :)
  enddo
  print*          

  print*, 'cham5'
  do i = 1, 5
    print*, cham5(i, :)
  enddo
  print*   

  print*, 'curv'
  do i = 1, 5
    print*, curv(i, :)
  enddo
  print* ;  read* 
endif 

!  ----- decrease the stepsize by half, this part has not been checked, to provide good guess for subroutine champ ------ 
  
! since all the computations are not very accurate, in this case we do the linear interpolation to add an additional 'middle' point, 
! but keep the number of available p.o.s, without doing the refinement, treat the new 'middile' point as a good one. 

ismin = 0  !flag to show if the stepsize is set to the minimum 

if (ds/2.d0 .le. ds_min .and. dabs( ds - ds_min) > 1.d-14 ) then 
  con_stop = 1 
  return 
    
elseif ( ds / 2.d0 > ds_min) then 
  ds = ds / 2.d0
  
else 
  ds = ds_min 
  ismin = 1
endif


if(debug == 1) then 
  print*, 'ds=', ds, ',	ismin=', ismin, ',	ds_min=', ds_min 
  if(ds < ds_min) read*
endif 
 
!  if ismin == 1 .or ipo ==2,  we only have 1 available good point  
!  else, we have at least 2 'good' p.o.s, we can add a middle point(linear interpolation) bewteen the last 2 good curvs
      
if (ipo == 2 .or. ismin == 1) then 
  print*, 'Only 1 available previous point!'      
  print*, 'check if ipo=2, nds=1', ipo, nds ;  ! read*
  
!  if(ipo == 2 .and. nds /= 1) ! read* 
  
! take the initial guess of the state and period(if needed) of the last available orbit 
! use the ds_min as the new smaller stepsize, so there is no available pervious points with step size as ds_min, nds = 1 

  cham(1,:)  = cham(min0(nds,4),:)
  cham5(1,:) = cham(1,:)
  curv(1,:) = curv(min0(nds,5),: )
                
  nds = 1 
  
else !  ipo > 2  
   
  if (debug == 1) then           
    print*, 'Devide ds by 2, ds=ds/2', ds    
    read*
  endif    
        
  v1 = cham( min0(nds,4) - 1, : )
  v2 = cham( min0(nds,4), : )
   
  if (debug == 1) then              
    print*, 'check the previous 2 vector field'
    print*, v1
    print*, v2
  endif 
              
  cham(1,:) =  v1
  cham(3,:) =  v2
                
  cham(2,:) = ( v1 + v2 ) / 2.d0
  cham(2,:) = cham(2,:) / dnrm2(nctr, cham(2,:), 1) ! normalize
    
  ! update curv to be the last 5 good 'p.o.'
  curv(1, :) = curv(min0(nds,5)-1, :)
  curv(3, :) = curv(min0(nds,5), :)                      
  curv(2,:) = ( curv(1,:) + curv(3,:) ) / 2.d0! replace the first row of curve by the middle point
            
  nds = 3
  cham5(1:3, : ) = cham(1:3, :) ! do not forget this one 
            
  ds = ds / 2.d0  ! decrease the step size by half   
    
  ! check the ds for the previous two points  
  ds2 = dnrm2 ( nctr,  ( curv(3,:)-curv(2,:) ), 1) 
  print*, 'ds for the middile point, and after the decrease:', ds2, ds 
!  read*       
            
endif ! ds/2 < ds_min or not
             
! update the prediction  
yctr = curv( min0(nds,5), : ) ! of dimension nctr, for gr_po case, 7(tp included)


if (debug == 1) then  !--ckd
  print*, 'After dsdecr, check!'
  print*, 'ds=',ds, 'nds=', nds

  print*, 'cham'
  do i = 1, 4
    print*, cham(i, :)
  enddo
  print*          

  print*, 'cham5'
  do i = 1, 5
    print*, cham5(i, :)
  enddo
  print*   

  print*, 'curv'
  do i = 1, 5
    print*, curv(i, :)
  enddo
  print*
  
  print*, 'After the decrease!, nds = ', nds, 'ds=', ds
  read* 
endif 

ipo = ipo -1 
return          
end subroutine dsdecr        
        

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
subroutine champ(np, g, cham, dir)
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
integer ::    i, nnp 
real(kind=dp) ::  a(nctr),  a2m, am_max, & ! the tangent vector on the characteristic curve of the current p.o. 
                  gi(ntar, ntar),  det, c  ! the submatrix of g by deleting the i-th column

! -- for bifurcation
real(kind=dp) ::  subG(ntar-1, ntar-1), c1(ntar-1), c2(ntar-1), & 
                  v1(ntar-1), v2(ntar-1), b1, b2, alpha 
                  
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
 
! ************* deal with the bifurcation *******************
! instead of producing meanlingless prediction when DG is not full-rank, that meas, detmat is equal to zero 

! we follow the strategy in Simo's book: 'On the Analytical and Numerical Approximation of Invariant Manifolds', P5 

!If we approach the bifurcation, the behaviour o fthe dominant A_j is of the type 0(d^k) with k odd, d being the distance at the bifurcation point. THe one should reverse the +- sign to continue along the new branch. In this case, we should do a few more steps of the continuation method before refining by using the modified Newton's method to prevent from bad conditioning of the matrix Df * Df^{-1}.
! This could be one approach. 

! I prefer the one in the Notes I took in his class, compute the two bases for the kernel in the degenerate case 
! Since the system is always indeterminant, and always has one more unknown than the equations. 
! The dimension is supposed to be n X n+1 
! Take the first (n-1) X (n-1) submatrix as subA, with the last two columns as c1, b1
! the kernel can be expressed as [v(n-1), b1, b2] ^ T, with b1,b2 as arbitrary numbers 
! so we have 
! subA*v  + c1*b1 + c2*b2 = 0 
! v = - subA^-1 * c1 * b1  -  subA^-1 * c2 * b2  
!  thus,  v1 = subA^-1 * c1 and  v2 = subA^-1 * c2 serve as the two bases of our kernel, any combinations will be the solution 

! And since v1 and v2 are to be normalized, we assume b1,b2 \in [0, 1] and rewrite v to be 
! v = cos(a)*v1 + sin(a)*v2, in such a way that we still get  a united vector will the direction only dependent on one parameter alpha 


! check for the input: 
print*, 'DG to compute the kernel:'
do i = 1, ntar, 1
  print*, g(i, :)
end do
print*; read*


! It is too hard to add the step size control -20160412 --- come on!!!
  
! The vector field for the current point, of which the Jacobi Matrix is G
  
! construct Gi, which is the submatrix of G discarding i-th column 
a2m = 0.d0
do i = 1, nctr

  if ( i== 1) then  ! g1 = g(:, 2:nctr)
    gi = g(:, 2 : nctr) 
       
  elseif (i < nctr) then  ! discard i-th column
    gi(:, 1 : i-1)    = g(:, 1   : i-1) 
    gi(:, i : nctr-1) = g(:, i+1 : nctr) 
     
  else ! i = nctr,  g_nctr = g(:, 1,:nctr-1)
    gi = g(:, 1 : nctr-1)  
  endif 

!  compute the determinant of the square matrix gi  
  call detmat(gi, ntar, det)
  
  if(debug == 1) then 
    print*, 'det=', det;  read*
  endif 
  
  a(i) = det * (-1)**(i+1)  
  ! ker ( DG ) =  (-1)**(i+1) * det (DG/i-th column)
!  if( dabs(a(i)) < 1.d-16) then 
!    print*, 'a(i) too small', a(i); print*; read*
!    a(i) = 0.d0 ! treat as zero if a(i) is too small 
!  endif   
enddo 


a2m = dnrm2(nctr, a, 1)
print*, '|a| = ', a2m, 'a =', a
print*; read*
  
if(a2m < 1.d-8) then 
  ! ---- this is when we encounter a bifurcation---- 
  print*, 'Small determinant of G, possible bifurcation. |G| < 1.d-8, a=', a2m
  print*, 'a = ', a
  
  print*, 'check if the system is indeterminant, nctr = ntar+1?', nctr, ntar 
  print*; read*
  
  ! --- compute the two bases, and  introduce a new angle alpha to control which branch to follow. 
  !  v1 = subA^-1 * c1 and  v2 = subA^-1 * c2  two bases 
  !  v = cos(alpha)*v1 + sin(alpha)*v2, a \in [0,2pi]
  !  subA is the first (n-1)*(n-1) submatrix 
  subG = g(1:ntar-1, 1:ntar-1)
  c1 = g(1:ntar-1, ntar); c2 = g(1:ntar-1, ntar+1); 
  v1 = matmul(subG, c1);  v2 = matmul(subG, c2);            ! two bases
  
  v1 = v1/dnrm2(ntar-1,v1,1);   v2 = v2/dnrm2(ntar-1,v2,1)  ! normalize 
  
  print*, 'Input the value of alpha: [0, 2pi]:  '
  read*, alpha 
  
  b1 = dcos(alpha); b2 = dsin(alpha)
  a(1:ntar-1) = b1*v1 + b2*v2
  a(ntar:ntar+1) = (/b1, b2/)

 
  ! check if a is unit one 
  print*, 'check for the degenerate case:'
  print*, 'c1, c2, ', c1, c2 
  print*, 'b1, b2', b1, b2 
  print*, 'subG', subG 
  print*, 'v1, v2', v1, v2 
  print*, 'final kernel: ', a 
  print*, '|a| = 1? ', dnrm2(a)
  print*; read*
  
else 

   !  good prediction    
   a = dble(dir) * a / a2m
   
endif  



! update the matrix cham, which contain the vector filed of the last four points
if(np .gt. 4)  then 
  cham(1:3, :) = cham(2:4, :)  
  cham(4, :)   = a 
      
else 
  cham(np, :) = a  
  if (np == 1) return
      
endif
 
nnp = min0(np, 4)

! The bifurcation happens when the sum of a(i)^2 = 0, a turning point 
if( a2m < 1.d-20) then 
  print*, 'Bifurcation happen!, Sigma a(i)^2 = ', a2m 
  cham(nnp, :) = 0.d0
  print*; read*;  return
endif 
  
  
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
print*
  
if(debug == 1) read*
  
return
end subroutine champ



