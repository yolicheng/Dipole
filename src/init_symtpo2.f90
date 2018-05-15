! Initialization for index of control and target variables.  symtpo + asymtpo 

! ****************************************************************************** 
! For symmetric p.o.:   idcor = (1,2,3), poincare-map, free-time, fixe-time
!  Specify the index of non-zero components for different types of symmetry for 
!  initial and final configurations :   idic,  idfc 
!     For simply-symmetric p.o.s, they are the same type, idic=idfc
!     for doubly-symmetric p.o.s,  they are not the same 

!  Because for the doubly-symmetry case, two options, it is possible to refine the 
!      1.  second crossing through the poincare section, which is of the same kind of mirror configuration as the initial point 
!      2.  the first crossing, with different type of mirror configuration 


!  ------------------- simply-symmetry --  T/2 --------------------------------- 
! Three kinds of planar symmetry, w.r.t.            I.C. and F.C.         non-zero components(control)      zero components(target) 
!    1. x=0 plane,  i.e., y-z plane        P1:  (0, y, z,    vx, 0,  0 )           (2, 3,   4 )                 (1,  5, 6 )
!    2. y=0 plane,  i.e., x-z plane        P2:  (x, 0, z,     0, vy, 0 )           (1, 3,   5 )                 (2,  4, 6 )           
!    3. z=0 plane,  i.e., x-y plane        P3:  (x, y, 0,     0, 0, vz )           (1, 2,   6 )                 (3,  4, 5 )

!  and also three kinds of axis symmetry, w.r.t 
!    4. x-axis                             A1: (x, 0, 0,    0, vy, vz)            (1,    5, 6)                  (2,3,   4)
!    5. y-axis                             A2: (0, y, 0,   vx,  0, vz)            (2,    4, 6)                  (1,3,   5) 
!    6. z-axis                             A3: (0, 0, z,   vx, vy,  0)            (3,    4, 5)                  (1,2,   6) 
 
! for simply-symmetric p.o., the problem of determination of p.o. reduces to find two mirror configuration of the same type at two different  epoches. 
! with the point on the symmetric element( in the axis for A type and on the plane for P type), and the velocity perpendicular to that element. 

! The control variables are the  non-zero components in the initial mirror configuration, with the left componnents are zeros 
! while the target variables are the zero components in the final mirror configuration 


! --- TODO: for the moment, we skip this case... 2017-02-16 10:13:55 
!  ------------------- doubly-symmetry --  T/4 --------------------------------- 
!  doubly symmetry, be careful with the index, which one to put first.... 
! -------- plane symmetry + axis symmetry------------ P for I.C. and A for F.C.
!    7. y-z plane   +   y-axis         P1 --> A2      P1:control  (2, 3,   4 )              A2: target   (1,3,   5)         
!    8. y-z plane   +   z-axis         P1 --> A3      P1:control  (2, 3,   4 )              A3: target   (1,2,   6)          

!    9. x-z plane   +   x-axis         P2 --> A1      P2:control  (1, 3,   5 )              A1: target   (2,3,   4)        
!   10. x-z plane   +   z-axis         P2 --> A3      P2:control  (1, 3,   5 )              A3: target   (1,2,   6)          


!   11. x-y plane   +   x-axis         P3 --> A1      P3:control  (1, 2,   6 )              A1: target   (2,3,   4)        
!   12. x-y plane   +   y-axis         P3 --> A2      P3:control  (1, 2,   6 )              A2: target   (1,3,   5)          

! How to save the index? By observing the value of  control variables and target variables for different cases, it is convevinet to introduce two 6-by-6 arrays
!where 1-3 rows refers to three P-type mirror configurations, 
!! while 4-6 rows refer  to  three A-type mirror configurations
  
!ctrl = [ 2,3, 4 ;  1,3, 5; 1,2, 6;   1, 5,6 ;  2, 4,6; 3, 4,5 ]
!tarl = [ 1, 5,6 ;  2, 4,6; 3, 4,5;   2,3, 4 ;  1,3, 5; 1,2, 6;]

! and we have the first three lines for planar-symmetry, and the last three ones for axis-symmetry

! for each case,  we pick the control varaibles from non-zero components 
! and the target variables from the  zero components 

! We introduce two flag for the type of mirror configuration of the I.C.(idic) and F.C.(idfc), the range of value of both are [1:6] 
!  i=1-6  ctrl(idic, :)         tarl(idfc, :) 

!       Input Varaibles 
!  idic          index of the type of initial configurarion 
!  idfc          index of the type of final configuration
!  idcor         index of the approach to compute  symmetric p.o. 
!                1 : poincare map approach to speicify the next epoch 
!                2 : free-time shooting
!                3 : fixed-time shooting? --TODO:  test
                 
                
! ****************************************************************************** 
  subroutine init_symtpo2(idic, idfc, idcor0, ind0)
  
  integer, intent(in) :: idic, idfc, idcor0, ind0   ! index of symmetry, used to assign the public values accordingly

! Local Varaibles
  integer ::  ctrl(6, 3), tarl(6,3) ! index of control and target varables for all six possible mirror configuration 
  integer ::  allocatestatus, i, j

! for symmetric p.o.   
  idcor = idcor0
  issymt = 1
  
  ! by default, period is not a control variale
  tfree = 0 
  
!---- initialization of the index array for control and  target varaibles----
  ctrl(1, :) = (/2,3, 4/)
  ctrl(2, :) = (/1,3, 5/)
  ctrl(3, :) = (/1,2, 6/)
  
  ctrl(4, :) = (/1, 5,6/)
  ctrl(5, :) = (/2, 4,6/)
  ctrl(6, :) = (/3, 4,5/)
  
  tarl(1:3, :) = ctrl(4:6, :)
  tarl(4:6, :) = ctrl(1:3, :)
!--------------------------------------------------------------------------------  
  
! ! --- allocate memory for ctr and tar for symmetric p.o.  
  if( idcor == 1 ) then 
  ! poincare map approach, by default, we search for the first crossing through the poincare section
    imax = 1
    ind = ind0 
    p0  = 0.d0 ! by default....
    
    ntar = 2
    nctr = 3
  elseif( idcor == 2 ) then 
  !  free time shooting
    tfree = 1
    nctr = 4 
    ntar = 3
  
  elseif( idcor == 3 ) then 
  !  fixed time shooting
    nctr = 3
    ntar = 3  
  end if
  

! allocate the memory for tar and ctr   
  allocate ( tar(ntar), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"

  allocate ( ctr(3), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"

! control variables are the same for different approach 
  ctr = ctrl(idic, :)
  
  
! target variables are the same for different approach  
  if( idcor /= 1 ) then 
    tar = tarl(idfc, :)
    
  else 
  
  ! poincare map appraoch, target variables should be tar(idfc,:) deleting the component that is equal to ind...
    j = 0
    do i= 1, 3, 1
     
      
      if(debug == 1) then 
        print*, 'ind=', ind, 'tarl(idfc, i)=', tarl(idfc, i)
        read*
      endif    
      
      if( tarl(idfc, i) /= ind ) then 
        
        j = j + 1
        tar(j) = tarl(idfc, i)
      endif 
      
    end do
    
  endif 


! --- the period ----  
!if( (idic-3) * (idfc-3) .ge. 0)  ! both <=3 or both > 3, no good to use this conditional statement, idic =3, idfc>3, cannot specify
!  if( idic.le. 3 .and. idfc .le. 3 .or.  idic .gt. 3 .and. idfc .gt. 3  ) then
!  ! -- planar  symmetry------           -- axis  symmetry------ 
!    ntp = 2 

!  else 
!  ! -- doubly  symmetry------ 
!    ntp = 4  
!   
!  end if   

! refine doubly symmetry, if idic /= idfc, we can see it is doubly symmetry,
! this is how we compute the near planar one that is both symmetric w.r.t y-z plane and x-z plane, where idic = 1, idfc =2  
  if( idic == idfc ) then
    ntp = 2
  else
    ntp = 4
  end if
  
  
  if(idcor == 1)   ntp = ntp / imax    

!  check the assignment 
  print*, 'idic=', idic, ' ,idfc=', idfc, ' ,imax=', imax, ' ,ntp=', ntp
  print*, 'ctr=', ctr, ' ,tar=', tar 
  read*
      
  
  if (debug == 1) then          
    print*, 'check init_po: ind, p0, index of control and target var'
    print*,  ind, p0, ',', ctr, ',', tar 
    read(*,*)
  endif   
  
  end subroutine init_symtpo2
  
  
  
  ! ****************************************************************************** 
! Idea: Fix the value a component of the initial state as the poincare section x = x0, 
!       and generally an orbit will pass any orbit on it twice, so we take the second intersection with the x = x0 
!       with the velocity in x component to be in the same direction with the initial one. And then use differential corrector method
!       to make sure the final point overlaps with the initial one . 

! A potential bug: what if the initial point with x=x0 is exactly the point that only has one intersection with the x=x0 plane.... 
!                  It rarely happens, but there is possibility.

! Solution: Take the most secure approach, use gr_po instead....

!  The energy h0 and x0 is fixed, look for the first return to x=x0 plane, we only need to ask four final component in state to be as the 
!  initial one, the fifth will be satisfied simultaneously.  And we have 5 initial components to modify(except x0)

! constraint funcion:   ds^2 = ( (yf-y0)^2 + (zf-y0)^2 + (vxf-vx0)^2 + (vyf-vy0)^2 ) to be minimum

! asym   		contol  		target 
! 1   	 (y0,z0, vx0,vy0,vz0) - (2,3, 4,5,6)   (yf, zf, vxf, vyf) 

!  To ensure that enough memory is available to allocate space for your array, make use of the STAT option of the ALLOCATE command:

!   ALLOCATE ( A(N), STAT = AllocateStatus)
!   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
!   Here, AllocateStatus is an integer variable. AllocateStatus takes the value 0 if allocation is successful or some other machine dependent value of there is insufficient memory.

!   An array can be released from memory by using the DEALLOCATE command
!   DEALLOCATE (A, STAT = DeAllocateStatus)
! ******************************************************************************  

  subroutine init_asymtpo(asymt)
  integer, intent(in) :: asymt  ! index of symmetry, used to assign the public values accordingly
  integer :: allocatestatus
 
! flag of symmetric p.o. 
  issymt = 0
  ntp = 1    
  imax = 2
  
! number of control and target variables
  ntar = 4 
  nctr = 5 
  
! fix the value of x(ind) = p0, and try to modify the other component to make the first return to x=x0 plane goes to the same point
  ind = asymt   
  
! allocate memory for ctr and tar   
  allocate ( ctr(nctr), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"
    
  allocate ( tar(ntar), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"

 
  if ( asymt == 1 ) then 
    tar = (/2,3,4,5/)   
    ctr = (/2,3,4,5,6/) 
  
  elseif (asymt == 2) then 
    tar = (/1,3,4,5/)   
    ctr = (/1,3,4,5,6/) 
    
  elseif (asymt == 3) then  
    tar = (/1,3,4,5/)   
    ctr = (/1,3,4,5,6/) 
       
  endif 
 
  if (debug == 1) then          
    print*, 'check init_asympo: asymt, ind,  control and target variables'
    print*, asymt, ',', ind, ',', ctr, ',', tar 
    read(*,*)
  endif
  
  end  subroutine init_asymtpo 

