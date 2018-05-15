! Initialization for index of control and target variables.  symtpo + asymtpo 

! ***********************************************************************
!           ---------  init_symtpo   -------

! 6- Symmetry - initial condition - final condtion
!    For the discussion of how to utilize the symmetry to compute the p.o., please refer to Rusell, global search for p.o. P6-symmetries
!    For x-z planar symmetry,       given initial condition(x0,0,z0,0 0,dy0,0) --> T/2 --> (x0,0,z0, 0,-dy0,0)
!    For x-axis symmetry(vt_lyap),  given initial condition(x0,0,z0,0 0,dy0,0) --> T/4 --> (x0,0,0,  0, dyf,dzy) 


! For symmetric p.o.:   idcor = (1,2,3), poincare-map, free-time, fixe-time
!  Specify the index of components for unknowns and equations (control and target )
!  non-zero components for control  --- initial configurations  
!  zero components to target        --- final configuration    

!  For simply-symmetric p.o.s, they are the same type, idic = idfc ( extracted from ctrl and tarl)
!  for doubly-symmetric p.o.s, they are not the same 


!  Because for the doubly-symmetry case, two options to refine,   ( TODO )
!      1.  second crossing through the poincare section, which is of the same kind of mirror configuration as the initial point 
!      2.  the first crossing, with different type of mirror configuration 


!  ------------------- simply-symmetry --  T/2 --------------------------------- 
! Three kinds of planar symmetry, w.r.t.           I.C. and F.C.         non-zero components(control)     zero components(target) 
!    1. x=0 plane,  i.e., y-z plane        P1:  (0, y, z,   vx,  0, 0)           (2, 3,   4 )               (1,  5, 6 )
!    2. y=0 plane,  i.e., x-z plane        P2:  (x, 0, z,    0, vy, 0)           (1, 3,   5 )               (2,  4, 6 )           
!    3. z=0 plane,  i.e., x-y plane        P3:  (x, y, 0,    0, 0, vz)           (1, 2,   6 )               (3,  4, 5 )

!  and also three kinds of axis symmetry, w.r.t 
!    4. x-axis                             A1: (x, 0, 0,    0, vy, vz)           (1,    5, 6)               (2,3,   4)
!    5. y-axis                             A2: (0, y, 0,   vx,  0, vz)           (2,    4, 6)               (1,3,   5) 
!    6. z-axis                             A3: (0, 0, z,   vx, vy,  0)           (3,    4, 5)               (1,2,   6) 
 
! For simply-symmetric p.o., the problem of determination of p.o. is reduced 
! to find two mirror configuration of the same type at two different epoches. 

! with the point on the symmetric element( in the axis for A type and on the plane for P type), 
! and the velocity perpendicular to that element. 

! The control variables are the  non-zero components in the initial mirror configuration,  the rest componnents are zeros 
! while the target variables are the zero components in the final mirror configuration 


! --- TODO: for the moment, we skip the doubly symmetric case (vertical lyapunov orbit)... 2017-02-16 10:13:55 
!  ------------------- doubly-symmetry --  T/4 --------------------------------- 

! start with A1(x,0,0) --- after T/4, goes to P2 (x,0,z), 

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
!                1 : poincare map approach to speicify the next epoch   ! the current used one 
!                2 : free-time shooting
!                3 : fixed-time shooting? --TODO:  test
                 
                
! ****************************************************************************** 
  subroutine init_symtpo(idic, idfc, idcor0, ind0)
  
  integer, intent(in) :: idic, idfc, idcor0, ind0   ! index of symmetry, used to assign the public values accordingly

! Local Varaibles
  integer ::  ctrl(6, 3), tarl(6,3) ! index of control and target varables for all six possible mirror configuration 
  integer ::  allocatestatus, i, j

! for symmetric p.o.   
  idcor  = idcor0
  issymt = 1
  
  ! by default, period is not a control variale
  tfree  = 0 
  
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
  ! the four-th control varaible TP will be updated within routines by check idcor == 2? 
    tfree = 1
    nctr = 4  
    ntar = 3
  
  elseif( idcor == 3 ) then 
  !  fixed time shooting
    nctr = 3
    ntar = 3  
  end if
  

! allocate memory for tar and ctr   
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

! refine doubly symmetry, if idic /= idfc, we can see it is doubly symmetry,
! this is how we compute the near planar one that is both symmetric w.r.t y-z plane and x-z plane, where idic = 1, idfc =2  
  if( idic == idfc ) then
    ntp = 2
!    if(idic > 3) then imax = 2 
  else
    ntp = 4
  end if
   
   
  if(idcor == 1)   ntp = ntp / imax    ! Poincare map method 

!  check the assignment 
  print*, 'idic=', idic, ', idfc=', idfc, ', imax=', imax, ', ntp=', ntp
  print*, 'ctr=', ctr, ', tar=', tar 
  read*
      
  
  if (debug == 1) then          
    print*, 'check init_po: ind, p0, index of control and target var'
    print*,  ind, p0, ',', ctr, ',', tar 
    read(*,*)
  endif   
  
  end subroutine init_symtpo
  
  
  
! ********************************************************************
!     ----- init_asymtpo -----
!  ** Idea ** :  
!  Fix the value of a component as the poincare section x = x0, 
!  in principle an orbit passing this initial point will cross the section twice, 
!  the second intersection with  x = x0 will return to the initial point, with the velocity
!  in x component  in the same direction as the initial one
!  this serves as the periodicity condtion. 
! 
!  ** A potential bug ** TODO  2017-02-21 10:08:43 
!  what if the initial point with x=x0 is exactly the point that only has one intersection with the x=x0 plane.... 
!  It rarely happens, but there is possibility.

!  Solution: Take the most secure approach, use gr_po instead....

!  The energy h0 and x0 is fixed, look for the first return to x=x0 plane,
!  we only need to ask four final components in state to be as the initial one, 
!  the fifth will be satisfied simultaneously by the energy conservation law. 
!  And we have 5 initial components to modify(except x0)

! constraint funcion:   ds^2 = ( (yf-y0)^2 + (zf-y0)^2 + (vxf-vx0)^2 + (vyf-vy0)^2 ) to be minimum

!   asym          contol            target 
!   1    (y0,z0, vx0,vy0,vz0) - (2,3, 4,5,6)   (yf, zf, vxf, vyf) 

!  To ensure that enough memory is available to allocate space for your array, make use of the STAT option of the ALLOCATE command:

!   ALLOCATE ( A(N), STAT = AllocateStatus)
!   IF (AllocateStatus /= 0) STOP  "*** Not enough memory ***"
!   DEALLOCATE (A, STAT = DeAllocateStatus)

! ******************************************************************************  

  subroutine init_asymtpo(asymt)
  integer, intent(in) :: asymt  ! index of symmetry, used to assign the public values accordingly
  integer :: allocatestatus
 
! flag of symmetric p.o. 
  issymt = 0
  ntp  = 1    
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

! ****************************************************************************** 
! The initializaiton for asymmetric p.o.:  idcor = (5,7) 
!  Initialize the index of the control variables and target variables as in the state vector

! Module based varaibles
!  tfree, idcor (5,7)

! ****************************************************************************** 
  subroutine init_asymtpo2(idcor0, ind0, x0_fxd0)
  
  ! this is try to add a constraint to fix a component of the state vector, y(ind0)= x0_fxd0
  ! so the control variables are still all the six components, but the number of target function increases by 1 
  
  integer,intent(in) :: idcor0, ind0
  real(kind=dp),intent(in) :: x0_fxd0
 
  
  integer ::  allocatestatus, nrow, ncol, i, j


!  module-based variables  
  idcor = idcor0 ! supposed to be 5 
  issymt = 0 
  tfree  = 0 ! by default, fixed-time shooting
 
 
! ! --- allocate memory for ctr and tar  for asymmetric p.o.  
  
  if( idcor == 5 ) then   !-- only finish this one 
  !-- doesn't work very well, the continuation only proceed for a few orbits 
  !  free time shooting + additional constraint
    ind = ind0 
    x0_fxd = x0_fxd0
    nctr = 7   ! the last one is to save period...
    ntar = 6
    
    tfree = 1
    nrow = ntar - 1
    ncol = nctr - 1 
  
  elseif (idcor == 7) then  ! improve idcor=5 strategy
  ! free time shooting, deliminate ind from tar to satisfy x(ind) = x0_fxd 
    ind = ind0
    nctr = 6  ! five state components + TP
    ntar = 5
    
    tfree = 1 
    nrow = ntar 
    ncol = nctr-1   
    
  ! ----------- TODO  ------------------
  elseif( idcor == 4 ) then !-- TODO  
    ! poincare map approach, by default, we search for the first crossing through the poincare section
    imax = 1
    ind = ind0
    p0 = x0_fxd0 
    
    ntar = 4
    nctr = 5
    
    nrow = ntar
    ncol = nctr 
    
  elseif( idcor == 6 ) then 
    ! free time shooting -- no additional constraint....
    tfree = 1
    nctr = 7
    ntar = 6 
     
    nrow = ntar
    ncol = nctr - 1 
    
  end if
  
    
  ntp = 1
 
 
 ! allocate the memory for tar and ctr 
  allocate ( tar(nrow), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"

  allocate ( ctr(ncol), stat = allocatestatus)
  if (allocatestatus /= 0) stop "*** not enough memory ***"
  
  
 ! assign the vaules for tar and ctr 
  if ( idcor == 5  ) then ! the current used one, additional constraint: y(ind) = x0_fxd
    tar = (/1,2,3,4,5/)   ! x_6^f is satisfied by conservation energy law
    ctr = (/1,2,3,4,5,6/) 
  
  elseif( idcor == 6  ) then ! --TODO with full state vector as the control variables    
    tar = (/1,2,3,4,5,6/) 
    ctr = (/1,2,3,4,5,6/) 
  
  elseif( idcor == 7)  then  ! -- to approve idcor = 5 
    tar = (/1, 2, 3, 4, 5 /)
     
    ! instead of an additional equation, fix x(ind)=x0_fxd as the initial condition,
    ! and eliminate ind-th component from ctr, as well as ind-th column in Phi
    j = 0
    do i= 1, 6, 1
!      print*, 'ind=', ind, 'i=', i 
!      read*
      
      if( i /= ind ) then 
        j = j + 1
        ctr(j) = i
      endif 
    end do
    
      
  end if  
  
  
  print*,' ind, x0_fxd', ind, x0_fxd
!  read*
  
  print*, 'idcor=', idcor, 'nctr=', nctr, 'ntar=', ntar
  print*, 'ctr=', ctr, 'tar=', tar 
  read* 
  
  end subroutine init_asymtpo2


! the problem with this approach lies in the continuation 
! for the computation of the kernel  v = DF, if we keep the ind-th component as one of the control variables, 
!  then the last row of DF = [1 0 0 0 0 0 0 ], question: the last element shoule be dx(ind) or 0? with 0 fails to compute the det, because it is equal to  0
