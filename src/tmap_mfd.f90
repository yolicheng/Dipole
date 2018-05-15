!***********************************************************************
!     ****   tmap_mfd   ****
! This routine is to compute the Tmap of a manifold associated to the eigenvector vep in the interval [t0: dt: tf], 
! just to save place for main_tmapmfd_plf.f90

!  ****** Input Variables ******
!  y0       I.C. of the associated p.o. 
!  tp       period of the p.o.
!  ndim     dimenison of the state vector 
!  nmf      number of orbits on the manifolds to compute 
!  vep      eigenvector (direction)
!  epsl     diviation along the eigendiretion
!  tdir     time sense, 1 : Unstable manifold;  -1: stable manifold 
!  t0,tf    starting and ending epoch for the Time-T map 
!  dt       time interval for two consecutive maps 
!  finit    file tag for I.C. of the manifold along the vep direction 
!  ftmap    file tag for the Time-T map 
!  ispl     1: plot the orbit of the manifold; others: no 
!  fob      file tag for the orbits on the manifold if ispl == 1 

!  Routine Used:
!     gr_mfdinst   plob_n

!  Finally Revised by Yu -- 20160

! TODO: check the plot of orbit during the computation of Time-t map
!***********************************************************************
subroutine tmap_mfd(y0, tp, ndim, nmf, vep, epsl, tdir, t0, tf, dt, finit, ftmap, &
                    ispl, fob,  deriv, deriv_n, gr_cj)

use dp_mod
implicit none
 
! Input  and Output Declaration   
integer, intent(in)            ::  ndim, finit, ftmap, fob, tdir, nmf, ispl
real(kind=dp), intent(in)      ::  y0(ndim), tp, vep(ndim), epsl,  & ! mfd 
                                   dt, tf, t0    
 
! Local Variable
integer        :: i, k !, ntmap
real(kind=dp)  :: ymfd(nmf, ndim), yn(nmf, ndim), ti0,  tif, thf, yi(ndim), cj, yfmfd(nmf, ndim)

external :: deriv, deriv_n,  gr_cj
  

  ! -- I.C.
  print*, 'check tmap_mfd'; print*; read*; !ck
  
  print*, 'y0, ndim, tp, vep, nmf, epsl, finit', y0, ndim, tp, vep, nmf, epsl, finit
!  print*, 'deriv, gr_cj', deriv, gr_cj
  print*; read*; !ck
  
  call gr_mfdinst(y0, ndim, tp, vep, nmf, epsl, finit, ymfd, deriv, gr_cj) 
  close(finit)
  ! -- Time-T map 
   
  ! -- Discard the I.C. of the manifold 
  ! Do not need to  save the initial curve at t=0  which is the initial condition saved in mfd_init**.dat
  ! do k = 1, nmf, 1
  !   yi = ymfd(k,:);   call gr_cjplf(yi, cj)
  !   write(234, *) 0.d0, yi, cj
  ! end do
  ! write(234, *) 
  
  yn =  ymfd 

  ! Check [t0: dt : tf] to see how the loops appear, but when it's close to (t0+tf)/2, use smaller stepsize 

  ! ntmap = idint( (tf - t0) / dt ) + 1
  
  ! The initial start time for the first epoch 
  tif = t0 - dt 
  thf = (t0 + tf) / 2.d0
  
  print*, 't0, tf, thf', t0, tf, thf 
  ! Evaluation at ntmap epoches
  i = 0
  
  do  
    i = i + 1
    
    ti0 = tif
    tif = ti0 + dt
    if(dabs( tif - thf) <= 2*dt ) tif = ti0 + dt / 2.d0
    
    
    if(i == 1) ti0 = 0.d0 
    if(tif > tf) exit 
    
    print*, 'ti0, tif, tdir', ti0, tif, tdir; print*; read*; !ck
    
    call plob_n(yn, ndim, nmf, ti0, tif, tdir, ispl, fob, yfmfd, deriv_n, gr_cj) 
  
    ! save the points at tf 
    do k = 1, nmf, 1
      yi = yfmfd(k,:) ;    call gr_cj(yi, cj)
      write(ftmap, *) tdir*tif, yi, cj
    end do
    write(ftmap, *);  write(ftmap, *); ! add two blank lines to seperate epoches to avoid extra line segments when plotting with lines 
  
    ! update the I.C. for next time at map
    yn = yfmfd
  end do

  close(ftmap) 
  if(ispl == 1) close(fob)

  return  
  
end subroutine tmap_mfd

