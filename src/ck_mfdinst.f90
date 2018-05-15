subroutine ck_mfdinst(ypo, n, tp,  vep, nmf, epsl, ftag, ymfd, deriv, gr_cj)

! check how the stable(or unstable) eigenvector evolves along the manifold
! 1 -- by variational matrix to transform the initial vep 
! 2 -- at the new point, compute the momoddromy and eigenvector, to compare with the first approach 

  
! Discretimize:  take points equally spaced in time along the periodic orbit 
!                dt = TP / nmf  
!                and take t =  (0 : nmf-1 ) * dt
  

! Compute initial state of the unstable manifold, be the dominant eigenvalue
! The monodromy matrix is phi, use power method (gr_mpower) to compute the dominant eigenvector
! Sometimes when the insteability is mild, power method may fail, so we use eigrg here 

! The initial state: 
! add a small perturbation on the PO of magnitude epsl, along the stable (unstable) eigenvector 

 
  !calculate the eigenvector of the initial point on p.o. 
  ! Reference:
  !    Dynamics and Mission Design Near Libration Points, Volume 3, P96, 
  !    4.2 local approximation of stable manifold.

  ! The eigenvector is scaled in such a way that the three components of these vectors 
  ! relative to the position have Euclidean norm equal to one. So that we write as:
  ! vep = vep / sqrt(vep(1)**2 + vep(1)**2 + vep(1)**2  )
  !  
  !  In this way if we want to take initial conditions at a selected physical distance from 
  !  point ymfd(i,:) which is in the p.o., we only have to multiply vep by the selected 
  !  distance in km with the desired sign.

  ! For velocity, the perturbation is defined as. 
  ! Update 20150602--- the perturbation is 1.d-6, both for position and velocity

  !! this idea is from Parker's paper, give up  
  ! ymfd(1,1:3) = y0(1:3) + epsl * vep(1:3) / norm2(vep(1:3)) !position
  ! ymfd(1,4:6) = y0(4:6) + epsl * vep(4:6) / norm2(vep(4:6))!velocity
 
! 20150409 -- Default as unstable manifold
!   The initial state y0 is for unstable manifold, it's approximated by
!   the state on the p.o. + a small perturbation(1d-6 is small enough) along the eigenvector
 
!   It's computed by dominant eigenvector, with the corresponding eigenvalue > 1 
! 
!    x0_wu = x0_po + epsl* vep/vepm

!  Note: 
!    The direction vep must be normalized to unit vector, to keep the perturbation 
!    on different orbits of the same magnitude

!  Input Variable:
!     ypo(n)   the initial point on the p.o. :  POS+VEL
!     n        dimension of the state vector
!     tp       the period of the p.o. 
!     nmf      number of orbits computed on the manifold   
!     epsl     the  displacement of the initial condition along the vector field to compute manifold 
!              1.d-6 as default 
!    ftag      the file to save all the initial conditions 
  
!  Output Variable
!   ymfd(nmf,n)   the initial state (n-dimension) for the orbits on mfd

!  idea-- the first way suggested by Gerard
!     use the STM to compute the final perturbation from the eigenvector, 
!     do not need to  compute the monodromy matrix for each point


! function used: gr_mpower
use dp_mod
use pi_mod
implicit none

! ------------ Input and Output---------------
integer, intent(in)        :: n, nmf, ftag  
real(kind=dp), intent(in)  :: ypo(n), tp,  epsl, vep(n) 
real(kind=dp), intent(out) :: ymfd(nmf, n)
  
! -- Local Variables
integer       ::  i,  nvar 
real(kind=dp) ::  r(13, n*(n+1)), b(n*(n+1)), f(n*(n+1)), yf( n*(n+1) ),  &
                  hmin, hmax, e1, t, h, & ! rk78  
                  dt, cj,  dh, tf,  tc, yc( n*(n+1) ), phi(n,n), yi(n), vepi(n)  
! ck vepi
real(kind=dp) ::  yc2 (n*(n+1) ), phi2(n,n), yf2(n*(n+1)), wr(n), wi(n), vr(n, n)


real(kind=dp)     ::  dnrm2  
external          ::  deriv, gr_cj ! gr_rtbp: rtbp;   gr_lf: lorentz force
character(LEN=64) ::  fmt1  


  ! Format to save the I.C.
  nvar = n*(n+1)  
  write(fmt1,  fmt='(a,i0,a)')  '(', n+2, 'e24.14)' ! time + state vector + energy
  
  write(ftag, *)  ! add one blank line to seperate different manifold
 
  ! Time interval between 2 consecutive points
  dt =  tp / nmf
    
  ! integeration parameters for  rk78 integrator
  hmin = dmin1(1.d-6, dt / 10.d0)
  hmax = 1.d0
  e1   = 1.d-14
  
  ymfd  = 0.d0 ! set all components to be 0 

! print*; print*, 'Check mfdinst: initial sate',  ypo  
 
  t = 0.d0
  h = dmin1(1.d-3, dt / 10.d0)


 ! the first orbit on the manifold, t = 0 
  yi =  ypo  + epsl * vep/ dnrm2(n, vep, 1) ! save as the last one 
  ymfd(1, :) = yi 
  
  call gr_cj( yi,  cj )
  write(ftag, fmt1)  0.d0,  yi,  cj
   
 
  ! To calculate the eigenvector at each discreted point along the p.o., 
  ! No need to compute the Monodromya matrix each time.  We can treat vep as a small perturbation, 
  ! the following vep can be obtained through the State Transit Matrix, which
  ! is computed along with integration of the state by setting nvar = n*(n+1)

  !***  CKD ***  the two approaches are very similar 
    
    
  ! ymfd(i, : ) here is initialized as 0.d0
  ! For every point, the STM phi should be updated
  t  = 0.d0
  yf = 0.d0
  yf(1:n) = ypo
  yf(n+1 : nvar : n+1) = 1.d0 
  
  
  lpmfd: do i = 2, nmf, 1  
   
    tf = (i-1)*dt
    
    do while(t .lt. tf )
      call gr_rk78(t, yf, nvar, h,hmin,hmax, e1,r,b,f, deriv)
    end do
 
    tc = t
    yc = yf
   
    ! integate backwards using a very small step size dh = tf - tc
    dh = tf - tc
    call gr_rk78(tc, yc, nvar, dh, hmin,hmax, e1,r,b,f, deriv )
!    print*, 'Check, t =', tc, '== tf?', tf ;  print*; read*
  
    phi = reshape(yc(n+1 : nvar), (/n, n/)) 
    vepi = matmul(phi, vep)
    
    print*, 'Vep by STM:', vepi; print*; read* 
    
    ! another approach to compute the small deviation is to do modulus on the position 
    ! and that has geometric meaning, we displace a small distance epsl away for the p.o. 
    
    yi =  yc(1:n) + epsl * vepi/ dnrm2( n, vepi, 1) 
    call gr_cj(yi,  cj)  
    write(ftag, fmt1)   tc, yi, cj
    
!    write(*, fmt1) yi, cj; read* !ck
    ymfd(i,:) = yi
    
    
    ! ---- check the eigenvector vep by approach 2 
    ! CKD!!! ----- both approaches are very similar 
    yc2 = 0.d0
    yc2(1:n) = yc(1:n)
    yc2(n+1 : nvar : n+1) = 1.d0 
  
    call plob(yc2, 0.d0, tp, nvar, 1, 6, 0, deriv, gr_cj,  yf2) 
    phi2 = reshape(yf2(n+1 : nvar), (/n, n/)) 
    print*, 'yc2, yf2 the same?'
    print*,  yc2(1:n);
    print*, yf2(1:n); print*; read*
    
    call eigrg(phi2, n, 1, wr,wi, vr)
    print*, 'Check vep by MM';  read* 
    
  end do lpmfd
  
  
end subroutine ck_mfdinst




