subroutine gr_mfdinst(ypo, n, tp,  vep, nmf, epsl, ftag, ymfd, deriv, gr_cj)

! Compute initial state of the unstable manifold, be the dominant eigenvalue
! The monodromy matrix is phi, use power method (gr_mpower) to compute the dominant eigenvector

! the initial state: 
! add a small perturbation on the PO with magnitude of epsl,  long the eigenvector at each point
! 
! 20161011 - modify for general use
! 20160108 - modify for use in lorentz force case 
! 20150512  -- to f90 format - save ymfd to file ./dat/mfdInst.dat

! Modified to a more complete subroutine, with the initial condition and period, and number of 
! orbits on the manifold, compute the initial conditions for each orbit
 
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

  ! http://ccar.colorado.edu/asen5050/projects/projects_2012/truesdale/problem.html
  ! The nature of ε has a large impact on creating the manifolds. The value must be small to correctly compute the manifolds, 
  ! but too small a value will result in   very lengthy computation times. While technically the correctness of the manifolds 
  ! is approached asymptotically as ε shrinks to zero, in practice ε need only be small relative to the mass parameter of the system.
  ! In practice, a value on the order of 1.d-4 is used for the Earth-Moon system, 
  !   1d-6 is used for the Sun-Earth system [13].

  !  Another consideration is the difference between position and velocity perturbations. 
  !  While ε is often given in standard distance units – i.e. 1000 km in the Sun-Earth system 
  !   – a perturbation of similar magnitude would be huge for the velocity terms. In order to correctly scale the position perturbation,
  !  we propose the following:
  !eq9
  ! where R and V are the position and velocity vectors. Since each has a time dependence, 
  ! it is important to scale the perturbation correctly at each point.
  ! This process is used to create perturbed initial conditions in the function manifold_mono.m.

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
!     do not need to discrete the p.o., and compute the monodromy matrix
!     for each point

! 20150525 --- thorough check of this subroutine  -ckd!
!     including phi update and perturbation modulus rm 

! function used: gr_mpower
use dp_mod
use pi_mod
implicit none

! ------------ Input and Output---------------
integer, intent(in)        :: n, nmf, ftag  
real(kind=dp), intent(in)  :: ypo(n), tp,  epsl, vep(n) 
real(kind=dp), intent(out) :: ymfd(nmf, n)
  
! -- Local Variables
integer       ::  i, nvar 
real(kind=dp) ::  r(13, n*(n+1)), b(n*(n+1)), f(n*(n+1)), yf( n*(n+1) ),  &
                  hmin, hmax, e1, t, h, & !rk78  
                  cj,  dh, tf,  tc, yc( n*(n+1) ), phi(n,n), yi(n), dt 

real(kind=dp) ::  dnrm2  
external deriv, gr_cj ! gr_rtbp: rtbp;   gr_lf: lorentz force
character(LEN=64) :: fmt1  


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
  h = dmin1(1.d-4, dt / 10.d0)


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
    yi = matmul(phi, vep)
    
    yi =  yc(1:n) + epsl * yi / dnrm2( n, yi, 1) 
    call gr_cj(yi,  cj)  
    
    write(ftag, fmt1)   tc, yi, cj
    
!    write(*, fmt1) tc, yi, cj; read* ! ck
    
    ymfd(i,:) = yi
  end do lpmfd
  
end subroutine gr_mfdinst




