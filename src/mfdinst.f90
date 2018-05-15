subroutine mfdinst(ypo, tp, mmat, nmf, epsl, ymfd)
! Compute initial state of the unstable manifold, be the dominant eigenvalue
! The monodromy matrix is phi, use power method (gr_mpower) to compute the dominant eigenvector

! the initial state: 
! add a small perturbation on the PO with magnitude of epsl,  long the eigenvector at each point
! 
! 20160108 - modify for use in lorentz force case 
    
! 20150512  -- to f90 format - save ymfd to file ./dat/mfdInst.dat
! Modified to a more complete subroutine, with the initial condition and period, and number of 
! orbits on the manifold, compute the initial conditions for each orbit
 
 
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
!     ypo(6)   the initial point on the p.o.(x,y,z,vx,vy,vz)
!     nmf      number of orbits computed on the manifold   
!     epsl     the  displacement of the initial condition along the vector field to compute manifold 
!              1.d-6 as default 
 
!  Output Variable
!       ymfd(nmf,6) the initial state (6-dimension) for the orbits on mfd

!  idea-- the first way suggested by Gerard
!     use the STM to compute the final perturbation from the eigenvector, 
!     do not need to discrete the p.o., and compute the monodromy matrix
!     for each point

! 20150525 --- thorough check of this subroutine  -ckd!
!     including phi update and perturbation modulus rm 

! function used: gr_mpower

!use emconst, only: pi, xmu, dp
use po_mod, only : gr_rk78, deriv
implicit none

! constant parameters 
integer, parameter:: n = 42, dp=kind(1.d0)
 character(LEN=64) :: fmt1 = "(6d20.10)"

real(kind=dp), parameter :: pi = datan(1.d0) ! little trick to get pi

integer :: imfd, i, j, k, itag 
real(kind=dp):: r(13,n),b(n),f(n), y0(n), hmin, hmax, e1, t, h, tf, & !rk78  
		vep(6), phi(6,6),yf(n), dat(6,6,1), vap,   & ! phi 
		tp, epvel, cj
  
integer, intent(in) :: nmf 
real(kind=dp), intent(in):: ypo(6), mmat(6,6), epsl
real(kind=dp), intent(out):: ymfd(nmf,6)
  
!external deriv ! gr_rtbp: rtbp;   gr_lf: lorentz force

! integeration parameters for  rk78 integrator
hmin = 1.d-10
hmax = 1.d0
e1 = 1.d-13
      
y0 = 0.d0  
y0(1:6) = ypo
y0(7:42:7) = 1.d0

ymfd = 0.d0 ! set all components to be 0 

print*, 'mfdinst: initial sate', y0(1:6)  
 
t = 0.d0
h = 1.d-3  !4 days=100h,  about 6min

itag = 10 
open(itag,file='./dat/mfdInSt.dat', access='append',status='replace')

 
dat(:,:,1) = mmat  !monodromy
 
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

!http://ccar.colorado.edu/asen5050/projects/projects_2012/truesdale/problem.html
!The nature of ε has a large impact on creating the manifolds. The value must be small to correctly compute the manifolds, but too small a value will result in very lengthy computation times. While technically the correctness of the manifolds is approached asymptotically as ε shrinks to zero, in practice ε need only be small relative to the mass parameter of the system.
! In practice, a value on the order of 1.d-4 is used for the Earth-Moon system, 
!   1d-6 is used for the Sun-Earth system [13].

!Another consideration is the difference between position and velocity perturbations. While ε is often given in standard distance units – i.e. 1000 km in the Sun-Earth system – a perturbation of similar magnitude would be huge for the velocity terms. In order to correctly scale the position perturbation, we propose the following:
!eq9
!where R and V are the position and velocity vectors. Since each has a time dependence, it is important to scale the perturbation correctly at each point. This process is used to create perturbed initial conditions in the function manifold_mono.m.

! For velocity, the perturbation is defined as. 
! Update 20150602--- the perturbation is 1.d-6, both for position and velocity

!! this idea is from Parker's paper, give up  
!ymfd(1,1:3) = y0(1:3) + epsl * vep(1:3) / norm2(vep(1:3)) !position
!ymfd(1,4:6) = y0(4:6) + epsl * vep(4:6) / norm2(vep(4:6))!velocity


call gr_mpower(dat,1,vep,vap,1)

ymfd(nmf,:) =  y0(1:6) + epsl * vep/ norm2(vep) ! save as the last one 


!call gr_cjrtbp(ymfd(1,:), xmu, cj)  !ck cj
  	 
   
 
! To calculate the eigenvector at each discreted point along the p.o., 
! No need to compute the Monodromya matrix each time. Compute one time 
! is enough, just at the initial point. Treat vep as a small perturbation, 
! the following vep can be obtained through the State Transit Matrix, which
! is computed along with integration of the state by setting n=42.

! ymfd(i, : ) here is initialized as 0.d0
! For every point, the STM phi should be updated
t = 0.d0
yf = y0

lpmfd: do i = 1, nmf-1, 1  !try to discretimize the whole peroid, keep the final point
  write(*,*) i, 'th orbit on the manifold'
 
  h = 1.d-2
  tf = tp/(nmf-1)*i
  
  do while(t .lt. tf )
    call gr_rk78(t,yf,n,h,hmin,hmax,e1,r,b,f )
  end do
 
  if (dabs(t-tf) .gt. 1.d-12) then !ckd with do while loop, results are same
    h = tf - t
    call gr_rk78(t, yf,n,h,hmin,hmax,e1,r,b,f )
  endif !to get the integration time: tf
  
  
  
  phi = reshape(yf(7:42),(/6,6/)) 
  
  do  j = 1,6
    do  k = 1,6
      ymfd(i,j) = ymfd(i,j) + phi(j,k)*vep(k) 
    end do
  end do
  
  ymfd(i,:) =  yf(1:6) + epsl * ymfd(i,:)/ norm2(ymfd(i,:)) 
 
! give up this idea to use different pertubation for position and velocity 
!  ymfd(i,1:3) =  yf(1:3) + epsl * ymfd(i+1,1:3) / norm2(ymfd(i+1,1:3))!position
!  ymfd(i,4:6) =  yf(4:6) + epsl * ymfd(i+1,4:6) / norm2(ymfd(i+1,4:6)) !velocity
 
  
!  call gr_cjrtbp(ymfd(i+1,:), xmu, cj)  !ck cj
  !write(itag,fmt1) (/ymfd(i,:), cj/)  
  write(itag,fmt1)  ymfd(i,:)
  
end do lpmfd

! add the last orbit on the manifold
write(itag, fmt1)   ymfd(nmf,:)

 close(itag) ! mfdInSt.dat is written inside this subroutine
  
end subroutine mfdinst




