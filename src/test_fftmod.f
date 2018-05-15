c   test fft_mod

      use sort_mod 
      use fft_mod
      
      implicit real*8(a-h,o-z)
      integer, parameter :: dp = kind(1.d0)
      character*14 fisort
      
      ! the number of points 
      parameter (m = 17,n = 2**m,nd2 = n/2,nd2p1 = nd2+1)
      
      dimension ain(n)  
      dimension a(n), b(n), fcout(100,5), ind_sort(10)
!      dimension fnois(300), ccnois(300), ssnois(300)
      real*8 t, lt !, maxdif 
     
      real(kind=dp), allocatable  ::  fc(:,:)! fcin(:,:) !, ff(:), fc(:,:)
      integer ::  AllocateStatus
      real*8, parameter ::  pi2 = 8.d0*datan(1.d0) 
     
      write(*,*) ' Maximum amplitude of the residual ?'
      read(*,*)  tolres
      
      write(*,*) ' Tolerance in correction for Newton Method 
     1 for refinement?'
      read(*,*) tol 
      
!      tol = 1.d-7
      
!      write(*,*)' Maximum magnitude of the amplitude of the residual ?'
!      read(*,*) tolres
      
      fisort = 'fur.out' 
      open(31,file=fisort)
      

! instead of map, replace ain by the data I want to do the fourier analysis
!      lt = 10.d0
!      np = 2**10
       lt = 64
       np = 64
      
      do i = 1, np  
      ! **NOTE**
      ! It seems we start with i=1 is better than i=0, tested by the same 
        t =  i * lt/ np  !* 2 * pi ! we have 2 cycles
!        t =  (i-1) * lt/ (np-1) !* 2 * pi ! we have 2 cycles
!        a1 = 1.8d0 ; ! 10
!        a2 = -2.d0 ! -2.d0;  ! 27
!        a3 = 11.25d0 ! 4.1
!  
!        a(i) =  10.32+a2*dsin(pi2*1.28d0*t) 
!     1   + a3 *dsin(pi2*4.584d0*t)+ 3*a3 *dsin(pi2*8.293d0*t)

     ! keep this example, and put as an example for the fft methodology report.
     
           a(i) = dcos(pi2*0.13d0*t) 
     1   -.5d0*dsin(pi2*0.27d0*t)  +.75d0*dsin(pi2*0.41d0*t)
  
!     1   + 5*a1 *dsin(pi2*20d0*t)  + 3*a2 *dsin(pi2*27d0*t) ! discard this big frequencies
!        +20 * dcos(pi2*0.136d0*t)  
!        write(111, *) t,  a(i)
      end do
      
      read*
      
      
!      call system('gnuplot dib.plt')  
      ain = a 
      
      nitmax = 15
      call init_fft( np, lt, 10, tol, tolres, nitmax)
      
!      subroutine furian( din, ffc, nf, fc_refn, sep, resmax)
      call furian( ain, 6, nf, fcout, npa, resmax)      
      print*, 'Minimum distance between 2 freqs: ', npa
      
      write(31, *) '# Number of frequencies detected:', nf  
      write(31, *), '# freq     cos     sin     amp      order' 
     
      print*, 'Number of frequencies detected:', nf; read*
      print*, 'freq     cos     sin     amp      order'
      
      ! the costant term + nf, we have nf+1 rows in fcout
      do i = 1, nf+1 
        write(31, '(4f16.10,f10.0)') fcout(i,:)
        write(*, '(4f16.10, f10.0)'), fcout(i,:)
      enddo 
      print*;  
      
 
! ----Discard Third step: Sort by amplitude (4-th col) in decreasing order 
!    omit the first row(the costant term) 

! - It is already in decreasing order due to the way we do the refinement 
!    start from the maximum amplitude 

      allocate( fc(nf,4), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
      
      fc = fcout(2:nf+1, :) 
      
      call sort_2d(fc, nf, 4, 4, -1, ind_sort) 
      
      print*; print*, 'Sorted by amplitude'
      do i = 1, nf, 1
        write(*,*) fc(i,:)
      end do   
      print*; read*
      
 ! ---- Fourth step: Detect the basic frequencies -----------  
      print*; 
      print*, '****Detect the basic frequencies(yes=1)****'
      read(*,*) isbas
      
!      call frebas( fc, nf, ffc,  nfbas, fbas, delmax)
      if(isbas == 1)  call frebas( fc, nf, 6, nfbas, freba, delmax)   
         
      close(31)
      stop
      end

