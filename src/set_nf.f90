! -- check the Maximum norm of csf and sif, and choose an appropriate value for nf 
! FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
!inf = nf_init

!do 
! inf =  inf / 2 
! ckm = dlange('m', n, inf, csfn(:, nf_init-inf+1 : nf_init), n, 0.d0)
! skm = dlange('m', n, inf, sifn(:, nf_init-inf+1 : nf_init), n, 0.d0)
! 
! do k = 1, n, 1
!   write(*,'(10e24.14)') csfn(k, nf_init-inf+1 : nf_init)
!   write(*,'(10e24.14)') sifn(k, nf_init-inf+1 : nf_init)    
! end do
! print*
! 
!print*, 'maximum norm of CK and SK for the last ', inf, 'Fourier harmonics'
!print*,  ckm, skm 

!print*; read*

! if(inf <= 2) exit  
!enddo 


!inf = nf_init/2

!do 
! inf = inf/2
! ckm = dlange('m', n, inf, csfn(:, nf_init/2-inf+1 : nf_init/2), n, 0.d0)
! skm = dlange('m', n, inf, sifn(:, nf_init/2-inf+1 : nf_init/2), n, 0.d0)
! 
! do k = 1, n, 1
!   write(*,'(10e24.14)') csfn(k, nf_init/2-inf+1 : nf_init/2)
!   write(*,'(10e24.14)') sifn(k, nf_init/2-inf+1 : nf_init/2)    
! end do
! print*
! 
!print*, 'maximum norm of CK and SK for the last', inf, 'Fourier harmonics before nf_init'
!print*,  ckm, skm 

!print*; read*

! if(inf <= 2) exit  
!enddo 
