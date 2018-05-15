!  **** My fft Module*******
! TODO : replace the sines and cosine by the recurrent procedure done by Gerard and me in 20160719 

!  Try to work the refinement independently from laskar.f, and work directly with the result by fftsc 

!  We deal with the frequency one-by-one, refine one frequency and the coefficients and substract from the orginal data, when a new frequency is detected, we refine together with all the previous detect frequencies again. 

! 1.  Given equally spaced input samples, first multiply the Hanning window
!     function, and then do simple FFT with the hanning-flitered data to
!     aviod leakage 

! 2.  Detect only one frequency with the maximal amplitude, and use Newton
!     Method to refine the only coefficients, better than simultaneously....

! ---discard---- this step proves uneccessary when the modulus of the amplitude 
!                of the residual is greater than tolerance\
!! 3.  Refine again with all the exiting frequencies using the orginal samples
!     if we have more than 1 frequencies.
!     
! 3.  Substract all the refined frequency from the original data. 

! 4.  Start again from step 1 with the residual instead of the original sample
 
! --- repeat this process until the modulus of the residual is smaller than tol--- 

! 5.  Do a final simultaneous refinemnt(freq+coef) of all the detected frequecies
!    . 

!/* --- from Josep Maria's comment
! * Computes Laskar's approximation, but using DFT instead of a
! * numerical quadrature formula (the DFT is a better approximation of the
! * continuous integral, due to the discrete Poisson's summation
! * formula.
! * - test: TEST_FURIAN_LAPROX
! */
 
! Question: Is is good to start with the first peak? 
!           or with the one has the maximum amplitude? ---better!!

! So we need to save the original input as constant variables, it is convenient to use module 

! --- finally tested, 20160614 
!  with example: 
!       a(i) =  10.32+1.8*dsin(pi2*1.28d0*t) 
!     1   + 11.25 *dsin(pi2*4.584d0*t)+ 33.75 *dsin(pi2*8.293d0*t)
!      
module fft_mod

implicit none
save 

! subroutines from libraries 
external  ::  deltx, trigrec_c, trigrec_s   ! -- libgr.a 


! --- parameters------
integer, parameter, private :: dp = kind(1.d0)
real(kind=dp), parameter, private ::  pi2 = 8.d0*datan(1.d0)

! --- Samples, Number of samples ---
integer, private :: np, npm, npd2, npd2p1  
!real(kind=dp), allocatable, private ::  din(:) ! input data -- discard

! --- For control 
integer, private       ::  nh, nitmax, sepmin, nfmax 
real(kind=dp), private ::  lt, tol, tolres 

! --- For save 
! to write the output of refinment process
integer, private :: fwork, fout 



! ---Module private variables--- must be initialization before any call of this module 
! tolres:      tolerance of the maximum amplitude of the residual for frequency detecton 
!            
! tol:     tolerance for the correction in Newton Method for the refinement   
!  np        number of samples
!  lt        length of the time interval of the samples
!            the epoch of samples:  t = i*lt/np, i = 1, np  

contains

!  !***********************************************************************
!     ****   init_fft  ****
!  Initialize the commonly shared variables
!-- by input
!   np, din, tol, tolres  
!     -- the last two suggestion value: 1.e-8, 1.e-3

!-- by default
!   nh  

!***********************************************************************
subroutine init_fft( np0, lt0, nfmax0, tol0, tolres0, nitmax0)
implicit none

! Input  and Output Declaration  
integer, intent(in)       ::  np0, nitmax0, nfmax0 
real(kind=dp), intent(in) ::  lt0,  tol0, tolres0   
 
! Local Variable
! integer :: AllocateStatus
  
! ----- default values      
  ! order of Hanning window function     
  nh = 2              

!  minimun difference in order of two detected frequecies 
  sepmin = nh+1 
  
! ----- by assignment in the callee
  lt     = lt0
  tol    = tol0
  nfmax  = nfmax0
  tolres = tolres0 
  nitmax = nitmax0
  
! ---- commonly used value for dimension declaration  
  np  = np0
  npm = int( dlog( dble(np) ) / dlog(2.d0) ) ! log_2(np)
  npd2 = np / 2 
  npd2p1 = npd2 + 1
  

!  for save
   fwork = 6
!  fwork = 222 
!  open(fwork, file='fur.work')
  
  
  return  
end subroutine init_fft

  
!***********************************************************************
!     ****   furian   ****
!  Main routine to do Furier Analysis

! extract the frequencies that has magnitude greater than tolra, the number is returned in nfre 

! To avoid frequencies to be too close,  we only saved the peak,  am_k-1 < am_k > am_k+1

! 1.  Given equally spaced input samples, first multiply the Hanning window
!     function, and then do simple FFT with the hanning-flitered data to
!     aviod leakage 

! 2.  Detect only one frequency with the maximal amplitude, and use Newton
!     Method to refine  the coefficients. To avoid leakage, the distance 
!     between 2 frequecies should be greater than sepmin=1+nh

! 3.  If we have more than 1 frequencies, refine again with all the exiting frequencies
!     using the orginal samples for freq+coef simultaneously 
!     
! 4.  Substract all the refined frequency from the original data. 

! 5.  Start again from step 1 with the residual instead of the original sample
 
! --- repeat this process until the modulus of the residual  
!     or the amplitude of the new detectd frequecy is smaller than tol         ----

!       Input Variables 
!  din       input data array, equi-spaced samples,  1 column     
!  ffc       file tag to write all the detected frequecies
    
!       Output Variables 
!  nf       number of frequencies detected and refined
!  fcout    dimension nf-by-4: freq-cc-ss-amp-order          

!  Module-based Varaibles
!   tol, tolf, dp  

!  Routine Used:
!     

!  Finally Revised by Yu -- 20160612
!***********************************************************************
subroutine furian( din, ffc, nf, fc_refn, sep, resmax)
      
implicit none

! Input  and Output Declaration 
integer, intent(in)      ::  ffc   
integer, intent(out)     ::  nf, sep
!real(kind=dp), intent(in)   ::  lt, tol, tolres! by int_fft  
real(kind=dp), intent(out)      ::  fc_refn(100, 5), resmax  
 real(kind=dp), intent(inout)   ::  din(np)
 
! Local Variable
real(kind=dp)  :: dwork(np), dwork0(np), fc_fft(npd2p1, 4)

! for fftsc  
! wk is not used if np is a power of 2
integer        :: iwk( npm) 
real(kind=dp)  :: st(npd2p1), ct(npd2p1),  wk(6*npd2+150) 
double complex :: cwk(npd2p1)  

! for refinement 
integer       ::  AllocateStatus, isref, nfmax, i, j, sepi, sepimin, feqnew, nf1, nff 
real(kind=dp) ::  maxferr,  fcin(2,4), fc_ava(npd2-1) , prom
real(kind=dp), allocatable :: fcinall(:,:)   

! compute the residual from the average  of ain
call promig(din,np, prom) ! TODO this has none-sense, at least I didn't find it
!print*, 'prom=', prom ; read*

! working array, not touch din, keep as constant
dwork= din 

nf = 0
sep = npd2 ! the minimum distance between 2 consecutive frequecies  

sepimin = npd2  ! the distance of the new detected frequncy from the previously refined ones.

! initialize the constant term 

fc_refn(1,:) = 0.d0 
do   
! 1.  Given equally spaced input samples, first multiply the Hanning window
!     function, and then do simple FFT with the hanning-flitered data to
!     aviod leakage   

! compute the mean and work with the residual(the deviation from the mean for FFT)
! or work directly with fft? -- choose this option, easier, but it works bad with non-zero constant term...
  feqnew = 1
  dwork0 = dwork 
  
  call hannin(dwork, np) 
  
  call fftsc(dwork, np, st, ct, iwk, wk, cwk) 
  call modulus_fftsc(st, ct, np, lt, fc_fft)
  
! 2.  Detect only one frequency with the maximal amplitude, and use Newton
!     Method to refine (freq+coef) simultaneously.

!     TODO---To avoid leakage,  when a new frequecy is detected, we check 
!            the distance  in order of the two frequecies should be >= sepmin
!            sepmin = 1+nh 

!     we then will work with fc_fft(nfmax+1, :) for the refinement
  
  if(nf == 0 ) then 
    fc_ava  = fc_fft(2:npd2, 4)
    nfmax = maxloc(fc_ava,1)  ! the location of the maximal amplitude
  
  else   ! check the distance...

    do 
      feqnew = 1 ! flag to show the detection of new frequecy
      
      ! the avaible feq+coef (delete the aliased one by setting amp=0)
      fc_ava  = fc_fft(2:npd2, 4)
      nfmax = maxloc(fc_ava, 1) ! the location of the maximal amplitude

      write(fwork,*) 'nfmax =',  nfmax 
      write(fwork,*)  fc_fft(nfmax+1, :);  write(fwork,*) 
      
      do i = 2, nf+1 
        sepi = iabs( nfmax - int(fc_refn(i,5)) )
!        print*, 'sepi, nfmax', sepi, nfmax; read* !ck
        
        if( sepi <= sepmin ) then 
        ! if the frequecy is already refined, but not so accurate, 
        ! we forcely reset the amplitude to be 0, and find the frequecy
        ! with the second maximal amplitude 
           
           write(fwork,*) 'Aliased Frequency!'
           write(fwork,*) 'Detected:',  fc_refn(i,:)
           write(fwork,*) 'Aliased: ',  fc_fft(nfmax+1,:) 
           write(fwork,*) ; read*
           
           fc_fft(nfmax+1, 4) = 0.d0 ! 
           feqnew = 0
           exit
        endif
       
        if(feqnew == 1 .and. sepi < sepimin )  sepimin = sepi ! the minimum distance between two frequencies 
      enddo 
      
      ! check the amplitude of the new frequecies
      if(fc_ava(nfmax) < tolres) then 
        write(fwork,*)  'Amplitude of new freq in residual is less than tolres', fc_ava(nfmax-1)
        feqnew = 0
        exit ! to terminate the refinement
      endif 
      
      if(feqnew == 0) cycle
      if(sepimin < sep)  sep = sepimin
      exit !succeed a new frequecy
    enddo 
    
  endif 
  
  if(feqnew == 0) exit 
! now we obtain the avaible new frequency, which has the maximal amplitude in the residual  
  write(fwork,*),'The maximum amplitude'
  write(fwork,*), fc_fft(nfmax+1, :)
  write(fwork,*)
      
  nf = nf + 1
  
! update the constant term 
!  fc_refn(1,2)   =   fc_fft(1,2) ! cc 
  fc_refn(1,2)  = fc_refn(1,2) + fc_fft(1,2) ! cc 
  fc_refn(1,4)   = dabs(fc_refn(1,2) ) ! amp
  
!  print*, 'nf, The costant term:', nf,  fc_refn(1,2); !read*!ckd
  
  fc_refn(nf+1, 5)   = nfmax  ! keep the order in 5-th column to check the distance 
  fc_refn(nf+1, 1:4) = fc_fft(nfmax+1, :)

! ---- for refinment of only 1 frequecy in coefficients ---- 
  fcin(1,:)   = fc_fft(1, :)     ! the costant term 
  fcin(2,:)   = fc_refn(nf+1, 1:4) 
   
! when refine only 1 frequecy, we refine only the coefficients and fix the frequecy
! refine approach: 2
  nf1 = 1
!  call refine( fcin, nf1,  dwork0, 2, maxferr, isref)  ! no big difference...
    call refine( fcin, nf1,  din, 2, maxferr, isref)    ! this is better? 
! call refine( fcin, 1,  din, 3, maxferr, isref) ! disaster, for only 1 frequecies...
  
  if( isref == 1 ) then
  ! update the refined  coef
    fc_refn(nf+1, 1:4) = fcin(2, :)
    write(fwork,*) nf, '-th new frequecy'
    write(fwork,*) fcin(2, :)
    write(fwork,*)
    
  else
    write(fwork,*) 'Fail to refine! nf=', nf
    exit
  end if
! ------------------------------------------------------------

! NOTE: tried simultaneous refinement here, NO good, because the frequecies is not accurate 
!       will produce bigger error in frequecies---- discard 

! 3.  If we have more than 1 frequencies, refine again with all the exiting frequencies
!     using the orginal samples for freq+coef simultaneously 
  if( nf >  1) then  
    ! for dgelsd solver, we have to pass the input array with exact dimension as the input(nr,nc)
    allocate( fcinall(nf+1,5), stat = AllocateStatus)
    if (AllocateStatus /= 0)  stop "*** Not enough memory ***"
    
    ! all the detected frequecies are to be refined
    fcinall   =  fc_refn(1:nf+1, 1:4) 
      
    ! refine approach: 3  --- freq+coef 
    call refine( fcinall, nf, din,  3, maxferr, isref)
    
    if( isref == 1 ) then
    ! update the refined freq+coef
      fc_refn(1:nf+1, 1:4)  = fcinall(1:nf+1,:) 
     
      write(fwork,*)
      write(fwork,*)  'refinement of all the',nf, 'detected frequecies'
      do j = 2, nf+1
         write(fwork,*) fc_refn(j, :)
      enddo 
      write(fwork,*)  
      
    else
      nf = nf-1
      write(fwork,*)  'Fail to refine! nf=', nf
      exit
    end if
    deallocate(fcinall)  
  endif
    
! 4.  Substract all the refined frequency from the original data. 
! Start again from step 1 with the residual instead of the original sample 
!  dwork here is updated by the residual
  dwork = din
  
  write(fwork,*) 'Before substaction'
  do j = 2, nf+1
    write(fwork,*)  fc_refn(j, :)
  enddo 
  write(fwork,*)  
  
! allocate array for the nf+1 frequecies   
  allocate( fcinall(nf+1,4), stat = AllocateStatus)
  if (AllocateStatus /= 0)  stop "*** Not enough memory ***"
    
  ! all the detected frequecies are to be refined
  fcinall   =  fc_refn(1:nf+1, 1:4)  
  call substr( dwork, np,  nf,  fcinall)
  
  deallocate(fcinall) 
  
  ! if the maximum of the residual is less than tol, stop the refinement
  resmax = maxval( dabs(dwork))
  write(fwork,*)   'Maximal residual', resmax;  write(fwork,*)
  read*
  
  if(resmax < tolres) then 
    write(fwork,*)  'Maximal residual is less than tolres, resmax=', resmax
    exit
  endif 
  
  if(nf > nfmax) then 
    write(fwork,*) 'Number of freqs exceeds the nfmax: ', nf
    exit
  endif  
enddo  

! update the constant term... in fact, no need to compute the average prom... 
! anyhow, do it this way.... 
fc_refn(1, 2) = fc_refn(1, 2) + prom 

! ---- write to the file the detected frequencies ....
write(ffc, *) ' #  detect ', nf, ' frequecies! '
write(ffc, *) ' #  Maximal modulus of residual = ', resmax
write(ffc, *) ' #  feq    cc    ss    amp      order'
write(ffc, *)

nff = min0(nf, nfmax) 
do i = 1, nff + 1 
  write(ffc, '(4f26.10, f10.0)') fc_refn(i, :)
enddo   
write(ffc, *)  ! add a blank line 

return  
end subroutine furian

!***********************************************************************
!     *****  Substr  *****

! Substract the refined frequencies from the orginal data din(or the residual from last the substaction of the previous detected frequecies) 
! TODO: how to express the signal by a given frequency???
 
!       Input Variables 
!  a        the input array to work with( original data or residual)            
!  np       number of samples 
!  nf       number of detected frequecies
!  fc       

!       Output Variables 
!  a        the updated residual after the substaction          

!  Routine Used:
!     None   

!  Finally Revised by Yu -- 20160815
!***********************************************************************
subroutine substr( a, np, nf, fc)

implicit none

!  Input  and Output Declaration   
integer, intent(in)     ::  np, nf  
real(kind=dp), intent(in)         ::  fc(nf+1, 4)  
real(kind=dp), intent(inout)      ::  a(np) 
 
!  Local Variable
integer :: i, j 
real(kind=dp)  :: t, freq, cc, ss, g, arg
   
  do i = 1, np  
  ! the epoch of  i-th sample, equispaced np points in the interval [0, lt] 
    t = i * lt / np 

  !  the constant term, cc0, saved as the second componen in the first line in fc   
    g = fc(1, 2) 
    
  !  Sigma for k = 1,...,nf of {  cc(k) * cos(2*pi*t * freq(k) ) +  ss(k) * sin( 2*pi*t * freq(k) )  }
    
    do j = 2, nf + 1   
      
      freq = fc(j, 1) ! the evaluation frequency 
      cc = fc(j, 2) 
      ss = fc(j, 3) 
      
      ! the term related to freq: cc * cos(2pi*t * freq) + ss * sin(2pi*t * freq)
      arg = pi2 * t   
      g    = g + cc * dcos(arg) + ss * dsin(arg)
    enddo  
    
    a(i) = a(i) - g
  end do  

  return  
end subroutine substr


subroutine refine( f, nf,  din, idref,  maxferr, isref)
! Thie routine is to refine the frequencieies and amplitude simultaneously detected by fai

! TODO: sines and cosine, the one by mod(t, 2Pi) is not accurate enough. 

!       Input Variables 
!  f            all the frequencieies to be defined, dimesion nf-by-4
!               format of the data structure, with the first row as the constant term: 
!               freq -  ccos -- ssine -- amp  
!  nf           number of  frequencieies
! idref         flag of refinement approach: 
!               1.  freq ;    2: coef(cc+ss);   3: freq+coef 
      
!       Output Variables 
!  f            only the first column is the updated by the refined frequencieies   
!  maxferr       maximum difference between the input samples and the computed
!               approximation, an estimate of the detection by FFTSC+REFINE
           

!     Module-base Varaibles
!  np           number of samples
!  lt           length of the time interval of the samples
!               the epoch of samples:  t = i*lt/np, i = 1, np 
!  tol          tolerance( in maximum error) for Newton's method for refinement of more than 1 frequecies, 
!               1.e-6 would be enough? the algorithm is stopped if maxferr is smaller than tol  

!  tolf        tolerance for the correction(in frequencies) of Newton
!              Method for refinement for only 1 frequecy, i.e., 2.d-5


!     Local Varaibles
!  dif         the difference between the approximated value and the real one, of dimension np-by-nrhs, in this case, it is 1 
!              so we fix the value nrhs == 1
!              Note: the value varies between different problems.
 
!  Routine Used:
!    dnrm2, deltx

! Finally Revised by Yu -- 20160607
!----------------------------------------

  implicit none

! Input  and Output Declaration   
  integer, intent(in)      ::  idref ! , np
  integer, intent(out)     ::  isref
  integer, intent(inout)   ::  nf 
  real(kind=dp), intent(in)       ::  din(np) !, tol, lt, tolf
  real(kind=dp), intent(out)      ::  maxferr
  real(kind=dp), intent(inout)    ::  f(nf+1, 4) 
  
! Local Variable
  real(kind=dp)  ::  g(np), dg(np, idref*nf), ferr(np,1),  df(idref*nf), &
                     dfm_pre, dfm 
  integer        ::  iter, i, j, info, inccor, nrhs
  
  real(kind=dp), external :: dnrm2 
  
  isref = 1
  
  nrhs = 1  ! only one column in right-hand side 
   
! refine iteratively 
  iter = 0        ! index of iteration 
  inccor = 0       ! number of time when the modulus correction increases

  do 
    iter = iter + 1 
    
    if(iter > nitmax) then 
      write(fwork,*) 'Too many iterations! ', iter
      isref = 0 
      return
    end if 
    
     write(fwork,*) '***********', iter, 't-th iteration ***************'
  ! the computed approximation g and  differential wrt the frequencies 
    call refine_gdg( f, nf, idref, g, dg)
    
    ! the difference between the original data and the computed approximation
    ferr(:,1) = din - g  
    
!   because of the large amount of samples 
!   we only need to check maximum the error  
    maxferr =  maxval( dabs(ferr) )
 
!-TODO: do we need to do the check on the modulus of the residual here? 
!       Maybe NO! because if all the domimant frequecies are not detected, this vaulue 
!       will be quite large, otherwise, we want to refine more accurately.
!    if(maxferr < tolres) then 
!       write(fwork,*) 'Maximum error in g is less than tolres',  maxferr 
!      return
!    endif
    if( iter > 1)  dfm_pre = dfm
    call deltx( np, idref*nf, nrhs,  dg,  ferr, df, info)

    if( info /= 0 ) then 
    
      write(fwork,*)  'SSL Solver Failed to converge!';  read*
      isref = 0
      return
      
    else   
    ! check the modulus of  correction, if it is smaller than tol terminate the iteration 
    ! TODO -- although it works for this example, I think it is better to just check the 
    ! correction in frequecies rather than freq+coef
 
 ! for the termination of the iteration, one criteria could be the change trend of the modulus of the correction, if it fluctuate, no good, stop
      dfm  =  dnrm2(idref*nf, df, 1)
      
      ! TODO --- for sure we should stop when |df| starts increasing
      !          As Gerard suggests, it is possible for the error to increase if the number of frequecies increase and we do the refinement simultaneously
      !          better to check the second increase of error.... 
      
      
      if(iter > 1 .and. dfm > dfm_pre) then 
      
        inccor = inccor + 1 
        
        if(inccor >= 3 ) then ! allow the correction to increase no more than twice... 
          isref = 0
          print*, 'No more correction! |df| starts increasing! '
          read*
          return 
        endif 
          
      endif  
      
!      print*, 'dfm=',dfm, 'tol=', tol; read* !ck
      if(dfm < tol) then 
        write(fwork,*) 'Modulus of correction less than tol!', dfm 
        return
      endif  
    
      ! print the correction to check
      write(fwork,*) 'Old:  freq     cos     sin     amp'
      do j = 1, nf 
        write(fwork,*) f(j+1, : )
      enddo
      
      write(fwork,*)  
      if(idref == 1 ) write(fwork,*) 'Correction: freq     cos     sin'
      if(idref == 2 ) write(fwork,*) 'Correction: cos      sin'
      if(idref == 3 ) write(fwork,*)'Correction: freq     cos     sin'
             
      do j = 1, nf 
         write(fwork,*) df(j: idref*nf: nf) 
      enddo 
        
      ! -- update the frequencies and coefficients
      do j = 1, nf 
        ! freq
        if(idref == 1 .or. idref == 3 ) f(j+1, 1) = f(j+1, 1) + df(j)
        if(idref == 1)  cycle
        
        ! coef + amp
        if(idref == 2) then 
          f(j+1, 2) = f(j+1, 2) + df(j)
          f(j+1, 3) = f(j+1, 3) + df(nf+j)
          
        elseif ( idref == 3 ) then
          f(j+1, 2) = f(j+1, 2) + df(nf+j)
          f(j+1, 3) = f(j+1, 3) + df(2*nf+j)
        endif  
        
        f(j+1, 4) = dsqrt(  f(j+1, 2)**2 + f(j+1, 3)**2 )
        
      enddo 
       
     !check the modulus of the last frequency, which has the smallest amplitude
      if( f(nf+1,4) < tolres .and. nf > 1) then 
        print*,  nf,'-th freq has amplitude less than tolres: ', f(nf+1,4)
        isref = 0
!        nf = nf-1
        return
      endif 
      
      write(fwork,*) 
      write(fwork,*) 'New:  freq   cos  sin  amp'
      do j = 1, nf 
        write(fwork,*), f(j+1, : )
      enddo
      write(fwork,*)   
  
    endif
   
  enddo    

  return  
end subroutine refine

!*********************************************************************** 
!  Compute the diffrential  dg, which is the derivative of the freq wrt (w,cc,ss)
!  and the difference between the sample data and the approximated value 
!  of function recoveried from  freq-cc-ss 

!       Input Variables 
!  f            all the frequencieies to be defined, nf-by-4
!               freq -  ccos -- ssine -- amp  
! idref         flag of refinement approach:  (see refine) 
!               1.  freq ;    2: coef(cc+ss);   3: freq+coef 
!  nf           number of  frequencieies

!       Output Variables 
!  g            the approximate value of dat 
! dg            the differential of g wrt freq            

!        Module-Based Varaibles
!  np           number of samples
!  lt           length of the time interval of the samples
!               the epoch of samples:  t = i*lt/np, i = 1, np  
! pi2 

!  Routine Used:
!     None

!  Finally Revised by Yu -- 20160607
!----------------------------------------
subroutine refine_gdg( f, nf, idref, g, dg)

implicit none

! Input  and Output Declaration   
integer, intent(in)          ::  nf, idref  
real(kind=dp), intent(in)    ::  f(nf+1, 5)  
real(kind=dp), intent(out)   ::  g(np), dg(np, idref*nf) 
 
! Local Variable
integer        :: i, j
real(kind=dp)  :: arg,  t, freq, cc, ss, dt, cn(np), sn(np), darg

  
! time step size for sampling  
  dt = pi2* lt/np
  
  do i = 1, np, 1
    g(i)  =  f(1, 2)
  end do

!  ----- deal with frequecy one by one  --------- 
  do j = 2, nf+1   
    
      freq = f(j, 1)
      cc = f(j, 2) 
      ss = f(j, 3) 
      
! the sum from 1 to nf of  cc*cos(2*pi*w * i*dt_step) +  ss*sin(2*pi*w* i*dt_step)  = cc*cn(i) + ss*sn(i) 
!  if we take darg = w * 2*pi* dt_step, where w is the current frequency 
      
   ! TODO:  mod approach is not reliable, we replace with recurrent approach of example 4  
   !      arg  = dmod( pi2* freq *t, pi2)
   !      g(i) = g(i) + cc * dcos(arg) + ss * dsin(arg)
      
      
   ! step size for sampling 
    darg =  freq * dt
  
  !  ---  The Fourier Series with the first np terms ----------
  !  ---Compute cos( n * darg ) and sin( n * darg )  
    call trigrec_c( darg, np, cn)  
    call trigrec_s( darg, np, sn)
  
    ! check trigrec_c. trigrec_s 
    
    ! approximate value g recoveried from the detected frequencieies and amplitudes
    ! TODO : Another approach is to add the constant term after all the other terms, advoid too many addtions 
    
    do i = 1, np  
      !  the constant term f(1, 2) ,  f:  freq-cs   
      g(i) =  g(i) +  cc * cn(i) + ss * sn(i)
      
      ! dg =  d g / d w (nf),  d g / d cc (nf), d g / d ss (nf)
      
      ! derivative of g w.r.t freq 
      !  dg / d w(i) = the sum from 1 to nf (cc * -sn(i) * 2pi*dt_step * i  + ss * cn(i) * 2pi*dt_step * i) 
      !                                     == -cc * sn(i) * dt * i + ss * cn(i) * dt * i
      
      if(idref == 1 .or. idref == 3 ) &
!         dg(i, j-1) = -cc*pi2*t * dsin(arg) + ss*pi2*t * dcos(arg)
         dg(i, j-1) = -cc* sn(i) * dt * i + ss * cn(i) * dt * i
       
      ! if we only need to refine frequecy  
      if(idref == 1) cycle 
      
      ! coef 
      if(idref == 2 ) then 
        dg(i, j-1)    = cn(i) ! w.r.t cc 
        dg(i, nf+j-1) = sn(i) ! w.r.t ss
       
      ! for both freq + coef, we have 1-nf: dg/dw, nf+1-2nf: dg/cc,  2nf+1-3nf: dg/ss  
      elseif(idref == 3) then
      !
        dg(i, nf+j-1)   =  cn(i) ! w.r.t cc 
        dg(i, 2*nf+j-1) =  sn(i) ! w.r.t ss
      endif 
          
    enddo  
  
  end do  

  return 
   
end subroutine refine_gdg

 
! *****************************************************************************
!     ****   Hanning Window function   ****

! To take in consistent with the time epoch of the samples, which are 
!   i * dt, for i = 1, ..., NP, with  dt = LT/NP
! the cos(n*x) here should also start with i = 1, instead of i=0 in the original function 

! The same as in JMM's paper:

!  H_nh = ( 1 - cos( 2*pi * t/LT) )**nh * qnh 
!       = ( 1 - cos( i * 2*pi/NP) )**nh * qnh
!  qnh = nh!/ (2*nh-1)!!

! and cos(n*x) is computed in a recurrent way:
!    cos(n*x)  = 2*cos(x) * cos( (n-1)*x ) - cos( (n-2)*x ), for n>=3
!    here, x = del = 2*pi / NP  

!  Module-based Varaibles: nh, pi2, dp
!*****************************************************************************

subroutine hannin(a, np)
implicit none 

integer, intent(in)     ::  np
real(kind=dp), intent(inout)        ::  a(np)

! Local Varaibles
integer ::  i
real(kind=dp) :: dt, cn(np), qnh !, & 
!                 a1(np), c0, c1, c2, z ! to check if cosnx works?
  
  if(nh == 0)   return 
  if(nh == 1)   qnh = 1
  if(nh == 2)   qnh = 2.d0/3
  if(nh == 3)   qnh = .4d0

  dt = pi2 / np  ! x
  
  call cosnx(dt, np, cn)
  
  ! ckd, the same as the recurrence 
  do i = 1, np 
    a(i) = a(i) * ( 1-cn(i) )**nh * qnh
  enddo
  
  return
end subroutine hannin

!***********************************************************************
!     ****   cosnx   ****   --- example 2 --- disard ---
!  Compute cos(n*x) in a recurrence way, and return all the cos(i*x) for 
!  i = 1, ..., n, we only care when n>=2 

! --Example 2, the worst one shown in the book
!    cos(n*x)  = 2*cos(x) * cos( (n-1)*x ) - cos( (n-2)*x ) 
!    cos(0) = 1

!       Input Variables 
!  x        the base stepsize         
!  n        number of points to be computed
 
!       Output Variables 
!  cn       dimension n,  cos(i*x) for i = 1, ..., n            

!  Routine Used:
!     None

!  Finally Revised by Yu -- 20160719 --TODO 
!***********************************************************************
subroutine cosnx( x, n, cn)

implicit none
integer, parameter :: dp = kind(1.d0)   

!  Input  and Output Declaration   
integer, intent(in)     ::  n
real(kind=dp), intent(in)      ::  x  
real(kind=dp), intent(out)     ::  cn(n) 
 
!  Local Variable
integer :: i
real(kind=dp)  :: c0, z 

  c0 = 1.d0          ! cos(0)
  cn(1) = dcos(x)    ! cos(x)
  z  = 2.d0 * cn(1)  ! 2cos(x) -- constantly used term
  
  do i = 2, n
    cn(i) = z*cn(i-1) - c0
    if(n == 2) return
    
    c0 = cn(i-1)
  enddo 
        
  return  
end subroutine cosnx

!***********************************************************************
!     ****   sinnx   **** --TODO  --- example 2 -- disard
!  Compute sin(n*x) in a recurrence way, and return all the sin(i*x) for 
!  i = 1, ..., n, we only care when n>=2 

!    sin(n*x)  = 2*cos(x) * sin( (n-1)*x ) - sin( (n-2)*x )
!    sin(0) = 0

!       Input Variables 
!  x        the base stepsize         
!  n        number of points to be computed
 
!       Output Variables 
!  sn       dimension n,  sin(i*x) for i = 1, ..., n            

!  Routine Used:    None

!  Finally Revised by Yu -- 20160607 --TODO: unckd...
!***********************************************************************
subroutine sinnx( x, n, sn)

implicit none
integer, parameter :: dp = kind(1.d0)   

!  Input  and Output Declaration   
integer, intent(in)     ::  n
real(kind=dp), intent(in)      ::  x  
real(kind=dp), intent(out)     ::  sn(n) 
 
!  Local Variable
integer :: i
real(kind=dp)  :: s0, z
  
  s0 = 0.d0            ! sin(0)
  sn(1) = dsin(x)      ! sin(x)
  z  = 2.d0 * dcos(x)  ! 2*cos(x) -- constantly used term
  
  do i = 2, n
    sn(i) = z*sn(i-1) - s0
    if(n == 2) return
    
    s0 = sn(i-1)
  enddo 
        
  return  
end subroutine sinnx




!*******************************************************************************
!  This subroutine is to compute the residual (deviation from the mean),
!  which is the difference between the observed data and its mean 
!  It has not been used....

!*******************************************************************************
subroutine promig(ain,np, prom)
implicit none 

integer, intent(in)     ::  np  
real(kind=dp), intent(inout)       ::  ain(np)   
real(kind=dp), intent(out)      ::  prom 

! Local Variables
real(kind=dp) ::  suma 
      
! the sum of all the components of the input array ain
  suma = sum(ain)
  
 ! the mean of sum(a)
  prom = suma / np
  write(*,*)' average =', prom
      
  ! deviation from the mean, also called residual
  ain = ain - prom 
  
  return
end subroutine promig


end module fft_mod
