!  **** My fft Module*******
! TODO : replace the sines and cosine by the recurrent procedure done by Gerard and me in 20160719 

! TODO: for the points, tak the interval (0, pi2], instead of [0, pi2) 
!       the first one seems better... 

! TODO: 2016-10-06 12:18:40 
!  substract the constant term at the very beginning to work with the residual
!  and remember to add it back on the constant term for deteted frequecies
  
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

! subroutines from private libraries claimed external 
external  ::   trigrec_c, trigrec_s   ! -- libgr.a  deltx,


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
integer, private :: fwork, fout, debug 



! ---Module private variables--- must be initialization before any call of this module 
! tolres:      tolerance of the maximum amplitude of the residual for frequency detecton 
!            
! tol        tolerance for the correction in Newton Method for the refinement   
!  np        number of samples
!  lt        length of the time interval of the samples
!            the epoch of samples:  t = i*lt/np, i = 1, np  

contains

!  !***********************************************************************
!     ****   init_fft  ****
!  Initialize the commonly shared variables
!-- by input
!   np, din, tol, tolres  
!     -- the last two suggestion value: 1.e-8, 1.e-4

!-- by default
!   nh   

!***********************************************************************
subroutine init_fft( np0, lt0, nfmax0, nitmax0)
implicit none

! Input  and Output Declaration  
integer, intent(in)       ::  np0, nitmax0, nfmax0 
real(kind=dp), intent(in) ::  lt0 !,  tol0, tolres0   
 
! Local Variable
integer ::  istol
  

! debug flag  
  debug = 0
  
! ----- default values      
  ! order of Hanning window function     
  nh = 2              

!  minimun difference in order of two detected frequecies 
  sepmin = nh+1 

! ----  precision control  
  print*, 'Maximum amplitude of the residual // Tolerance in correction for refinement? '
  print*, '( 1 : Default: 1.d-10, 1.d-10;  Otherwisze input new values)'
  read(*,*) istol

  if(istol == 1) then 
    tolres = 1.d-10; tol = 1.d-10
  else 
    read(*,*)  tolres,  tol
  endif 
  
! ----- by assignment in the callee
  lt     = lt0
  nfmax  = nfmax0
  nitmax = nitmax0
  
!    tolres = tolres0 
!    tol    = tol0
! ---- commonly used value for dimension declaration  
  np  = np0
  npm = int( dlog( dble(np) ) / dlog(2.d0) ) ! log_2(np)
  npd2 = np / 2 
  npd2p1 = npd2 + 1
  
  if(debug == 1) then 
    print*, 'Check init_fft,  lt, nfmax, np: ',  lt, nfmax, np; read*
  endif 
   
   ! save(/=6) or display(6) the intermediate results, for check purpose 

!   fwork = 6     ! print on screen
   
  fwork = 222   ! for save 
  open(fwork, file='fur.work')
  
  
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
!   tol, tolf, dp, nfmax  

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
integer       ::  AllocateStatus, isref, ind_max, i, j, sepi, sepimin, feqnew, nf1, nff 
real(kind=dp) ::  maxferr,  fcin(2,4), fc_ava(npd2-1) , prom
real(kind=dp), allocatable :: fcinall(:,:)   

!  print*, 'lt, nfmax, np: ',  lt, nfmax, np; read*  !ckd

  ! compute the residual from the average  of din, Note that we work with the residual all the time 
  ! TODO this has none-sense, at least I didn't find     any point in computing the average, since we do not need to substract it. ???
  call promig(din, np, prom) 

!  print*, 'prom=', prom ; read*

  ! working array dwork, do not touch din, keep din as constant
  dwork= din 

  !--- initialize the constant term 

  ! number of extracted frequecies, set to 0
  nf = 0

  ! the distance of the new detected frequncy from the previously refined ones
  sep = npd2 
 
  ! the minimum distance between 2 consecutive frequecies 
  sepimin = npd2  

  fc_refn(1,:) = 0.d0 
  
  do   
! 1.  Given equispaced input samples, first multiply the Hanning window
!     function, and then do simple FFT with the hanning-flitered data to
!     aviod leakage   

! compute the mean and work with the residual(the deviation from the mean for FFT)

! or work directly with fft? -- choose this option, easier, but it performs bad  with non-zero constant term? why? TODO 20161006
  
    feqnew = 1
    dwork0 = dwork 
  

!      print*, 'check input data: a'
!      write(*,'(8f18.12)') dwork(1:32)
!      print*; read*

     ! TODO 2016-10-09 20:29:37 :
     ! hannin--- recover the starting point as 1 instead of 0  
     
    call hannin(dwork, np) 
          
!      print*, 'check hanning-flitered data:'
!      write(*,'(8f18.12)') dwork(1:32) 
!      print*; read*
      
    call fftsc(dwork, np, st, ct, iwk, wk, cwk) 
    call modulus_fftsc(st, ct, np, lt, fc_fft)
  
      ! do sinine and cosine Fast Fourier analysis, compare with other furian routines
!      print*, 'check fftsc: st'
!      write(*,'(8f18.12)')  st(1:32)
!      print*; read*
!      
!      print*, 'check fftsc: ct'
!      write(*,'(8f18.12)')  ct(1:32)
!      print*; read*
!      
! 2.  Detect only one frequency with the maximal amplitude, and use Newton
!     Method to refine (freq+coef) simultaneously.

!     To avoid leakage,  when a new frequecy is detected, we check 
!     the distance  in order of the two frequecies should be >= sepmin
!     sepmin = 1+nh(order of Hanning function), see JM paper for more detail
! Work with the one has maximum magnitude, fc_fft(ind_max+1, :) for refinement 
  
    if(nf == 0 ) then 
    ! for the first freq, start with the one with maximal amplitude
  
      fc_ava  = fc_fft(2:npd2, 4) ! the amplitude 
      ind_max = maxloc(fc_ava,1)  ! the location of the maximal amplitude
  
    else    
  
    ! For the second freq and on, check  the distance  to avoid aliasing
    
      do 
        feqnew = 1 ! flag to show the detection of new frequecy, by default yes
      
        ! the available feq+coef (delete the aliased one by setting amp=0)
        fc_ava  = fc_fft(2:npd2, 4) 
        ind_max = maxloc(fc_ava, 1) ! the location of the maximal amplitude

        if(debug == 1) then 
          write(fwork,*) 'ind_max =',  ind_max 
          write(fwork,*)  fc_fft(ind_max+1, :);  write(fwork,*) 
        endif 
        
        
        do i = 2, nf+1 
          sepi = iabs( ind_max - int(fc_refn(i,5)) )
!         print*, 'sepi, ind_max', sepi, ind_max; read* !ck
        
          if( sepi <= sepmin ) then 
          ! if the frequecy is already refined, but still not accurate enough, 
          ! we forcely reset the amplitude to be 0, and find the frequecy
          ! with the second maximal amplitude 
            
             if( debug == 1 ) then
               write(fwork,*) 'Aliased Frequency!'
               write(fwork,*) 'Detected:',  fc_refn(i,:)
               write(fwork,*) 'Aliased: ',  fc_fft(ind_max+1,:) 
               write(fwork,*) ; read*
             end if
           
             fc_fft(ind_max+1, 4) = 0.d0  
             feqnew = 0
             exit
          endif
          
          ! the minimum distance between two frequencies 
          if(feqnew == 1 .and. sepi < sepimin )  sepimin = sepi
           
        enddo 
      
        ! check the amplitude of the new frequecies, terminate the refinement
        if(fc_ava(ind_max) < tolres) then 
          if(debug==1) write(fwork,*)  'Amplitude of new freq in residual is less than tolres', fc_ava(ind_max-1)
          feqnew = 0
          exit   
        endif 
      
        if(feqnew == 0) cycle
        
        ! succeed a new frequecy
        if(sepimin < sep)  sep = sepimin
        exit 
      enddo 
    
    endif 
    
    if(feqnew == 0)   exit 
    
    
! now we obtain the avaible new frequency, which has the maximal amplitude in the residual  
    if( debug==1 ) then
      write(fwork,*),'The maximum amplitude'
      write(fwork,*), fc_fft(ind_max+1, :)
      write(fwork,*)
    end if
    
      
    nf = nf + 1
  
!   update the constant term 
!    fc_refn(1,2)   =   fc_fft(1,2) ! cc 
    fc_refn(1,2)  = fc_refn(1, 2) + fc_fft(1, 2) ! cc 
    fc_refn(1,4)  = dabs( fc_refn(1,2) ) ! amp
  
    if(debug == 1) print*, 'nf, The costant term:', nf,  fc_refn(1,2); !read*!ckd
  
    ! save the order in 5-th column to check the distance 
    fc_refn(nf+1, 5)   = ind_max  
    fc_refn(nf+1, 1:4) = fc_fft(ind_max+1, :)

! ---- for refinment of only 1 frequecy in coefficients ---- 
    fcin(1,:)   = fc_fft(1, :)     ! the costant term 
    fcin(2,:)   = fc_refn(nf+1, 1:4) 
   
! For only 1 frequecy, we refine only the coefficients and fix the frequecy
! refine approach: 2
! ** NOTE ** ---- discard approach 3 for 1 frequency 
! We have tried simultaneous refinement(freq+coef) for 1 frequecy, Disaster  ! 
! because the frequecies is not accurate, the correction will produce larger error in frequecy

    nf1 = 1
!  call refine( fcin, nf1,  dwork0, 2, maxferr, isref)  ! no big difference...
    call refine( fcin, nf1,  din, 2, maxferr, isref)    ! this is better? 
    
    ! Approach 3,  disaste for only 1 frequecy !!!
!  call refine( fcin, 1,  din, 3, maxferr, isref) 
  
    if( isref /= 0 ) then
    ! update the refined coef
      fc_refn(nf+1, 1:4) = fcin(2, :)
      
      if(debug == 1 ) then 
        write(fwork,*) nf, '-th new frequecy'
        write(fwork,*) fcin(2, :)
        write(fwork,*)
      endif 
      
    else
      
      if(debug == 1)  write(fwork,*) 'Fail to refine! nf=', nf
        
      exit
    end if
! ------------------------------------------------------------


! 3.  If we have more than 1 frequencies, refine again with all the exiting frequencies
!     using the orginal samples din for freq+coef simultaneously 
    if( nf >  1) then  
  
      ! Allocate input array of exact dimension (nr X nc) for dgelsd solver 
      allocate( fcinall(nf+1,5), stat = AllocateStatus)
      if (AllocateStatus /= 0)  stop "*** Not enough memory ***"
    
      ! all the detected frequecies are to be refined
      fcinall   =  fc_refn(1:nf+1, 1:4) 

      ! refine approach: 3  --- freq+coef 
      call refine( fcinall, nf, din,  3, maxferr, isref)
      
      ! TODO: check the modulus of the last freq? or put isref=2? 
      
      if( isref /= 0 ) then
      ! update the refined freq+coef
        fc_refn(1:nf+1, 1:4)  = fcinall(1:nf+1,:) 
     
        if( debug == 1 ) then
          write(fwork,*)
          write(fwork,*)  'refinement of all the',nf, 'detected frequecies'
          do j = 2, nf+1
             write(fwork,*) fc_refn(j, :)
          enddo 
          write(fwork,*) 
        end if
         
      
      else
        nf = nf-1
        write(fwork,*)  'Fail to refine! nf=', nf
        exit
      end if
      deallocate(fcinall)  
    endif
    
! -- 4.  Substract all the refined frequencies from the original data. 
!   Start again from step 1 with the residual instead of the original sample 
!   dwork here is updated by the residual
    dwork = din

  ! allocate array for the nf+1 frequecies   
    allocate( fcinall(nf+1,4), stat = AllocateStatus)
    if (AllocateStatus /= 0)  stop "*** Not enough memory ***"
    
    ! all the detected frequecies are to be refined
    fcinall   =  fc_refn(1:nf+1, 1:4)  
  
    ! substract the refined frequecies 
    call substr( dwork, np,  nf,  fcinall)
  
    deallocate(fcinall) 
  
  ! if the maximum of the residual is less than tol, stop the refinement
    resmax = maxval( dabs(dwork))
    write(fwork,*)   'Maximal residual', resmax;  write(fwork,*)
    
    ! Good place to pause for check 
    if(debug == 1)     read*
  
    if(resmax < tolres) then 
      write(fwork,*)  'Maximal residual is less than tolres, resmax=', resmax
      exit
    endif 
  
    if(nf > nfmax) then 
      write(fwork,*) 'Number of freqs exceeds the nfmax: ', nf
      read* 
      exit
    endif  
  enddo  

! Update the constant term, which is save in the first row column 2
! in fact, I think no need to compute the average prom... 
! anyhow, do it this way... TODO 
  fc_refn(1, 2) = fc_refn(1, 2) + prom 
  fc_refn(1, 4) = dabs( fc_refn(1, 2) )
! ---- write to the file the detected frequencies ....
  write(ffc, *) ' #  detect ', nf, ' frequecies! '
  write(ffc, *) ' #  Maximal modulus of residual = ', resmax
  write(ffc, *) ' #  feq    cc    ss    amp      order'
  write(ffc, *)

!  only keep at most nfmax frequencies
  nff = min0(nf, nfmax) 
  do i = 1, nff + 1 
    write(ffc, '(4f26.16, f10.0)') fc_refn(i, :)
  enddo   
  write(ffc, *)  ! add a blank line 

  return  
end subroutine furian

!***********************************************************************
!     *****  fun_four  *****
!  Evaluate the fourier representation at np equispaced points in the time interval [0, LT]
!  given the nf frequencies and the coorsponding coefficients, usually for quasi-periodic functions  

!  f(t) = cc(1) + Sigma for k=1,...,nf of { cc(k+1) * cos( 2pi*freq(k) * t) +  ss(k+1) * sin( 2pi*freq(k) * t) }
!         where t = 0, LT/NP, ..., LT*(NP-1)/NP

! TODO: this subroutine is not used yet, may replace 'substr' in the near future -20160816

!  **NOTE**  
!  -- 1. the idea is not the same as in Fourier series for periodic curve 
!  -- 2. compute trigometric functions cos(nx) and sin(nx) recurrent by trigrec_c and trigrec_s 

!       Input Variables                                                                                          
!  np       number of samples 
!  nf       number of detected frequecies
!  fc       the known freq + coefficients, dimension nf+1-by-4: freq - cc - ss - amp 
!           the first line is the constant term 
!           use this array is more convenient for the fft_mod functions 

!       Output Variables 
!  fun      dimension np, the evaluated values at np equispaced points          

!  Routine Used:
!     None   

!  Finally Revised by Yu -- 20160816
!***********************************************************************
subroutine fun_four( fc, nf, lt, np, fun)

implicit none

!  Input  and Output Declaration   
integer, intent(in)     ::  nf, np  
real(kind=dp), intent(in)       ::  fc(nf+1, 4), lt 
real(kind=dp), intent(out)      ::  fun(np) 
 
!  Local Variable
integer :: i, k 
real(kind=dp)  ::  freq, cc, ss,  arg, cn(np), sn(np)
  
  ! the constant term, which is saved as the second component in the first line in fc   
  fun = fc(1, 2)  ! array operation, works on every component 
  
  do k = 2, nf + 1  
    ! -------  the evaluated frequence ---------
    freq = fc(k, 1) 
    cc   = fc(k, 2) 
    ss   = fc(k, 3) 
    
    ! basic angle for each frequecy
    arg = pi2 * lt / np * freq  ! arg = 2pi*t *freq, where t = i*lt/np 
   
    ! the epoch of  i-th sample:  t = (i-1) * lt / np,  equispaced np points in the interval [0, lt] 
    
    ! for the first point P1, t=0, we have cos(0) = 1.d0, sin(0) = 0.d0
    ! but cn(sn) start from  i = 1, so we need to handle i=1 seperately
    fun(1) = fun(1) -  cc 
    
    ! the trigometric function for each point, i = 2, np -th point
    call trigrec_c( arg, np, cn)  
    call trigrec_s( arg, np, sn)
    
    ! Be careful: cos(i*arg) refers to i+1 th point,  since for i+1 -th point, we have time t = (i+1 -1)*lt/np 
    do i = 1, np-1   

      ! the term related to freq: cc * cos(2pi*t * freq) + ss * sin(2pi*t * freq)
      fun(i+1) = fun(i+1) - cc * cn(i) - ss * sn(i)
      
    enddo  
  end do  

  return  
end subroutine fun_four

! TODO: check the mod_fftsc, and see how we define the first frequency
! current debug position!!!!!!!!

!***********************************************************************
!     *****  Substr  *****

! Substract the refined frequencies from the orginal data din (or the residual from last the substaction of the previous detected frequecies) 
! Using recurrent equation to compute trigometric functions cos(nx) and sin(nx)

!  The Fourier series can be written as a function of t, given the refined frequecies :

!    f(t) = cc(1) + Sigma for k=1,...,nf of { cc(k+1) * cos( 2pi*freq(k) * t) +  ss(k+1) * sin( 2pi*freq(k) * t) }

!     where t is the epoch of the sample point within the interval [0, LT]

! ** NOTE ** be default, we have  t = 0, LT/NP, ..., LT*(NP-1)/NP

!  SO if we want to express Fourier series in terms of the argument of the evaluated point, of which the value range is [0,2pi)
!    we have, at time t the argument theta = 2*pi*t, and, 
 
!    f(theta) = cc(1) + Sigma for k=1,...,nf of { cc(k+1) * cos( theta * freq(k) ) +  ss(k+1) * sin( theta * freq(k)) }        
 
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
integer :: i, k 
real(kind=dp)  :: darg, freq, cc, ss, arg, cn(np), sn(np) ! t,  g,
  
  ! first, substract the constant term, which is saved as the second componen in the first line in fc   
  
  ! TODO: be careful the sample start with t
  a   =  a -  fc(1, 2)  ! array operation, works on every component 
  darg = pi2 * lt / np
  
  do k = 2, nf + 1  
    ! -------  the evaluated frequency  ---------
    freq = fc(k, 1) 
    cc   = fc(k, 2) 
    ss   = fc(k, 3) 
    
    ! basic angle for each frequecy
    arg =  darg  *  freq  ! arg = 2pi*t *freq, where t = i*lt/np 
    ! the trigometric function for each point, i = 2, np -th point
    call trigrec_c( arg, np, cn)  
    call trigrec_s( arg, np, sn)
    
    ! the epoch of  i-th sample:  t = (i-1) * lt / np,  equispaced np points in the interval [0, lt)
    ! so arg = 2pi* (i-1) * lt / np * freq = (i-1) *  2pi*lt/np *freq 
    
    ! for the first point P1, t=0, we have cos(0) = 1.d0, sin(0) = 0.d0
    ! but cn(sn) start from  i = 1, so we need to handle i=1 seperately
!    a(1) = a(1) -  cc 
!    
!    ! Be careful: cos(i*arg) refers to i+1 th point,  since for i+1 -th point, we have time t = (i+1 -1)*lt/np 
!    do i = 1, np-1   

!      ! the term related to freq: cc * cos(2pi*t * freq) + ss * sin(2pi*t * freq)
!      a(i+1) = a(i+1) - cc * cn(i) - ss * sn(i)
!      
!    enddo  

! instead of starting from 0, we start from 1*darg
    
    do i = 1, np   
      a(i) = a(i) - cc * cn(i) - ss * sn(i)
    enddo  
    
  end do  
 
!  check the substract 
  if(debug == 1) then 
    print*, 'Check the substraction:'
    write(*,'(10e20.10)') a(1:40)
    print*;read*
  endif 
   
  return  
end subroutine substr


subroutine refine( f, nf,  din, idref,  maxferr, isref)
! Refine the frequencieies and amplitude according to idref, we have three approaches:
! 1. only freq 
! 2. only coef 
! 3. freq + coef simultaneously

!       Input Variables 
!  f            all the frequencieies to be defined, dimesion nf-by-4
!               format of the array structure, with the first row as the constant term: 
!               freq -  ccos -- ssine -- amp  
!  nf           number of  frequencieies
! idref         flag of refinement approach: 
!               1.  freq ;    2: coef(cc+ss);   3: freq+coef 
      
!       Output Variables 
!  f            updated  frequencies+coefs
!  maxferr      maximum difference between the input samples and the computed
!               approximation, an estimate of the detection by FFTSC+REFINE
           

!     Module-base Varaibles
!  np           number of samples
!  lt           length of the time interval of the samples
!               the epoch of samples:  t = (i-1)*lt/np, i = 1, np 
!  tol          tolerance( in maximum error) for Newton's method 
!               1.e-6 would be enough? the process will be stopped if maxferr is smaller than tol  

!  tolf        tolerance for the correction(in frequencies) of Newton
!              Method for refinement of only 1 frequecy, i.e., 2.d-5


!     Local Varaibles
!  ferr        the difference between the approximated value and the real one, 
!              of dimension np-by-nrhs, we fix the value nrhs == 1
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
  real(kind=dp)  ::  g(np), dg(np, idref*nf), ferr(np),  df(idref*nf), &
                     dfm_pre, dfm 
  integer        ::  iter, i, j, info, inccor, nrhs
  
  real(kind=dp), external :: dnrm2 
  
  isref = 1
  
  nrhs = 1  ! only one column in right-hand side 
   
! refine iteratively 
  iter   = 0        ! index of iteration 
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
    
    ! check g and dg 
    if(debug == 1) then   
      print*, 'check g:'
      write(*, '(10e20.10)') g(1:20)
      print*; read*
    
      print*, 'check dg:'
      do i = 1, 2
        write(*, '(10e20.10)') dg(i, :) 
      enddo  
      print*; read*
    endif 

    ! the difference between the original data and the computed approximation
    ferr = din - g  
    
!   because of the large amount of samples 
!   we only need to check maximum the error  
!    maxferr =  maxval( dabs(ferr) )  ! discard... use Euclidean norm instead
    maxferr =   dnrm2(np, ferr, 1) / dsqrt( dble(np) ) 
    
    
!-TODO: do we need to do the check on the modulus of the residual here? 
!       Maybe NO! because if all the domimant frequecies are not detected, this vaulue 
!       will be quite large, otherwise, we want to refine more accurately.
!    if(maxferr < tolres) then 
!       write(fwork,*) 'Maximum error in g is less than tolres',  maxferr 
!      return
!    endif

    if( iter > 1)  dfm_pre = dfm
    
    
    !subroutine deltx( nr, nc, a, b, x, info)
    call deltx( np, idref*nf, dg,  ferr, df, info)
    
    dfm  =  dnrm2(idref*nf, df, 1) / dsqrt( dble(idref*nf) )
    
    ! TODO: this we may want to keep
    if(debug == 1) write(fwork, *), '|dx|, |ferr|', dfm, maxferr
      
   
!    if(debug == 1)  then 
!      print*, df ;
!    endif  
    ! TODO: if we donot want terminate here, comment next line    
!    read*
    
    if( info /= 0 ) then 
    
      write(fwork,*)  'SSL Solver Failed to converge!';  read*
      isref = 0
      return
      
    else   
    ! check the modulus of  correction, if it is smaller than tol terminate the iteration 
    ! TODO -- although it works for this example, I think it is better to just check the 
    ! correction in frequecies rather than freq+coef
 
 ! for the termination of the iteration, one criteria could be the change trend of the modulus of the correction,
 ! if it fluctuates, no good, stop. 
 ! TODO: check the maximum norm, or the Euclidean norm?
 !   dfm  =  dnrm2(idref*nf, df, 1) / dsqrt( dble(idref*nf) )
 !   print*, '|dx|, |ferr|', dfm, maxferr
      
      ! TODO --- for sure we should stop when |df| starts increasing
      !          As Gerard suggests, it is possible for the error to increase if the number of frequecies increase and we do the refinement simultaneously
      !          better to check the second increase of error.... 
      
      
      if(iter > 1 .and. dfm > dfm_pre) then 
      
        inccor = inccor + 1 
        
        if(inccor >= 4 ) then ! allow the correction to increase no more than twice... 
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
      if(debug == 1) then 
        write(fwork,*) 'Old:  freq     cos     sin     amp'
        do j = 1, nf 
          write(fwork,*) f(j+1, : )
        enddo
      
        write(fwork,*)  
        if(idref == 1 ) write(fwork,*) 'Correction: freq     cos     sin'
        if(idref == 2 ) write(fwork,*) 'Correction: cos      sin'
        if(idref == 3 ) write(fwork,*) 'Correction: freq     cos     sin'
             
        do j = 1, nf 
           write(fwork,*) df(j: idref*nf: nf) 
        enddo 
      endif  
      
      
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
     ! but the one with this modulus improves the previous freqs a lot
     ! so keep the refined result from an extra frequecy
      if( f(nf+1,4) < tolres/10 .and. nf > 1) then 
        print*,  nf,'-th freq has amplitude less than tolres: ', f(nf+1,4)
        isref = 2
        nf = nf-1
        return
      endif 
      
      ! print the updated freq+coef to check
      if(debug == 1) then 
        write(fwork,*) 
        write(fwork,*) 'New:  freq   cos  sin  amp'
        do j = 1, nf 
          write(fwork,*), f(j+1, : )
        enddo
        write(fwork,*)   
      endif 
      
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
!               the epoch of samples:  t = (i-1)*lt/np, i = 1, np  
! pi2 

!  Routine Used:
!     None

!  Finally Revised by Yu -- 20160607
!----------------------------------------
subroutine refine_gdg( f, nf, idref, g, dg)

implicit none

! Input  and Output Declaration   
integer, intent(in)          ::  nf, idref  
real(kind=dp), intent(in)    ::  f(nf+1, 4)  
real(kind=dp), intent(out)   ::  g(np), dg(np, idref*nf) 
 
! Local Variable
integer        :: i, j
real(kind=dp)  :: freq, cc, ss, dt, cn(np), sn(np), darg  ! arg,  t,
real(kind=dp), dimension(0:np) :: cnb0, snb0

  
! time step size for sampling  
  dt = pi2* lt / np
  
  g = f(1, 2) ! array opertation 

!  ----- deal with frequecy one by one  --------- 
  do j = 2, nf+1   
    
      freq = f(j, 1)
      cc   = f(j, 2) 
      ss   = f(j, 3) 
      
! the sum for i from 0 to nf-1 of  { cc*cos(2*pi*freq * i*dt_step) +  ss*sin(2*pi*freq* i*dt_step)  = cc*cn(i) + ss*sn(i) }
!  if we take arg = 2*pi*freq * dt_step, where freq is the current frequency 
      
   ! TODO:  mod(arg, 2pi) approach is not reliable, we replace with recurrent approach of example 4  
   !      arg  = dmod( pi2* freq *t, pi2)
   !      g(i) = g(i) + cc * dcos(arg) + ss * dsin(arg)
      
      
   ! step size for sampling 
    darg =  freq * dt
  
  !  ---  The Fourier Series with the first np terms ----------
  !  ---Compute cos( n * darg ) and sin( n * darg )  
    call trigrec_c( darg, np, cn)  
    call trigrec_s( darg, np, sn)
    
    cnb0(0) = 1.d0; snb0(0) = 0.d0
    cnb0(1:np) = cn; snb0(1:np) = sn
    
    ! check trigrec_c. trigrec_s 
    
    ! approximate value g recoveried from the detected frequencieies and amplitudes
    ! TODO : Another approach is to add the constant term after all the other terms, advoid too many addtions 
   
!    print*, 'check refine_gdg, nf=', nf !ckd
!    print*, 'ipt=1, ifreq, g', 1, j-1, g(1); read* !ckd
    
!    do i = 0, np-1 
!      !  the constant term f(1, 2) ,  f:  freq-cs  
!       g(i+1) =  g(i+1) +  cc * cnb0(i) + ss * snb0(i) ! start from 0

!!      print*,'ipt, ifreq, g', i+1, j-1, g(i+1); read* !ckd 
!      
!      ! dg =  d g / d w (nf),  d g / d cc (nf), d g / d ss (nf)
!      
!      ! derivative of g w.r.t freq 
!      !  dg / d w(i) = the sum from 1 to nf (cc * -sn(i) * 2pi*dt_step * i  + ss * cn(i) * 2pi*dt_step * i) 
!      !                                     == -cc * sn(i) * dt * i + ss * cn(i) * dt * i
!      
!      if(idref == 1 .or. idref == 3 ) &
!!         dg(i, j-1) = -cc*pi2*t * dsin(arg) + ss*pi2*t * dcos(arg)
!         dg(i+1, j-1) = -cc* snb0(i) * dt * i + ss * cnb0(i) * dt * i
!       
!      ! if we only need to refine frequecy  
!      if(idref == 1) cycle 
!      
!      ! coef 
!      if(idref == 2 ) then 
!        dg(i+1, j-1)    = cnb0(i) ! w.r.t cc 
!        dg(i+1, nf+j-1) = snb0(i) ! w.r.t ss
!       
!      ! for both freq + coef, we have 1-nf: dg/dw, nf+1-2nf: dg/cc,  2nf+1-3nf: dg/ss  
!      elseif(idref == 3) then
!      !
!        dg(i+1, nf+j-1)   =  cnb0(i) ! w.r.t cc 
!        dg(i+1, 2*nf+j-1) =  snb0(i) ! w.r.t ss
!      endif 
!          
!    enddo  
  
  
     do i =1, np 
       g(i) =  g(i) +  cc * cn(i) + ss * sn(i)
       if(idref == 1 .or. idref == 3 ) &
         dg(i, j-1) = -cc* sn(i) * dt * i + ss * cn(i) * dt * i
    
      if(idref == 1) cycle   ! if we only need to refine frequecy  
      
      ! coef 
      if(idref == 2 ) then 
        dg(i, j-1)    = cnb0(i) ! w.r.t cc 
        dg(i, nf+j-1) = snb0(i) ! w.r.t ss
       
      ! for both freq + coef, we have 1-nf: dg/dw, nf+1-2nf: dg/cc,  2nf+1-3nf: dg/ss  
      elseif(idref == 3) then
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
  
  ! ckd, the same as the recurrence  ! TODO:  recovery this part
  do i = 1, np 
    a(i) = a(i) * ( 1-cn(i) )**nh * qnh
  enddo
  
  ! ck if we start from t=0 ?
!  a(1) = 0.d0 
!  do i = 2, np 
!    a(i) = a(i) * ( 1-cn(i-1) )**nh * qnh
!  enddo
  
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
!  write(*,*)' average =', prom
      
  ! deviation from the mean, also called residual
  ain = ain - prom 
  
  return
end subroutine promig


!***********************************************************************
!     ****   deltx   ****
!  Linear equation solver from Lapack: deglsd
!        A* X = B, min || A * X - B||

!       Input Variables 
! nr     number of row in A               
! nc     number of columns in A, also the dimension of dx
! a      the coefficient matrix, for Newton Method, it is usually the Jacobi matix 
! b      the right hand array (always the error) 

!       Output Variables 
! x       the solution of A*X = B         

!  Routine Used:
!    deglsd  

!  Finally Revised by Yu -- 20160
!***********************************************************************
subroutine deltx( nr, nc, a, b, x, info)
implicit none 

! Input  and Output Declaration   
integer, intent(in)         ::  nr, nc 
integer, intent(out)        ::  info
real(kind=dp), intent(in)   ::  a(nr, nc), b(nr) 
real(kind=dp), intent(out)  ::  x(nc)
 
 
!  --- explation of the variables usde in dgelsd -------------
! dgelsd (M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO)
!  M       number of rows 
!  N       number of columns 
! NRHS     number of right hand sizes, i.e., number of column of X and B 
! A        LDA-by-N matrix, the coefficient matrix 
! LDA      leading dimension of A,  >=max(1,M)

! B        LDB-by-NRHS matrix, the right hand side matrix is DOUBLE PRECISION array,
!          dimension (min(M,N)). On entry, the M-by-NRHS right hand side matrix B.
!          On exit, B is overwritten by the N-by-NRHS solution matrix X.
! LDB      leading dimension of B,  >=max( 1, max(M,N) )

! S        The singular values of A in decreasing order. dimension (min(M,N))
!          The condition number of A in the 2-norm = S(1)/S(min(m,n)). 

!RCOND     RCOND is DOUBLE PRECISION
!          RCOND is used to determine the effective rank of A.
!          Singular values S(i) <= RCOND*S(1) are treated as zero.
!          If RCOND < 0, machine precision is used instead.
! RANK     The effective rank of A, i.e., the number of singular values
!          which are greater than RCOND*S(1).
! WORK     WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
! LWORK    For good performance, LWORK should generally be larger.
!           12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2
!            NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
!            SNLSIZ ~= 25

! IWORK     dimension (MAX(1,LIWORK))
!           LIWORK >= max(1, 3 * MINMN * NLVL + 11 * MINMN),MINMN = MIN( M,N ).
!           On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK.           

! INFO     = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  the algorithm for computing the SVD failed to converge;
!                if INFO = i, i off-diagonal elements of an intermediate
!                bidiagonal form did not converge to zero

! -------- for deglsd  ---------- 
external  ::  dgelsd  ! lapack lls solver
      
! length of work,  3*n + 64*(m+1), use a big value to allocate enough memory
integer, parameter ::  lwork = 2**16 
integer            ::  nrhs, rank, iwork(lwork) 
real(kind=dp)      ::  s(np), rcond, work(lwork) , acopy(nr, nc) 
! -----------------------------------------------------------
      
! Local Variable
real(kind=dp)  :: bwork(nr)
 
! LLS  slover initialization 
  nrhs = 1  ! only one column in B 
  acopy = a
! Choose RCOND to reflect the relative accuracy of the input data, for double precision, 
! use a small one, a good option is one smaller than the cpu precision. 
  RCOND = 1.d-14
       
  bwork = b   ! working array, not to change b... 
  
  CALL DGELSD(nr, nc, nrhs, acopy, nr, bwork, nr, S, RCOND, RANK, WORK, LWORK, IWORK, INFO)
  x = bwork(1:nc)
  return  
end subroutine deltx
 
 

end module fft_mod
