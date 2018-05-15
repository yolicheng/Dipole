      !***********************************************************************
      !     ****   frebas   ****
      !  Detect the basic frequecies, by linear combiation of order less than odrmax, if we are not satisfied with the result, it provides possibility to redetect with new values for error control (ordmax, tolf) 
      !  
      !  to avoid allocate and deallocate all the time, we declare array fc of dimesion  100-by-5, just make sure the callee pass the array with the same dimension 
      
      ! ** NOTE ** 
      !  the basic frequecies are not uniquely defined, they do not have to have some physical meanings
      !  we have to deal with the pecific problems 
      
      !       Input Variables 
      !  fc       input freq+coef+amp, dimension 100-by-5
      !  nf       number of frequecies to be detected for basic freqs
      !  ffbas    file tag to save all the result of detection process
      !  nfbas    number of basic frequecies detected 
      !  fcbas    file tag to save the basic frequecies and error info
               
      
      !       Output Variables 
      !  delmax   max error of linear combination, estimate the accuracy in basic frequency detection           
      
      !  Routine Used:
      !     resfre, genlex,  ( sort_mod )
      
      !  Finally Revised by Yu -- 20160620
      !***********************************************************************
      subroutine frebas(fcin, nf, ffbas, ordmax, tolf,  
     .                  nfbas, fcbas, delmax)
      use sort_mod
!     use fft_mod ! to pass ordmax and tolf 
      
      implicit real*8(a-h,o-z)
      integer, parameter :: dp = kind(1.d0)   
      
      ! Input and Output Declaration   
      integer, intent(in)            ::  nf, ffbas  
      integer, intent(inout)         ::  ordmax
      real(kind=dp), intent(inout)   ::  tolf
      real(kind=dp), intent(in)      ::  fcin(100, 5) 
      integer, intent(out)           ::  nfbas        
      real(kind=dp), intent(out)     ::  fcbas(20, 5), delmax 
!       
      ! Local Variable
      integer        :: isnew ! ind_sort(100),ordmax
!      real(kind=dp)  :: tolf 
        
     
!     character*14 fitin, fitout
      dimension fren(100), amp(100), cc(100), ss(100), deltav(100)
      dimension freba(20), kas(20), kascom(20,100), fc(100,5), tmp(5)
      
      fc = fcin
      print*, nf, 'Frequencies to be detected :'
      do i = 1, nf 
         write(*, '(I5, 4e24.14, f10.0)') i, fc(i, :)
      enddo   
      write(*, *)  ! add a blank line 
!      read*
       
      ! Initialize 
111   delmax = 0.d0
      print* 
      print*, 'Wanna dirty trick to set manually the basic frequecy?  
     .   (1 = Yes)'
      read*, isetbas
      
      if(isetbas==1) then 
        print*, 'Input the index of the frequency to put as first!'
        read*, nr_fbas  
        print*, fc(nr_fbas, :); print*
        tmp   =  fc(nr_fbas, :)
        fc(2:nr_fbas, :) = fc(1:nr_fbas-1, :)
        fc(1, :) = tmp
      endif 
      
      
      fren(1:nf) = fc(1:nf, 1)
      cc(1:nf)   = fc(1:nf, 2)
      ss(1:nf)   = fc(1:nf, 3)
      amp(1:nf)  = fc(1:nf, 4)
      
      deltav = 0.d0 
      kascom = 0
      
      ! Error control, better to put to the module, so we do not need 
      ! to update each time we call it, at least for one orbit 
      
!      ordmax = 5
!      tolf   = 5.d-6
      
      print*;
      write(*,*) 'Default vaule of Maximum order (sum(abs(coefs)) and
     1 Tolerace for the freq identification?'
      write(*,*) ordmax, tolf 
      
      print*; write(*,*) 'Do you want to use new vaules (1 = Yes)'
      read(*,*)  isnew 
      
      if(isnew == 1) then 
        write(*,*) 'Please Input New values: '
        read(*,*)  ordmax,tolf
      endif
        
!  sort all the detected frequecies by amplitude in decreasing order 
!  since we detect the frequecies by amplitude, not necessary
!      call ordfre(nf,fren,cc,ss,amp)
     
      ! subroutine sort_2d( a, nr, nc, ic, order, ind_sorted)
      ! dimension nf-by-5, 4-th column: amp
!      call sort_2d(fc, nf, 5, 4, 1, ind_sorted(1:nf)) 
      
      
! the first one has the maximal amplitude, treat as the first basic freq by default      
      nfbas = 1
      jj = 1
      
6     freba(nfbas)   = fren(jj)
      fcbas(nfbas, :) = fc(jj, :) !for output
      
      write(*,*)' I have one basic freq=', fren(jj)
      kascom(nfbas,jj) = 1
      
      do 4 j = nfbas+1, nf
         if(amp(j).lt.1.d-14) go to 4
         
         fre = fren(j)
         call resfre(fre, nfbas, freba, ordmax,tolf, kas,iexit,delta)
         if(iexit.eq.0) go to 4
         
         amp(j)=0
         deltav(j) = delta
         do 41 k = 1, nfbas
41       kascom(k,j) = kas(k)
4     continue

      do 5 j = nfbas+1, nf
         if(amp(j).lt.1.d-14) go to 5
         
         nfbas = nfbas + 1
         jj = j
         go to 6
5     continue

      write(*,*)  'Total number of basic freqs:', nfbas
      
      do j = 1, nfbas
         write(*,*)     freba(j)
      enddo
      write(*,*)

      write(*,*)' Frequencies, linear combinations & amplitudes'
      
      do 8 j = 1, nf
      write(*,100)     j,fren(j),deltav(j),dsqrt(cc(j)**2+ss(j)**2)
100   format(' ',i3,f15.10, 2d15.7)

      kasum=0
      delmax = dmax1(delmax,deltav(j))
      
      do 9 k = 1, nfbas
9     kasum = kasum + iabs(kascom(k,j))

8     write(*,101)     kasum,(kascom(k,j),k=1,nfbas)
!8     write(ffbas,101) kasum,(kascom(k,j),k=1,nfbas)
101   format(' ',i6,'   ',20i3)

      write(*,*)   'Max error of the linear combinations =', delmax
      
      ! if we are not satisfied with this detection, restart with new values of tolf
      print*; write(*,*) 'Redetect the basic frequecies (1 = Yes)?'
      read(*,*) isnew 
      if(isnew == 1) goto  111
      
      ! write the output in file with tag ffbas 
      write(ffbas,*) '# Total number of basic freqs:', nfbas
      write(ffbas,*) '# Max error of the linear combinations= ', delmax
      do j = 1, nfbas
         write(ffbas,'(4e24.14, f10.0)') fcbas(j, :)
      enddo
      write(ffbas,*) 
      
      return
      end subroutine frebas
      
      
      
!**********************************************************************
!
!**********************************************************************
      subroutine resfre(fre,nfreba,freba,maxor,tolfre,kas,iexit,delta)
      implicit real*8(a-h,o-z)
      dimension freba(20),kas(20),iexp(20)
      do 100 ior=1,maxor
         iexp(1)=ior
         do 101 jj=2,nfreba
101      iexp(jj)=0
102      do 108 jj=1,nfreba
108      kas(jj)=iexp(jj)
         go to 109
106      jneg=nfreba
110      if(kas(jneg).gt.0)then
            kas(jneg)=-kas(jneg)
            go to 109
         else
            kas(jneg)=-kas(jneg)
            jneg=jneg-1
            if(jneg.eq.0)go to 107
            go to 110
         endif
109      fretes=0
         do 105 ll=1,nfreba
105      fretes=fretes+kas(ll)*freba(ll)
         delta=fretes-fre
         if(dabs(delta).lt.tolfre)go to 103
         go to 106
107      if(iexp(nfreba).eq.ior)go to 100
         call genlex(iexp,nfreba)
         go to 102
100   continue
      iexit=0
      go to 1
103   iexit=1
1     continue
      return
      end   

! --- discard this one, just to sort by the amplitude, bubble sorting---
! One reason is we detect the frequecies by amplitude, they are already sorted 
! if needed, we have more efficient sorting method in sort_mod

      subroutine ordfre(nf,fren,cc,ss,amp)
      implicit real*8(a-h,o-z)
      dimension fren(100),amp(100),cc(100),ss(100)
      do 1 j=1,nf
1     amp(j)=dsqrt(cc(j)*cc(j)+ss(j)*ss(j))

2     ic=0
      do 3 j=1,nf-1
         if(amp(j).lt.amp(j+1))then
            a=amp(j)
            amp(j)=amp(j+1)
            amp(j+1)=a
            
            a=fren(j)
            fren(j)=fren(j+1)
            fren(j+1)=a
            
            a=cc(j)
            cc(j)=cc(j+1)
            cc(j+1)=a
            
            a=ss(j)
            ss(j)=ss(j+1)
            ss(j+1)=a
            
            ic=ic+1
         endif
3     continue
      if(ic.gt.0)go to 2
      return
      end

c
c   subroutine genlex: donat un ordre lexicogr…fic d'un cert grau, genera
c   la segment tira. Quan arriba a (0,0,...,grau total) genera de nou la
c   inicial. Per crides successives puc generar totes les tires.
c
      subroutine genlex(iexp, nvar)
      dimension iexp(nvar)
      is=iexp(nvar)
      do 1 nv=nvar-1,1,-1
         is=is+iexp(nv)
         if(iexp(nv).ne.0)then
            iexp(nv)=iexp(nv)-1
            iexp(nv+1)=is-iexp(nv)
            return
         else
            iexp(nv+1)=0
         endif
1     continue
      iexp(1)=is
      return
      end
