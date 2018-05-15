      implicit real*8(a-h,o-z)
      
      character*14 fitin,fitout
      
      dimension fren(990), amp(990),cc(990),ss(990),deltav(990)
      dimension freba(20), kas(20),kascom(20,990)
      
!      write (*,*) ' Input file'
!      read (*,*) fitin

!      fitin = '111.dat'
      
      fitin = 'fcdat.dat'
      fitout = 'fcout.dat'
      
!      write (*,*) ' Output file'
!      read (*,*) fitout
!      
      open(1, file=fitin)
      open(2, file=fitout)
      
      
      delmax=0
      
      nfre=0

1     read(1,*,end=2)  x,y,z, a, a2
      write (*,*) x,y, z,a  
 
!      read*
 
      nfre=nfre+1
      fren(nfre)=x
      cc(nfre)=y
      ss(nfre)=z
      go to 1
      
2     write(*,*)' I have',nfre,'  frequencies'
      write(2,*)' I have',nfre,'  frequencies'
      
      
      do 3 j=1,400
         deltav(j)=0
         do 3 k=1,20
3     kascom(k,j)=0

      write(*,*)' Maximum order (sum(abs(coefs)) and
     1 Tolerace for the freq identification?'
     
      read(*,*)  maxor,tolfre
      write(2,*) maxor,tolfre
      call ordfre(nfre,fren,cc,ss,amp)
      
      
      nfreba=1
      jj=1
6     freba(nfreba)=fren(jj)
      write(*,*)' I have one basic freq=',fren(jj)
      kascom(nfreba,jj)=1
      do 4 j=nfreba+1,nfre
         if(amp(j).lt.1.d-14)go to 4
         fre=fren(j)
         call resfre(fre,nfreba,freba,maxor,tolfre,kas,iexit,delta)
         if(iexit.eq.0)go to 4
         amp(j)=0
         deltav(j)=delta
         do 41 k=1,nfreba
41       kascom(k,j)=kas(k)
4     continue
      do 5 j=nfreba+1,nfre
         if(amp(j).lt.1.d-14)go to 5
         nfreba=nfreba+1
         jj=j
         go to 6
5     continue
      write(*,*)' Total number of basic freqs:',nfreba
      write(2,*)' Total number of basic freqs:',nfreba
      do 7 j=1,nfreba
         write(*,*)freba(j)
         write(2,*)freba(j)
7     continue
      write(*,*)' Frequencies, linear combinations & amplitudes'
      do 8 j=1,nfre
      write(*,100)j,fren(j),deltav(j),dsqrt(cc(j)**2+ss(j)**2)
      write(2,100)j,fren(j),deltav(j),dsqrt(cc(j)**2+ss(j)**2)
100   format(' ',i3,f15.10,2d15.7)
      kasum=0
      delmax=dmax1(delmax,deltav(j))
      do 9 k=1,nfreba
9     kasum=kasum+iabs(kascom(k,j))
      write(*,101)kasum,(kascom(k,j),k=1,nfreba)
8     write(2,101)kasum,(kascom(k,j),k=1,nfreba)
101   format(' ',i6,'   ',20i3)
      write(*,*)' Max error of the linear combinations =',delmax
      write(2,*)' Max error of the linear combinations=',delmax
      stop
      end

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

      subroutine ordfre(nfre,fren,cc,ss,amp)
      implicit real*8(a-h,o-z)
      dimension fren(990),amp(990),cc(990),ss(990)
      do 1 j=1,nfre
1     amp(j)=dsqrt(cc(j)*cc(j)+ss(j)*ss(j))
2     ic=0
      do 3 j=1,nfre-1
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
c   subroutine genlex: donat un ordre lexicografic d'un cert grau, genera
c   la segment tira. Quan arriba a (0,0,...,grau total) genera de nou la
c   inicial. Per crides successives puc generar totes les tires.
c
      subroutine genlex(iexp,nvar)
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
