      PROGRAM dxtwofft
C     driver for routine dtwofft
      INTEGER N,N2
      REAL*8 PER,PI
      PARAMETER(N=32,N2=2*N,PER=8.0,PI=4*datan(1.d0))
      INTEGER i,isign
      REAL*8 data1(N),data2(N),fft1(N2),fft2(N2), x, dx
      
!      dx = PI/3
!      do 11 i=1,N
!        x=2.0*PI*i/PER
!        data1(i)= 3.d0/4.d0*dcos(2*pi*2* x ) 
!        data2(i)= 0.2*dsin(2*pi*7* x ) 
!11    continue

      do 11 i=1,N
        x=2.0*PI*i/PER
        data1(i)= dcos(x ) 
        data2(i)= dsin(x ) 
11    continue

!      do 11 i=1,N
!        x=2.0*PI*i/PER
!        data1(i)=nint(cos(x ))
!        data2(i)=nint(sin(x ))
!11    continue

c 	check the initial data
	print*, 'data1' 
      do 12 i=1,N
        print*, data1(i) 
12    continue
	print*, 'data2' 
      do 13 i=1,N
        print*, data2(i) 
13    continue

      call dtwofft(data1,data2,fft1,fft2,N)
      write(*,*) 'Fourier transform of first function:'
      call prntft1(fft1,N2)
      write(*,*) 'Fourier transform of second function:'
      call prntft1(fft2,N2)
      
C     invert transform
      isign=-1
      call dfour1(fft1,N,isign)
      write(*,*) 'Inverted transform = first function:'
      call prntft1(fft1,N2)
      
      call dfour1(fft2,N,isign)
      write(*,*) 'Inverted transform = second function:'
      call prntft1(fft2,N2)
      END



      SUBROUTINE prntft1(data,n2 )
      INTEGER i,n2,nn2,m
      REAL*8 data(n2)
      
      write(*,'(1x,t7,a,t13,a,t24,a,t35,a,t47,a)')
     *     'n','Real(n)','Imag.(n)','Real(N-n)','Imag.(N-n)'
      write(*,'(1x,i6,4f12.6)') 0,data(1),data(2),data(1),data(2)
      
      do 11 i=3,(n2/2)+1,2
        m=(i-1)/2
        nn2=n2+2-i
        write(*,'(1x,i6,4f12.6)') m,data(i),data(i+1),
     *       data(nn2),data(nn2+1)
11    continue

      write(*,'(/1x,a)') ' press RETURN to continue ...'
      read(*,*)
      return
      END
