! ***** fast Fourier transform of two real functions ******
      SUBROUTINE dtwofft(data1,data2,fft1,fft2, n)
      INTEGER n
      REAL*8 data1(n),data2(n)
      DOUBLE COMPLEX fft1(n),fft2(n)
      
CU    USES four1
      INTEGER j,n2
      DOUBLE COMPLEX h1,h2,c1,c2
      
      c1=dcmplx(0.5,0.0)
      c2=dcmplx(0.0,-0.5)
      
      do 11 j=1,n
        fft1(j)=dcmplx(data1(j),data2(j))
11    continue

      call dfour1(fft1,n,1)
      fft2(1)=dcmplx(dimag(fft1(1)),0.0)
      fft1(1)=dcmplx(dble(fft1(1)),0.0)
      n2=n+2
      
      do 12 j=2,n/2+1
        h1=c1*(fft1(j) + dconjg(fft1(n2-j)))
        h2=c2*(fft1(j) - dconjg(fft1(n2-j)))
        fft1(j)=h1
        fft1(n2-j) = dconjg(h1)
        
        fft2(j)=h2
        fft2(n2-j) = dconjg(h2)
12    continue

      return
      END
