      SUBROUTINE prntft(data, n2, fdat)
      INTEGER i,n2,nn2,m, fdat
! 	Input
!  	data  	dimesnion n*n, to be printed, 
! 	fdat  	the file to be write into, 6: the screen
      
      REAL*8 data(n2)
      
      write(*,'(1x,t7,a,t13,a,t24,a,t35,a,t47,a)')
     *     'n','Real(n)','Imag.(n)','Real(N-n)','Imag.(N-n)'
     
      write(fdat,'(1x,i6,4f20.10)') 0,  data(1),data(2),data(1),data(2)
      
      do 11 i=3,(n2/2)+1,2
        m=(i-1)/2
        nn2=n2+2-i
        write(fdat,'(1x,i6,4e20.10)') m,data(i),data(i+1),
     *       data(nn2),data(nn2+1)
11    continue

      write(fdat,*) ! add one blank line
      write(*,'(/1x,a)') ' press RETURN to continue ...'
      read(*,*)
      return
      END
