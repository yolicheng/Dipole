      subroutine prt_eigvec( n, ftag, wi, v )
      integer          n, ftag
      double precision wi( * ), v( n, * )

      double precision zero
      parameter        ( zero = 0.0 )
      integer          i, j

!     write(ftag,*)
      do i = 1, n
         j = 1
         do while( j.le.n )
            if( wi( j ).eq.zero ) then
               write(ftag,9998,advance='no') v( i, j )
               j = j + 1
            else
               write(ftag,9999,advance='no') v( i, j ), v( i, j+1 )
               write(ftag,9999,advance='no') v( i, j ), -v( i, j+1 )
               j = j + 2
            end if
         end do
         write(ftag,*)
      end do
      write(ftag,*)
      

! better without the parenthesis.

 9998 format( 11(f12.8, 1x) )
 9999 format( 11(f12.8, f12.8, 2x) )
      return
      end

