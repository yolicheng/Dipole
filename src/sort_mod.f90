! This module inclues two of the most popular sorting algorithms

! Optimized quicksort and insertion sort. 

! Idea:  Pack the index together with the 1 dimensional array according to which to sort. 
!        And work with the packed type (group: value + ind) during the swap.

!        in this way,  the sort of a 2 dimensinal array is simpified to sorting a 1d array 
!        plus a swap of the array according to the ordered index (which is group%ind). 

! 2016-08-05 18:30:55 
!   Return also the ordered index of the chosen column, for possible use or just to check. 

!subroutine sort_1d( a, na, order)
!subroutine sort_2d( a, nr, nc, ic, order)

module sort_mod

implicit none

integer, parameter, private :: dp = kind(1.d0)

! pack the value to be sorted together with its original order of the index 
 
type group 
  integer :: ind             ! original index of unsorted data 
  real(kind=dp) ::  value    ! value to be sorted 
end type group 


contains


!***********************************************************************
!   Sort 1d array 
!   by default, it is in increasing order, is order = -1, we do a flipover 
!***********************************************************************
subroutine sort_1d( a, na, order, ind_sorted)

!       Input Variables 
!  A         the 1d array to be sorted         
! nA         size of to array to be sorted
! order      1: sort in increasing order;  -1: sort in decreasing order 

!       Output Variables
! A          the sorted array              
!ind_sorted  the sorted index of the 1d array  according to which to sort 

!  Routine Used:
!    quicksort, insertionsort  

!  Finally Revised by Yu -- 20160709
!----------------------------------------
implicit none
!integer, parameter :: dp = kind(1.d0)
   
! DUMMY ARGUMENTS
integer, intent(in)              :: nA, order  
real(kind=dp), intent(inout)     :: A(na)
integer, intent(out)             :: ind_sorted(na)

! LOCAL VARIABLES
integer :: i 
type(group) :: awork(na)   ! working array: order + value 

! Initialize awork
do i = 1, na, 1
  awork(i)%value = a(i)
  awork(i)%ind = i
end do

call quicksort(awork, na)   
    
! reassign the array and the index array, try this array operation, if it fails try element-by-element assignment
a            = awork%value
ind_sorted   = awork%ind 

! if we need to sort in decreasing order, flipover the sorted array and the index array.
if(order == -1) then 
  a         = a(na:1:-1)
 ind_sorted = ind_sorted(na:1:-1)
endif   

return  
end subroutine sort_1d


!***********************************************************************
!   Sort 2d array
!***********************************************************************
subroutine sort_2d( a, nr, nc, ic, order, ind_sorted)
!       Input Variables 
! a          the array to be sorted         
! nr         number of rows of the array to be sorted
! nc         number of columns of the array
! ic         index of column to be sorted according to
!            TODO: for the moment, only sort by column, consider no row, according to my habbit of saving data
! order      1: sort in increasing order;  -1: sort in decreasing order 

!       Output Variables
! a          the sorted array              
!ind_sorted  the sorted index of the 1d array  according to which to sort 

!  Routine Used:
!    quicksort, insertionsort  

!  Finally Revised by Yu -- 20160709
!----------------------------------------
implicit none
!integer, parameter :: dp = kind(1.d0)
   
! DUMMY ARGUMENTS
integer, intent(in)          :: nr, nc, ic, order  
real(kind=dp), intent(inout) :: a(nr,nc)
integer, intent(out)         :: ind_sorted(nr)



! LOCAL VARIABLES
integer     :: i   
type(group) :: awork(nr)   ! working array: ind + value 

! Initialize awork using the ic-th column, 

do i = 1, nr, 1
  awork(i)%value = a(i, ic)
  awork(i)%ind   = i
end do

call quicksort(awork, nr)   

! reassign the array and the index array using array operation  
ind_sorted = awork%ind
a          = a(ind_sorted, :)


if(order == -1)  then 
   a          = a(nr:1:-1, :)
   ind_sorted = ind_sorted(nr:1:-1)
endif    

return  
end subroutine sort_2d


!***********************************************************************
!   Quicksort
!   Main routine is quicksort, combined with an INSERTION sort for smaller sized partitions, (na <= 20) 
!   where QuickSort is less efficient. 
!***********************************************************************
recursive subroutine QuickSort(A, nA)

! DUMMY ARGUMENTS
integer, intent(in) :: nA
type(group), intent(in out) :: A(na)

! We only call quicksort when the size of A is bigger than limit
! otherwisze, do insertion sort, which is more efficient for small array 
! suggestion value : 20, set as the default value 

integer, parameter :: limit = 20

! LOCAL VARIABLES
integer :: left, right, i, eighth

real(kind=dp)   :: pivot, temp
integer         :: marker
real(kind=dp)   :: sample(9)

! Rroblem with the size of sample, 

! working array: ind + value 
type(group) ::  tempwork
    eighth = nA/8
!     print*, 'nA =', nA, 'eighth = nA/8 =', eighth; read*

    if (nA > 1) then
        if (nA > limit) then ! Do quicksort for large groups
            ! ************************
            ! 9-SAMPLE PIVOT METHOD
            eighth = (nA-1)/8;  !nA/8-1
           
!            print*, 'nA =', nA, 'eighth = ', eighth; read*
!            print*, a(1:nA:eighth)
!            read*
!            sample = a(1:nA:eighth)%value
            
            do i = 1, 9, 1
              sample(i) = a(1+(i-1)*eighth)%value
            end do 


            ! Sort Network for N=9, using Batcher's Merge-Exchange. Skip some steps because I only care about the median (5)
            if (sample(1) > sample(9)) then; temp = sample(1); sample(1) = sample(9); sample(9) = temp; end if
            if (sample(1) > sample(5)) then; temp = sample(1); sample(1) = sample(5); sample(5) = temp; end if
            if (sample(2) > sample(6)) then; temp = sample(2); sample(2) = sample(6); sample(6) = temp; end if
            if (sample(3) > sample(7)) then; temp = sample(3); sample(3) = sample(7); sample(7) = temp; end if
            if (sample(4) > sample(8)) then; temp = sample(4); sample(4) = sample(8); sample(8) = temp; end if
            if (sample(5) > sample(9)) then; temp = sample(5); sample(5) = sample(9); sample(9) = temp; end if
            if (sample(1) > sample(3)) then; temp = sample(1); sample(1) = sample(3); sample(3) = temp; end if
            if (sample(2) > sample(4)) then; temp = sample(2); sample(2) = sample(4); sample(4) = temp; end if
            if (sample(5) > sample(7)) then; temp = sample(5); sample(5) = sample(7); sample(7) = temp; end if
            if (sample(6) > sample(8)) then; temp = sample(6); sample(6) = sample(8); sample(8) = temp; end if
            if (sample(3) > sample(9)) then; temp = sample(3); sample(3) = sample(9); sample(9) = temp; end if
            if (sample(3) > sample(5)) then; temp = sample(3); sample(3) = sample(5); sample(5) = temp; end if
            if (sample(4) > sample(6)) then; temp = sample(4); sample(4) = sample(6); sample(6) = temp; end if
            if (sample(7) > sample(9)) then; temp = sample(7); sample(7) = sample(9); sample(9) = temp; end if
            if (sample(1) > sample(2)) then; temp = sample(1); sample(1) = sample(2); sample(2) = temp; end if
            if (sample(3) > sample(4)) then; temp = sample(3); sample(3) = sample(4); sample(4) = temp; end if
            if (sample(5) > sample(6)) then; temp = sample(5); sample(5) = sample(6); sample(6) = temp; end if
            if (sample(7) > sample(8)) then; temp = sample(7); sample(7) = sample(8); sample(8) = temp; end if
            if (sample(2) > sample(9)) then; temp = sample(2); sample(2) = sample(9); sample(9) = temp; end if
            if (sample(2) > sample(5)) then; temp = sample(2); sample(2) = sample(5); sample(5) = temp; end if
            if (sample(4) > sample(7)) then; temp = sample(4); sample(4) = sample(7); sample(7) = temp; end if
           !if (sample(6) > sample(9)) then; temp = sample(6); sample(6) = sample(9); sample(9) = temp; end if ! skipped
           !if (sample(2) > sample(3)) then; temp = sample(2); sample(2) = sample(3); sample(3) = temp; end if ! skipped
            if (sample(4) > sample(5)) then; temp = sample(4); sample(4) = sample(5); sample(5) = temp; end if
           !if (sample(6) > sample(7)) then; temp = sample(6); sample(6) = sample(7); sample(7) = temp; end if ! skipped
           !if (sample(8) > sample(9)) then; temp = sample(8); sample(8) = sample(9); sample(9) = temp; end if ! skipped
            pivot = sample(5)
            ! ************************
            
            left = 0
            right = nA + 1
            
            do while (left < right)
            
                right = right - 1
                do while (a(right)%value > pivot)
                  right = right - 1
                end do
                
                left = left + 1
                do while (a(left)%value < pivot)
                  left = left + 1
                end do
                
                if (left < right) then
                    tempwork = a(left) 
                    a(left) = a(right)
                    a(right) = tempwork
                end if
            end do

            if (left == right) then
                marker = left + 1
            else
                marker = left
            end if

            call quicksort(a(:marker-1), marker-1 )
            call quicksort(a(marker:), nA-marker+1 )

        else
        
            call InsertionSort(a, na)    ! Insertion sort for small groups is faster than Quicksort
            
        end if
    end if


return
end subroutine QuickSort


!**********************************************************************
!  Insertion sort 
!  Sort the 1 dimensional array A of type group(ind+value) using insertion sort method
!**********************************************************************

subroutine InsertionSort(a, na)

! DUMMY ARGUMENTS
integer, intent(in)        :: na 
type(group), intent(inout) :: a(na)

! LOCAL VARIABLES
type(group)    :: temp 
integer        :: i, j

! if A is multi-dimensional array, deal with the ic-th column 
! otherwisze, ic and ind is not used

  do i = 2, na
    j = i - 1
    temp = a(i)
        
    do
      if (j == 0) exit  
      if (a(j)%value <= temp%value) exit 
           
      a(j+1) = a(j) 
      j = j - 1
    enddo  
    
    ! place the i-th element in A to the right place j+1  
    a(j+1)  = temp
  end do  

  return
end subroutine InsertionSort

end module sort_mod

