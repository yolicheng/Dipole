subroutine matnorm(a, row, col, nm)
! 20160222 
! check the matrix norm associated to the vector norm
! 	Input Varaibles
!  a 		the matrix 
!  row,col 	the size of the input matrix row-by-col

! 	Output Varaibles
!  nm 		the matrix norm, which is the maximum absolute row sum of the matrix, the infinity norm 


implicit none
integer, parameter ::  dp = kind(1.d0) 

integer, intent(in)::  row, col
real(kind=dp), intent(in) ::  a(row, col)
real(kind=dp), intent(out) ::  a(row, col)

! local variables
integer :: i, j, ifam, dir, imax, id, fpo, fpoinst, fmmat, fegp, fegv, tdir

real(kind=dp) ::  dlf3(n,n), dlfsswap(n),swap2(n), & ! differential of vector field of lorentz force
                  wr(n),wi(n), vr(n,n), y0(neq), &   ! stability of equilibrium points
                  poinst(6), ds,tol, vf(neq), ynew(npo,8), epsl_po, cj,  & !pofam
                  po0(6), tpo, pof(6), & ! plpo
                  mmat(n,n), wr_mm(n,n), wi_mm(n,n), vr_mm(n,n), & ! MM 
                  ymfd(nmf,6), epsl, prsc ! mfd 
                  
external :: gr_lf, gr_cjlf , &
	    gr_rtbp, gr_cjrtbp
