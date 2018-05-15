!***********************************************************************
!     ****   dcj_da_debug   ****
!  check dcj_da 

!  Finally Revised by Yu -- 20170528 ! ckd --- seems ok.
!***********************************************************************
subroutine dcj_da_debug(xi, step, n, m, nf, c0, ck, sk, gr_cj, dcj_da)
use dp_mod
use curve_time_mod, only : varphi 
implicit none
 
! Input  and Output Declaration  
integer, intent(in)     ::  n, m,  nf  
real(kind=dp), intent(in)      ::   xi, step, c0(n*m), ck(n*m,nf), sk(n*m, nf) 
real(kind=dp), intent(out)     ::   dcj_da(n*(2*nf+1))  
 
external :: gr_cj
 
! Local Variable
integer        :: nf2p1, ir, i, j
real(kind=dp)  :: cjp, cjm,  pvp(n), pvm(n), pv(n), cj , &
                  c0p(n*m), ckp(n*m,nf), skp(n*m, nf), c0m(n*m), ckm(n*m,nf), skm(n*m, nf) 
                  
                  
  
  call varphi(xi, n, 0, nf, c0, ck, sk, pv)
  call gr_cj(pv, cj )
  
!  print*, 'Original cj:'
!  print*, pv,  cj; print*; read* 
  
  nf2p1 = 2*nf + 1
  
  
  do i = 1, n
    c0p = c0; c0m = c0  
    
    ir = nf2p1*(i-1)
    
    ! dcj/ dc0 
    c0p(i) = c0(i) + step; 
    c0m(i) = c0(i) - step  
    call varphi(xi, n, 0, nf, c0p, ck, sk, pvp)
    call varphi(xi, n, 0, nf, c0m, ck, sk, pvm)
    
    call gr_cj(pvp, cjp)
    call gr_cj(pvm, cjm)
    
    dcj_da(ir+1) = (cjp - cjm) / step / 2.d0 
!   
    ! dcj/ dck  
    do j = 1, nf 
      
      ckp = ck; ckm = ck  
      
      ckp(i,j) = ck(i,j) + step; 
      ckm(i,j) = ck(i,j) - step  
      call varphi(xi, n, 0, nf, c0, ckp, sk, pvp)
      call varphi(xi, n, 0, nf, c0, ckm, sk, pvm)
    
      call gr_cj(pvp, cjp)
      call gr_cj(pvm, cjm)
    
      dcj_da(ir+1+j) = (cjp - cjm) / step / 2.d0 
      
      
      skp = sk;  skm = sk
      
      skp(i,j) = sk(i,j) + step; 
      skm(i,j) = sk(i,j) - step 
       
      call varphi(xi, n, 0, nf, c0, ck, skp, pvp)
      call varphi(xi, n, 0, nf, c0, ck, skm, pvm)
    
      call gr_cj(pvp, cjp)
      call gr_cj(pvm, cjm)
      
      dcj_da(ir+1+nf+j) = (cjp - cjm) / step / 2.d0
      
    enddo  
   
  enddo     
    
  return  
end subroutine dcj_da_debug


















