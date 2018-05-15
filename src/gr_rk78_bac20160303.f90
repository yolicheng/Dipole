! try to put gr_rk78 inside this module to see if it works
!********************************************************************
!
!  THIS ROUTINE IS AN IMPLEMENTATION OF A RUNGE-KUTTA-FEHLBERG
!  METHOD OF ORDERS 7 AND 8. USING A TOTAL OF 13 STEPS (AND
!  EVALUATIONS OF THE VECTORFIELD) IT COMPUTES TWO DIFFERENT
!  ESTIMATIONS OF THE NEXT POINT. THE DIFFERENCE BETWEEN BOTH
!  ESTIMATIONS (WITH LOCAL ERRORS OF ORDER 8 AND 9) IS COMPUTED
!  AND THE L1 NORM IS OBTAINED. THIS NORM IS DIVIDED BY N (THE
!  NUMBER OF EQUATIONS). THE NUMBER OBTAINED IN THIS WAY IS REQUIRED
!  TO BE LESS THAN A GIVEN TOLERANCE E1 TIMES (1+0.01*DD) WHERE DD
!  IS THE L1 NORM OF THE POINT COMPUTED TO ORDER 8. IF THIS
!  REQUIREMENT IS SATISFIED THE ORDER 8 ESTIMATION IS TAKEN AS THE
!  NEXT POINT. IF NOT, A SUITABLE VALUE OF THE STEP H IS OBTAINED
!  AND THE COMPUTATION IS STARTED AGAIN.
!  IN ANY CASE, WHEN THE NEXT POINT IS COMPUTED, A PREDICTION OF
!  THE STEP H, TO BE USED IN THE NEXT CALL OF THE ROUTINE, IS
!  DONE.
!
!  INPUT DATA:
!
!       X  CURRENT VALUE OF THE INDEPENDENT VARIABLE.
!    Y(i) i=1,N  THE CURRENT VALUE OF THE DEPENDENT VARIABLE.
!       N  THE DIMENSION OF THE DEPENDENT VARIABLE.
!       H  THE TIME STEP TO BE USED.
!     HMIN  THE MINIMUM ALLOWED VALUE FOR THE ABSOLUTE VALUE OF H.
!    HMAX  THE MAXIMUM ALLOWED VALUE FOR THE ABSOLUTE VALUE OF H.
!      E1  A TOLERANCE.
!   DERIV  THE NAME OF THE ROUTINE COMPUTING THE VECTOR FIELD (TO
!          BE DECLARED EXTERNAL IN THE CALLING PROGRAM).
!
!  OUTPUT DATA:
!
!       X  THE NEXT VALUE OF THE INDEPENDENT VARIABLE.
!     Y(i) i=1,N  THE ESTIMATED NEXT VALUE FOR THE DEPENDENT
!          VARIABLE.
!       H  THE TIME STEP TO BE USED IN THE NEXT CALL OF THIS
!          ROUTINE.
!
!  AUXILIARY PARAMETERS:
!
!  R,B,F   A MATRIX OF DIMENSIONe
!  ROUTINES USED: DERIV
!
!********************************************************************
  SUBROUTINE GR_RK78 (X,Y,N,H,HMIN0,HMAX0,E0,R,B,F,DERIV)
  implicit none
  integer, parameter :: dp = kind(1.d0) 
  
  integer, intent(in) :: n 
  real(kind=dp), intent(in) ::  hmin0, hmax0, e0
  real(kind=dp), intent(inout) :: x, y(n), h
  external deriv
  
  ! Local Varaibles
  integer :: ii, j, j1, jk, k, l 
  DOUBLE PRECISION  A,BET,ALFA(13), R(13,N),B(N),F(N), &
            BETARK(79),C(11),CP(13),D,DD,E3,E4, FAC
            
  DATA II/0/
  SAVE II,ALFA,BETARK,C,CP
         
  print*, 'ii, h, hmin, hmax, e1', ii, h, hmin0, hmax0, e0 !, tmax
  read(*,*) 
         
! IF (II.NE.0) GOTO 12
  if ( ii == 0) then 
    II=1
    ALFA(1)=0.D0
    ALFA(2)=2.D0/27.D0
    ALFA(3)=1.D0/9.D0
    ALFA(4)=1.D0/6.D0
    ALFA(5)=5.D0/12.D0
    ALFA(6)=.5D0
    ALFA(7)=5.D0/6.D0
    ALFA(8)=1.D0/6.D0
    ALFA(9)=2.D0/3.D0
    ALFA(10)=1.D0/3.D0
    ALFA(11)=1.D0
    ALFA(12)=0.D0
    ALFA(13)=1.D0
    
    BETARK(1)=0.D0
    BETARK(2)=2.D0/27.D0
    BETARK(3)=1.D0/36.D0
    BETARK(4)=1.D0/12.D0
    BETARK(5)=1.D0/24.D0
    BETARK(6)=0.D0
    BETARK(7)=1.D0/8.D0
    BETARK(8)=5.D0/12.D0
    BETARK(9)=0.D0
    BETARK(10)=-25.D0/16.D0
    BETARK(11)=-BETARK(10)
    BETARK(12)=.5D-1
    BETARK(13)=0.D0
    BETARK(14)=0.D0
    BETARK(15)=.25D0
    BETARK(16)=.2D0
    BETARK(17)=-25.D0/108.D0
    BETARK(18)=0.D0
    BETARK(19)=0.D0
    BETARK(20)=125.D0/108.D0
    BETARK(21)=-65.D0/27.D0
    BETARK(22)=2.D0*BETARK(20)
    BETARK(23)=31.D0/300.D0
    BETARK(24)=0.D0
    BETARK(25)=0.D0
    BETARK(26)=0.D0
    BETARK(27)=61.D0/225.D0
    BETARK(28)=-2.D0/9.D0
    BETARK(29)=13.D0/900.D0
    BETARK(30)=2.D0
    BETARK(31)=0.D0
    BETARK(32)=0.D0
    BETARK(33)=-53.D0/6.D0
    BETARK(34)=704.D0/45.D0
    BETARK(35)=-107.D0/9.D0
    BETARK(36)=67.D0/90.D0
    BETARK(37)=3.D0
    BETARK(38)=-91.D0/108.D0
    BETARK(39)=0.D0
    BETARK(40)=0.D0
    BETARK(41)=23.D0/108.D0
    BETARK(42)=-976.D0/135.D0
    BETARK(43)=311.D0/54.D0
    BETARK(44)=-19.D0/60.D0
    BETARK(45)=17.D0/6.D0
    BETARK(46)=-1.D0/12.D0
    BETARK(47)=2383.D0/4100.D0
    BETARK(48)=0.D0
    BETARK(49)=0.D0
    BETARK(50)=-341.D0/164.D0
    BETARK(51)=4496.D0/1025.D0
    BETARK(52)=-301.D0/82.D0
    BETARK(53)=2133.D0/4100.D0
    BETARK(54)=45.D0/82.D0
    BETARK(55)=45.D0/164.D0
    BETARK(56)=18.D0/41.D0
    BETARK(57)=3.D0/205.D0
    BETARK(58)=0.D0
    BETARK(59)=0.D0
    BETARK(60)=0.D0
    BETARK(61)=0.D0
    BETARK(62)=-6.D0/41.D0
    BETARK(63)=-3.D0/205.D0
    BETARK(64)=-3.D0/41.D0
    BETARK(65)=-BETARK(64)
    BETARK(66)=-BETARK(62)
    BETARK(67)=0.D0
    BETARK(68)=-1777.D0/4100.D0
    BETARK(69)=0.D0
    BETARK(70)=0.D0
    BETARK(71)=BETARK(50)
    BETARK(72)=BETARK(51)
    BETARK(73)=-289.D0/82.D0
    BETARK(74)=2193.D0/4100.D0
    BETARK(75)=51.D0/82.D0
    BETARK(76)=33.D0/164.D0
    BETARK(77)=12.D0/41.D0
    BETARK(78)=0.D0
    BETARK(79)=1.D0
    
    C(1)=41.D0/840.D0
    C(2)=0.D0
    C(3)=0.D0
    C(4)=0.D0
    C(5)=0.D0
    C(6)=34.D0/105.D0
    C(7)=9.D0/35.D0
    C(8)=C(7)
    C(9)=9.D0/280.D0
    
    C(10)=C(9)
    C(11)=C(1)
    CP(1)=0.D0
    CP(2)=0.D0
    CP(3)=0.D0
    CP(4)=0.D0
    CP(5)=0.D0
    CP(6)=C(6)
    CP(7)=C(7)
    CP(8)=C(8)
    CP(9)=C(9)
    CP(10)=C(10)
    CP(11)=0.D0
    CP(12)=C(1)
    CP(13)=C(1)
    
!9        CONTINUE

!12       JK=1
  else 
    
    do 
     
      JK=1
!         DO 3 J=1,13
!         DO 6 L=1,N
!6        B(L)=Y(L)
!         A=X+ALFA(J)*H
!         IF(J.EQ.1)GO TO 13
!         J1=J-1
!         DO 4 K=1,J1,1
!         JK=JK+1
!         BET=BETA(JK)*H
!         DO 4 L=1,N
!4        B(L)=B(L)+BET*R(K,L)
!13       CONTINUE
!         CALL DERIV (A,B,N,F)
!         DO 3 L=1,N
!3        R(J,L)=F(L)

      DO   J=1,13 ! 3
        DO   L=1,N !6
          B(L)=Y(L)
        enddo 

        A=X+ALFA(J)*H
         
        IF(J .ne. 1) then !
          J1=J-1
          DO  K=1,J1,1
            JK=JK+1
            BET=BETARK(JK)*H
            DO   L=1,N
              B(L)=B(L)+BET*R(K,L)
            enddo 
          enddo  
        endif

        CALL DERIV (A,B,N,F)
        DO   L=1,N
          R(J,L)=F(L)
        enddo 
      enddo  ! 3     

      D=0
      DD=0
         
      DO L=1,N !- 1
      
        B(L)=Y(L)
        F(L)=Y(L)
        
        DO K=1,11 !--5 
          BET=H*R(K,L)
          B(L)=B(L)+BET*C(K)
          F(L)=F(L)+BET*CP(K)
        enddo  ! -5
         
        F(L)=F(L)+H*(CP(12)*R(12,L)+CP(13)*R(13,L))
        D=D+DABS(F(L)-B(L))
        DD=DD+DABS(F(L))
      enddo ! -1

      D=D/N 
      FAC=1.+DD*1.D-2
      E3=E0*FAC
         
      IF (DABS(H).LT.HMIN0.OR.D.LT.E3)  exit !go to 7
      H=H*0.9D0*(E3/D)**0.125D0
      IF(DABS(H).LT.HMIN0)  H=HMIN0*H/DABS(H)
    enddo 
         
    X=X+H ! -- 7/
    IF(D.LT.E3) D=DMAX1(D,E3/256)
    H=H*0.9D0*(E3/D)**0.125D0
         
    IF(DABS(H).GT.HMAX0)  H=HMAX0*H/DABS(H)
    IF(DABS(H).LT.HMIN0)  H=HMIN0*H/DABS(H)
   
  
    DO L=1,N
      Y(L)=F(L)
    enddo 
         
    B(1)=D
  endif 
  
  RETURN
  END subroutine gr_rk78
         
