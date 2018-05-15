C********************************************************************
C
C  THIS ROUTINE IS AN IMPLEMENTATION OF A RUNGE-KUTTA-FEHLBERG
C  METHOD OF ORDERS 7 AND 8. USING A TOTAL OF 13 STEPS (AND
C  EVALUATIONS OF THE VECTORFIELD) IT COMPUTES TWO DIFFERENT
C  ESTIMATIONS OF THE NEXT POINT. THE DIFFERENCE BETWEEN BOTH
C  ESTIMATIONS (WITH LOCAL ERRORS OF ORDER 8 AND 9) IS COMPUTED
C  AND THE L1 NORM IS OBTAINED. THIS NORM IS DIVIDED BY N (THE
C  NUMBER OF EQUATIONS). THE NUMBER OBTAINED IN THIS WAY IS REQUIRED
C  TO BE LESS THAN A GIVEN TOLERANCE E1 TIMES (1+0.01*DD) WHERE DD
C  IS THE L1 NORM OF THE POINT COMPUTED TO ORDER 8. IF THIS
C  REQUIREMENT IS SATISFIED THE ORDER 8 ESTIMATION IS TAKEN AS THE
C  NEXT POINT. IF NOT, A SUITABLE VALUE OF THE STEP H IS OBTAINED
C  AND THE COMPUTATION IS STARTED AGAIN.
C  IN ANY CASE, WHEN THE NEXT POINT IS COMPUTED, A PREDICTION OF
C  THE STEP H, TO BE USED IN THE NEXT CALL OF THE ROUTINE, IS
C  DONE.
C
c  TODO: 
c     there is only one more thing to consider, when the velocity is high
c     at some point , for instance, close to the origin.
c     We will get stuck in the dead loop to compute available h for this step 
c     as well as a prediction for the next one ....

c  Ask Alex about this problem..... and debug gr_rk78.f   plob.f90  and plob_n.f90 
c  please alway keep in mind the time sense of integration....

C  INPUT DATA:
C
C       X  CURRENT VALUE OF THE INDEPENDENT VARIABLE.
C    Y(i) i=1,N  THE CURRENT VALUE OF THE DEPENDENT VARIABLE.
C       N  THE DIMENSION OF THE DEPENDENT VARIABLE.
C       H  THE TIME STEP TO BE USED
C          Comment by Yu: but cannot alway be achieved, if it violates the tolerance
c          use new and suitable value of h instead, the real value of h is x(output) - x(input)

C     hmin  THE MINIMUM ALLOWED VALUE FOR THE ABSOLUTE VALUE OF H.
C    HMAX  THE MAXIMUM ALLOWED VALUE FOR THE ABSOLUTE VALUE OF H.
C      E1  A TOLERANCE.
C   DERIV  THE NAME OF THE ROUTINE COMPUTING THE VECTOR FIELD (TO
C          BE DECLARED EXTERNAL IN THE CALLING PROGRAM).
C
C  OUTPUT DATA:
C
C       X  THE NEXT VALUE OF THE INDEPENDENT VARIABLE.
C     Y(i) i=1,N  THE ESTIMATED NEXT VALUE FOR THE DEPENDENT
C          VARIABLE.
C       H  THE TIME STEP TO BE USED IN THE NEXT CALL OF THIS
C          ROUTINE
C
C  AUXILIARY PARAMETERS:
C
C  R,B,F   A MATRIX OF DIMENSION (13,N), AND TWO VECTORS OF
C          DIMENSION N TO BE USED AS A WORKING SPACE.
C
C  ROUTINES USED: DERIV

C TODO: pay attention to the extreme case of dead loop 

C********************************************************************
         SUBROUTINE GR_RK78 (X,Y,N,H,hmin,HMAX,E1,R,B,F,DERIV)
         IMPLICIT REAL*8 (A-H, O-Z)
         DOUBLE PRECISION X,Y(N),H,E1,R(13,N),B(N),F(N),A,BET,ALFA(13),
     &   BETA(79),C(11),CP(13),D,DD,E3,hmin,HMAX,FAC
     
         DATA II/0/
         SAVE II,ALFA,BETA,C,CP
         
         ibug = 0
         nincr = 0 
         iter = 0
         
         ! ** NOTE **  Debug 2016-09-24 09:48:28 
         !  never try to fix hmin = hmax, this will violate the error tolerance
         !  and stuck the routine into dead loop for smaller step for line 230... 
          
         if(ibug == 1 .and. dabs(hmin) > 1.d-3) then  
          
          ! this may happen in the last step for Newton Method, we use dabs(dh) as hmin
           print*, '---- start gr_rk78 ----- '
           print*, 'hmin, hmax, e1', hmin, hmax, e1
           print*; read*; !ck 
         endif 
         
         IF (II.NE.0) GOTO 12
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
         BETA(1)=0.D0
         BETA(2)=2.D0/27.D0
         BETA(3)=1.D0/36.D0
         BETA(4)=1.D0/12.D0
         BETA(5)=1.D0/24.D0
         BETA(6)=0.D0
         BETA(7)=1.D0/8.D0
         BETA(8)=5.D0/12.D0
         BETA(9)=0.D0
         BETA(10)=-25.D0/16.D0
         BETA(11)=-BETA(10)
         BETA(12)=.5D-1
         BETA(13)=0.D0
         BETA(14)=0.D0
         BETA(15)=.25D0
         BETA(16)=.2D0
         BETA(17)=-25.D0/108.D0
         BETA(18)=0.D0
         BETA(19)=0.D0
         BETA(20)=125.D0/108.D0
         BETA(21)=-65.D0/27.D0
         BETA(22)=2.D0*BETA(20)
         BETA(23)=31.D0/300.D0
         BETA(24)=0.D0
         BETA(25)=0.D0
         BETA(26)=0.D0
         BETA(27)=61.D0/225.D0
         BETA(28)=-2.D0/9.D0
         BETA(29)=13.D0/900.D0
         BETA(30)=2.D0
         BETA(31)=0.D0
         BETA(32)=0.D0
         BETA(33)=-53.D0/6.D0
         BETA(34)=704.D0/45.D0
         BETA(35)=-107.D0/9.D0
         BETA(36)=67.D0/90.D0
         BETA(37)=3.D0
         BETA(38)=-91.D0/108.D0
         BETA(39)=0.D0
         BETA(40)=0.D0
         BETA(41)=23.D0/108.D0
         BETA(42)=-976.D0/135.D0
         BETA(43)=311.D0/54.D0
         BETA(44)=-19.D0/60.D0
         BETA(45)=17.D0/6.D0
         BETA(46)=-1.D0/12.D0
         BETA(47)=2383.D0/4100.D0
         BETA(48)=0.D0
         BETA(49)=0.D0
         BETA(50)=-341.D0/164.D0
         BETA(51)=4496.D0/1025.D0
         BETA(52)=-301.D0/82.D0
         BETA(53)=2133.D0/4100.D0
         BETA(54)=45.D0/82.D0
         BETA(55)=45.D0/164.D0
         BETA(56)=18.D0/41.D0
         BETA(57)=3.D0/205.D0
         BETA(58)=0.D0
         BETA(59)=0.D0
         BETA(60)=0.D0
         BETA(61)=0.D0
         BETA(62)=-6.D0/41.D0
         BETA(63)=-3.D0/205.D0
         BETA(64)=-3.D0/41.D0
         BETA(65)=-BETA(64)
         BETA(66)=-BETA(62)
         BETA(67)=0.D0
         BETA(68)=-1777.D0/4100.D0
         BETA(69)=0.D0
         BETA(70)=0.D0
         BETA(71)=BETA(50)
         BETA(72)=BETA(51)
         BETA(73)=-289.D0/82.D0
         BETA(74)=2193.D0/4100.D0
         BETA(75)=51.D0/82.D0
         BETA(76)=33.D0/164.D0
         BETA(77)=12.D0/41.D0
         BETA(78)=0.D0
         BETA(79)=1.D0
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
         
9        CONTINUE

         iter = iter + 1   
         
12       JK=1
         DO 3 J=1,13
         DO 6 L=1,N
6        B(L)=Y(L)
         A=X+ALFA(J)*H
         
         IF(J.EQ.1)GO TO 13
         J1=J-1
         DO 4 K=1,J1,1
         JK=JK+1
         BET=BETA(JK)*H
         
         DO 4 L=1,N
4        B(L)=B(L)+BET*R(K,L)

13       CONTINUE

         CALL DERIV (A,B,N,F)
         
         DO 3 L=1,N
3        R(J,L)=F(L)

         D=0
         DD=0
         
         DO 1 L=1,N
         B(L)=Y(L)
         F(L)=Y(L)
         DO 5 K=1,11
         BET=H*R(K,L)
         B(L)=B(L)+BET*C(K)
5        F(L)=F(L)+BET*CP(K)
         F(L)=F(L)+H*(CP(12)*R(12,L)+CP(13)*R(13,L))
         D=D+DABS(F(L)-B(L))
1        DD=DD+DABS(F(L))

         D=D/N
         FAC=1.d0+DD*1.D-2
         E3=E1*FAC
         

         if(ibug == 1) then 
           print*, 'f: iter, e3, d,h', iter, e3, d,h;   read*  !ck
         endif
         
         !TODO 2016-12-01 09:57:52 :
         !  one bug is the dead loop due to the needed step h violates 
         !  the error tolerance, that is, h < hmin
         !  but at the same time we ask alway dabs(h) >= hmin 
         
         ! ** NOTE **  
         ! --1 Never update hmin, hmax and e, it is not allowed !!!!
         
         ! --2 A failed trial:-- 
         ! allow temporarily the lower bound of hmin to be 1.d-16
         ! (the limit of double precision)
         ! and the new value hmin will remain for the whole integration process
         ! until we start a new orbit, since all the values of hmin, hmax 
         ! will be reassigned for a new orbit
         
         ! -- Too small stepsize will lead to too big accumulated integration error
         ! so never allow h to be smaller than 1.d-6 
         
         ! --3 One option is: 
         ! check the number of iteration to obtainan appropriate value of h
         ! when it is bigger than 6, we set h to be smaller than hmin, 
         ! and do one more step check in the callee to decide if this call of gr_rk78
         ! is failed or not.
         
         
         ! To avoid too many iterations for an appropriate h, we ask no more than 6
         ! and we check the value of h outside this routine  
         if(iter >= 6 .and. D .GT. E3 ) then 
           print*, 'Integration failed! h<hmin required'; print*; read*
           
           ! return and let the callee to check the minimum step 
           h = hmin/10.d0
           print*, 't, e3, d, hmin, new h', x, e3, d, hmin, h ;  !ck
!           print*; read*
           return  
         endif  
         
         IF (DABS(H).LT.hmin.OR.D.LT.E3) GO TO 7
         H=H*0.9D0*(E3/D)**0.125D0
         
         IF(DABS(H).LT.hmin) H=hmin*H/DABS(H)
         
         GOTO 9
         
 7       X=X+H

         IF(D.LT.E3)D=DMAX1(D,E3/256)
         H=H*0.9D0*(E3/D)**0.125D0
         IF(DABS(H).GT.HMAX) H=HMAX*H/DABS(H)
         IF(DABS(H).LT.hmin) H=hmin*H/DABS(H)
         
         DO 10 L=1,N
10       Y(L)=F(L)

         B(1)=D
           
         RETURN
         END
