      SUBROUTINE GR_MPOWER(DAT,M,VEP,VAP,IS)
C**********************************************************************
C  POWER METHOD IN A PRODUCT OF M MATRICES IN DIRECT OR INVERSE SENSE. 
C  GIVEN MATRICES 6*6 STORED IN DAT(*,*,k), k=1,NIT, IT COMPUTES THE 
C  LOG10 OF THE DOMINANT EIGENVALUE, VAP, AND THE ASSOCIATED 
C  EIGENVECTOR, VEP, OF THE MATRIX A WHERE A IS:
C   A=DAT(*,*,M)*DAT(*,*,M-1)*...*DAT(*,*,1) IF IS=1, or,
C   A=DAT(*,*,1)*DAT(*,*,2)*...*DAT(*,*,M)   IF IS=-1.
C  USING THE POWER METHOD.
C
C  ROUTINES USED: NONE
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C PRECISION AND MAX NUMBER OF ITERATIONS
      PARAMETER (PREC=1.D-14,ITM=20)
      DIMENSION DAT(6,6,*),VEP(6), BUX(6),CUX(6)
      KI=1
      KF=M
      INC=1
      IF (IS.EQ.-1) THEN
      KI=M
      KF=1
      INC=-1
      ENDIF
      IT=0
      DO 7 I=1,6
      VEP(I)=1.D0/DSQRT(6.D0)
7     CONTINUE
6     VAP=0.D0
      DO 1 K=KI,KF,INC
      VN=0.D0
      DO 10 I=1,6
      BUX(I)=0.D0
      DO 12 J=1,6
      BUX(I)=BUX(I)+DAT(I,J,K)*VEP(J)
12    CONTINUE
      VN=VN+BUX(I)*BUX(I)
10    CONTINUE
      VN=DSQRT(VN)
      VAP=VAP+DLOG10(VN)
      DO 2 I=1,6
      VEP(I)=BUX(I)/VN
2     CONTINUE
1     CONTINUE
      IF (IT.EQ.0) THEN
      VAPV=VAP
      DO 3 I=1,6
      CUX(I)=VEP(I)
3     CONTINUE
      ELSE
      DVAP=DABS((VAP-VAPV)/VAP)
      VN=0.D0
      DO 4 I=1,6
      CUX(I)=VEP(I)-CUX(I)
      VN=VN+CUX(I)*CUX(I)
4     CONTINUE
      VN=DSQRT(VN)
      IF (IT.GT.ITM-5) THEN
      WRITE (*,*) 'GR_MPOWER. REL. ERR IN VAP: ',DVAP
      WRITE (*,*) 'GR_MPOWER. ABS. ERR IN VECTOR: ',VN
      ENDIF
      IF (DVAP.LT.PREC.AND.VN.LT.PREC) RETURN
      VAPV=VAP
      DO 5 I=1,6
      CUX(I)=VEP(I)
5     CONTINUE
      ENDIF
      IT=IT+1
      IF (IT.GT.ITM) THEN
      WRITE (*,*) 'GR_MPOWER. POWER METHOD SEEMS NOT CONVERGENT'
      WRITE (*,*) '        VALUE OF IS: ',IS
      STOP
      ENDIF
      
      GOTO 6
      END

