	SUBROUTINE ADAMS(NP,X,HH,CAM)
c
C++++++++++++++++++++++++++  
 
C***********************************************************************
C
C     ADAMS-BASHFORTH FORMULAS
C          INPUT PARAMETERS:
C     NP         NUMBER OF COMPUTED PERIODIC ORBITS
C     X(*)       INITIAL CONDITIONS OF THE LAST P.O.
C     HH         INTEGRATION STEP ALONG ARC LENGTH OF SOLUTION LOCUS
C     CAM(*,*)   (SEE SUBROUTINE CAMP)
C          OUTPUT PARAMETERS:
C     X(*)       NEW APROXIMATED INITIAL CONDITIONS FOR P.O.
C
C***********************************************************************

        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION X(3),CAM(4,3)
        
        DO 10 I=1,3
        IF(NP.GT.4)GO TO 4
        GO TO (1,2,3,4),NP
1       X(I)=X(I)+HH*CAM(1,I)
        GO TO 10
2       X(I)=X(I)+.5*HH*(3.*CAM(2,I)-CAM(1,I))
        GO TO 10
3       X(I)=X(I)+HH*(23.*CAM(3,I)-16.*CAM(2,I)+5.*CAM(1,I))/12
        GO TO 10
4       X(I)=X(I)+HH*(55.*CAM(4,I)-59.*CAM(3,I)+37.*CAM(2,I)-9.*CAM
     -  (1,I))/24.
10      CONTINUE
        RETURN
        END  
