C****************************************************************************  
C
C  SUBROUTINE GR_RTBP 
C                                                                              
C  VECTORFIELD OF THE SPATIAL RTBP AND ITS VARIATIONAL EQUATIONS IF DESIRED.
C
C  IN THIS RTBP THE BIG PRIMARY IS LOCATED AT (XMU,0,0) WITH MASS 1-XMU AND
C  THE SMALL ONE AT (XMU-1,0,0) WITH MASS XMU.
C
C  INPUT PARAMETERS:                                                   
C
C    T     RTBP TIME
C    X(*)  RTBP-VARIATIONAL COORDINATES (X(1),X(2),...,X(42))
C    N     NUMBER OF EQUATIONS (6 or 42 ACCORDING IF WE WANT JUST THE
C          VECTORFIELD OR MOREOVER THE VARIATIONAL FLOW).

C    XMU   RTBP MASS PARAMETER THROUGH Module massrat_rtbp_mod 
C
C  OUTPUT PARAMETERS:
C
C    F(*)  VECTORFIELD. THE FIRST 6 COMPONENTS CORRESPOND TO THE RTBP
C          EQUATIONS, AND THE REMAINING ONES TO THE VARIATIONAL EQUATIONS
C          STORED BY COLUMNS. 
C
C  NOTE: IF N EQUALS 6 THE VARIATIONAL EQUATIONS ARE SKIPPED 
C
C***************************************************************************** 
        SUBROUTINE GR_RTBP(T,X,N,F)
        
        use massrat_rtbp_mod, only : mu => xmu  
        
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION X(N),F(N)
        
!        COMMON /RTBPC/XMU
        
        F(1)=X(4)
        F(2)=X(5)
        F(3)=X(6)
        
        
        Y1=X(1)-XMU
        Y12=Y1*Y1
        Y2=X(2)*X(2)
        Y3=X(3)*X(3)
        R1=Y12+Y2+Y3
        R1A=R1*DSQRT(R1)
        R15=R1A*R1
        R1=R1A
        R2=(Y1+1.D0)*(Y1+1.D0)+Y2+Y3
        R2A=R2*DSQRT(R2)
        R25=R2A*R2
        R2=R2A
        P1=(1.D0-XMU)/R1
        P2=XMU/R2
        Q=-1.D0*(P1+P2)
        F(4)=2.D0*X(5)+X(1)-(X(1)-XMU)*P1-(X(1)-XMU+1.D0)*P2
        F(5)=-2.D0*X(4)+X(2)*(1.D0+Q)
        F(6)=X(3)*Q
        IF (N.EQ.6) RETURN
        RR1=3.D0*(1.D0-XMU)/R15
        RR2=3.D0*XMU/R25
        QP=RR1*Y1+RR2*(Y1+1.D0)
        QPP=(RR1+RR2)*X(2)
        O11=1.D0+Q+RR1*Y12+RR2*(Y1+1.D0)*(Y1+1.D0)
        O12=QP*X(2)
        O13=QP*X(3)
        O22=1.D0+Q+(RR1+RR2)*Y2
        O23=QPP*X(3)
        O33=Q+(RR1+RR2)*Y3
        DO 1 J=1,6
        J6=6*J
        DO 2 I=1,3
        F(J6+I)=X(J6+I+3)
2       CONTINUE
        F(J6+4)=O11*X(J6+1)+O12*X(J6+2)+O13*X(J6+3)+2.D0*X(J6+5)
        F(J6+5)=O12*X(J6+1)+O22*X(J6+2)+O23*X(J6+3)-2.D0*X(J6+4)
        F(J6+6)=O13*X(J6+1)+O23*X(J6+2)+O33*X(J6+3)
1       CONTINUE
        RETURN
        END

