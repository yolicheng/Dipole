        SUBROUTINE GR_CJRTBP(X,XMU,CJA)
C******************************************************
C COMPUTATION OF THE JACOBI CONSTANT, CJA, OF A POINT X
C IN THE SPATIAL RTBP.
C IN THIS RTBP THE BIG PRIMARY IS AT (XMU,0,0) WITH
C MASS 1-XMU, AND THE SMALL ONE AT (XMU-1,0,0) WITH
C MASS XMU.
C******************************************************
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION X(6)
        X1=X(1)*X(1)
        X2=X(2)*X(2)
        X3=X(3)*X(3)
        R1=DSQRT((X(1)-XMU)**2+X2+X3)
        R2=DSQRT((X(1)-XMU+1)**2+X2+X3)
        OME=(X1+X2)*.5+(1-XMU)/R1+XMU/R2+.5*XMU*(1-XMU)
        CJA=2*OME-X(4)*X(4)-X(5)*X(5)-X(6)*X(6)
        RETURN
        END

