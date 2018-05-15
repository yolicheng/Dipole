        SUBROUTINE GR_FOUN(F, N, M, CSF,SIF)
! ** NOTE ** 
!  Modify the fourier series to the form that I think it is right.         
C****************************************************************
C Computes the Fourier series of a function F evaluated at N
C equally spaced points in a period, F(k+1)=F(Xk) where
C Xk=k*(PERI/N), k=0..N-1 (N must be odd).
C M is the number of harmonics to be computed.

C It gives:
C  CSF(i) i=0,M containig the cosinus coefficients.
C  SIF(i) i=0,M containig the sinus coefficients.

!C        both of the term (SIN COS)(2*PI*i/ W), W is the freq. !--- I think it is wrong 

!C This is:   !
!C  F(t)= csf(0)+SUM i=1_M (CSF(i)* cos(i* 2pi/prd ) + SIF(i)* sin(2PIi/) ) ! --- I think this is right


! Update as the right one 

C        both of the term (SIN COS)(2*PI*i/ M), M is the number of harmonics 

C This is:
C  F(t)= csf(0) + SUM i = 1_M (CSF(i)*cos(2PIi/M) + SIF(i)*sin(2PIi/M))


C****************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (MMAX=500000)
        DIMENSION F(N),CSF(0:M),SIF(0:M), CO(MMAX),SI(MMAX)
        DATA CO(1)/1.D0/,SI(1)/0.D0/
        SAVE CO, SI
        
        PI2=8.D0*DATAN(1.D0)
        CC=DCOS(PI2/N)
        CS=DSIN(PI2/N)
        
        CO(2)=CO(1)
        SI(2)=SI(1)
        
        DO 10 I=1,M
        CA=CO(2)
        SA=SI(2)
        CO(2)=CA*CC-SA*CS
        SI(2)=CA*CS+CC*SA
        
        DO 1 K=3,N
        CO(K)=CO(K-1)*CO(2)-SI(K-1)*SI(2)
        SI(K)=CO(K-1)*SI(2)+SI(K-1)*CO(2)
1       CONTINUE

        PC=0.D0
        PS=0.D0
        DO 12 J=1,N
        PC=PC+F(J)*CO(J)
        PS=PS+F(J)*SI(J)
12      CONTINUE

        CSF(I)=2.0D0*PC/N
        SIF(I)=2.0D0*PS/N
10      CONTINUE

        CSF(0)=0.D0
        DO 14 I=1,N
        CSF(0)=CSF(0)+F(I)
14      CONTINUE

        CSF(0)=CSF(0)/N
        SIF(0)=0.D0
        RETURN
        END
        
