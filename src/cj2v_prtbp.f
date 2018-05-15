        SUBROUTINE cj2v_prtbp(X, CJA,  indv, isvy )
C******************************************************
C Given the Jacobi Constant, and the three values of (x,y, vx(or vy)),  compute 
c the positive value of the missing velocity 
c To make it general, introduce indv as the index of the missing components 
c we also need to make sure the value of the missing component as input is zero
c 

C MASS ratio XMU.
C******************************************************
        use rtbpconst_mod, only : mu 
        
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION X(4)
        
        ! default value set as 1, succeed to compute vy 
        isvy = 1 
        
        X1=X(1)*X(1)
        X2=X(2)*X(2)
        
        R1  = DSQRT((X(1)-MU)**2+X2)
        R2  = DSQRT((X(1)-MU+1)**2+X2)
        OME = (X1+X2)*.5+(1-MU)/R1+MU/R2+.5*MU*(1-MU)     
         
!       CJA = 2*OME- X(3)*X(3) - X(4)*X(4)   

        vy2 = 2*OME - cja - X(3)*X(3) - X(4)*X(4) 
        
        if(vy2 < 0.d0) then 
          isvy = 0 
          print*, 'vy**2 < 0!', vy2;
!           read*
!          x(indv) = 0.d0 
          return 
        endif 
        
        vy = dsqrt(vy2)
        x(indv) = vy 
        
!        X3=X(3)*X(3)
!        R1=DSQRT((X(1)-XMU)**2+X2+X3)
!        R2=DSQRT((X(1)-XMU+1)**2+X2+X3)
!        OME=(X1+X2)*.5+(1-XMU)/R1+XMU/R2+.5*XMU*(1-XMU)
!        CJA=2*OME-X(4)*X(4)-X(5)*X(5)-X(6)*X(6)
        RETURN
        END

