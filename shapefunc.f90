C
C *********************************************************
C *                                                       *
C * CALCULATE THE SHAPE FUNCTIONS AND THEIR DERIVATIVES   *
C * FOR 6 NODE TRIANGULAR ELEMENTS                        *
C *                                                       *
C *********************************************************     
C
        SUBROUTINE SPFUNC6(SI,ETA,SF,DSF)
      	IMPLICIT  NONE

	  REAL*8,INTENT(IN) :: SI,ETA
	  REAL*8,INTENT(OUT):: SF(8),DSF(2,8)
C
        SF(1)=(1.0D0-SI-ETA)*(1.0D0-2.0D0*SI-2.0D0*ETA)
        SF(2)=SI *(2.0D0*SI -1.0D0)
        SF(3)=ETA*(2.0D0*ETA-1.0D0)
        SF(4)=4.0D0*SI*(1.0D0-SI-ETA)
        SF(5)=4.0D0*SI*ETA
        SF(6)=4.0D0*ETA*(1.0D0-SI-ETA)
c
        DSF(1,1)=4.0D0*SI+4.0*ETA-3.0D0
        DSF(1,2)=4.0D0*SI-1.0D0
        DSF(1,3)=0.0D0
        DSF(1,4)=4.0D0-8.0D0*SI-4.0D0*ETA
        DSF(1,5)=4.0D0*ETA
        DSF(1,6)=-4.0*ETA
C
        DSF(2,1)=4.0D0*SI+4.0D0*ETA-3.0D0
        DSF(2,2)=0.0D0
        DSF(2,3)=4.0D0*ETA-1.0D0
        DSF(2,4)=-4.0*SI
        DSF(2,5)=4.0*SI
        DSF(2,6)=4.0-8.0D0*ETA-4.0*SI
 
C
	RETURN
	END SUBROUTINE SPFUNC6
C
C
C
C
C *********************************************************
C *                                                       *
C * CALCULATE THE SHAPE FUNCTIONS AND THEIR DERIVATIVES   *
C * FOR 8 NODE QUADRILATERIAL ELEMENTS                    *
C *                                                       *
C *********************************************************     
C
        SUBROUTINE SPFUNC8(SI,ETA,SF,DSF)
        IMPLICIT  NONE

	  REAL*8,INTENT(IN) :: SI,ETA
	  REAL*8,INTENT(OUT):: SF(8),DSF(2,8)
C
        SF(1)=-0.25D0*(1.0D0-SI)*(1.0D0-ETA)*(1.0D0+SI+ETA)
        SF(2)= 0.5D0 *(1.0D0-SI*SI)*(1.0-ETA)
        SF(3)= 0.25D0*(1.0D0+SI)*(1.0D0-ETA)*(SI-ETA-1.0D0)
        SF(4)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0D0+SI)
        SF(5)= 0.25D0*(1.0D0+SI)*(1.0D0+ETA)*(SI+ETA-1.0D0)
        SF(6)= 0.5D0 *(1.0-SI*SI)*(1.0D0+ETA)
        SF(7)=-0.25D0*(1.0D0-SI)*(1.0D0+ETA)*(SI-ETA+1.0D0)
        SF(8)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0-SI)            
      
        DSF(1,1)= 0.25D0*(2.0D0*SI+ETA)*(1.0D0-ETA)
        DSF(1,2)=-SI*(1.0D0-ETA)     	
        DSF(1,3)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0-ETA)
        DSF(1,4)= 0.5D0*(1.0D0-ETA*ETA)
        DSF(1,5)= 0.25D0*(2.0*SI+ETA)*(1.0D0+ETA)
        DSF(1,6)=-SI*(1.0D0+ETA)     
        DSF(1,7)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0+ETA)
        DSF(1,8)=-0.5D0*(1.0-ETA*ETA)
C
        DSF(2,1)= 0.25D0*(SI+2.0D0*ETA)*(1.0D0-SI)
        DSF(2,2)=-0.5D0*(1.0-SI*SI)
        DSF(2,3)= 0.25D0*(1.0D0+SI)*(2.0D0*ETA-SI)
        DSF(2,4)=-(1.0D0+SI)*ETA
        DSF(2,5)= 0.25*(1.0D0+SI)*(SI+2.0D0*ETA)
        DSF(2,6)= 0.5D0*(1.0D0-SI*SI)
        DSF(2,7)=-0.25D0*(1.0D0-SI)*(SI-2.0D0*ETA)
        DSF(2,8)=-(1.0D0-SI)*ETA   
C
	RETURN
	END SUBROUTINE SPFUNC8


C
C *********************************************************
C *                                                       *
C * CALCULATE THE SHAPE FUNCTIONS AND THEIR DERIVATIVES   *
C * FOR 6 NODE TRIANGULAR ELEMENTS                        *
C *                                                       *
C *********************************************************     
C
        SUBROUTINE SPFUNC6_1(SI,ETA,SF,DSF,DDSF)
      	IMPLICIT  NONE

	  REAL*8,INTENT(IN) :: SI,ETA
	  REAL*8,INTENT(OUT):: SF(8),DSF(2,8),DDSF(3,8)
C
        SF(1)=(1.0D0-SI-ETA)*(1.0D0-2.0D0*SI-2.0D0*ETA)
        SF(2)=SI *(2.0D0*SI -1.0D0)
        SF(3)=ETA*(2.0D0*ETA-1.0D0)
        SF(4)=4.0D0*SI*(1.0D0-SI-ETA)
        SF(5)=4.0D0*SI*ETA
        SF(6)=4.0D0*ETA*(1.0D0-SI-ETA)
! ------------------------------------------
!
        DSF(1,1)=4.0D0*SI+4.0*ETA-3.0D0
        DSF(1,2)=4.0D0*SI-1.0D0
        DSF(1,3)=0.0D0
        DSF(1,4)=4.0D0-8.0D0*SI-4.0D0*ETA
        DSF(1,5)=4.0D0*ETA
        DSF(1,6)=-4.0*ETA
C
        DSF(2,1)=4.0D0*SI+4.0D0*ETA-3.0D0
        DSF(2,2)=0.0D0
        DSF(2,3)=4.0D0*ETA-1.0D0
        DSF(2,4)=-4.0*SI
        DSF(2,5)=4.0*SI
        DSF(2,6)=4.0-8.0D0*ETA-4.0*SI

!  ---------------------------
        DDSF(1,1)= 4.D0
        DDSF(1,2)= 4.0D0     	
        DDSF(1,3)= 0.0D0
        DDSF(1,4)=-8.0D0
        DDSF(1,5)= 0.0D0
        DDSF(1,6)= 0.0D0     

        DDSF(2,1)= 4.0D0
        DDSF(2,2)= 0.0D0
        DDSF(2,3)= 4.0D0
        DDSF(2,4)= 0.0D0
        DDSF(2,5)= 0.0D0
        DDSF(2,6)=-8.0D0

        DDSF(3,1)= 4.0D0
        DDSF(3,2)= 0.0D0   	
        DDSF(3,3)= 0.0D0
        DDSF(3,4)=-4.0d0
        DDSF(3,5)= 4.0D0
        DDSF(3,6)=-4.0d0     

	  RETURN
	 END SUBROUTINE SPFUNC6_1
C



        SUBROUTINE SPFUNC8_1(SI,ETA,SF,DSF,DDSF)
        IMPLICIT  NONE

        REAL*8,INTENT(IN) :: SI,ETA
        REAL*8,INTENT(OUT):: SF(8),DSF(2,8),DDSF(3,8)
!
        SF(1)=-0.25D0*(1.0D0-SI)*(1.0D0-ETA)*(1.0D0+SI+ETA)
        SF(2)= 0.5D0 *(1.0D0-SI*SI)*(1.0-ETA)
        SF(3)= 0.25D0*(1.0D0+SI)*(1.0D0-ETA)*(SI-ETA-1.0D0)
        SF(4)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0D0+SI)
        SF(5)= 0.25D0*(1.0D0+SI)*(1.0D0+ETA)*(SI+ETA-1.0D0)
        SF(6)= 0.5D0 *(1.0-SI*SI)*(1.0D0+ETA)
        SF(7)=-0.25D0*(1.0D0-SI)*(1.0D0+ETA)*(SI-ETA+1.0D0)
        SF(8)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0-SI)            
! ------------------------------------------------------------
      
        DSF(1,1)= 0.25D0*(2.0D0*SI+ETA)*(1.0D0-ETA)
        DSF(1,2)=-SI*(1.0D0-ETA)     	
        DSF(1,3)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0-ETA)
        DSF(1,4)= 0.5D0*(1.0D0-ETA*ETA)
        DSF(1,5)= 0.25D0*(2.0*SI+ETA)*(1.0D0+ETA)
        DSF(1,6)=-SI*(1.0D0+ETA)     
        DSF(1,7)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0+ETA)
        DSF(1,8)=-0.5D0*(1.0-ETA*ETA)
!
        DSF(2,1)= 0.25D0*(SI+2.0D0*ETA)*(1.0D0-SI)
        DSF(2,2)=-0.5D0*(1.0-SI*SI)
        DSF(2,3)= 0.25D0*(1.0D0+SI)*(2.0D0*ETA-SI)
        DSF(2,4)=-(1.0D0+SI)*ETA
        DSF(2,5)= 0.25*(1.0D0+SI)*(SI+2.0D0*ETA)
        DSF(2,6)= 0.5D0*(1.0D0-SI*SI)
        DSF(2,7)=-0.25D0*(1.0D0-SI)*(SI-2.0D0*ETA)
        DSF(2,8)=-(1.0D0-SI)*ETA   
! --------------------------------------------------------------

        DDSF(1,1)=0.5D0*(1.0D0-ETA)
        DDSF(1,2)=-(1.0D0-ETA)     	
        DDSF(1,3)= 0.5D0*(1.0D0-ETA)
        DDSF(1,4)= 0.0D0
        DDSF(1,5)= 0.5D0*(1.0D0+ETA)
        DDSF(1,6)=-(1.0D0+ETA)     
        DDSF(1,7)= 0.5D0*(1.0D0+ETA)
        DDSF(1,8)= 0.0D0

        DDSF(2,1)= 0.50D0*(1.0D0-SI)
        DDSF(2,2)= 0.0D0
        DDSF(2,3)= 0.50D0*(1.0D0+SI)
        DDSF(2,4)=-(1.0D0+SI)
        DDSF(2,5)= 0.50D0*(1.0D0+SI)
        DDSF(2,6)= 0.0D0
        DDSF(2,7)= 0.50D0*(1.0D0-SI)
        DDSF(2,8)=-(1.0D0-SI) 

        DDSF(3,1)= 0.25D0*(-2.0D0*SI-ETA*2.0D0+1.0D0)
        DDSF(3,2)= SI     	
        DDSF(3,3)= 0.25D0*(-2.0D0*SI+ETA*2.0D0-1.0D0)
        DDSF(3,4)= -ETA
        DDSF(3,5)= 0.25D0*(2.0D0*SI+ETA*2.0D0+1.0D0)
        DDSF(3,6)=-SI     
        DDSF(3,7)= 0.25D0*(2.0D0*SI-ETA*2.0D0-1.0D0)
        DDSF(3,8)= ETA

        RETURN
	  END SUBROUTINE SPFUNC8_1
