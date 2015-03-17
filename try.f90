	module test
	
	INTEGER :: NTIME,ITIME,ISYS,NSYS,IORDER
	INTEGER :: NELEM,NELEMB,NELEMF,NNODE,NNODED,NNB,NNBD,NNF,NNTCH
	INTEGER Node_Singular,IPOL !IPOL added March.17

	REAL*8 G,RHO,PI,PI4
        REAL*8 H,AMP,BETA,W1,V
        REAL*8 WK,TPER
	REAL*8 Tstep,TIME,TimeRK,RAMPF,RAMPV
        real*8 Line_sum,Area_sum,Area0_sum,Area1_sum
	real*8 Area2_sum,SUMF_1,SUMF_2


   	INTEGER, ALLOCATABLE:: NCN(:),IETYPE(:),NCON(:,:),NCOND(:,:)
!	ncn node number in each element
! IETYPE: type of the element; =1, on body surface; =2, on free surface
!	NCON second part of mesh data
!	NCOND third part of mesh data
	INTEGER, ALLOCATABLE:: NNORMC(:)
        INTEGER, ALLOCATABLE:: NODELE(:,:),NODNOE(:)
	INTEGER, ALLOCATABLE:: NODELJ(:,:),NODQUA(:)
   	INTEGER, ALLOCATABLE:: NCONB(:,:),NCONDB(:,:)
	

	REAL*8,ALLOCATABLE:: XYZE(:,:,:),DXYZE(:,:,:),XYZ(:,:),DXYZ(:,:)
	!XYZ, node coordinate 3 row x col
	REAL*8,ALLOCATABLE:: DAMPE(:,:),DAMPF(:)
	REAL*8,ALLOCATABLE:: XYZB(:,:),DXYZB(:,:)

        REAL*8,ALLOCATABLE:: TXYZE(:,:,:)


	REAL*8,ALLOCATABLE:: BKN(:,:),BKN_O(:,:),UNKN(:,:),UNKN_O(:,:)
	REAL*8,ALLOCATABLE:: ET(:,:),ET_O(:,:),DPDT(:,:)

	

        REAL*8,ALLOCATABLE:: SAMB(:,:,:),SAMBXY(:,:,:)
	REAL*8,ALLOCATABLE:: DSAMB(:,:,:),ANGLE(:)
!
! SAMB:
! SAMBXY: Coordinates of Gaussin points
! DSAMB:  Normal direvatives at Gaussian points
! ANGLE:  Solid angle
!
	REAL*8,ALLOCATABLE:: DH(:,:,:),DP(:,:,:),Dposi(:,:)
!
! DH:
! DP:
! Dposi:
!
!  For linear equations
!
	INTEGER, ALLOCATABLE:: INDX(:,:)
	REAL*8,ALLOCATABLE::   AMATA(:,:,:),CMATA(:,:,:),BMATA(:,:)
!       REAL*8,ALLOCATABLE::   LEFT(:,:,:),RIGHT(:,:,:),RIGHT1(:,:)

        DATA G /9.807d0/
	DATA PI /3.141592653589793238d0/
	DATA RHO /1023.0d0/

! ------------------------------------
! Temporary arrayes
!
   	INTEGER, ALLOCATABLE:: NNORMN(:)
	REAL*8,  ALLOCATABLE:: XYZTP(:,:),DXYZTP(:,:),DAMPTP(:)
!
	DATA AL/1.2/
	DATA MAX_LEVEL/ 8/
	end module

!
! Variable for potential, force and body motion
! =============================================
!
        MODULE PVar_mod
!
	INTEGER  NNT
	PARAMETER (NNT=2000)    

        REAL*8 XC,YC,ZC,XTC,YTC,ZTC
        REAL*8 ARE,XF,YF,XK2,YK2,XCF
        REAL*8 VOLM,XB,YB,ZB

        REAL*8 AMAS(6,6),BDMP(6,6)  
        REAL*8 RMAS(6,6),CDMP(6,6),CRS(6,6),STKM(6,6),XIA(3,3)
!
        REAL*8 FORCEW(6),FORCE0(6),FORSCD(6),AMPJ(6)
	REAL*8 FORCE(6),FORCE_O(6)
!
        REAL*8 TRMAS(6,6),VISC(6,6)

        REAL*8  DISP(6),DSDT(6),DISP_O(6)
	REAL*8  DSDT_O(6),DSDDTL(6),TRXYZ(3,3)   

        REAL*8  RESPR(6),RESPI(6)  


        END MODULE PVar_mod

      MODULE SEBSM_MOD
	REAL*8,ALLOCATABLE:: A_SEBSM(:,:),B_SEBSM(:)
      INTEGER N_SEBSM
      END MODULE SEBSM_MOD

!
! =========================================
!  Varibles used in the Tri-pole transform 
!							  
       MODULE TRVar_mod

	 INTEGER NOSAMP
	 REAL*8 XYNOD(3,50),DXYNOD(6,50),SAMNOD(50,0:8)

       END MODULE TRVar_mod

