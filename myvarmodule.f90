C  
C **************************************************************** 
C *                                                              * 
C *  Programme for the first order diffraction and radiation     * 
C *  wave analysis of 3-dimensional bodies                       * 
C *															   *
C * Revised 											           *
C *  4 into NSYS	 Sept.16, 2000                                 *
C *                                                              * 
C **************************************************************** 
C 
        SUBROUTINE main
!        
	USE test
	
!
        IMPLICIT  NONE  
	  INTEGER I,IFWKO,NDRFA,IS,IND,M,IP,K,NDMAX !IPOL removed March.17
	  INTEGER NTnum,IFLAG_T
	  !COMMON /DATA/ IPOL
!
	  REAL*8  WL,AMFJ,R,EX(4),Alpha,WAVES  
	  REAL*8  FAMPR(6),FAMPI(6),FORCER(6),FORCEI(6)
	  REAL*8  PL_AMP(6),FORAMP
	  REAL*8  FCD_AMR,  FCD_AMI
!
! IFLAG_T,FCD_AMR, FCD_AMI: not used in this program
!
        DATA EX /  1.,  1., -1., -1./ 
!
! FAMPR, FAMPI:   real and imaginary parts of body motion 
!                 from frequency domain calculation
! FORCER, FORCEI: real and imaginary parts of wave force 
!                 from frequency domain calculation
! PL_AMP:         amplitude of plotting variable
C
C ----------------------------------------
C Input data files
        OPEN(1, FILE='INPUT/DATIN.txt',      STATUS='OLD') 
        OPEN(2, FILE='INPUT/DATBDMS.txt',    STATUS='OLD') 
        OPEN(3, FILE='INPUT/DATFSMS.txt',    STATUS='OLD') 
        OPEN(7, FILE='INPUT/DATMASS.txt',    STATUS='OLD')  
        OPEN(4, FILE='INPUT/SNODE.txt',      STATUS='OLD') 
!       This is for input of src point
! ---------------------------------------------
!	  OPEN(8, FILE='INPUT\SOLIDANGLE.txt',    STATUS='OLD') 
!	  open(67, FILE='INPUT\M.txt',    STATUS='OLD')
!	read(67,*) wl
!	stop
!
!
! -----------------------------------
!  Output data files
!
        OPEN(6,  FILE='OUTPUT/OUTScreen.txt',    STATUS='UNKNOWN')
        OPEN(9,  FILE='OUTPUT/OUTPUT1.txt',    STATUS='UNKNOWN')
	  OPEN(10, FILE='OUTPUT/OUTPUT.txt' ,    STATUS='UNKNOWN')
	  OPEN(101,FILE='OUTPUT/OUTAMT.txt' ,    STATUS='UNKNOWN')
	  
!	  OPEN(1022,FILE='OUTPUT\OUTBMT.txt' ,    STATUS='UNKNOWN')
	  
	  OPEN(13,FILE='OUTPUT/OUTELA.txt' ,    STATUS='UNKNOWN')
	  OPEN(17,FILE='OUTPUT/OUTPUT7.txt' ,    STATUS='UNKNOWN')  
!~~
          OPEN(14,FILE='OUTPUT/ETI.txt' ,    STATUS='UNKNOWN')
	  open(15,FILE='OUTPUT/INETI.txt' ,    STATUS='UNKNOWN')

!        OPEN(1041,FILE='OUTPUT\ETI2.txt' ,    STATUS='UNKNOWN')
!	  open(1051,FILE='OUTPUT\INETI2.txt' ,    STATUS='UNKNOWN')

C                                                                              
C	   
C
C ----------------------------------------------- 
C
C
!       MFREE=1
!	 NBETA=0
	 
        WRITE(*, *) ' Test on huper-singular integration'

C
C ---------------------------------------- 
C
C  IFWKO=0 : input wave number, otherwise input wave frequency
C  H<0: infinite water depth
C  Nwave: Number of waves to simulate
C
        READ(1,*)      IFWKO

        IF (IFWKO .EQ. 0)  THEN
          READ(1,*)      H, AMP, WK, BETA
          IF (H .LE. 0.0D0) THEN
            W1=DSQRT(G*WK)
          ELSE
            W1=DSQRT(G*WK*DTANH(WK*H))
          END IF
        ELSE
          READ(1,*)      H, AMP, W1, BETA
          IF(H .LE. 0.0D0) THEN
            WK=W1**2/G
          ELSE
!--------------------To be updated------------------------------------
!            CALL WAVECK(W1,H,WK)				 ! Compute wave number
          END IF
        END IF
!
	  READ(1,*)  FCD_AMR,  FCD_AMI
	  READ(1,*)  IFLAG_T,WAVES,  NTnum
!	  READ(1,*)  NPLOUT, NDRFA    
!
!  Waves: number of waves to simulate (may be 10 or 0.5)
!  NTnum: number of time steps in one wave
!  NPLOUT=1-5, Output the wave profile with step intervals 
!  NDRFA=0, plot force; =1, plot response
!

        READ(4,*)     Node_Singular


	  PI4=4.0*PI
        IORDER=1
   
        WRITE(*,*)  '  After 10'
!  

!  PI4£¬IORDER µÄ¶šÒåÔÚmvarÖÐ perturbation order of the problem 
! 
        W1=DSQRT(G*WK*DTANH(WK*H)) 
        TPER=2.*PI/W1 
        BETA=BETA*PI/180.0D0 
C C  beta»¯³É»¡¶È
        V=W1*W1/G 
	  WL=2.0D0*PI/WK
C 
	  TStep=TPer/NTnum
	  NTIME=INT(WAVES*TPER/TStep)
!  NTIME ×ÜÊ±Œä²œÊý
        WRITE(10,*) 
        WRITE(9,*)
	  WRITE(10,*) '                   ================='
        WRITE(9,*) '                   ================='
C
        WRITE(10,1111)  H,AMP,WK,V,WL,W1,TPER,BETA*180./PI 
        WRITE(9,1111)  H,AMP,WK,V,WL,W1,TPER,BETA*180./PI 

!

!  ------------------------------------
!bodmassÔÚmass.fÖÐ
!--------------------------------------------------------------------------------------
!        CALL BODMASS 				            ! Read in data of body mass, etc.

        WRITE(10,*),'  After BODMASS' 

!------------------moved to read mesh by S.Gao March 17------------------------
!    mvarÖÐISYS: number of symmetric planes
!    NELEMB: number of elements on body surface
!    NNB: number of nodes on the body surface according to coordinate
!    NNBD: number of nodes on the body surface according to directives 
        READ(2,*) ISYS !number of symmetric planes
        READ(2,*)   NELEMB, NNB, NNBD, IPOL

        NNODE=NNB
        NNODED=NNBD
        NELEM=NELEMB
!---------------------------to here-------------------------------------------


!    µ¥ÔªÊý£¬œÚµãÊý£¬œÚµãµŒÊý£¬×ø±êžöÊý
!  ? 8ÊÇËÄ±ßÐÎµ¥Ôª£¬3ÊÇÈýžö·œÏò?
   	 ALLOCATE (NCONB(NELEMB,8),NCONDB(NELEMB,8))
       ALLOCATE (XYZ(3,NNB),DXYZ(6,NNBD))

!	  INTEGER NELEM,NELEMB,NELEMF,NNODE,NNODED,NNB,NNF,NNTCH

        IF(ISYS.EQ.0) NSYS=1
        IF(ISYS.EQ.1) NSYS=2
        IF(ISYS.EQ.2) NSYS=4
C
!    mvarÖÐNELEMF:number of elements on the free surface
         nelemf=0
!        READ(3,*)   NELEMF
!
	  WRITE(10,*) ' ISYS=',ISYS,' NSYS=',NSYS
	  WRITE(10,*) ' NELEMB=',NELEMB,' NELEMF=',NELEMF
!    mvarÖÐNELEM: number of total elements
	  NELEM=NELEMB+NELEMF
!
	  WRITE(10,*) ' NELEM=',NELEM,'  IOPL=',IPOL



        ALLOCATE(SAMB(NELEM,16,0:8),SAMBXY(NELEM,16,3),
	1		   DSAMB(NELEM,16,6),NCN(NELEM),NCON(NELEM,8),
     1		   NCOND(NELEM,8),IETYPE(NELEM),NNORMN(8*NELEM) )
        ALLOCATE( XYZE(3,8,NELEM),DXYZE(3,8,NELEM),DAMPE(8,NELEM))
        ALLOCATE( TXYZE(3,8,NELEM))
	  ALLOCATE( XYZTP(3,8*NELEM),DXYZTP(3,8*NELEM),DAMPTP(8*NELEM))
! mvarÖÐµÄ¶šÒå NCN: number of nodes in the element
! IETYPE: type of the element; =1, on body surface; =2, on free surface
! SAMBXY: Coordinates of Gaussin points
! DSAMB:  Normal direvatives at Gaussian points
! XYZE  : Initial Coordinates of nodes of body mesh
! TXYZE : Coordinates of nodes of body mesh at the simulation time 
! --------------------------------------------
! MESHFS4ºÍMESHBDÔÚmeshda4.fÖÐ
!        CALL MESHFS4   		        ! Read in data on free surface mesh
        
        WRITE(10,*),'  After MESHFS4' 

        CALL MESHBD(IPOL) 		    ! Read in data on body mesh,IPOL 
        WRITE(10,*),'  After MESHBD4' 
!
	  CLOSE(2)
       OPEN(50, FILE='OUTPUT/DATBDMS.txt',    STATUS='UNKNOWN')

!	  CALL CONVSB
!        WRITE(10,*),'  After CONVSB' 
        
! mvrÖÐ NNODE:  total number of nodes according to the coordinate
! NNODED: total number of nodes according to the normals
	  WRITE(10,*) ' NNODE=',NNODE
	  WRITE(10,*) ' NNODED=',NNODED

	  NDMAX=MAX(NNODE,NNODED)
       ALLOCATE(NODELE(NNODE,264),NODNOE(NNODE),NODELJ(NNODE,264),
     1         NODQUA(NNODE),NNORMC(NNODED),
     2		   ANGLE(NNODE),DAMPF(NNF))
!
! ---------------------------------------------------
!
!
!  ---------------------------------
!
!        CALL CABLE 		! Read in data of cable curve of force-displacement
!	                    ! Expanding by B-spline function
!	  PRINT *,' AFTER CABLE'
!
!	  CALL FENDER		! Read in data of Fender curve of force-displacement
                          ! Expanding by B-spline function
!
!	  CLOSE(11)
!	  CLOSE(12)
!
!  --------------------------------------------

       ALLOCATE(AMATA(NNODE,NNODE,NSYS),CMATA(NNODE,NNODED,NSYS),
	1		  BMATA(NNODE,NSYS), INDX(NNODE,NSYS))
!       write(*,*)  sizeof(amata),sizeof(cmata)
!       pause
       
!	 ALLOCATE(
!	1		   INDX(NNODE,NSYS))
       ALLOCATE(UNKN(NNODE,NSYS),  BKN(NNODED,NSYS),
	1		  UNKN_O(NNODE,NSYS),BKN_O(NNODED,NSYS),
	1	      ET(NNF,NSYS),    ET_O(NNF,NSYS),DPDT(NNODE,NSYS))
!	 ALLOCATE(HEIGHT(4,NNF,NSYS),PFREEN(4,NNF,NSYS))
	 ALLOCATE(DH(4,NNF,NSYS),DP(4,NNF,NSYS),Dposi(4,6))
!	 ALLOCATE(LEFT(NNODE,NNODE,NSYS),RIGHT(NNODE,NNODED,NSYS),
!     1         RIGHT1(NNODEd,NSYS) )
!
!      
! Identify the Gaussian sampling points and evaluate the    
! Corresponding values for the shape functions and Jacobian matrices       
!
!        CALL BODYFD                  
        WRITE(10,*)  'AFTER BODYFD' 
!  ---------------------------------
!
!   
! Assembling matrix and computing diffraction and radiation potentials
!
        CALL TASSB0   
	  WRITE(10,*) 'AFTER ASSEMB0' 

! 
	  ITIME=0
        TIME=0.0d0
!	  
!	  ET(:,:)  =0.0
!	  BKN(:,:) =0.0
!	  UNKN(:,:)=0.0
!	  UNKN_O(:,:)=0.0
!	  FORCE(:) =0.0
!	  DISP(:)  =0.0
!	  DSDT(:)  =0.0
!
!	  CALL PLOTOUT8
!
!  ==================================================

        
 
!       DEALLOCATE(INDX)
	 DEALLOCATE(AMATA,CMATA,BMATA,INDX)
!       DEALLOCATE(LEFT,RIGHT,RIGHT1)
C                              
1010    FORMAT(F7.3,1x,F7.3,1x,6E14.5) 
C 
1111    FORMAT(//,'  WATER DEPTH=',F9.3,'    WAVE AMPLITUDE=', F6.2,/,
     1    '  WAVE NUMBER=',F9.5,'  K0=',F9.5,'  WAVE LENGTH=',F9.4,/, 
     3    '  ANGULAR FREQU.=',F9.5,'   WAVE PERIOD=',F7.3,/,      
     2    '  WAVE DIRECTION:',F7.3,'  Degree',/)
C 
1113    FORMAT(/,15x,' FIRST  ORDER PROBLEM')  
1114    FORMAT(/,15x,' SECOND ORDER PROBLEM') 
1115	  FORMAT(' I_time=',I5,'    at the time:',F10.5,' s')
1200	FORMAT(2x,I3,3F12.5,1x,2F13.5)
C 


        END    
