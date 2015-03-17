! 
! ********************************************* 
! * SOLUTION OF LINEAR EQUATIONS  [A]{X}= {B} * 
! *           LU DECOMPOSITION                * 
! * FROM 'NUMERICAL RECIPES'     pp. 35-37    * 
! ********************************************* 
!
!
! ---------------------------------------- 
! NSYS: \BE\D8\D5\F3A[:,:]\B5ĸ\F6\CA\FD\A3\AC\D3\C3\D3ڿ\AA\CA\FD\D7\E9 
! IP  : \B1\BE\B4μ\C6\CB\E3\B5ľ\D8\D5\F3\D0\F2\BA\C5
! A   :
! N   :
! NP  :
! LI  :
! NMOD:
! INDX:
! B   :
! ---------------------------------------- 
! 

        SUBROUTINE RLUDCMP(IP,A,N,NP,NSYS,INDX,D)           
        IMPLICIT REAL*8(A-H,O-Z)  
        PARAMETER (NMAX=5000, TINY=1.0E-20) 
        DIMENSION INDX(NP,NSYS) 
        REAL*8 SUM,DUM,AAMAX 
        REAL*8 A(NP,NP,NSYS),VV(NMAX) 
C 
        D=1. 
        DO 12 I=1, N 
        AAMAX=(0., 0.) 
          DO 11 J=1, N 
          IF ( DABS(A(I,J,IP)).GT. DABS(AAMAX) )  
     1                        AAMAX=A(I,J,IP) 
11      CONTINUE 
!
        IF (DABS(AAMAX) .EQ. 0.0)  THEN	  
	    Print  *, ' SINGULAR MATRIX   inside RLUDCMP' 
	    Print  *, ' IP=',IP,' I=',I
		PAUSE
        ENDIF
!
	  VV(I)=1./AAMAX 
12      CONTINUE 
        DO 19 J=1, N 
          DO 14 I=1, J-1 
            SUM=A(I,J,IP) 
            DO 13 K=1, I-1 
              SUM=SUM-A(I,K,IP)*A(K,J,IP) 
13            CONTINUE 
            A(I,J,IP)=SUM 
14          CONTINUE 
          AAMAX=(0., 0.) 
          DO 16 I=J, N                          
            SUM=A(I,J,IP) 
            DO 15 K=1, J-1 
              SUM=SUM-A(I,K,IP)*A(K,J,IP) 
15            CONTINUE 
            A(I,J,IP)=SUM 
            DUM=VV(I)*SUM 
            IF (DABS(DUM) .GE. DABS(AAMAX)) THEN 
            IMAX=I 
            AAMAX=DUM 
            END IF 
16          CONTINUE 
         IF(J .NE. IMAX) THEN 
         DO 17 K=1, N 
           DUM=A(IMAX,K,IP) 
           A(IMAX,K,IP)=A(J,K,IP) 
           A(J,K,IP)=DUM 
17         CONTINUE 
          D=-D 
          VV(IMAX)=VV(J) 
          END IF 
          INDX(J,IP)=IMAX 
          IF(A(J,J,IP).EQ.0.) A(J,J,IP)=TINY 
          IF(J.NE.N) THEN 
          DUM=1./A(J,J,IP) 
          DO 18 I=J+1,N 
            A(I,J,IP)=A(I,J,IP)*DUM 
18          CONTINUE 
          END IF 
19        CONTINUE 
        RETURN 
	 END SUBROUTINE RLUDCMP





!    
! ---------------------------------------- 
! A   :
! N   :
! NP  :
! IP  :
! NSYS:
!
! B   :
! LI  :
! NMOD:
!
! INDX:
! ----------------------------------------
!
        SUBROUTINE RLUBKSB(IP,A,N,NP,LI,NSYS,NMOD,INDX,B)  
        IMPLICIT REAL*8(A-H,O-Z)  
        DIMENSION INDX(NP,NSYS) 
        REAL*8 SUM,A(NP,NP,NSYS),B(NP,NMOD,NSYS) 
C                                     
        II=0 
        DO 12 I=1,N 
         LL=INDX(I,IP) 
        SUM=B(LL,LI,IP) 
        B(LL,LI,IP)=B(I,LI,IP) 
        IF(II.NE.0) THEN 
        DO 11 J=II, I-1 
        SUM=SUM-A(I,J,IP)*B(J,LI,IP) 
11      CONTINUE 
        ELSE IF (DABS(SUM).NE. 0.) THEN 
        II=I 
        END IF 
        B(I,LI,IP)=SUM 
12      CONTINUE 
        DO 14 I=N, 1, -1 
        SUM=B(I,LI,IP) 
        DO 13 J=I+1, N 
        SUM=SUM-A(I,J,IP)*B(J,LI,IP) 
13      CONTINUE 
        B(I,LI,IP)=SUM/A(I,I,IP) 
14      CONTINUE 
        RETURN 
	  END SUBROUTINE RLUBKSB

!
!
!
C 
!
C 
C ********************************************* 
C * Inverse of a matric  [A] by               * 
C *           LU decomposition                * 
C * From 'NUMERICAL RECIPES'     pp. 35-37    * 
C ********************************************* 
C 
       SUBROUTINE INVERSE(a,y,n,np)
       INTEGER n,np,indx(n)
       REAL*8 a(np,np),b(np),y(np,np),d

	 CALL ludcmp(a,n,np,indx,d)
	  y=0.0d0
       do 10 i=1, n
	  y(i,i)=1.0d0
10	 continue
c
	 do 100 j=1, n
	  do 20 i=1, n
	   b(i)=y(j,i) !i,j
20	  continue
c
        call  lubksb(a,n,np,indx,b)
c
	  do 30 i=1, n
	   y(j,i)=b(i) !j,i
30	  continue

100	 continue


      END	 SUBROUTINE INVERSE

!
C 
C ********************************************* 
C * SOLUTION OF LINEAR EQUATIONS  [A]{X}= {B} * 
C *           LU DECOMPOSITION                * 
C * FROM 'NUMERICAL RECIPES'     pp. 35-37    * 
C ********************************************* 
C 
      SUBROUTINE ludcmp(a,n,np,indx,d)

      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (DAbs(a(i,j)).gt.aamax) aamax=DAbs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*DAbs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax) then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n) then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END	 SUBROUTINE ludcmp




      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0) then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12     continue
!       
	 do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14     continue
       return
       END	 SUBROUTINE lubksb

