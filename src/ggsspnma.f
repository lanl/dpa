c$$$c%---------------------------------------------------------------%
c                                                                |
C     Perturbed Normal Mode Analysis                             |
c        using Elastic Network Potential                         |
c        for Calclating Dynamics of                              |
c        of PROTEIN/LIGAND Complex                               |
c        with Alternative Binding Sites                          |
c                                                                |
c       By  Dengming Ming                                        |
c       Los Alamos National Lab                                  |
c       July 2004                                                |
c       Jan  2005  Revised by Ming                               |
c       May  2005  add Backbone Enhance by Ming                  | 
c       Feb  2006  developed Perturbation Calculation by Ming    |
c                                                                | 
c%---------------------------------------------------------------%
	Program PerturbedNormalModeAnalysiswithElasticNetWorkMethod
c
c
c       wcc: scaling interactions between CA -- CA
c       wlg: scaling interactions between CA -- Ligand
c       wll: scaling interactions between Ligand -- Ligand
c       cutcc: cutoff for interactions between CA -- CA
c       cutlg: cutoff for interactions between CA -- Ligand
c       cutll: cutoff for interactions between Ligand -- Ligand
c
c	write(*,*)
c	write(*,*) 'Program for system with less than 500 CA'
c	write(*,*) 'This restriction can be removed by changing: Nat0,DD'
c	write(*,*)
c
	implicit none
	character*100  file_cacrd,file_surfp,file_dpa,infocalc
	integer   CALC_INFO
	integer   nca,nsurf
	real*8 wcc,wctc,wlg,wll,cutcc,cutlg,cutll
	real*8 ev_threshold	!a value to estimate the calculation Success or Failure
	real   time1,time2
c
	read(*,*) file_cacrd,file_surfp,file_dpa
	read(*,*) wcc,wctc,wlg,wll
	read(*,*) cutcc,cutlg,cutll	
	read(*,*) ev_threshold
c       
	call initial(file_cacrd,nca)
	call initial(file_surfp,nsurf)
c
	call gene_surfdpa(calc_info,infocalc,ev_threshold,file_cacrd,
     &       file_surfp,file_dpa,nca,nsurf,cutcc,cutlg,cutll,wcc,
     &       wctc,wlg,wll,time1,time2)

	if(calc_info.eq.0) then
	   write(*,'(1x,I5,3x,F10.3)') CALC_INFO,time2-time1
	else
	   write(*,'(1x,I5,3x,A)') CALC_INFO, infocalc
	endif
	stop
c       
	end
 	subroutine gene_surfdpa(calc_info,infocalc,ev_threshold,
     &       file_cacrd,file_surfp,file_dpa,nca,nsurf,cutcc,cutlg,
     &       cutll,wcc,wctc,wlg,wll,time1,time2)
	implicit none
	character*100 infocalc,file_cacrd,file_surfp,file_dpa
	integer   calc_info,nca,nsurf,nca3,ndimca,ndimhet,i,k
	INTEGER   nat0,nat30,nhet,nhet3,one
	PARAMETER (nat0=1200,nat30=3*nat0,nhet=1,nhet3=3*nhet,one=1)
	real*8    cacrd(3*nca),surfcrd(3*nsurf),hetcrd(nhet3)
	real*8    wcc,wctc,wlg,wll,cutcc,cutlg,cutll,pert
	real*8  ev_threshold,ev0(nat30),vec0(nat30*nat30)
	common  /MODE0/ ev0,vec0
	real      time1,time2
c
	nca3=3*nca
	ndimca=(nca3*(nca3+1))/2
	ndimhet=(nhet3*(nhet3+1))/2
c       
	call SECOND(time1)
	open(unit=20,file=file_cacrd,STATUS='OLD',Err=5000) !Protein Backbone Coordinates
	do i=1,nca
	   read(20,*,Err=5010) cacrd(3*i-2),cacrd(3*i-1),cacrd(3*i)
	enddo
	close(20)
	open(unit=20,file=file_surfp,STATUS='OLD',Err=5500) !Protein Backbone Coordinates
	do i=1,nsurf
	   read(20,*,Err=5510) surfcrd(3*i-2),surfcrd(3*i-1),surfcrd(3*i)
	enddo
	close(20)
c
	call hessian_apoprotein(calc_info,infocalc,cacrd,nca,nca3,
     &       ndimca,cutcc,wcc,wctc,ev_threshold)
	if(calc_info.ne.0) return
c
	open(20,file=file_dpa,status='unknown')
	do i=1,nsurf
	   hetcrd(1)=surfcrd(3*i-2)
	   hetcrd(2)=surfcrd(3*i-1)
	   hetcrd(3)=surfcrd(3*i)
	   call DPA_1th_perturb(calc_info,infocalc,cacrd,hetcrd,
     &          nca+nhet,nca3,nhet3,ndimca,ndimhet,cutlg,cutll,wlg,
     &          wll,pert)
	   if(calc_info.ne.0) return
	   write(20,'(1x,3F9.3,2x,F12.6,I10)') (hetcrd(k),k=1,3),pert,i
	enddo
	close(20)
c
	call SECOND(time2)
c
	calc_info=0
	return
c
 5000	calc_info=-9999
	infocalc='open cacrd file Error'
	return
 5010	calc_info=-9998
	infocalc='read cacrd Error'
	return
 5500	calc_info=-9997
	infocalc='open surfcrd file Error'
	return
 5510	calc_info=-9996
	infocalc='read surfcrd Error'
	return
c
	return
	end

c
c------------------------------------------------------------------------------
c* Desk :Find the eignvalues and eignvactors: i.e. Modes:  
c !Modes start from NDD
c
      subroutine CharmmDiagq(N,NDD,mode_need,ev,vec)
      implicit none
      INTEGER   N,mode_need,NDD,NADD,k !N=3*npart
      REAL*8 VEC(*),EV(*)
      REAL*8 A(N+1),B(N+1),P(N+1),W(N+1),TA(N+1),TB(N+1),Y(N+1)
c
      NADD=NDD-1                !NDD is the first mode to be found.
      call DIAGQ(N,mode_need,VEC,A,B,P,W,EV,TA,TB,Y,NADD)
c
      return
      end
c
c-----------------------------------------------------------------------
c
c
      SUBROUTINE DIAGQ(NX,NFRQX,VEC,A,B,P,W,EV,TA,TB,Y,NADD)
C
C     THIS ROUTINE IS A CONGLOMERATION OF GIVEN, HOUSEC, AND EIGEN
C     WHERE THE BEST FEATURES OF EACH WERE KEPT AND SEVERAL OTHER
C     MODIFICATIONS HAVE BEEN MADE TO INCREASE EFFICIENCY AND ACCURACY.
C
C   By Bernard R. Brooks   1981
C
C   NX      - ORDER OF MATRIX
C   NFRQX   - NUMBER OF ROOTS DESIRED
C   DD      - SECOND DERIVATIVE MATRIX IN UPPER TRIANGULAR FORM
C   VEC     - EIGENVECTORS RETURNED (NX,NFRQX)
C   EV      - EIGENVALUES RETURNED (NX)
C   A,B,P,W,TA,TB,Y - ALL SCRATCH VECTORS (NX+1)
C   NADD    - NUMBER OF LOWEST ROOTS TO SKIP
C
C
      implicit none
! This is from  include 'number.fcm'
CHARMM Element source/fcm/number.fcm 1.1
C
C This file contains floating point numbers.
C
C positive numbers
      REAL*8     ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX,
     &           SEVEN, EIGHT, NINE, TEN, ELEVEN, TWELVE, THIRTN,
     &           FIFTN, NINETN, TWENTY, THIRTY
      PARAMETER (ZERO   =  0.D0, ONE    =  1.D0, TWO    =  2.D0,
     &           THREE  =  3.D0, FOUR   =  4.D0, FIVE   =  5.D0,
     &           SIX    =  6.D0, SEVEN  =  7.D0, EIGHT  =  8.D0,
     &           NINE   =  9.D0, TEN    = 10.D0, ELEVEN = 11.D0,
     &           TWELVE = 12.D0, THIRTN = 13.D0, FIFTN  = 15.D0,
     &           NINETN = 19.D0, TWENTY = 20.D0, THIRTY = 30.D0)
C
      REAL*8     FIFTY, SIXTY, SVNTY2, EIGHTY, NINETY, HUNDRD,
     &           ONE2TY, ONE8TY, THRHUN, THR6TY, NINE99, FIFHUN, THOSND,
     &           FTHSND,MEGA
      PARAMETER (FIFTY  = 50.D0,  SIXTY  =  60.D0,  SVNTY2 =   72.D0,
     &           EIGHTY = 80.D0,  NINETY =  90.D0,  HUNDRD =  100.D0,
     &           ONE2TY = 120.D0, ONE8TY = 180.D0,  THRHUN =  300.D0,
     &           THR6TY=360.D0,   NINE99 = 999.D0,  FIFHUN = 1500.D0,
     &           THOSND = 1000.D0,FTHSND = 5000.D0, MEGA   =   1.0D6)
C
C negative numbers
      REAL*8     MINONE, MINTWO, MINSIX
      PARAMETER (MINONE = -1.D0,  MINTWO = -2.D0,  MINSIX = -6.D0)
C
C common fractions
      REAL*8     TENM20,TENM14,TENM8,TENM5,PT0001,PT0005,PT001,PT005,
     &           PT01, PT02, PT05, PTONE, PT125, PT25, SIXTH, THIRD,
     &           PTFOUR, PTSIX, HALF, PT75, PT9999, ONEPT5, TWOPT4
      PARAMETER (TENM20 = 1.0D-20,  TENM14 = 1.0D-14,  TENM8  = 1.0D-8,
     &           TENM5  = 1.0D-5,   PT0001 = 1.0D-4, PT0005 = 5.0D-4,
     &           PT001  = 1.0D-3,   PT005  = 5.0D-3, PT01   = 0.01D0,
     &           PT02   = 0.02D0,   PT05   = 0.05D0, PTONE  = 0.1D0,
     &           PT125  = 0.125D0,  SIXTH  = ONE/SIX,PT25   = 0.25D0,
     &           THIRD  = ONE/THREE,PTFOUR = 0.4D0,  HALF   = 0.5D0,
     &           PTSIX  = 0.6D0,    PT75   = 0.75D0, PT9999 = 0.9999D0,
     &           ONEPT5 = 1.5D0,    TWOPT4 = 2.4D0)
C
C others
      REAL*8 ANUM,FMARK
      REAL*8 RSMALL,RBIG
      PARAMETER (ANUM=9999.0D0, FMARK=-999.0D0)
      PARAMETER (RSMALL=1.0D-10,RBIG=1.0D20)
C
C Machine constants (these are very machine dependent).
C
C RPRECI should be the smallest number you can add to 1.0 and get a number
C that is different than 1.0.  Actually: the following code must pass for
C for all real numbers A (where no overflow or underflow conditions exist).
C 
C         B = A * RPRECI
C         C = A + B
C         IF(C.EQ.A) STOP 'precision variable is too small'
C 
C The RBIGST value should be the smaller of:
C 
C         - The largest real value               
C         - The reciprocal of the smallest real value 
C
C If there is doubt, be conservative.
C 
C!!!!! NOTE:  Some of the values have not been checked....
C!!!!! Please fix these values - BRB .......
C
      REAL*8 RPRECI,RBIGST
C
      PARAMETER (RPRECI = 2.22045D-16, RBIGST = 4.49423D+307)
C
c End of 'number.fcm'
c
c This is from:  include 'stream.fcm'
CHARMM Element source/fcm/stream.fcm 1.1
C
C     This is the STREAM data block.
C     It contains information abount the current runstream.
C
C     MXSTRM - Maximum number of active stream files.
C     POUTU  - Default output unit number.
C     NSTRM  - Number of active input streams.
C     ISTRM = JSTRM(NSTRM) - Current input unit number
C     JSTRM(*) - stack of input streams numbers.
C     OUTU   - Unit number for all standard CHARMM output.
C     PRNLEV - Print level control for all writing to OUTU
C     IOLEV  - -1 to 1  -1=write no files.  1= write all files.
C     WRNLEV - -5 TO 10  0=SEVERE ONLY, 10=LIST ALL WARNINGS
C     LOWER  - if .true. all files with names not in double quotes
C              will be opened in lower case for write. For read
C              UPPER case will be tried first and if not succesful
C              file name will be converted to lower case.
C     QLONGL - Use long lines in the output where appropriate.
C              (Otherwise, restrict output lines to 80 characters)
C
      LOGICAL LOWER,QLONGL
      INTEGER MXSTRM,POUTU
      PARAMETER (MXSTRM=20,POUTU=6)
      INTEGER   NSTRM,ISTRM,JSTRM,OUTU,PRNLEV,WRNLEV,IOLEV
C
      COMMON /CASE/   LOWER, QLONGL
      COMMON /STREAM/ NSTRM,ISTRM,JSTRM(MXSTRM),OUTU,PRNLEV,WRNLEV,IOLEV
C
C End of 'stream.fcm'
c
      INTEGER   NX,NFRQX,NADD
      REAL*8 VEC(*)
      REAL*8 A(NX),B(NX),P(NX),W(NX),EV(NX),TA(NX),TB(NX),Y(NX)

      REAL*8 ETA,THETA,DEL1,DELTA,SMALL,DELBIG,THETA1,TOLER,ETAR
      REAL*8 RPOWER,RPOW1,RAND1,FACTOR,ANORM,U,ANORMR
      REAL*8 SUM1,BX,S,SGN,TEMP,XKAP,EXPR,ALIMIT,ROOTL,ROOTX,TRIAL,F0
      REAL*8 AROOT,ELIM1,ELIM2,T,EPR,XNORM,XNORM1,EVDIFF
      INTEGER   N,NEV,NEVADD,NTOT,I,IPT,J,IJ,NN,MI,MI1,JI,JI2,II
      INTEGER   ML,ML1,L,M,K,MJ,MJ1,NOMTCH,NOM,IA,ITER
      INTEGER   J1,MK,MK1,KK
C
CCC
      real*8 anumx
      INTEGER   NAT0,N_DIM0
      PARAMETER (NAT0=1200, N_DIM0=(9*NAT0*NAT0+3*NAT0)/2)
      Real*8  DD(N_DIM0)
      common /DIAG/ DD
CCC

!	do i=1,nx
!	   j1=nx*(i-1)-(i-1)*(i-2)/2
!	   write(*,'(A, 20F8.3)') 'DD>',(DD (j1+j-i+1),j=i,nx)
!	enddo
      anumx=zero
      do i=1,nx
         A(I)=anumx
         B(I)=anumx
         P(I)=anumx 
         W(I)=anumx
         EV(I)=anumx
         TA(I)=anumx
         TB(I)=anumx
         Y(I)=anumx
      enddo
CCC
C
      ETA=RPRECI
      THETA=RBIGST
C
      N=NX
      NEV=NFRQX
      NEVADD=NEV+NADD
C
      DEL1=ETA/100.0
      DELTA=ETA**2*100.0
      SMALL=ETA**2/100.0
      DELBIG=THETA*DELTA/1000.0
      THETA1=1000.0/THETA
      TOLER=100.0*ETA
      ETAR=1.0/ETA
      RPOWER=8388608.0
      RPOW1=RPOWER*0.50
      RAND1=RPOWER-3.0
C
C Find largest element.
      FACTOR=ZERO
      NTOT=(N*(N+1))/2
      DO I=1,NTOT
         FACTOR=MAX(FACTOR,ABS(DD(I)))
      ENDDO
C
C Check for zero matrix.
      IF(FACTOR.LE.THETA1) THEN
         IF(WRNLEV.GE.2) WRITE(OUTU,811)
 811     FORMAT(' WARNING FROM <DIAGQ>. Zero matrix passed.',
     1     ' Identity matrix returned.')
         DO I=1,NEV
            EV(I)=ZERO
            IPT=(I-1)*N
            DO J=1,N
               IPT=IPT+1
               VEC(IPT)=ZERO
               IF(I+NADD.EQ.J) VEC(IPT)=ONE
            ENDDO
         ENDDO
         RETURN
      ENDIF
C
C Compute norm of matrix
      FACTOR=ONE/FACTOR
      IJ=0
      ANORM=ZERO
      DO I=1,N
         DO J=I,N
            IJ=IJ+1
            U=(DD(IJ)*FACTOR)**2
            IF(I.EQ.J) U=U*HALF
            ANORM=ANORM+U
         ENDDO
      ENDDO
C
C Scale the matrix
      ANORM=SQRT(ANORM+ANORM)/FACTOR
      ANORMR=ONE/ANORM
      DO I=1,NTOT
         DD(I)=DD(I)*ANORMR
      ENDDO
C
      NN=N-1
      MI=0
      MI1=N-1
C
C Perform trigiagonalization
      DO I=1,NN
         SUM1=ZERO
         B(I)=ZERO
         JI=I+1
         IPT=MI+I
         A(I)=DD(IPT)
         IPT=IPT+1
         BX=DD(IPT)
         JI2=JI+1
         DO J=JI2,N
            IPT=IPT+1
            SUM1=SUM1+DD(IPT)*DD(IPT)
         ENDDO
         IF(SUM1.LT.SMALL) THEN
            B(I)=BX
            DD(MI+JI)=ZERO
         ELSE
            S=SQRT(SUM1+BX**2)
            SGN=SIGN(ONE,BX)
            TEMP=ABS(BX)
            W(JI)=SQRT(HALF*(ONE+(TEMP/S)))
            IPT=MI+JI
            DD(IPT)=W(JI)
            II=I+2
            IF(II.LE.N) THEN
               TEMP=SGN/(TWO*W(JI)*S)
               DO J=II,N
                  IPT=IPT+1
                  W(J)=TEMP * DD(IPT)
                  DD(IPT)=W(J)
               ENDDO
            ENDIF
            B(I)=-SGN*S
C
            DO J=JI,N
               P(J)=ZERO
            ENDDO
            ML=MI + MI1
            ML1=MI1-1
            DO L=JI,N
               IPT=ML+L
               DO M=L,N
                  BX=DD(IPT)
                  P(L)=P(L)+BX*W(M)
                  IF(L.NE.M) P(M)=P(M)+BX*W(L)
                  IPT=IPT+1
               ENDDO
               ML=ML +ML1
               ML1=ML1-1
            ENDDO
C
C
            XKAP=ZERO
            DO K=JI,N
               XKAP=XKAP+W(K)*P(K)
            ENDDO
            DO L=JI,N
               P(L)=P(L)-XKAP*W(L)
            ENDDO
            MJ=MI+MI1
            MJ1=MI1-1
            DO J=JI,N
               DO K=J,N
                  EXPR=(P(J)*W(K))+(P(K)*W(J))
                  DD(MJ+K)=DD(MJ+K)-EXPR-EXPR
               ENDDO
               MJ=MJ+MJ1
               MJ1=MJ1-1
            ENDDO
         ENDIF
         MI=MI+MI1
         MI1=MI1-1
      ENDDO
C
C Begin sturm bisection method.
C
      A(N)=DD(MI+N)
      B(N)=ZERO
C
      ALIMIT=ONE
      DO I=1,N
         W(I)=B(I)
         B(I)=B(I)*B(I)
      ENDDO
      DO I=1,NEVADD
         EV(I)=ALIMIT
      ENDDO
      ROOTL=-ALIMIT
C
      DO I=1,NEVADD
         ROOTX=ALIMIT
         DO J=I,NEVADD
            ROOTX=MIN(ROOTX,EV(J))
         ENDDO
         EV(I)=ROOTX
C
 130     CONTINUE
            TRIAL=(ROOTL+EV(I))*HALF
c##IF CRAY (oldcode_for_cray)
c            IF(TRIAL.EQ.ROOTL.OR.TRIAL.EQ.EV(I)) GOTO 200
c##ELSE (oldcode_for_cray)
            EVDIFF=ABS(ROOTL-EV(I))
            IF(EVDIFF.LT.THETA1) GOTO 200
            IF(EVDIFF*ETAR.LT.ABS(TRIAL)) GOTO 200
c##ENDIF (oldcode_for_cray)
            NOMTCH=N
            J=1
 150     CONTINUE
            F0=A(J)-TRIAL
 160     CONTINUE
            IF(ABS(F0).LT.THETA1) GOTO 170
            IF(F0.GE.ZERO) NOMTCH=NOMTCH-1
            J=J+1
            IF(J.GT.N) GOTO 180
            F0=A(J)-TRIAL-B(J-1)/F0
            GOTO160
 170     CONTINUE
            J=J+2
            NOMTCH=NOMTCH-1
            IF(J.LE.N) GOTO 150
 180     CONTINUE
            IF(NOMTCH.GE.I) GOTO 190
            ROOTL=TRIAL
            GOTO 130
 190     CONTINUE
            EV(I)=TRIAL
            NOM=MIN(NEVADD,NOMTCH)
C            NOM=MIN0(NEVADD,NOMTCH)   		!original 
            EV(NOM)=TRIAL
            GOTO 130
 200     CONTINUE
      ENDDO
C
C Finished computing requested eigenvalues
      DO I=1,NEV
         EV(I)=EV(I+NADD)
      ENDDO
C
C Compute eigenvectors (backtransformation)
      DO I=1,NEV
         AROOT=EV(I)
         DO J=1,N
            Y(J)=ONE
         ENDDO
         IA=IA+1
         IF(I.EQ.1) THEN
            IA=0
         ELSE
            IF(ABS(EV(I-1)-AROOT).GE.TOLER) IA=0
         ENDIF
         ELIM1=A(1)-AROOT
         ELIM2=W(1)
         DO J=1,NN
            IF(ABS(ELIM1).LE.ABS(W(J))) THEN
               TA(J)=W(J)
               TB(J)=A(J+1)-AROOT
               P(J)=W(J+1)
               TEMP=ONE
               IF(ABS(W(J)).GT.THETA1) TEMP=ELIM1/W(J)
               ELIM1=ELIM2-TEMP*TB(J)
               ELIM2=-TEMP*W(J+1)
            ELSE
               TA(J)=ELIM1
               TB(J)=ELIM2
               P(J)=ZERO
               TEMP=W(J)/ELIM1
               ELIM1=A(J+1)-AROOT-TEMP*ELIM2
               ELIM2=W(J+1)
            ENDIF
            B(J)=TEMP
         ENDDO
C
         TA(N)=ELIM1
         TB(N)=ZERO
         P(N)=ZERO
         P(NN)=ZERO
         ITER=1
         IF(IA.NE.0) GOTO 460
C
 320     L=N+1
         DO J=1,N
            L=L-1
 330        CONTINUE
            IF(L.EQ.N) THEN
               ELIM1=Y(L)
            ELSE IF(L.EQ.N-1) THEN
               ELIM1=Y(L)-Y(L+1)*TB(L)
            ELSE
               ELIM1=Y(L)-Y(L+1)*TB(L)-Y(L+2)*P(L)
            ENDIF
C
C Overflow check
            IF(ABS(ELIM1).GT.DELBIG) THEN
               DO K=1,N
                  Y(K)=Y(K)/DELBIG
               ENDDO
               GOTO 330
            ENDIF
            TEMP=TA(L)
            IF(ABS(TEMP).LT.DELTA) TEMP=DELTA
            Y(L)=ELIM1/TEMP
         ENDDO
C
         IF(ITER.EQ.2) GOTO 500
         ITER=ITER+1
C
 420     CONTINUE
         ELIM1=Y(1)
         DO J=1,NN
            IF(TA(J).EQ.W(J)) THEN
               Y(J)=Y(J+1)
               ELIM1=ELIM1-Y(J+1)*B(J)
            ELSE
               Y(J)=ELIM1
               ELIM1=Y(J+1)-ELIM1*B(J)
            ENDIF
         ENDDO
         Y(N)=ELIM1
         GOTO 320
C
 460     CONTINUE
         DO J=1,N
            RAND1=MOD(4099.0*RAND1,RPOWER)
            Y(J)=RAND1/RPOW1-ONE
         ENDDO
         GOTO 320
C
C Orthog to previous
 500     IF(IA.EQ.0) GOTO 550
         DO J1=1,IA
            K=I-J1
            TEMP=ZERO
            IPT=(K-1)*N
            DO J=1,N
               IPT=IPT+1
               TEMP=TEMP+Y(J)*VEC(IPT)
            ENDDO
            IPT=(K-1)*N
            DO J=1,N
               IPT=IPT+1
               Y(J)=Y(J)-TEMP*VEC(IPT)
            ENDDO
          ENDDO
 550      CONTINUE
          IF(ITER.EQ.1) GOTO 420
C
C Normalize
 560     CONTINUE
         ELIM1=ZERO
         DO J=1,N
            ELIM1=MAX(ELIM1,ABS(Y(J)))
         ENDDO
         TEMP=ZERO
         DO J=1,N
            ELIM2=Y(J)/ELIM1
            TEMP=TEMP+ELIM2*ELIM2
         ENDDO
         TEMP=ONE/(SQRT(TEMP)*ELIM1)
         DO J=1,N
            Y(J)=Y(J)*TEMP
            IF(ABS(Y(J)).LT.DEL1) Y(J)=ZERO
         ENDDO
         IPT=(I-1)*N
         DO J=1,N
            IPT=IPT+1
            VEC(IPT)=Y(J)
         ENDDO
      ENDDO
C
      DO I=1,NEV
         IPT=(I-1)*N
         DO J=1,N
            IPT=IPT+1
            Y(J)=VEC(IPT)
         ENDDO
C
         L=N-2
         MK=(N*(N-1))/2-3
         MK1=3
C
         DO J=1,L
            T=ZERO
            K=N-J-1
            M=K+1
            DO KK=M,N
               T=T+DD(MK+KK)*Y(KK)
            ENDDO
            DO KK=M,N
               EPR=T*DD(MK+KK)
               Y(KK)=Y(KK)-EPR-EPR
            ENDDO
            MK=MK-MK1
            MK1=MK1+1
         ENDDO 
C
         T=ZERO
         DO J=1,N
            T=T+Y(J)*Y(J)
         ENDDO
         XNORM=SQRT(T)
         XNORM1=ONE/XNORM
         DO J=1,N
            Y(J)=Y(J)*XNORM1
         ENDDO
C
         IPT=(I-1)*N
         DO J=1,N
            IPT=IPT+1
            VEC(IPT)=Y(J)
         ENDDO
      ENDDO
C
      DO I=1,N
         EV(I)=EV(I)*ANORM
      ENDDO
C
      RETURN
      END
C*Desk :CharmDiagq: End!
	subroutine initial(file_input,n_point)
c       
c       return the number of point in file: file_input
c       
	implicit none
	integer   n_point
	character*100 file_input
c       
	open(unit=20,file=file_input,status='old',err=100)
	n_point=0
 10	read(20,*,err=100,end=200) 
	n_point=n_point+1
	goto 10
 50	write(*,*) -1, '  error occur in OPEN input crd file:',file_input
	stop 
 100	write(*,*) -1, '  error occur in read input crd file:',file_input
	stop 
 200	close(20)
c       
	return
	end
c*Desk :Hessian Matrix :Begin!
	subroutine hessian_apoprotein(calc_info,infocalc,cacrd,nca,
     &       nca3,ndimca,cutcc,wcc,wctc,ev_threshold)
	implicit none
	character infocalc*100
	integer   calc_info,nca,nca3,ndimca,i,j,k,i1,j1,ki0,kj0,kh
	real*8    cacrd(*),x(nca),y(nca),z(nca)
!	real*8 hessian(ndimca)
	real*8 hessian1,hessian2,hessian3
	real*8 hessian4,hessian5,hessian6
	real*8 xx,yy,zz,sij,dsqrt,rr2ij
	real*8 cutoff,cutcc,wcc,wctc,wiwj
	integer   nat0,nat30,ione,ndim0
	parameter (nat0=1200,nat30=3*nat0,ione=1,ndim0=(nat30**2+nat30)/2)
	real*8  ev_threshold,ev0(nat30),vec0(nat30*nat30)
	real*8  hessian(ndim0)
	common  /DIAG/  hessian
	common  /MODE0/ ev0,vec0
c
	if(nca.gt.nat0) then
	   calc_info=-1000
	   infocalc='  The Declared DIMESION of Hessian is too small'
	   return
	endif
	do i=1,ndimca
	   hessian(i)=0.0
	enddo
c   
c--     Reading CA 
	do i=1,nca
	   x(i)=cacrd(3*i-2)
	   y(i)=cacrd(3*i-1)
	   z(i)=cacrd(3*i)
	enddo
c       
c--  GENERATE HESSIAN
c	!Off-diagonalize elements for lower-triangle off-diagonal blocks.
c       i1=4 !(4,1) is the 1st element of the 2nd block of 1st column. 
c       j1=1
c       |1 2 3|
c       |  4 5| main diagonal block
c       |    6|
c       
c
	do i=1,nca		!atom i
	   do j=i+1,nca	                 !atom j
	      i1=3*(i-1)+1
	      ki0=(i1-1)*nca3-((i1-1)*(i1-2))/2
	      j1=3*(j-1)+1
	      xx=x(j)-x(i)
	      yy=y(j)-y(i)
	      zz=z(j)-z(i)
	      sij=sqrt(xx*xx+yy*yy+zz*zz)    !distance between  atom I and atom J
	      rr2ij=1.0/sij/sij
	      if(sij.lt.0.3) then
		 calc_info=-5000
		 infocalc='  Input Protein CRD error: Pair distance < 0.3'
		 return
	      endif
		              	        !determin CUTOFF & INTERACTION strengh
	      if(abs(abs(i-j)-1.).le.1.D-5) then
		 wiwj=wctc
				! write(*,*) i,j,wctc,wcc
	      else
		 wiwj=wcc
	      endif
	      cutoff=cutcc
c       
	      if(sij.le.cutoff) then
c
		 hessian1=-xx*xx*rr2ij*wiwj
		 hessian2=-yy*xx*rr2ij*wiwj
		 hessian3=-zz*xx*rr2ij*wiwj
		 hessian4=-yy*yy*rr2ij*wiwj
		 hessian5=-zz*yy*rr2ij*wiwj
		 hessian6=-zz*zz*rr2ij*wiwj
		 i1=3*(i-1)+1
		 j1=3*(j-1)+1
		 ki0=(i1-1)*NCA3-((i1-1)*(i1-2))/2
		 kh=ki0+j1-i1+1 !kh for "line i1, column j1"
		 hessian(kh)=hessian1 !line i1, column j1
		 hessian(kh+1)=hessian2 !line i1, column j1+1
		 hessian(kh+2)=hessian3 !line i1, column j1+2
c       
		 kh=ki0+NCA3-i1+1+j1-i1 !kh for "line i1+1, column j1"
		 hessian(kh)=hessian2 !line i1+1, column j1
		 hessian(kh+1)=hessian4 !line i1+1, column j1+1
		 hessian(kh+2)=hessian5 !line i1+1, column j1
c       
		 kh=ki0+NCA3-i1+1+NCA3-i1+j1-i1-1 !kh for "line i1+2, column j1"
		 hessian(kh)=hessian3 !line i1+2, column j1
		 hessian(kh+1)=hessian5 !line i1+2, column j1
		 hessian(kh+2)=hessian6 !line i1+2, column j1   
	      endif
	   enddo
	enddo

ccccccccccccccccccccccc
! --Diagonal-Block of Hessian
	do i=1,nca
	   i1=(i-1)*3+1 
	   ki0=(i1-1)*NCA3-((i1-1)*(i1-2))/2
	   hessian1=0.0
	   hessian2=0.0
	   hessian3=0.0
	   hessian4=0.0
	   hessian5=0.0
	   hessian6=0.0
	   do j=i+1,nca
	      j1=3*(j-1)+1
	      kh=ki0+j1-i1+1	!kh for "line i1, column j1"
	      hessian1=hessian1+hessian(kh) !line i1, column j1
	      hessian2=hessian2+hessian(kh+1) !line i1, column j1+1
	      hessian3=hessian3+hessian(kh+2) !line i1, column j1+2

	      kh=ki0+NCA3-i1+1+j1-i1 !kh for "line i1+1, column j1"
	      hessian4=hessian4+hessian(kh+1) !line i1+1, column j1+1
	      hessian5=hessian5+hessian(kh+2) !line i1+1, column j1+2
c       
	      kh=ki0+NCA3-i1+1+NCA3-i1+j1-i1-1 !kh for "line i1+2, column j1"
	      hessian6=hessian6+hessian(kh+2) !line i1+2, column j1+2
	   enddo	 
	   do j=1,i-1		!set j as row index now.
	      j1=3*(j-1)+1
	      kj0=(j1-1)*NCA3-((j1-1)*(j1-2))/2
c       
	      kh=kj0+i1-j1+1	!kh for "line j1, column i1"
	      hessian1=hessian1+hessian(kh) !line j1, column i1
	      hessian2=hessian2+hessian(kh+1) !line j1, column i1+1
	      hessian3=hessian3+hessian(kh+2) !line j1, column i1+2
c       
	      kh=kj0+NCA3-j1+1+i1-j1 !kh for "line j1+1, column i1"
	      hessian4=hessian4+hessian(kh+1) !line j1, column i1+1
	      hessian5=hessian5+hessian(kh+2) !line j1, column i1+2
c       
	      kh=kj0+NCA3-j1+1+NCA3-j1+i1-j1-1 !kh for "line j1+2, column i1"
	      hessian6=hessian6+hessian(kh+2) !line j1, column i1+2  
c     
	   enddo  
	   kh=ki0+1
	   hessian(kh)=-hessian1
	   hessian(kh+1)=-hessian2
	   hessian(kh+2)=-hessian3
	   kh=ki0+NCA3-i1+1
	   hessian(kh+1)=-hessian4
	   hessian(kh+2)=-hessian5
	   kh=ki0+NCA3-i1+1+NCA3-i1
	   hessian(kh+1)=-hessian6    
	enddo
c
	call CharmmDiagq(nca3,ione,nca3,ev0,vec0)
	if(dabs(ev0(6)).gt.ev_threshold.or.ev0(7).lt.ev_threshold) then 
	   calc_info=-1
	   infocalc='  CUTCC for protein is too small'
	   return
	endif
c
	calc_info=0
c$$$	open(unit=20,file='test.ev_s0000',status='unknown') 
c$$$	write(20,'(A,2x,I5,8x,A27)') '#',CALC_INFO, '(0:success; else: failure)' 
c$$$	write(20,90) '#',nca,0,cutcc,cutlg,cutll,wcc,wctc,wlg,wll,time2-time1
c$$$	write(20,100) '#',ndd,mode_need+ndd-1     
c$$$	do j=1,mode_need
c$$$	   write(20,'(F15.8,2x,I5)') ev0(j),j
c$$$	enddo
c$$$	close(20)
c
	return
	end
c*Desk :Hessian Matrix :Begin!
	subroutine DPA_1th_perturb(calc_info,infocalc,cacrd,hetcrd,
     &       nsys,nca3,nhet3,ndimca,ndimhet,cutlg,cutll,wlg,wll,
     &       lgpert)
	implicit none
	character*100 infocalc
	integer   calc_info,nsys,nca,nhet,nca3,nhet3,ndimca,ndimhet
	real*8    cacrd(*),hetcrd(*),x(nsys),y(nsys),z(nsys)
	real*8    xx,yy,zz,sij,dsqrt,rr2ij
	real*8    cutoff,cutlg,cutll,wlg,wll,wiwj
	real*8    hessian1,hessian2,hessian3,hessian4,hessian5,hessian6,tmp
        real*8    kev(nhet3),kvec(nhet3*nhet3)
	real*8    kessian(ndimhet),gessian(nca3,nhet3),tmpw(nhet3),g_sqrt_invk(nca3,nhet3)
	integer   nat0,nat30,ndim0,nat1,nat31,ione !nat1: maximum number of atoms that ligand binds
	parameter (nat0=1200,nat30=3*nat0,nat1=200,nat31=3*nat1)
	parameter (ione=1,ndim0=(nat30**2+nat30)/2)
	real*8    hessian(ndim0),ev0(nat30),vec0(nat30*nat30),ev(nca3)
	integer   rp(nca3),idgnz(nca3)
	integer   col_id(nat30*nat31),h_id(nat30*nat31)
	real*8    pessian(nat30*nat31),hv(nca3),lgpert
	common  /MODE0/ ev0,vec0
	common  /DIAG/  hessian

        integer   i,j,k,l,i0,j0,i1,j1,k1,ki0,kj0,kh
c
        nca=nca3/3
        nhet=nhet3/3
c   
c--     Reading CA & HETATM coordinates
	do i=1,nca
	   x(i)=cacrd(3*i-2)
	   y(i)=cacrd(3*i-1)
	   z(i)=cacrd(3*i)
	enddo
	do i=1,nhet
	   x(nca+i)=hetcrd(3*i-2)
	   y(nca+i)=hetcrd(3*i-1)
	   z(nca+i)=hetcrd(3*i)
	enddo
c       
c--  GENERATE HESSIAN
c	!Off-diagonalize elements for lower-triangle off-diagonal blocks.
c       i1=4 !(4,1) is the 1st element of the 2nd block of 1st column. 
c       j1=1
c       |1 2 3|
c       |  4 5| main diagonal block
c       |    6|
c       
	do i=1,nca3
	   idgnz(i)=0		!idex for NONZERO elements in GESSIAN
	enddo
	do i=1,nca3
	   do j=1,nhet3
	      gessian(i,j)=0.0
	   enddo
	enddo
	do i=1,ndimhet
	   kessian(i)=0.0
	enddo
c
	do i=1,nsys-1		!i<=nca: Protein residues; i>nca: Ligand residues
	   do j=nca+1,nsys	!ligand residues
	      xx=x(j)-x(i)
	      yy=y(j)-y(i)
	      zz=z(j)-z(i)
	      sij=sqrt(xx*xx+yy*yy+zz*zz) !distance between  atom I and atom J
	      rr2ij=1.0/sij/sij
	      if(sij.lt.1.) then
		 calc_info=-5100
		 infocalc='  Input CRD data Error: Protein/Ligand distance < 1' 
		 return
	      endif
				!determin CUTOFF & INTERACTION strengh
	      if(i.le.nca) then !i: protein/ligand interaction :: protein atom; j: ligand
		 wiwj= wlg
		 cutoff=cutlg
	      else		!i,j: interaction within ligand atoms
		 wiwj=wll
		 cutoff=cutll
	      endif
c            
	      if(sij.le.cutoff) then
		 hessian1=-xx*xx*rr2ij*wiwj
		 hessian2=-yy*xx*rr2ij*wiwj
		 hessian3=-zz*xx*rr2ij*wiwj
		 hessian4=-yy*yy*rr2ij*wiwj
		 hessian5=-zz*yy*rr2ij*wiwj
		 hessian6=-zz*zz*rr2ij*wiwj
		 if (i.le.nca) then ! For GESSIAN matrix -- Protein/Ligand Interaction portion
		    i1=3*(i-1)+1
		    j1=3*(j-nca-1)+1
		    gessian(i1,j1)=hessian1     !line i1,   column j1
		    gessian(i1,j1+1)=hessian2   !line i1,   column j1+1
		    gessian(i1,j1+2)=hessian3   !line i1,   column j1+2
		    gessian(i1+1,j1)=hessian2   !line i1+1, column j1
		    gessian(i1+1,j1+1)=hessian4 !line i1+1, column j1+1
		    gessian(i1+1,j1+2)=hessian5 !line i1+1, column j1+2
		    gessian(i1+2,j1)=hessian3   !line i1+2, column j1
		    gessian(i1+2,j1+1)=hessian5 !line i1+2, column j1+1
		    gessian(i1+2,j1+2)=hessian6 !line i1+2, column j1+2  
		    idgnz(i1)=1
		    idgnz(i1+1)=1
		    idgnz(i1+2)=1
		 else		! For KESSIAN matrix -- Ligand portion
		    i1=3*(i-nca-1)+1
		    ki0=(i1-1)*nhet3-((i1-1)*(i1-2))/2
		    j1=3*(j-nca-1)+1		    

		    kh=ki0+j1-i1+1 !kh for "line i1, column j1"
		    kessian(kh)=hessian1 !line i1, column j1
		    kessian(kh+1)=hessian2 !line i1, column j1+1
		    kessian(kh+2)=hessian3 !line i1, column j1+2
c       
		    kh=ki0+nhet3-i1+1+j1-i1 !kh for "line i1+1, column j1"
		    kessian(kh)=hessian2 !line i1+1, column j1
		    kessian(kh+1)=hessian4 !line i1+1, column j1+1
		    kessian(kh+2)=hessian5 !line i1+1, column j1
c       
		    kh=ki0+nhet3-i1+1+nhet3-i1+j1-i1-1 !kh for "line i1+2, column j1"
		    kessian(kh)=hessian3 !line i1+2, column j1
		    kessian(kh+1)=hessian5 !line i1+2, column j1
		    kessian(kh+2)=hessian6 !line i1+2, column j1   
		 endif
	      endif
	   enddo
	enddo
c
!diagonal-Block of Kessian
	do i=1,nhet
	   i1=(i-1)*3+1 
	   ki0=(i1-1)*nhet3-((i1-1)*(i1-2))/2
	   hessian1=0.0
	   hessian2=0.0
	   hessian3=0.0
	   hessian4=0.0
	   hessian5=0.0
	   hessian6=0.0
	   do j=i+1,nhet
	      j1=3*(j-1)+1
	      kh=ki0+j1-i1+1	!kh for "line i1, column j1"
	      hessian1=hessian1+Kessian(kh) !line i1, column j1
	      hessian2=hessian2+Kessian(kh+1) !line i1, column j1+1
	      hessian3=hessian3+Kessian(kh+2) !line i1, column j1+2
c
	      kh=ki0+nhet3-i1+1+j1-i1 !kh for "line i1+1, column j1"
	      hessian4=hessian4+Kessian(kh+1) !line i1+1, column j1+1
	      hessian5=hessian5+Kessian(kh+2) !line i1+1, column j1+2
c       
	      kh=ki0+nhet3-i1+1+nhet3-i1+j1-i1-1 !kh for "line i1+2, column j1"
	      hessian6=hessian6+Kessian(kh+2) !line i1+2, column j1+2
	   enddo	 
	   do j=1,i-1		!set j as row index now.
	      j1=3*(j-1)+1
	      kj0=(j1-1)*nhet3-((j1-1)*(j1-2))/2
c       
	      kh=kj0+i1-j1+1	!kh for "line j1, column i1"
	      hessian1=hessian1+Kessian(kh) !line j1, column i1
	      hessian2=hessian2+Kessian(kh+1) !line j1, column i1+1
	      hessian3=hessian3+Kessian(kh+2) !line j1, column i1+2
c       
	      kh=kj0+nhet3-j1+1+i1-j1 !kh for "line j1+1, column i1"
	      hessian4=hessian4+Kessian(kh+1) !line j1, column i1+1
	      hessian5=hessian5+Kessian(kh+2) !line j1, column i1+2
c       
	      kh=kj0+nhet3-j1+1+nhet3-j1+i1-j1-1 !kh for "line j1+2, column i1"
	      hessian6=hessian6+Kessian(kh+2) !line j1, column i1+2  
	   enddo  
	   do j=1,nca
	      j1=3*(j-1)+1
	      hessian1=hessian1+gessian(j1,i1)
	      hessian2=hessian2+gessian(j1,i1+1)
	      hessian3=hessian3+gessian(j1,i1+2)
	      hessian4=hessian4+gessian(j1+1,i1+1)
	      hessian5=hessian5+gessian(j1+1,i1+2)
	      hessian6=hessian6+gessian(j1+2,i1+2)
	   enddo
	   kh=ki0+1
	   kessian(kh)=-hessian1
	   kessian(kh+1)=-hessian2
	   kessian(kh+2)=-hessian3
	   kh=ki0+nhet3-i1+1
	   kessian(kh+1)=-hessian4
	   kessian(kh+2)=-hessian5
	   kh=ki0+nhet3-i1+1+nhet3-i1
	   kessian(kh+1)=-hessian6
	enddo
c       
!Diagonalize Kessian
	do i=1,ndimhet
	   hessian(i)=kessian(i)
	enddo
	call CharmmDiagq(nhet3,ione,nhet3,kev,kvec)
	do i=1,nhet3		!check inver of  Kessian        
	   if(Dabs(kev(i)) .le. 1.d-8) then
              calc_info=-9
	      infocalc=' : Failure: SINGULAR MATRIX FOR KESSIAN'
	      return 
	   endif
	enddo	
c
c-- 2)generate quasi inverse of Kessian
	do i=1,nhet3
	   tmpw(i)=1./dsqrt(kev(i))
	enddo	
	j0=0
	do j=1,nhet3
	   do i=1,nhet3
	      kvec(j0+i)=kvec(j0+i)*tmpw(j)
	   enddo
	   j0=j0+nhet3
	enddo
c
c--     3) generate qausi-perturbed Matrix (QpM): Gessian x SQRT(Kessian^-1) = GESS x KESS^(-1/2)
	do i=1,nca
	   do j=1,nhet3
	      g_sqrt_invk(i,j)=0.
	   enddo
	enddo
        do i=1,nca3
	   if(idgnz(i).ne.0) then
	      j0=0
	      do j=1,nhet3
		 tmp=0.
		 do k=1,nhet3
		    tmp=tmp+gessian(i,k)*kvec(j0+k)
		 enddo
		 g_sqrt_invk(i,j)=tmp
		 j0=j0+nhet3
	      enddo
	   endif
	enddo
c--
	do i=1,nca3
	   rp(i)=0		!column index of NONZERO element in ith row of  Perturbation Matrix 
	enddo
c
	kh=0
	i1=0
        do i=1,nca3		!pessian(i,i)=0 if ith atom do NOT inter with Ligands
	   j1=i1
           do j=i,nca3
	      if(idgnz(j).ne.0.and.idgnz(i).ne.0) then
c
c-- CALC. PESSIAN
		 tmp=0.
		 do k=1,nhet3
		    tmp=tmp+g_sqrt_invk(i,k)*g_sqrt_invk(j,k) !from  "G x (K^-1) x G"
		 enddo    
		 if(j.le.(i+mod(3-mod(i,3),3))) then !Diagonal Term as required from sum_i_k=0
		    i0=mod(j,3)
		    if(i0.eq.0) i0=3
		    do j0=i0,nhet3,3 
		       tmp=tmp+gessian(i,j0)
		    enddo 
		 endif
		 kh=kh+1
		 rp(i)=rp(i)+1
		 if(rp(i).gt.nat31) then
		    calc_info=-2000
		    infocalc=' Declared Dimension of RP (neighboring contacts) is too small'
		    return
		 endif
		 pessian(kh)=-tmp 
c
c-- INDEX for PESSIAN
c		 col_id(i,rp(i))=j !in row i, column index of NONZERO element is j in origianl Hessian
c		 h_id(i,rp(i))=kh !in row i, nonzero element's value is the kh element of New Hessian
		 col_id(i1+rp(i))=j
		 h_id(i1+rp(i))=kh
c
		 if(j.ne.i)  then
		    rp(j)=rp(j)+1
		    if(rp(j).gt.nat31) then
		       calc_info=-2000
		       infocalc=' Declared Dimension of RP (neighboring contacts) is too small'
		       return
		    endif
c		    h_id(j,rp(j))=kh !in row j, nonzero element's value is kh_th element of New Hessian
c		    col_id(j,rp(j))=i !in row j,column index of NONZERO element is i_th in origianl Hessian
		    col_id(j1+rp(j))=i
		    h_id(j1+rp(j))=kh
c
		 endif
	      endif
	      j1=j1+nat31
           enddo
	   i1=i1+nat31
        enddo
c
	i0=0
	do i=1,nca3
	   k1=0
	   do k=1,nca3
	      tmp=0.
	      do j=1,rp(k)
c		 tmp=tmp+pessian(h_id(k,j))*vec0(i0+col_id(k,j)) !Sum_j {P(k,j)*V(j,i)}
		 tmp=tmp+pessian(h_id(k1+j))*vec0(i0+col_id(k1+j)) !Sum_j {P(k,j)*V(j,i)}
	      enddo
	      hv(k)=tmp		!H(k,i)
	      k1=k1+nat31
	   enddo
	   tmp=0.
	   do k=1,nca3
	      tmp=tmp+vec0(k+i0)*hv(k) !Sum_k {V(k,i)*H(k,i)}
	   enddo
	   ev(i)=ev0(i)+tmp	!write(11,'(1x,I8,2x,3F12.5)') i,tmp,ev(i),ev0(i)
	   i0=i0+nca3
	enddo
c	
	lgpert=0.0
	do i=7,nca3
	   if(ev(i).le.0.001) then
	      lgpert=-9999.
	      return
	   else
	      lgpert=lgpert+dlog(ev(i)/ev0(i))
	   endif
	enddo
	calc_info=0
c
	return
	end
