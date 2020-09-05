! TIDAL PREDICTION PROGRAM   FOR G GODIN   OCEANOGRAPHY   DEC/72
! THE MOST RECENT REVISION WAS DONE IN OCTOBER 1983 BY M FOREMAN.
!  SEE PMSR 77-10
!
! INPUT
! =====
!  1) THE STANDARD CONSTITUENT DATA PACKAGE IS READ BY ROUTINE ASTRO.
!     SEE THIS SUBROUTINE FOR A DESCRIPTION OF THE INPUT REQUIRED. IN
!     ORDER TO INSURE A CONSISTENT INTERPRETATION OF THE AMPLITUDES AND
!     PHASES, THE DATA PACKAGE SHOULD BE IDENTICAL TO ONE USED BY THE
!     TIDAL HEIGHTS ANALYSIS PROGRAM.
!
!  2) ONE RECORD WITH THE TIDAL STATION INFORMATION ISTN,(NA(J),J=1,4),
!     ITZONE,LAD,LAM,LOD,LOM IN THE FORMAT
!     (5X,I4,1X,3A6,A4,A3,1X,I2,1X,I2,2X,I3,1X,I2), WHERE
!          ISTN    = STATION NUMBER
!          NA      = STATION NAME
!          ITZONE  = TIME ZONE
!          LAD,LAM = LATITUDE IN DEGREES AND MINUTES OF THE STATION
!          LOD,LOM = LONGITUDE IN DEGREES AND MINUTES OF THE STATION.
!     FOLLOWED BY ONE RECORD FOR EACH CONSTITUENT TO BE INCLUDED IN THE
!     PREDICTION, WITH THE CONSTITUENT NAME(KON), AMPLITUDE(AMP), AND
!     PHASE LAG(G) IN THE FORMAT (5X,A5,28X,F8.4,F7.2). THE PHASE LAG
!     UNITS SHOULD BE DEGREES WHILE THE UNITS OF THE PREDICTED TIDAL
!     HEIGHTS WILL BE THE SAME AS THOSE OF THE INPUT AMPLITUDES.
!     THE LAST CONSTITUENT RECORD SHOULD BE FOLLOWED BY A BLANK RECORD.
!
!  3) ONE RECORD CONTAINING THE FOLLOWING INFORMATION ON THE PERIOD AND
!     TYPE OF PREDICTION DESIRED, THE FORMAT IS (3I3,1X,3I3,1X,A4,F9.5)
!		 ICE0=INITIAL CENTURY
!        ICEE= END CENTURY
!        IDY0,IMO0,IYR0 = DAY, MONTH, AND YEAR OF THE BEGINNING OF THE
!                         PREDICTION.
!        IDYE,IMOE,IYRE = DAY, MONTH, AND YEAR OF THE END OF THE
!                         PREDICTION.
!        ITYPE          = 'EQUI' IF EQUALLY SPACED PREDICTIONS ARE
!                         DESIRED,
!                         'EXTR' IF DAILY HIGH LOWS ARE DESIRED.
!        DT             = TIME SPACING OF THE PREDICTED VALUES IF
!                         ITYPE='EQUI',
!                         TIME STEP INCREMENT USED TO INITIALLY BRACKET
!                         A HIGH OR LOW IF ITYPE='EXTR'.
!     RECOMMENDED VALUES FOR DT WHEN ITYPE='EXTR' ARE GIVEN IN THE USER
!     MANUAL. EQUALLY SPACED PREDICTIONS BEGIN AT DT HOURS ON THE FIRST
!     DAY AND END AT 2400 HOURS OF THE LAST DAY.
!
!     TYPE 3 DATA MAY BE REPEATED ANY NUMBER OF TIMES . A BLANK RECORD
!     CAUSES THE PROGRAM TO GO BACK TO READ MORE TYPE 2 DATA . TWO
!     BLANKS TERMINATE THE JOB.
!
! OUTPUT
! ======
!     EQUI-SPACED PREDICTIONS ARE PRINTED 8 VALUES PER LINE . HI-LO
!     PREDICTIONS ARE PRINTED IN THE FORMAT OF THE CANADIAN TIDES AND
!     WATER LEVELS DAILY HI-LO CARD.
!     AS WELL THE PREDICTIONS ARE PUT ONTO LOGICAL UNIT 10 IN SIMILAR
!     CARD IMAGE FORMAT , WITH A BLANK CARD AFTER EACH SET OF PREDS.
!     A LIST OF THE CONSTITUENTS USED IS PRINTED AFTER EACH SET OF
!     PREDICTIONS.
!
! METHOD OF COMPUTING HIGH-LOW PREDICTIONS
! ========================================
!     THE FIRST DERIVATIVE OF THE TIDAL HGT IS CALCULATED AS THE SUM OF
!     A NUMBER OF TRIGONOMETRIC TERMS (THE CONSTITUENTS) . THE DERIV-
!     ATIVE IS FOUND AT SUCCESIVE TIMES IN STEPS OF DT UNTIL A CHANGE IN
!     SIGN IS NOTED. THE LAST TIME INTERVAL IS THEN BISECTED REPEATEDLY
!     UNTIL THE ROOT IS BRACKETED TO WITHIN 6 MINUTES AND THEN LINEAR
!     INTERPOLATION IS USED TO GET THE ROOT .  THE CHEBYSHEV ITERATION
!     FORMULA ... T(N+1) = 2*COS(SIGMA*DT)*T(N)-T(N-1) ...
!     IS USED TO UPDATE EACH CONSTITUENT IN BOTH THE STEPPING AND
!     BISECTION PROCESSES .  FINALLY THE TIDAL HGT AT EXTREMUM IS
!     FOUND USING TABLE LOOK-UP FOR THE COSINES NEEDED.
!     FOR EQUI-SPACED PREDICTIONS THE CHEBYSHEV FORMULA IS APPLIED
!     DIRECTLY.  THE NODAL CORRECTIONS AND THE CHEBYSHEV  ITERATION ARE
!     RENEWED AT THE START OF EACH MONTH .
!
! ROUTINES NEEDED
! ===============
!     ASTRO ... GETS FREQUENCIES,ASTRO ARGS AND NODAL CORRECTIONS
!     PUT ... BUFFERS OUT HI-LO PREDICTIONS
!     HPUT ... BUFFERS OUT EQUI-SPACE PREDICTIONS
!     GDAY ... A DAY CALENDAR PROGRAM
!
! TIME OF EXECUTION ON CDC 6400 COMPUTER
! ======================================
!     68 SEC FOR 1 YR OF HI-LOS FOR QUEBEC (STN 3250) 174 CONS,DT=3 HRS
!     54 SEC FOR 1 YR OF HOURLY HEIGHTS FOR QUEBEC
!     44 SEC FOR 1YR OF HI-LOS FOR VICTORIA (STN 7120) 62 CONS,DT=.5HRS
!
!
!
! DESCRIPTION OF IMPORTANT VARIABLES OR ARRAYS
! ============================================
! MTAB                  IS THE NUMBER OF CONSTITUENTS INCLUDED IN THE
!                       STANDARD CONSTITUENT DATA PACKAGE. AT PRESENT
!                       THIS IS 146.
! M                     IS THE NUMBER OF CONSTITUENTS TO BE INCLUDED IN
!                       THE PREDICTION.
! KONTAB(I),SIGTAB(I)   ARE ARRAYS CONTAINING THE CONSTITUENT NAMES AND
!                       FREQUENCIES AS THEY ARE READ IN AND CALCULATED
!                       FROM THE STANDARD CONSTITUENT DATA PACKAGE
! V(I),U(I),F(I)        ARE ARRAYS CONTAINING THE ASTRONOMICAL ARGUMENT,
!                       SATELLITE(NODAL MODULATION) PHASE CORRECTION AND
!                       SATELLITE AMPLITUDE CORRECTION RESPECTIVELY, FOR
!                       CONSTITUENT KONTAB(I).
! KON(I),SIG(I)         ARE ARRAYS CONTAINING THE NAMES AND FREQUENCIES
!                       OF CONSTITUENTS TO BE INCLUDED IN THE
!                       PREDICTION.
! AMP(I),G(I)           ARE ARRAYS CONTAINING THE AMPLITUDE AND PHASE
!                       LAG RESPECTIVELY FOR CONSTITUENT KON(I).
! INDX(I)               IS AN ARRAY WHOSE VALUE IS THE LOCATION OF
!                       CONSTITUENT KON(I) IN THE LIST KONTAB(K).
! TWOC(I),BTWOC(I,K)    ARE ARRAYS FOR STORING THE VALUES
!                       2*COS(TWOPI*SIG(I)*DT) AND
!                       2*COS(TWOPI*SIG(I)*DT/2**K) RESPECTIVELY.
! CHP(I),CH(I),CHA(I),  ARE ARRAYS USED FOR STORING TIDAL HEIGHTS IN
! CHB(I),CHM(I)         THE CHEBYSHEV ITERATION OR THE BISECTION
!                       SEARCH METHOD.
! ANGO(I),AMPNC(I)      ARE ARRAYS USED FOR STORING THE NODALLY
!                       CORRECTED CONSTITUENT ARGUMENTS AND AMPLITUDES
!                       AT THE BEGINNING OF A CHEBYSHEV ITERATION.
! COSINE(I)             IS AN ARRAY CONTAINING 2002 COSINE VALUES IN
!                       THE RANGE OF ZERO TO TWOPI.
! KONTAB,SIGTAB,V,U, AND F SHOULD HAVE MINIMUM DIMENSION MTAB, WHILE
! KON,SIG,AMP,G,INDX,TWOC,CH,CHP,CHA,CHB,CHM,ANGO, AND AMPNC SHOULD
! HAVE MINIMUM DIMENSION M.  BTWOC SHOULD HAVE A MINIMUM DIMENSION OF M
! BY THE MAXIMUM NUMBER OF ITERATIONS REQUIRED TO LOCATE A HIGH OR LOW
! (PRESENTLY THIS IS SET TO 10).
!
!
!integer(4) function WinMain (hInstance, hPrevInstance, lpszCmdLine, nCmdShow)
!MS$ATTRIBUTES STDCALL, ALIAS : '_WinMain@16' :: WinMain
!USE:: MSFLIB
use Subproc
      TYPE HARMONIC_DATA
                character*8 NameH
                real*8      G
                real*8      A
      END TYPE HARMONIC_DATA

      TYPE DATE
                integer*4 Year
                integer*2 Month
                integer*2 Day
                integer*2 Hour
                integer*2 Minute
      END TYPE DATE

      TYPE COORDINATE
                integer*2 Degree
                integer*2 Minutes
      END TYPE COORDINATE


      TYPE OUTPUT
                integer*4 Year
                integer*2 Month
                integer*2 Day
                integer*2 Hour
                integer*2 Minute
                real*8    Value
      END TYPE OUTPUT

!       quantity of harmonics
      integer*4 Count
      TYPE (HARMONIC_DATA) HARMS(300)
      TYPE (OUTPUT) OUT(10000)
      TYPE (DATE) BDate, EDate
      TYPE (COORDINATE) Lat, Lon
      real*4    DT
!       Astronomic data file
        Character(256) AFName
!       output data file
        Character(256) RFName
!       temporary data file
        Character(256) TFName

      CHARACTER*100 FICHERO,SALIDA
      character*5 kon,kblank,kontab,konk
      REAL*8  NA(4),HOURM,HOUR0,TM,T,TTOT,H
      REAL*8  DATP(100000),HGT
	  REAL*8  TWOC,CHP,CH,DTPY,DD,PI,TWOPI
      REAL*8  BTWOC,CHA,CHB,CHM,ANG,DRAD
	  REAL*8  SIGTAB,V,U,F,AMP,G,SIG
      integer*2  di(100000),mi(100000),ho(100000),me(100000)
      integer*2  hi,mini,hf,minf
      integer*4  anoini,anofin,an(100000)
	  integer    indx(150)
      dimension KONTAB(150),kon(150)
	  dimension sigtab(150),V(150),U(150),F(150)
	  dimension AMP(150),G(150),SIG(150),H(8)
      dimension TWOC(150),CHP(150),CH(150)
      dimension BTWOC(150,10),CHA(150),CHB(150),CHM(150)
      DIMENSION ANG0(150),AMPNC(150),COSINE(2002)
      DATA KBLANK/'     '/,KEQUI/'EQUI'/
!    CHARACTER*5 DTCHAR
!      REAL*4 DT
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
      CALL get_command_argument(1, AFName)
      CALL get_command_argument(2, RFName)
      CALL get_command_argument(3, TFName)

      OPEN(9,FILE=TFName,STATUS='UNKNOWN')

      READ(9,*) Lat%Degree,Lat%Minutes
      READ(9,*) Lon%Degree,Lon%Minutes
      READ(9,*) BDate%Hour,BDate%Minute,BDate%Day,BDate%Month,BDate%Year
      READ(9,*) EDate%Hour,EDate%Minute,EDate%Day,EDate%Month,EDate%Year
      READ(9,*) DT
      DO  Count=1,1000
        READ(9,*,end=111) Harms(Count)%NameH,Harms(Count)%A,Harms(Count)%G
	  end do
111   close(9)

      OPEN(8,file=AFName,STATUS='OLD')
	  LAD=Lat%Degree
	  LAM=Lat%Minutes
	  LOD=Lon%Degree
	  LOM=Lon%Minutes

      PI=3.14159265358979
      TWOPI=2.*PI
      KR=8
!      LP=6
      DO I=1,2002
      ANGLE=TWOPI*(I-1.0)/2000.
      COSINE(I)=COS(ANGLE)
      END DO
!
! READ ASTRO ARG , DOODSON AND NODAL CORRECTION DATA .
! ======================================================================
!
!      CALL OPNAST(KONTAB,SIGTAB,V,U,F,MTAB,HOURM,XLAT)
      CALL ASTRO(KONTAB,SIGTAB,V,U,F,MTAB,HOURM,XLAT)

      IF(MTAB.EQ.0) then
      WRITE(6,*) "ERROR 2100"
      GO TO 2100
      end if
!
! READ STATION NUMBER, NAME, AND LOCATION.
! ======================================================================
!

!      write(*,*)'initial date for predictions:hour, minute,
!     1 day,month and year'
!      READ(*,'(5(I2))')HI,MINI,IDY0,IMO0,IYR0

        HI = BDate%Hour
        MINI=BDate%Minute
        IDY0=BDate%Day
        IMO0=BDate%Month
        IYR0=BDate%Year

!     write(*,*)'end date for predictions:hour, minute,
!     1 day,month and year'
!      READ(*,'(5(I2))')HF,MINF,IDYE,IMOE,IYRE

        HF = EDate%Hour
        MINF=EDate%Minute
        IDYE=EDate%Day
        IMOE=EDate%Month
        IYRE=EDate%Year
!	    iyr0=iyr0+1900
!	    iyre=iyre+1900
		ice0=int((iyr0*1.)/100.)
		icee=int((iyre*1.)/100.)
		anoini=iyr0
		iyr0=anoini-1900
	    if(ice0.gt.19)iyr0=anoini-2000
		anofin=iyre
	    iyre=anofin-1900
 	    if(icee.gt.19)iyre=anofin-2000

!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************

      ITYPE='EQUI'

      ! DT=0.25
 !     CALL get_command_argument(1, DTCHAR)
 !     READ(DTCHAR,*)DT
        dif=mini-dt*60.
        cp=int(dif/(dt*60))
		if(cp.eq.0.and.dif.lt.0)cp=-1
        qp=(hi*(1./dt))+cp
		cf=int(minf/(dt*60))
        qf=(24-hf)*(1./dt)-cf

!       dt=1.0
!     se establece el numero de horas por encima de 1:00 h para
!     el primer d�a y el n�mero de horas por debajo de las 24:00
!     para el �ltimo d�a

      DO 1900 IBB=1,100
!      READ(9,50) ISTN,CARAC,ITZONE,LAD,LAM,LOD,LOM
!     READ(KR,50) ISTN,(NA(K),K=1,4),ITZONE,LAD,LAM,LOD,LOM
! 50   FORMAT(5X,I4,1X,A22,1X,A3,I3,I2,I3,I2)
      IF(ibb.gt.1) GO TO 2100
!
! UNLESS SPECIFIED OTHERWISE, A STATION LATITUDE OF 50 DEGREES IS
! ASSUMED FOR THE PURPOSE OF CALCULATING SATELLITE TO MAIN CONSTITUENT
! AMPLITUDE RATIOS.
!=======================================================================
!
      XLAT=LAD+ISIGN(LAM,LAD)/60.
      IF(ABS(XLAT).LT.5.) XLAT=SIGN(5.,XLAT)
!
! READ THE NAME, AMPLITUDE, AND GREENWICH PHASE LAG OF THE CONSTITUENTS
! TO BE USED IN THE PREDICTION.  THE CORRESPONDING FREQUENCIES ARE
! FOUND BY SCANNING THE LIST SIGTAB.  THE SEARCH METHOD WILL BE MOST
! EFFICIENT IF ON INPUT, THESE CONSTITUENTS ARE LISTED IN ORDER OF
! INCREASING FREQUENCY.
!=======================================================================
!


!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
      KT1=1
      DO 110 K=1,Count
        KON(k)=Harms(k)%NameH
        AMP(K)=Harms(k)%A
        G(K)=Harms(k)%G
!      READ(9,*,end=120) KON(K),AMP(K),G(K)
!      write(*,*) KON(K),AMP(K),G(K)
!70    FORMAT(5X,A5,28X,F8.4,F7.2)
!      IF(KON(K).EQ.KBLANK) GO TO 120
      G(K)=G(K)/360.
      KONK=KON(K)
      KTE=KT1-1+MTAB
      DO 80 KKT=KT1,KTE
      KT=KKT
      IF(KT.GT.MTAB) KT=KT-MTAB

      IF(KONTAB(KT).EQ.KONK) GO TO 100
 80   CONTINUE
!      WRITE(*,90) KONK
! 90   FORMAT(/' ROUTINE ASTRO CANNOT FIND CONSTITUENT ',A5 /)
!      STOP
 100  INDX(K)=KT
      SIG(K)=SIGTAB(KT)
!	  WRITE(11,*)SIG(K),kt
      KT1=KT+1
 110  CONTINUE

 120  M=Count

!      CLOSE(9)
!      OPEN(8,FILE=FICHERO,STATUS='OLD')
!     aqui se corrigen las fases del fichero de ctes arm�nicas
!     en funcion de la cantidad dif=mini/60
!       read(KR,'(A100)')LEEME
!      do i=1,m
!      READ(KR,'(11X,F11.8)')FREQU(I-1)
!       WRITE(*,*)FREQU(I-1)
!      END DO
!      IF(HC.EQ.0)CORR=0.25

! READ START AND END DATES , TYPE OF PREDICTIONS AND TIME INCREMENT .
! ======================================================================
!       CORRECCION DIA 23-3-88 PEMA
!     DO 1800 ICC=1,100
!      READ(KR,130) IDY0,IMO0,IYR0,IDYE,IMOE,IYRE,ITYPE,DT
!130  FORMAT(3I3,1X,3I3,1X,A4,F9.5)
!     IF(IDY0.EQ.0) GO TO 1900
!C    DT=AMIN1(6.0,DT)
!     FIN DE LA CORRECCION SE PUSIERON C DESDE PRIN HASTA AQUI
! READ START AND END DATES , TYPE OF PREDICTIONS AND TIME INCREMENT .
! ======================================================================
!       CORRECCION DIA 23-3-88 PEMA
!     DO 1800 I=1,100

!      READ(KR,130) IDY0,IMO0,IYR0,IDYE,IMOE,IYRE,ITYPE,DT
!      WRITE(*,130) IDY0,IMO0,IYR0,IDYE,IMOE,IYRE,ITYPE,DT
!130   FORMAT(1X,I2,1X,I2,1X,I2,2X,I2,1X,I2,1X,I2,1X,A4,F9.5)
!      READ(KR,(1X,I2)) IDY0
      IF(IDY0.EQ.0) GO TO 1900
      DT=AMIN1(6.0,DT)
	  IMOP=0
      IMO=IMO0
      IYR=IYR0
	  ice=ice0
      T=0.0
      CALL GDAY(IDY0,IMO0,IYR0,ICE0,KD0)
      CALL GDAY(IDYE,IMOE,IYRE,ICEE,KDE)
!	  PRINT*,IYR0
      HOUR0=KD0*24.0
	  HOURM=KDE*24.0
	  HOURM=0.5*(HOUR0+HOURM)
      CALL ASTRO(KONTAB,SIGTAB,V,U,F,MTAB,HOURM,XLAT)
	  do 666 k=1,m
666	  sig(k)=sigtab(indx(k))

      if(cont.eq.1)dif=-dif
      DO I=1,count
      G(I)=G(I)+sig(I)*((DIF/60.)+corr)
      END DO
	  CORR=0.


!     FIN DE LA CORRECCION SE PUSIERON C DESDE PRIN HASTA AQUI

!      IDY0=1
!      IMO0=1
!      print *,' PARA QUE ANIO ES EL CALCULO   DOS CIFRAS FINALES'
!      READ (*,913)IYR0
!      IDYE=31
!      IMOE=12
!      IYRE=IYR0
! 913  FORMAT(I2)
!      ITYPE='EXTR'
!      DT=3
! SET UP COEFFICIENTS FOR CHEBYSHEV ITERATION AND INITIALIZE .
! ======================================================================
!
 135  CONTINUE
      LMAX=0.99+ALOG(DT/0.1)*1.442695
      DO 150 K=1,M
      DRAD=TWOPI*SIG(K)*DT
      TWOC(K)=2.0*COS(DRAD)
      IF(ITYPE.EQ.KEQUI)  GO TO 150
      DO 140 L=1,LMAX
      DRAD=DRAD/2.0
      TPY=2.0*COS(DRAD)
      BTWOC(K,L)=TPY
      IF(ABS(TPY).GT. 0.01) GO TO 140
      DT=DT*.99
      GO TO 135
 140  CONTINUE
 150  CONTINUE
!
!
!
! HERE CALCULATIONS AND OUTPUT OF EQUI-SPACED PREDICTIONS ARE
! PERFORMED. THEY BEGIN AT DT HOURS INTO THE FIRST DAY AND ARE
! SPECIFIED EVERY DT HOURS THEREAFTER UP TO AND INCLUDING THE LAST
! DAY.  THE ASTRONOMICAL ARGUMENT V, AND NODAL MODULATION FACTORS,
! F AND U, AND RECALCULATED AT THE 16-TH OF EVERY MONTH ARE ASSUMED
! TO BE CONSTANT THROUGHOUT THE MONTH.
!=======================================================================
!
!WRITE(*,805) ISTN,(NA(K),K=1,4)
! 805  FORMAT('1', 8X,'EQUALLY SPACED PREDICTIONS FOR STATION',I5,' , ',3
!     1A6,A4)
!     WRITE(*,171) LAD,LAM,LOD,LOM,ITZONE,DT
!
!        CALL OPNHP(10)
      NN=1
      DO 900 IDD=1,30000
      IF(IMO.EQ.IMOP) GO TO 830
      CALL GDAY(16,IMO,IYR,ICE,KDM)
      HOURM=KDM*24.0
      CALL ASTRO(KONTAB,SIGTAB,V,U,F,MTAB,HOURM,XLAT)
      TM=HOURM-HOUR0
      DO 810 K=1,M
      INDK=INDX(K)
      DTPY=F(INDK)*AMP(K)
      ANG=V(INDK)-(TM-T)*SIG(K)+U(INDK)-G(K)
      ANG=ANG-AINT(ANG)
      CHP(K)=DTPY*COS((ANG-DT*SIG(K))*TWOPI)
 810  CH(K)=DTPY*COS(ANG*TWOPI)
!
!
 830  T=T+DT
      HGT=0.0
      DO 840 K=1,M
      DTPY=CH(K)
      CH(K)=TWOC(K)*DTPY-CHP(K)
      CHP(K)=DTPY
 840  HGT=HGT+CH(K)
      TTOT=HOUR0+T
      KD=TTOT/24.0
      HR=TTOT-KD*24.0
      IMOP=IMO
      CALL DMY(IDY,IMO,IYR,ICE,KD)
!     SE HA MODIFICADO EL ULTIMO PAR�METRO DAT DE LA SUBRUTINA EN LA
!     SIGUIENTE LINEA
      CALL HPUT(ISTN,HR,IDY,IMO,IYR,HGT,DT)
      datP(NN)=HGT

      IF(KD.GT.KDE) GO TO 910
      NN=NN+1
 900  CONTINUE
! 910  CALL CLSHP(D,D,D,D,D,D,D)
! 910  CALL CLSHP
!
! WRITE OUT THE LIST OF CONSTITUENTS USED.
! ======================================================================
!
 910 CONTINUE
!      DO 1010 K=1,M
! 1010 G(K)=G(K)*360.
!      WRITE(*,1020)
! 1020 FORMAT(/// ' NAME , AMPLITUDE AND GREENWICH PHASE LAG OF CONSTITUE
!     $NTS USED TO GET THE ABOVE PREDICTIONS ' /)
!      WRITE(*,1030) (K,KON(K),AMP(K),G(K),K=1,M)
! 1030 FORMAT(4(I5,1X,A4,F8.4,F7.2,5X))
!      DO 1040 K=1,M
! 1040 G(K)=G(K)/360.
!
!
! 1800 CONTINUE
 1900 CONTINUE
 2100 CONTINUE
      k=0
      do i=1+qp,nn-qf
      k=k+1
      datp(k)=datp(i)
      end do

!     aqu� se le quiere introducir la salida de las predicciones
!     en formato columna y con horas y fechas

!     ***** a partir de aqui se ha  insertado lo que viene 10/3/98***********

!*********************************************************************

!       CALCULO DE FECHAS A PARTIR DE LA FECHA INICIAL              *
!*********************************************************************
      mi(1)=mini
      ho(1)=hi
      di(1)=idy0
      me(1)=imo0
      an(1)=iyr0
	  if(ice0.gt.19)an(1)=iyr0+2000
	  if(ice0.eq.19)an(1)=iyr0+1900

      ifila=nn-qp-qf
      incrt=dt*60.
9     minu=mi(1)
      HORA=HO(1)
      DIA=DI(1)
      MES=ME(1)
      ANO=AN(1)
8     CONTINUE
      DO I=1,IFILA
      IF(MINU.EQ.60)THEN

      HORA=HORA+1
      MINU=0
       END IF
      IF(MINU.GT.60)THEN
      HORA=HORA+1
      MINU=MINU-60
        END IF
       IF(HORA.EQ.24)THEN
          DIA=DIA+1
          HORA=0
        END IF

      if(mes.eq.2)then
      idd=dia
      imm=mes
	  siglo=int((ano*1.)/100.)


      iyy=ano-1900
	  if(siglo.gt.19)iyy=ano-2000

	feb=28
	let1=((iyy/4)*4-iyy)
	  let2=((siglo/4)*4-siglo)

	  IF(imm.EQ.2.AND.let1.EQ.0)then
	  if(IYY.NE.0.OR.let2.EQ.0)  feb=29
	  end if


!

      end if

      IF(DIA.EQ.feb+1.AND.MES.EQ.2)THEN
          MES=MES+1
          DIA=1
          END IF
      IF(DIA.EQ.31.AND.(MES.EQ.4.OR.MES.EQ.6.OR.MES.EQ.9.OR.MES.EQ.11))THEN
          MES=MES+1
          DIA=1
        END IF
      IF(DIA.EQ.32.AND.(MES.EQ.1.OR.MES.EQ.3.OR.MES.EQ.5.OR.MES.EQ.7.OR.MES.EQ.8.OR.MES.EQ.10.OR.MES.EQ.12))THEN
          MES=MES+1
          DIA=1
      END IF
          IF(MES.EQ.13)THEN
          ANO=ANO+1
         MES=1
      END IF

          MI(I)=MINU
          HO(I)=HORA
          DI(I)=DIA
          ME(I)=mes
          AN(I)=ANO
          MINU=MINU+INCRT
      END DO

10    CONTINUE
 !     print*,'output file for predictions'
 !     READ(*,'(A100)')SALIDA
 !     OPEN(11,FILE=SALIDA,STATUS='UNKNOWN')

!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************



OPEN(11,FILE=RFName,STATUS='UNKNOWN')

write(11,*) ifila
20 format(5(i5),f7.3)
do i=1,ifila
  Out(i)%Hour=ho(i)
  Out(i)%Minute=mi(i)
  Out(i)%Day=di(i)
  Out(i)%Month=me(i)
  Out(i)%Year=an(i)
  Out(i)%Value=datP(i)
  write(11,20) Out(i)%Hour, Out(i)%Minute, Out(i)%Day, Out(i)%Month, Out(i)%Year, Out(i)%Value
end do
close(11)

STOP
END
