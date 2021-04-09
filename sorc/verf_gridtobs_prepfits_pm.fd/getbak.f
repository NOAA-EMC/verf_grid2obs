      SUBROUTINE getbak(lugb,lugi,numlev,vdate,fhour,src,istat)
C                .      .    .                                       .
C SUBPROGRAM:    GETBAK      RETRIEVE Background Gridded FIELDS
C   PRGMMR: Rogers/DiMego    ORG: W/NP22     DATE: 97-06-23
C
C ABSTRACT: RETRIEVES surface & UPPER AIR Background Gridded FIELDS
C
C PROGRAM HISTORY LOG:
C   97-06-14  Geoff DiMego   Modified Morone's GETETA for PREPFITS
C   97-06-20  Eric Rogers    Modified and TESTED and Re-named GETBAK
C   99-02-25  Keith Brill    Modified to make elon1 & elonv consistently
C                            0 -> 360 degrees.
C   99-05-18  Perry Shafran  Added CAPE, CIN, and LI
C   02-07-12  Eric Rogers    Initialized masks to zero instead of one
C                            to ensure that data outside grid domain is
C                            not used
C
C USAGE:    CALL GETBAK( LUGB, LUGI, NUMLEV, IBAK, ISTAT)
C
C   INPUT FILES:
C     LUGB     - INTEGER unit for GRIB data file
C     LUGI     - INTEGER unit for GRIB Index file
C
C   INPUT ARGUMENT LIST:
C     NUMLEV   - INTEGER NUMBER OF desired LEVELS
C
C   OUTPUT ARGUMENT LIST:
C     VDATE    - VERIFYING DATE OF BACKGROUND FIELDS
C     FHR      - LENGTH OF FORECAST 
C     ISTAT    - INTEGER =0  MEANS SUCCESSFUL COMPLETION
C                        =-1 MEANS GRID COULD NOT BE RETURNED
C   OUTPUT VIA COMMON BLOCK /GRID /
C     COMMON /GRID /PGD(ILIM,JLIM,MAXLEV),
C    1 Z(ILIM,JLIM,MAXLEV), MASKZ(ILIM,JLIM,MAXLEV),
C    2 T(ILIM,JLIM,MAXLEV), MASKT(ILIM,JLIM,MAXLEV),
C    3 U(ILIM,JLIM,MAXLEV), MASKU(ILIM,JLIM,MAXLEV),
C    4 V(ILIM,JLIM,MAXLEV), MASKV(ILIM,JLIM,MAXLEV),
C    5 Q(ILIM,JLIM,MAXLEV), MASKQ(ILIM,JLIM,MAXLEV),
C    6 ALNQ(ILIM,JLIM,MAXLEV)
C
C             REAL & INTEGER ARRAYS dimensioned (ILIM,JLIM,MAXLEV) 
C             Containing 3-D fields on arrays of (IMAX,JMAX,NUMLEV)
C             And its associated MASKx INTEGER ARRAY(IMAX,JMAX,NUMLEV)
C             MASK=0 INDICATES A GRIDPOINT WITH MISSING DATA 
C                PRESENT CONTENTS:
C                PGD PRESSURE (PA) no associated mask
C                Z   GEOPOTENTIAL HEIGHT (M)
C                T   TEMPERATURE (K)
C                U   WIND U-COMPONENT (M/S)
C                V   WIND V-COMPONENT (M/S)
C                Q   Specific HUMIDITY (G/G)
C                ALNQ Natural logarithm of Q no associated mask
C
C   OUTPUT VIA COMMON BLOCK /SURFCE/
C
C     COMMON /SURFCE/ 
C    1 PS(ILIM,JLIM),MSKPS(ILIM,JLIM),ZS(ILIM,JLIM),MSKZS(ILIM,JLIM),
C    2 TS(ILIM,JLIM),MSKTS(ILIM,JLIM),QS(ILIM,JLIM),MSKQS(ILIM,JLIM), 
C    3 US(ILIM,JLIM),MSKUS(ILIM,JLIM),VS(ILIM,JLIM),MSKVS(ILIM,JLIM),
C    4 PM(ILIM,JLIM),MSKPM(ILIM,JLIM),
C    5 CAPE(ILIM,JLIM),MSKCP(ILIM,JLIM),
C    6 CIN(ILIM,JLIM),MSKCN(ILIM,JLIM),
C    7 PLI(ILIM,JLIM),MSKLI(ILIM,JLIM)
C
C                REAL & INTEGER ARRAYS dimensioned (ILIM,JLIM) 
C                Containing surface fields on arrays of (IMAX,JMAX)
C                And its associated MSKxx INTEGER ARRAY(IMAX,JMAX)
C                IF =0 INDICATES A GRIDPOINT WITH MISSING DATA 
C                PRESENT CONTENTS:
C                PS SURFACE PRESSURE (PA)
C                ZS SURFACE GEOPOTENTIAL HEIGHT (M)
C                TS SURFACE TEMPERATURE (K)
C                   (ACTUALLY TEMP AT 2 METERS)
C                US SURFACE WIND U-COMPONENT (M/S)
C                   (ACTUALLY U AT 10 METERS)
C                VS SURFACE WIND V-COMPONENT (M/S)
C                   (ACTUALLY V AT 10 METERS)
C                QS SURFACE Specific HUMIDITY (G/G)
C                   (ACTUALLY Q AT 2 METERS)
C                CAPE SURFACE based CAPE (J/kg)
C                CIN SURFACE based Conv. inhibition (J/kg)
C                PLI  SURFACE to 500-MB lifted index (K)
C
C REMARKS:

C ATTRIBUTES:
C   LANGUAGE: FORTRAN-77
C   MACHINE:  CRAY C-90
C$$$
c GRIB2
      use grib_mod     
c GRIB2
c     use gridef
c     use grid3d
c     use surfce
c     use itotal
c     use vdates
c     use guser

      INCLUDE 'parm.inc'

      LOGICAL latlong, lambert, polarstereo
      COMMON /gridef/ imax, jmax, kmax, alat1, elon1, dxx, dyy, elonv,
     +          alatan, alatan1, alatan2, latlong, lambert, polarstereo
C
      real*8 d90
      real*8 D001,D0001,D00001
      real*8 D622
      real*8 D378
      real*8 D10
      real*8 D1
      DIMENSION RICDAT(5),IFDATE(8),IVDATE(8)
      INTEGER jpds(200), jgds(200), kpds(200), kgds(200)
c      INTEGER levs(40,6), ivar(6)
      INTEGER levs(40,10), ivar(6), idisp(6)
      integer fhour
      CHARACTER*8 modesc(110)
      character*10 src
      REAL*8 vdate

      COMMON /grid/ pgd(ilim,jlim,maxlev), z(ilim,jlim,maxlev),
     +            maskz(ilim,jlim,maxlev), t(ilim,jlim,maxlev),
     +            maskt(ilim,jlim,maxlev), u(ilim,jlim,maxlev),
     +            masku(ilim,jlim,maxlev), v(ilim,jlim,maxlev),
     +            maskv(ilim,jlim,maxlev), q(ilim,jlim,maxlev),
     +            maskq(ilim,jlim,maxlev), alnq(ilim,jlim,maxlev)

      COMMON /surfce/ ps(ilim,jlim), mskps(ilim,jlim), zs(ilim,jlim),
     +            mskzs(ilim,jlim), ts(ilim,jlim), mskts(ilim,jlim),
     +            qs(ilim,jlim), mskqs(ilim,jlim), us(ilim,jlim),
     +            mskus(ilim,jlim), vs(ilim,jlim), mskvs(ilim,jlim),
     +            pm(ilim,jlim), mskpm(ilim,jlim), cape(ilim,jlim),
     +            mskcp(ilim,jlim), cin(ilim,jlim), mskcn(ilim,jlim),
     +            pli(ilim,jlim), mskli(ilim,jlim),
     +            oz1s(ilim,jlim), mskoz1(ilim,jlim),
     +            oz8s(ilim,jlim), mskoz8(ilim,jlim)

C GRIB2
      type(gribfield) :: gfld
      integer ifile,j,jdisc,jpdtn,jgdtn,iret
      integer,dimension(200) :: jids,jpdt,jgdt
      real fldmin,fldmax,firstval,lastval
      logical :: unpack=.true.
c     character(255) :: cin
C GRIB2

C


      DATA levs /1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750,
     +            725, 700, 675, 650, 625, 600, 575, 550, 525, 500, 475,
     +            450, 425, 400, 375, 350, 325, 300, 275, 250, 225, 200,
     +            175, 150, 125, 100, 75, 50, 25, 1000, 975, 950, 925, 
     +            900, 875, 850, 825, 800, 775, 750, 725, 700, 675, 650,
     +            625, 600, 575, 550, 525, 500, 475, 450, 17 * 0, 1000,
     +            950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450,
     +            400, 350, 300, 250, 200, 150, 100, 50, 20 * 0, 1000, 
     +            925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70,
     +            50, 27 * 0, 1000, 850, 700, 500, 400, 300, 250, 200,
     +            150, 100, 30 * 0, 1000, 850, 700, 500, 400, 300, 250,
     +            200, 32 * 0, 1000, 975, 950, 925, 900, 875, 850, 825, 
     *            800, 775, 750, 725, 700, 650, 600, 550, 500, 450, 400, 
     *            350, 300, 275, 250, 225, 200, 175, 150, 125, 100, 
     *            11*0,1000, 850, 700, 500, 300, 250,
     +            34 * 0, 1000, 850, 700, 500, 250, 35 * 0,
     *            850, 550, 250, 37 * 0/
C
      DATA ivar /5, 0, 2, 3, 0, 1/
      DATA idisp /3, 0, 2, 2, 1, 1/
      DATA tzero /273.15/
      DATA bmiss /10E10/
      DATA D001 /0.01/
      DATA D0001 /0.001/
      DATA D00001 /0.000001/
      DATA D1 /1.0/
      DATA D10 /10.0/
      DATA D622 /0.622/
      DATA D378 /0.378/
      data d90 /-90.0/

C
      DATA modesc /38 * 'xxxxxxxx', '80KM NGM', 37 * 'xxxxxxxx', 
     +            'T126 AVN', 'T126 MRF', 'xxxxxxxx', 'T62 MRF ', 2 * 
     +            'xxxxxxxx', '80KM ETA', '32KM ETA', '29KM ETA', 
     +            '60KM RUC', 2 * 'xxxxxxxx', '48KM ETA', 15 * 
     +            'xxxxxxxx', '40KM RUC', 3 * 'xxxxxxxx', 'RTMA    ',
     *            '10KM ETA'/
C
      CHARACTER cbuf(mbuf)
      LOGICAL*1 lb(jf)
      REAL f(jf)
      INTEGER jens(200), kens(200)

      istat = 0

C     
C GRIB2
! Set GRIB2 field identification values to search for
      j=0              ! search from 0
      jdisc=0          ! for met field:0 hydro: 1, land: 2
!-- set id section
      jids=-9999
!-- set product def template, using template 4.0
      jpdtn=8
!-- set product def array
      jpdt=-9999
C GRIB2
C     READ INDEX FILE TO GET GRID SPECS 
C     
      irgi = 1
      irgs = 1
      kmax = 0
      ifdate = 0
      ricdat = 0
      print*,'imax,jmax=',imax,jmax
c     print*,'src=',src
C     
!for pdt, define catogory, parameter and level
!eg: tmp at 500 hpa
      j=0
      jpdt(1)=13
      jpdt(2)=193
      jpdt(10)=104
c     jpdt(20)=fhour
      if(src(:5).eq.'pmmax') then
        jpdt(1)=14
        jpdt(2)=202
      endif
!-- set grid def template
      jgdtn=-1
!-- set product def array
      jgdt=-9999
        call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *         unpack,j,gfld,iret)
       if(iret.ne.0) then
           print*,'iret=',iret
           print*,'No successful call to getgb2'
           stop
        endif
        ifdate(1)=gfld%idsect(6)
        if(ifdate(1).lt.2000) ifdate(1)=ifdate(1)+2000
        ifdate(2)=gfld%idsect(7)
        ifdate(3)=gfld%idsect(8)
        ifdate(5)=gfld%idsect(9)
        print*,'fhour=',fhour
        if(src(:5).eq.'pmave'.or.src(:5).eq.'pmmax') then
         if(fhour.eq.0) then
          ricdat(2)=0
         elseif(fhour.eq.24) then
          ricdat(2)=0
         elseif(fhour.eq.30) then
          ricdat(2)=6
         elseif(fhour.eq.54) then
          ricdat(2)=30
         elseif(fhour.eq.48) then
          ricdat(2)=24
         elseif(fhour.eq.72) then
          ricdat(2)=48
         elseif(fhour.eq.78) then
          ricdat(2)=54
         endif
        else
         ricdat(2)=gfld%ipdtmpl(9) +1 
        endif
c       print*,'gfld%ipdtmpl(10)=',gfld%ipdtmpl(10)
        print*,'ifdate(1)=',ifdate(1)
        print*,'ifdate(2)=',ifdate(2)
        print*,'ifdate(3)=',ifdate(3)
        print*,'ifdate(5)=',ifdate(5)
        print*,'ricdat(2)=',ricdat(2)
      
      itmpl = gfld%igdtnum
      print*,'itmpl=',itmpl
C     
C     USE KGDS(1) TO DETERMINE GRID PROJECTION
C     
C     KGDS(1) = 0 ----> LATITUDE/LONGITUDE
C     KGDS(1) = 1 ----> MERCATOR (NOT YET UESD)
C     KGDS(1) = 3 ----> POLAR STEREOGRAPHIC
C     KGDS(1) = 5 ----> LAMBERT CONFORMAL
C     
      if(itmpl.eq.0) then
        imax = gfld%igdtmpl(8)
        jmax = gfld%igdtmpl(9)
        latlong = .true.
        lambert = .false.
        polarstereo = .false.
      else if (itmpl.eq.30) then
        imax = gfld%igdtmpl(8)
        jmax = gfld%igdtmpl(9)
        print*,'imax,jmax=',imax,jmax
        latlong = .false.
        lambert = .true.
        polarstereo = .false.
      else if (itmpl.eq.20) then
        imax = gfld%igdtmpl(8)
        jmax = gfld%igdtmpl(9)
        latlong = .false.
        lambert = .false.
        polarstereo = .true.
      ELSE
        iret = 99
        WRITE (6,*) ' KGDS(1) = ', kgds(1)
        WRITE (6,*) ' GRID CAN NOT BE USED IN THIS CODE IRET= ', iret
        istat = -1
        RETURN
C     STOP777
      END IF

      jj1 = 1
      jjinc = 1

C     
C     SET THE REST OF THE GRID PARAMETERS BASED ON PROJECTION TYPE
C     
      IF (latlong) THEN
        alat1 = gfld%igdtmpl(12) * D00001
        elon1 = gfld%igdtmpl(13) * D00001
        alatend = gfld%igdtmpl(15) * D00001 
        IF (elon1.lt.0.0) elon1 = elon1 + 360.
        elonv = 0.0
        alatan1 = 0.0
        alatan2 = 0.0
        dxx = gfld%igdtmpl(17) * D00001
        dyy = gfld%igdtmpl(18) * D00001
C       IF we have a damn global grid with pt 1,1 at the NP, flip it
        IF (alat1.ge.alatend) THEN
          alat1 = alat1 - dyy * (jmax-1)
          alat1 = max(-90.,alat1)
          jj1 = jmax
          jjinc = -1
        END IF
      END IF
C     
      IF (lambert) THEN
        alat1 = gfld%igdtmpl(10) * 0.000001
        elon1 = gfld%igdtmpl(11) * 0.000001        
        IF (elon1.lt.0.0) elon1 = elon1 + 360.
        elonv = gfld%igdtmpl(14) * 0.000001
        IF (elonv.lt.0.0) elonv = elonv + 360.
        alatan1 = gfld%igdtmpl(19) * 0.000001
        alatan2 = gfld%igdtmpl(20) * 0.000001
        dxx = gfld%igdtmpl(15) * 0.000001
        dyy = gfld%igdtmpl(16) * 0.000001
      END IF
C     
      IF (polarstereo) THEN
        alat1 = gfld%igdtmpl(10) * 0.000001
        elon1 = gfld%igdtmpl(11) * 0.000001
        IF (elon1.lt.0.0) elon1 = elon1 + 360.
        elonv = gfld%igdtmpl(14) * 0.000001
        IF (elonv.lt.0.0) elonv = elonv + 360.
        alatan1 = 0.0
        alatan2 = 0.0
        dxx = gfld%igdtmpl(15) * 0.000001
        dyy = gfld%igdtmpl(16) * 0.000001
      END IF
C     
C     ADD IGF FORECAST HOURS TO START DATE TO GET VALID DATE
C     ADDDATE IS MODIFIED VERSION OF W3ADDATE
C     
      call w3movdat(RICDAT,IFDATE,IVDATE)
      print*,'ivdate=',ivdate
      vdate = (ivdate(1)*1000000+ivdate(2)*10000+ivdate(3)*100+
     1          ivdate(5))
      print*,'vdate=',vdate
C     
C     GET GRID NUMBER FROM PDS
C     
c     igdnum = kpds(3)
C     
      WRITE (6,*) ' WELCOME TO THE PREPFITS GRID DECODER '
c     WRITE (6,*) ' THE GRID YOU HAVE CHOSEN IS NUMBER ', igdnum
      IF (latlong) THEN
        WRITE (6,*) ' A LAT/LON GRID WITH RES= ', dxx, ' BY ', dyy, 
     +              ' DEG'
        WRITE (6,*) ' AND ORIGIN (1,1) @ ', alat1, ' LAT & ', elon1, 
     +              ' LONG'
      ELSE IF (polarstereo) THEN
        WRITE (6,*) ' A POLAR STEREO GRID CENTERED AT ', elonv, ' DEG E'
        WRITE (6,*) ' AND A HORIZONTAL RESOLUTION OF ', dxx, ' KM'
      ELSE IF (lambert) THEN
        WRITE (6,*) ' A LAMBERT CONFORMAL GRID CENTERED AT ', elonv, 
     +              ' DEG E'
        WRITE (6,*) ' AND A HORIZONTAL RESOLUTION OF ', dxx, ' KM'
      END IF
      WRITE (6,*) ' HORIZONTAL DIMENSIONS OF ', imax, ' X', jmax
c     modgen = kpds(2)
c     WRITE (6,*) ' THE MODEL GENERATION NUMBER IS ', modgen
c     WRITE (6,*) ' WHICH IS THE ', modesc(modgen)
C     
C     PROCESS THE GRIB FILE
C     
      numval = imax * jmax
C     initialize mask arrays
      mskoz1 = 0
C     default thin mandatory level set NUMLEV<=10
      levset = 5
      kmax = min(10,numlev)
c     VSREF set
      if (numlev.ge.3) then
        levset = 10
        kmax = min(4,numlev)
      endif
C     HREF set
      IF (numlev.ge.5) THEN
        levset = 9
        kmax = min(6,numlev)
      END IF
C     SREF para set
      IF (numlev.ge.6) THEN
        levset = 8
        kmax = min(6,numlev)
      END IF
C     full mandatory level set
      IF (numlev.ge.11) THEN
        levset = 4
        kmax = min(13,numlev)
      END IF
C     full mandatory level set but no 150 or 100 mb
      IF (numlev.ge.8) THEN
        levset = 6
        kmax = min(10,numlev)
      END IF
C     full mandatory level set
      IF (numlev.ge.11) THEN
        levset = 4
        kmax = min(13,numlev)
      END IF
C     every 50mb
      IF (numlev.ge.19) THEN
        levset = 3
        kmax = min(20,numlev)
      END IF
C     every 25mb but low levels only
      IF (numlev.ge.23) THEN
        levset = 2
        kmax = min(23,numlev)
      END IF
c RR set
      if (numlev.eq.29) then
       levset=7
       kmax = min(29,numlev)
      endif 
C     every 25mb
      IF (numlev.ge.37) THEN
        levset = 1
        kmax = min(40,numlev)
      END IF
C     
C     Fill pressure array
C     
      DO l = 1, kmax
        DO j = 1, jmax
          DO i = 1, imax
            pgd(i,j,l) = levs(l,levset)
          END DO
        END DO
      END DO
C     
      WRITE (6,1200) imax, jmax, numlev, kmax, levset
 1200 FORMAT (' IMAX,JMAX,NUMLEV,KMAX,LEVSET ',5I4)

C     -== GET SURFACE FIELDS ==-
      l = 0
      iv = 0

C     SURFACE PARTICULATE MATTER

      j = 0
      jids=-9999
c     jpdtn=0
      jpdt=-9999
      jpdt(1)=13
      jpdt(2)=193
      jpdt(10)=104
      if(src(:5).eq.'pmmax') then
        jpdt(1)=14
        jpdt(2)=202
      endif
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'PM iret=',iret
      if(iret.eq.0) then
       ii = 1
       jj = jj1
       do kk = 1, itot
        oz1s(ii,jj) = gfld%fld(kk)
c       if(oz1s(ii,jj).gt.0) print*,'ii,jj,oz1s=',ii,jj,oz1s(ii,jj)
         IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
          mskoz1(ii,jj) =1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 257
          END IF
        END DO
  257   CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
 1300   FORMAT (' IV, IVAR, L, IRET:  ',4I5)
      END IF
 

C         PRINT VALUES AT POINT IN MIDDLE OF GRID
          jer = jmax / 2
          ier = imax / 2
          WRITE (6,1400) ier, jer, pm(ier,jer), ps(ier,jer), 
     +                zs(ier,jer), ts(ier,jer), us(ier,jer), 
     +                vs(ier,jer), qs(ier,jer)
 1400     FORMAT (2I4,3X,F7.1,/,3X,F7.1,2X,F7.1,1X,3F8.2,2X,F8.6)
          print*,'numlev=',numlev
          DO ler = 1, numlev
            WRITE (6,1500) ler, pgd(ier,jer,ler), z(ier,jer,ler), 
     +                  t(ier,jer,ler), u(ier,jer,ler), v(ier,jer,ler),
     +                  q(ier,jer,ler), alnq(ier,jer,ler)
 1500       FORMAT (1X,I2,2X,F5.0,2X,F7.1,1X,3F8.2,2X,F8.6,2X,F7.3)
          END DO

          call gf_free(gfld)

          RETURN
          END
