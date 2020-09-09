      SUBROUTINE getbak(lugb,lugi,numlev,fhr,iearth,
     *    ismart,ipcp,istat,fhour)
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
      use gridef
      use grid3d
      use surfce
      use itotal
      use vdates
      use guser
      INCLUDE 'parm.inc'
C
c     real*8 z,t,u,v,q,alnq
C
      real*8 d90
      real*8 D001,D0001
      real*8 D622
      real*8 D378
      real*8 D10
      real*8 D1
      real*8 vap, vaps
      real*8 rh
      real*8 dgrid(itot)

C GRIB2
      type(gribfield) :: gfld
      integer ifile,j,jdisc,jpdtn,jgdtn,iret
      integer,dimension(200) :: jids,jpdt,jgdt
      real fldmin,fldmax,firstval,lastval
      logical :: unpack=.true.
c     character(255) :: cin
C GRIB2

c     double precision D001, gridd, pmd
      DIMENSION RICDAT(5),IFDATE(8),IVDATE(8)
      INTEGER jpds(200), jgds(200), kpds(200), kgds(200)
c      INTEGER levs(40,6), ivar(6)
      INTEGER levs(40,10), ivar(6), idisp(6)
      integer fhour
      CHARACTER*8 modesc(110)
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
      DATA ivar /7, 11, 33, 34, 51, 52/
      DATA idisp /3, 0, 2, 2, 1, 1/
      DATA tzero /273.15/
      DATA bmiss /10E10/
      DATA D001 /0.01/
      DATA D0001 /0.001/
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
      print*,'top of routine'

      allocate(z(ilim,jlim,maxlev))
      allocate(t(ilim,jlim,maxlev))
      allocate(u(ilim,jlim,maxlev))
      allocate(v(ilim,jlim,maxlev))
      allocate(q(ilim,jlim,maxlev))
      allocate(alnq(ilim,jlim,maxlev))

      allocate(pgd(ilim,jlim,maxlev))
      allocate(maskz(ilim,jlim,maxlev))
      allocate(maskt(ilim,jlim,maxlev))
      allocate(masku(ilim,jlim,maxlev))
      allocate(maskv(ilim,jlim,maxlev))
      allocate(maskq(ilim,jlim,maxlev))

      allocate(ps(ilim,jlim))
      allocate(zs(ilim,jlim))
      allocate(ts(ilim,jlim))
      allocate(qs(ilim,jlim))
      allocate(us(ilim,jlim))
      allocate(vs(ilim,jlim))
      allocate(pm(ilim,jlim))

      allocate(cape(ilim,jlim))
      allocate(cin(ilim,jlim))
      allocate(bcape(ilim,jlim))
      allocate(pli(ilim,jlim))

      allocate(vis(ilim,jlim))
      allocate(tmax(ilim,jlim))
      allocate(tmin(ilim,jlim))
      allocate(dpt(ilim,jlim))
      allocate(tocc(ilim,jlim))
      allocate(gust(ilim,jlim))
      allocate(pbls(ilim,jlim))

      allocate(pwo(ilim,jlim))
      allocate(trop(ilim,jlim))
      allocate(cdtz(ilim,jlim))
      allocate(pblris(ilim,jlim))
      allocate(pblmxs(ilim,jlim)) 

      allocate(aod(ilim,jlim))

      allocate(mskps(ilim,jlim))
      allocate(mskzs(ilim,jlim))
      allocate(mskts(ilim,jlim))
      allocate(mskqs(ilim,jlim))
      allocate(mskus(ilim,jlim))
      allocate(mskvs(ilim,jlim))
      allocate(mskpm(ilim,jlim))

      allocate(mskcp(ilim,jlim))
      allocate(mskcn(ilim,jlim))
      allocate(mskbc(ilim,jlim))
      allocate(mskli(ilim,jlim))

      allocate(mskvis(ilim,jlim))
      allocate(msktmax(ilim,jlim))
      allocate(msktmin(ilim,jlim))
      allocate(mskdpt(ilim,jlim))
      allocate(msktocc(ilim,jlim))
      allocate(mskgust(ilim,jlim))
      allocate(mskpbl(ilim,jlim))

      allocate(mskpw(ilim,jlim))
      allocate(msktrp(ilim,jlim))
      allocate(mskcdtz(ilim,jlim))
      allocate(mskpri(ilim,jlim))
      allocate(mskpmx(ilim,jlim))

      allocate(mskaod(ilim,jlim))

      allocate(grid(itot))
      allocate(grid2(itot))
      allocate(grid3(itot))
      allocate(grid4(itot))

      allocate(mask(itot))

C     
C     READ INDEX FILE TO GET GRID SPECS 
C     
C GRIB2
! Set GRIB2 field identification values to search for
      j=0              ! search from 0
      jdisc=0          ! for met field:0 hydro: 1, land: 2
!-- set id section
      jids=-9999
!-- set product def template, using template 4.0
      jpdtn=0
!-- set product def array
      jpdt=-9999
C GRIB2
      irgi = 1
      irgs = 1
      kmax = 0
      kpds = 0
      kgds = 0
      ifdate = 0
      ricdat = 0
c     jr = 0
c     kskip = 0
c     CALL getgi(lugi,kskip,mbuf,cbuf,nlen,nnum,irgi)
c     REWIND lugi
C     
C     NNUM IS THE NUMBER OF GRIB MESSAGES IN THE FILE. READ THE
C     FIRST PDS ONLY TO GET THE GRID NUMBER AND DATE
C     
C     NOTE: SINCE NGM HAS MANY GRIDS COMBINED INTO ONE GRIB FILE
C     SEARCH THROUGH THE INDEX FILE UNTIL WE GET TO A GRID 104
C     RECORD, THEN JUMP OUT OF LOOP AND SAVE GRID 104 INFO
C     IN COMMON BLOCK
C     
c     DO k = 1, nnum
c       jr = k - 1
c       jpds = -1
c       jgds = -1
c       CALL getgb1s(cbuf,nlen,nnum,jr,jpds,jgds,jens,kr,kpds,kgds,kens,
c    +              lskip,lgrib,irgs)
!for pdt, define catogory, parameter and level
!eg: tmp at 500 hpa
      jpdtn=48
      jpdt(1)=20      ! table 4.1
      jpdt(2)=102      ! table 4.2-0-0
      jpdt(3)=62000
      print*,'fhour=',fhour
      jpdt(20)=fhour
!-- set grid def template
      jgdtn=-1
!-- set product def array
      jgdt=-9999
c       print*,'I am here before getgb2'
        print*,'lugb=',lugb
        print*,'lugi=',lugi
        call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *         unpack,j,gfld,iret)
        if(iret.ne.0) then
          print*,'iret=',iret
          print*,'No successful call to getgb2'
          stop
        endif
        ifdate(1)=gfld%idsect(6)
        ifdate(2)=gfld%idsect(7)
        ifdate(3)=gfld%idsect(8)
        ifdate(5)=gfld%idsect(9)
        ricdat(2)=fhour
C       
C       USE  F I R S T  PDS TO GET START DATE & FORECAST LENGTH
C       
c       IF (k.eq.1) THEN
c         if(kpds(9).eq.1) kpds(8)=kpds(8)+1
c         iys = kpds(8)
c         if(kpds(9).eq.1) iys=iys+1
c         fhr(ibak) = kpds(14)
c         print*,'fhr=',fhr(ibak)
c         ifdate(2) = kpds(9)
c         ifdate(3) = kpds(10)
c         ifdate(5) = kpds(11)
c         ricdat(2) = kpds(14)
c         icent = kpds(21)
c         if(ifdate(2).lt.100) then
c           ifdate(1) = (icent-1)*100 + iys
c         else
c           ifdate(1) = icent * 100
c         endif
c       END IF
C       IF  N O T  NGM, THEN THIS PDS & GDS ARE FINE
c       IF (kpds(2).ne.39) THEN
c         WRITE (6,1000) (kgds(ker),ker = 1,14)
c1000     FORMAT (1X,14I8)
c         WRITE (6,1100) (kpds(ker),ker = 1,22)
c1100     FORMAT (1X,22I5)
c         GO TO 10
c       ELSE
c         IF (kpds(3).eq.104) THEN
c           WRITE (6,1000) (kgds(ker),ker = 1,14)
c           WRITE (6,1100) (kpds(ker),ker = 1,22)
c           GO TO 10
c         END IF
c       END IF
c     END DO
c  10 CONTINUE
C     
C     FILL IN GRIDEF COMMON BLOCK
C     THE FOLLOWING DEFINED REGARDLESS OF GRID PROJECTION
C     
c     print*,'imax,jmax'
c     imax = kgds(2)
c     jmax = kgds(3)
      itmpl = gfld%igdtnum
C     
C     USE KGDS(1) TO DETERMINE GRID PROJECTION
C     
C     KGDS(1) = 0 ----> LATITUDE/LONGITUDE
C     KGDS(1) = 1 ----> MERCATOR (NOT YET UESD)
C     KGDS(1) = 3 ----> POLAR STEREOGRAPHIC
C     KGDS(1) = 5 ----> LAMBERT CONFORMAL
C     
c     IF (kgds(1).eq.0) THEN
      print*,'itmpl=',itmpl
      if(itmpl.eq.0) then
        imax = gfld%igdtmpl(8)
        jmax = gfld%igdtmpl(9)
        print*,'imax,jmax=',imax,jmax
        latlong = .true.
        lambert = .false.
        polarstereo = .false.
c     ELSE IF (kgds(1).eq.3) THEN
      else if (itmpl.eq.30) then
        imax = gfld%igdtmpl(8)
        jmax = gfld%igdtmpl(9)
        latlong = .false.
        lambert = .true.
        polarstereo = .false.
      ELSE IF (kgds(1).eq.5) THEN
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
c       alat1 = kgds(4) * D0001
c       elon1 = kgds(5) * D0001
        alat1 = gfld%igdtmpl(10) * 0.000001
        elon1 = gfld%igdtmpl(11) * 0.000001
c       alatend = kgds(7) * D0001
        alatend = gfld%igdtmpl(15) * 0.001
        IF (elon1.lt.0.0) elon1 = elon1 + 360.
        if (nint(elon1).eq.360) elon1 = 0.
        elonv = 0.0
        alatan1 = 0.0
        alatan2 = 0.0
c       alatan = 0.0
c       dxx = kgds(9) * D0001
c       dyy = kgds(10) * D0001
        dxx = gfld%igdtmpl(17) * 0.000001
        dyy = gfld%igdtmpl(18) * 0.000001
C       IF we have a damn global grid with pt 1,1 at the NP, flip it
        IF (alat1.ge.alatend) THEN
          alat1 = alat1 - dyy * (jmax-1)
c         alat1 = dmax1(-90.,alat1)
          alat1 = dmax1(d90,alat1)
          jj1 = jmax
          jjinc = -1
        END IF
      END IF
C     
      IF (lambert) THEN
c       alat1 = kgds(4) * 0.001
c       elon1 = kgds(5) * 0.001
        alat1 = gfld%igdtmpl(10) * 0.000001
        elon1 = gfld%igdtmpl(11) * 0.000001
        IF (elon1.lt.0.0) elon1 = elon1 + 360.
c       elonv = kgds(7) * 0.001
        elonv = gfld%igdtmpl(14) * 0.000001
        IF (elonv.lt.0.0) elonv = elonv + 360.
c       alatan = kgds(12) * 0.001
c       alatan1 = kgds(12) * 0.001
c       alatan2 = kgds(13) * 0.001
        alatan1 = gfld%igdtmpl(19) * 0.000001
        alatan2 = gfld%igdtmpl(20) * 0.000001
c       dxx = kgds(8) * 0.001
c       dyy = kgds(9) * 0.001
        dxx = gfld%igdtmpl(15) * 0.000001
        dyy = gfld%igdtmpl(16) * 0.000001
      END IF
C     
      IF (polarstereo) THEN
        alat1 = kgds(4) * 0.001
        elon1 = kgds(5) * 0.001
        IF (elon1.lt.0.0) elon1 = elon1 + 360.
        elonv = kgds(7) * 0.001
        IF (elonv.lt.0.0) elonv = elonv + 360.
c       alatan = 0.0
        alatan1 = 0.0
        alatan2 = 0.0
        dxx = kgds(8) * 0.001
        dyy = kgds(9) * 0.001
      END IF
C     
C     ADD IGF FORECAST HOURS TO START DATE TO GET VALID DATE
C     ADDDATE IS MODIFIED VERSION OF W3ADDATE
C     
      call w3movdat(RICDAT,IFDATE,IVDATE)
      vdate(ibak) = (ivdate(1)*1000000+ivdate(2)*10000+ivdate(3)*100+
     1          ivdate(5))
C     
C     GET GRID NUMBER FROM PDS
C     
      igdnum = kpds(3)
c     ismart=0
c     if(igdnum.eq.196.or.igdnum.eq.197.or.igdnum.eq.198.or.
c    *    igdnum.eq.181) 
c    *    ismart=1
C     
      WRITE (6,*) ' WELCOME TO THE PREPFITS GRID DECODER '
      WRITE (6,*) ' THE GRID YOU HAVE CHOSEN IS NUMBER ', igdnum
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
      modgen = kpds(2)
      WRITE (6,*) ' THE MODEL GENERATION NUMBER IS ', modgen
c     WRITE (6,*) ' WHICH IS THE ', modesc(modgen)
C     
C     PROCESS THE GRIB FILE
C     
      numval = imax * jmax
C     initialize mask arrays
      maskz = 0
      maskt = 0
      masku = 0
      maskv = 0
      maskq = 0
      mskzs = 0
      mskts = 0
      mskus = 0
      mskvs = 0
      mskqs = 0
      mskps = 0
      mskpm = 0
      mskcp = 0
      mskcn = 0
      mskli = 0 
      mskpbl = 0
      mskpw = 0
      msktrp = 0
      mskcdtz = 0
      mskpri=0
      mskpmx=0
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
C     SURFACE PRESSURE (MB)

c     to start each new file with its index, set J=-1 for sfc pressure
c     j = -1
c     jpds = -1
c     jpds(3) = igdnum
c     jpds(5) = 1
c     jpds(6) = 1
c     jpds(13) = 1
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,
c    +            iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=0
      jpdt(10)=1
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, maxpts
c         ps(ii,jj) = grid(kk) * 0.01
          ps(ii,jj) = gfld%fld(kk)*0.01
          IF (mask(kk)) THEN
            mskps(ii,jj) = 1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 20
          END IF
        END DO
   20   CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
 1300   FORMAT (' IV, IVAR, L, IRET:  ',4I5)
      END IF

C     SEA-LEVEL PRESSURE (MB)

c     j = 0
c     jpds = -1
c     jpds(3) = igdnum
c     SHUELL MSL PRESSURE IS 2 (the default)
c     jpds(5) = 2
c     THE RUC'S MSL PRESSURE IS 129
c     IF (modgen.eq.86.or.modgen.eq.105) jpds(5) = 129
c     FEDOR'S MSL PRESSURE IS 130
c     IF (modgen.eq.83.or.modgen.eq.84.or.modgen.eq.85.or.modgen.eq.89
c    +            .or.modgen.eq.110) jpds(5) = 130
c     jpds(6) = 102
c     jpds(13) = 1
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,
c    +            iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=1
      jpdt(10)=101
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'pressure iret=',iret
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, maxpts
c         gridd=dble(grid(kk))
c         pmd=gridd*D001
c         pm(ii,jj) = grid(kk) * 0.01
          pm(ii,jj) = gfld%fld(kk) * 0.01
c         pm(ii,jj) = real(pmd)
c         print*,'ii,jj,pm=',ii,jj,pm(ii,jj),gridd,pmd,D001
c         print*,'ii,jj,pm=',ii,jj,pm(ii,jj)
          IF (mask(kk)) THEN
            mskpm(ii,jj) = 1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 30
          END IF
        END DO
   30   CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF

C     SURFACE GEOPOTENTIAL HEIGHT (M)

c     j = 0
c     jpds = -1
c     jpds(3) = igdnum
c     jpds(5) = 7
c     jpds(6) = 1
c     jpds(13) = 1
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,
c    +            iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=5
      jpdt(10)=1
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'sfc ht iret=',iret
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, maxpts
c         zs(ii,jj) = grid(kk)
          zs(ii,jj) = gfld%fld(kk)
          IF (mask(kk)) THEN
            mskzs(ii,jj) = 1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 40
          END IF
        END DO
   40   CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF

C     SURFACE TEMPERATURE (K) (ACTUALLY TEMP AT 2 METERS)

c     j = 0
c     jpds = -1
c     jpds(3) = igdnum
c     jpds(5) = 11
c     jpds(6) = 105
c     jpds(7) = 2
c     if(ismart.eq.1) then
c       print*,'ismart eq 1'
c       jpds(6) = 1
c       jpds(7) = 0
c     endif
c     jpds(13) = 1
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,
c    +            iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=0
      jpdt(2)=0
      jpdt(10)=103
      jpdt(12)=2
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'ts iret=',iret
c     print*,'ts kpds=',kpds
c     print*,'ts kgds=',kgds
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, maxpts
c         ts(ii,jj) = grid(kk)
          ts(ii,jj) = gfld%fld(kk)
c         print*,'ii,jj,ts=',ii,jj,ts(ii,jj),grid(kk)
          IF (mask(kk)) THEN
            mskts(ii,jj) = 1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 50
          END IF
        END DO
   50   CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF

C     SURFACE RELATIVE OR SPECIFIC HUMIDITY (ACTUALLY AT 2 METERS)

c     j = 0
c     jpds = -1
c     jpds(3) = igdnum
c     jpds(5) = 51
C     ETA OUTPUTS SPECIFIC HUMIDITY
c     IF (modgen.eq.83.or.modgen.eq.84.or.modgen.eq.85.or.modgen.eq.89
c    +            .or.modgen.eq.110) jpds(5) = 51
c     jpds(6) = 105
c     jpds(7) = 2
c     if(ismart.eq.1) then
c       print*,'ismart=1'
c       jpds(6) = 1
c       jpds(7) = 0
c     endif
c     jpds(13) = 1
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,
c    +            iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=1
      jpdt(2)=0
      jpdt(10)=103
      jpdt(12)=2
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'Q iret=',iret
C
C  1/31/05 : look for RH in Eta output if specific humidity is not present
C
      IF (iret.ne.0) THEN
c      jpds(5) = 52
c      jpds(6) = 105
c      jpds(7) = 2
c      jpds(13) = 1
       jpdt(2)=1
c      CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +            grid,iret)
       call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
       maxpts=gfld%ndpts
       print*,'RH iret=',iret
      ENDIF
c
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, maxpts
c         IF (jpds(5).eq.52) THEN
          IF (jpdt(2).eq.1) THEN
            IF (mskps(ii,jj).eq.1.and.mskts(ii,jj).eq.1) THEN
c             dgrid(kk)=grid(kk)
              dgrid(kk)=gfld%fld(kk)
c             rh = max(0.001,grid(kk)*.01)
c             rh = min(1.00,rh)
              rh = dmax1(D0001,dgrid(kk)*D001)
              rh = dmin1(D1,rh)
c             print*,'ii,jj,rh,ts=',ii,jj,rh,ts(ii,jj)
c             vaps = w3fa09(ts(ii,jj)) * 10.
              vaps = w3fa09a(ts(ii,jj)) * D10
              print*,'vaps=',vaps
              vap = rh * vaps
c             qs(ii,jj) = .622 * vap / (ps(ii,jj)-.378*vap)
              qs(ii,jj) = D622 * vap / (ps(ii,jj)- D378*vap)
c             print*,'ii,jj,qs(ii,jj),rh,dgrid,grid=',ii,jj,qs(ii,jj),
c    *          rh,dgrid(kk),grid(kk),ts(ii,jj)
            ELSE
              qs(ii,jj) = 0.0
              mask(kk) = .false.
            END IF
          ELSE
c           qs(ii,jj) = grid(kk)
            qs(ii,jj) = gfld%fld(kk)
c           print*,'ii,jj,qs(ii,jj)=',ii,jj,qs(ii,jj)
          END IF
          IF (mask(kk)) THEN
            mskqs(ii,jj) = 1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 60
          END IF
        END DO
   60   CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF

C     SURFACE WIND COMPONENTS (M/S) (ACTUALLY AT 10 METERS)

c     j = 0
c     jpds = -1
c     jpds(3) = igdnum
c     jpds(5) = 33
c     jpds(6) = 105
c     jpds(7) = 10
c     if(ismart.eq.1) then
c       print*,'ismart=1'
c       jpds(6) = 1
c       jpds(7) = 0
c     endif
c     jpds(13) = 1
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,
c    +            iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=2
      jpdt(2)=2
      jpdt(10)=103
      jpdt(12)=10
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'u iret=',iret
c
c Check to see if the winds are grid-relative or earth-relative
c
      if(btest(kgds(6),3) ) then
       iearth=0  !  grid relative
      else
       iearth=1  !  earth relative
      endif
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, maxpts
c         us(ii,jj) = grid(kk)
          us(ii,jj) = gfld%fld(kk)
          IF (mask(kk)) THEN
            mskus(ii,jj) = 1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 70
          END IF
        END DO
   70   CONTINUE
C       
        j = 0
c       jpds = kpds
c       jpds(5) = 34
        jpdt(2) = 3
c       CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
        IF (iret.eq.0) THEN
          ii = 1
          jj = jj1
          DO kk = 1, maxpts
c           vs(ii,jj) = grid(kk)
            vs(ii,jj) = gfld%fld(kk)
            IF (mask(kk)) THEN
              mskvs(ii,jj) = 1
            END IF
            iprev = ii
            ii = ii + 1
            IF (mod(iprev,imax).eq.0) THEN
              ii = 1
              jj = jj + jjinc
              IF (jj.gt.jmax.or.jj.lt.1) GO TO 80
            END IF
          END DO
   80     CONTINUE
        ELSE
          WRITE (6,1300) iv, jpds(5), l, iret
        END IF
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF

C  SURFACE-BASED CAPE
c     j=0
c     jpds=-1
c     jpds(5)=157
c     jpds(6)=1
c     jpds(7)=0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=7
      jpdt(2)=6
      jpdt(10)=1
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
c     jpds(6)=116
c     jpds(7)=46080
c           CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid2,iret)
      print*,'cape iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,maxpts
c         cape(ii,jj)=grid(kk)
          cape(ii,jj)=gfld%fld(kk)
          IF (mask(kk)) THEN
            mskcp(ii,jj)=1
          ENDIF
          iprev=ii
          ii=ii+1
          IF (mod(iprev,imax).eq.0) THEN
            ii=1
            jj=jj+jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 152
           ENDIF
         ENDDO
152     CONTINUE
        ELSE
         WRITE(6,1300) iv,jpds(5),l,iret
        ENDIF

C  BEST CAPE
c     j=0
c     jpds=-1
c     jpds(5)=157
c     jpds(6)=116
c     jpds(7)=46080
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=7
      jpdt(2)=6
      jpdt(10)=108
      jpdt(12)=18000
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'best cape iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,maxpts
c         bcape(ii,jj)=grid(kk)
          bcape(ii,jj)=gfld%fld(kk)
c         print*,'ii,jj,bcape(ii,jj)=',ii,jj,bcape(ii,jj),mask(kk)
          IF (mask(kk)) THEN
            mskbc(ii,jj)=1
          ENDIF
          iprev=ii
          ii=ii+1
          IF (mod(iprev,imax).eq.0) THEN
            ii=1
            jj=jj+jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 192
           ENDIF
         ENDDO
192     CONTINUE
        ELSE
         WRITE(6,1300) iv,jpds(5),l,iret
        ENDIF

C  SURFACE-BASED CIN

c     j=0
c     jpds=-1
c     jpds(5)=156
c     jpds(6)=1
c     jpds(7)=0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=7
      jpdt(2)=7
      jpdt(10)=1
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,maxpts
c         cin(ii,jj)=grid(kk)
          cin(ii,jj)=gfld%fld(kk)
          IF (mask(kk)) THEN
            mskcn(ii,jj)=1
          ENDIF
          iprev=ii
          ii=ii+1
          IF (mod(iprev,imax).eq.0) THEN
            ii=1
            jj=jj+jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 153
           ENDIF
         ENDDO
153     CONTINUE
        ELSE
         WRITE(6,1300) iv,jpds(5),l,iret
        ENDIF

C SURFACE-TO-500-MB LIFTED INDEX

c     j=0
c     jpds=-1
c     jpds(5)=24
c     jpds(6)=116
c     jpds(7)=7680
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=7
      jpdt(2)=10
      jpdt(10)=1
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'pli iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,maxpts
c         pli(ii,jj)=grid(kk)
          pli(ii,jj)=gfld%fld(kk)
          IF (mask(kk)) THEN
            mskli(ii,jj)=1
          ENDIF
          iprev=ii
          ii=ii+1
          IF (mod(iprev,imax).eq.0) THEN
            ii=1
            jj=jj+jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 154
           ENDIF
         ENDDO
154     CONTINUE
        ELSE
         WRITE(6,1300) iv,jpds(5),l,iret
        ENDIF

C MAX TEMPERATURE
          
c     j=0
c     jpds=-1
c     jpds(5)=15
c     jpds(6)=105
c     jpds(7)=2
c     if(ismart.eq.1) then
c       print*,'ismart=1'
c       jpds(6) = 1
c       jpds(7) = 0
c     endif
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=0
      jpdt(2)=4
      jpdt(10)=1
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'tmax iret=',iret
      print*,'tmax kpds=',kpds
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,maxpts
c        tmax(ii,jj)=grid(kk)
         tmax(ii,jj)=gfld%fld(kk)
         IF (mask(kk)) THEN
            msktmax(ii,jj)=1
         ENDIF
          iprev=ii
          ii=ii+1
          IF (mod(iprev,imax).eq.0) THEN
            ii=1
            jj=jj+jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 155
           ENDIF
         ENDDO
155     CONTINUE
        ELSE
         WRITE(6,1300) iv,jpds(5),l,iret
        ENDIF

C MIN TEMPERATURE
      
c     j=0
c     jpds=-1
c     jpds(5)=16
c     jpds(6)=105
c     jpds(7)=2
c     if(ismart.eq.1) then
c       print*,'ismart=1'
c       jpds(6) = 1
c       jpds(7) = 0
c     endif
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=0
      jpdt(2)=5
      jpdt(10)=1
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'tmin iret=',iret
      print*,'tmin kpds=',kpds
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,maxpts
c        tmin(ii,jj)=grid(kk)
         tmin(ii,jj)=gfld%fld(kk)
         IF (mask(kk)) THEN
            msktmin(ii,jj)=1
         ENDIF
          iprev=ii
          ii=ii+1
          IF (mod(iprev,imax).eq.0) THEN
            ii=1
            jj=jj+jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 158
           ENDIF
         ENDDO
158     CONTINUE
        ELSE
         WRITE(6,1300) iv,jpds(5),l,iret
        ENDIF

C VISIBILITY
      
c     j=0
c     jpds=-1
c     jpds(3)=igdnum
c     jpds(5)=20
c     jpds(6)=1
c     jpds(7)=0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=19
      jpdt(2)=0
      jpdt(10)=1
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'vis iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,maxpts
c        vis(ii,jj)=grid(kk)
         vis(ii,jj)=gfld%fld(kk)
         if(vis(ii,jj).gt.16090.0.and.vis(ii,jj).lt.bmiss) then
           vis(ii,jj)=16090.0
         endif
         IF (mask(kk)) THEN
            mskvis(ii,jj)=1
         ENDIF
          iprev=ii
          ii=ii+1
          IF (mod(iprev,imax).eq.0) THEN
            ii=1
            jj=jj+jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 156
           ENDIF
         ENDDO
156     CONTINUE
        ELSE
         WRITE(6,1300) iv,jpds(5),l,iret
        ENDIF

C DEW POINT
   
c     j=0
c     jpds=-1
c     jpds(5)=17
c     jpds(6)=105
c     jpds(7)=2
c     if(ismart.eq.1) then
c       print*,'ismart=1'
c       jpds(6) = 1
c       jpds(7) = 0
c     endif
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=0
      jpdt(2)=6
      jpdt(10)=103
      jpdt(12)=2
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'dpt iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,maxpts
c        dpt(ii,jj)=grid(kk)-tzero
         dpt(ii,jj)=gfld%fld(kk)-tzero
         IF (mask(kk)) THEN
            mskdpt(ii,jj)=1
         ENDIF
          iprev=ii
          ii=ii+1
          IF (mod(iprev,imax).eq.0) THEN
            ii=1
            jj=jj+jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 157
           ENDIF
         ENDDO
157     CONTINUE
        ELSE
         WRITE(6,1300) iv,jpds(5),l,iret
        ENDIF

C TOTAL CLOUD COVER
      
c     j=0
c     jpds=-1
c     jpds(5)=71
c     jpds(6)=200
c     jpds(7)=0
c     if(ismart.eq.1) then
c       print*,'ismart=1'
c       jpds(6) = 1
c       jpds(7) = 0
c     endif
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=6
      jpdt(2)=1
      jpdt(10)=200
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'tocc iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,maxpts
c        tocc(ii,jj)=grid(kk)
         tocc(ii,jj)=gfld%fld(kk)
         IF (mask(kk)) THEN
            msktocc(ii,jj)=1
         ENDIF
          iprev=ii
          ii=ii+1
          IF (mod(iprev,imax).eq.0) THEN
            ii=1
            jj=jj+jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 159
           ENDIF
         ENDDO
159     CONTINUE
        ELSE
         WRITE(6,1300) iv,jpds(5),l,iret
        ENDIF

C  CLOUD BASE HEIGHT
      
c     j=0
c     jpds=-1
c     jpds(3)=221
c     jpds(5)=7
c     jpds(6)=2
c     jpds(7)=0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=6
      jpdt(2)=11
c     jpdt(10)=10
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'cldb iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,maxpts
c         cdtz(ii,jj)=grid(kk)
          cdtz(ii,jj)=gfld%fld(kk)
          if(grid(kk).lt.0) cdtz(ii,jj)=bmiss
          IF (mask(kk)) THEN
            mskcdtz(ii,jj)=1
          ENDIF
          iprev=ii
          ii=ii+1
          IF (mod(iprev,imax).eq.0) THEN
            ii=1
            jj=jj+jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 169
           ENDIF
         ENDDO
169     CONTINUE
        ELSE
         WRITE(6,1300) iv,jpds(5),l,iret
        ENDIF

C SURFACE WIND GUST

c     j=0
c     jpds=-1
c     jpds(5)=180
c     jpds(6)=1
c     jpds(7)=0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=2
      jpdt(2)=22
      jpdt(10)=1
c     jpdt(12)=10
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'gust iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,maxpts
c        gust(ii,jj)=grid(kk)
         gust(ii,jj)=gfld%fld(kk)
         IF (mask(kk)) THEN
            mskgust(ii,jj)=1
         ENDIF
          iprev=ii
          ii=ii+1
          IF (mod(iprev,imax).eq.0) THEN
            ii=1
            jj=jj+jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 160
           ENDIF
         ENDDO
160     CONTINUE
        ELSE
         WRITE(6,1300) iv,jpds(5),l,iret
        endif

C PRECIPITABLE WATER

c     j=0
c     jpds=-1
c     jpds(5)=54
c     jpds(6)=200
c     jpds(7)=0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=1
      jpdt(2)=3
      jpdt(10)=200
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'pw iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,maxpts
c        pwo(ii,jj)=grid(kk)
         pwo(ii,jj)=gfld%fld(kk)
         IF (mask(kk)) THEN
            mskpw(ii,jj)=1
         ENDIF
          iprev=ii
          ii=ii+1
          IF (mod(iprev,imax).eq.0) THEN
            ii=1
            jj=jj+jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 161
           ENDIF
         ENDDO
 161     CONTINUE
        ELSE
         WRITE(6,1300) iv,jpds(5),l,iret
        ENDIF

C  TROPOPAUSE LEVEL
                                                                                            
c     j=0
c     jpds=-1
c     jpds(5)=1
c     jpds(6)=7
c     jpds(7)=0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=0
      jpdt(10)=7
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'trop iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,maxpts
c         trop(ii,jj)=grid(kk)/100.0
          trop(ii,jj)=gfld%fld(kk)/100.0
          IF (mask(kk)) THEN
            msktrp(ii,jj)=1
          ENDIF
          iprev=ii
          ii=ii+1
          IF (mod(iprev,imax).eq.0) THEN
            ii=1
            jj=jj+jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 162
           ENDIF
         ENDDO
162     CONTINUE
        ELSE
         WRITE(6,1300) iv,jpds(5),l,iret
        ENDIF

C Categorical snow
      
c     grid=0
c     grid2=0 
c     grid3=0
c     grid4=0
c     j=0
c     jpds=-1
c     jpds(5)=143
c     jpds(5)=65
c     jpds(6)=1
c     jpds(7)=0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
c     print*,'snow iret=',iret
c     print*,'snow kpds=',kpds
c     print*,'kpds(14),kpds(15)=',kpds(14),kpds(15)
c     ipcp=kpds(15)-kpds(14)
c     print*,'ipcp=',ipcp
c     jpds=-1
c     jpds(5)=142
c     jpds(6)=1
c     jpds(7)=0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid2,iret)
c     jpds=-1
c     jpds(5)=141
c     jpds(6)=1
c     jpds(7)=0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid3,iret)
c     jpds=-1
c     jpds(5)=140
c     jpds(6)=1
c     jpds(7)=0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid4,iret)
c     IF (iret.eq.0) THEN
c       ii=1
c       jj=jj1
c       DO kk=1,itot
c        csnow(ii,jj)=grid(kk)
c        if(csnow(ii,jj).gt.0.0) then
c        if(grid(kk).gt.0.0) then
c        print*,'jj,ii,csnow(ii,jj),snow=',ii,jj,csnow(ii,jj),grid2(kk)
c    *     ,grid3(kk),grid4(kk)
c        endif
c        if(grid(kk).gt.50) then
c          csnow(ii,jj)=100.0
c        else
c          csnow(ii,jj)=0.0
c        endif
c        IF (mask(kk)) THEN
c           msksnow(ii,jj)=1
c        ENDIF
c         iprev=ii
c         ii=ii+1
c         IF (mod(iprev,imax).eq.0) THEN
c           ii=1
c           jj=jj+jjinc
c           IF (jj.gt.jmax.or.jj.lt.1) GO TO 162
c          ENDIF
c        ENDDO
c162     CONTINUE
c       ELSE
c        WRITE(6,1300) iv,jpds(5),l,iret
c       ENDIF

c     jpds=-1
c     jpds(5)=65
c     jpds(6)=1
c     jpds(7)=0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
c     IF (iret.eq.0) THEN
c       ii=1
c       jj=jj1
c       DO kk=1,itot
c        csnow(ii,jj)=grid(kk)
c        if(csnow(ii,jj).gt.0.0) then
c        print*,'jj,ii,csnow(ii,jj)=',ii,jj,csnow(ii,jj)
c        endif
c         iprev=ii
c         ii=ii+1
c         IF (mod(iprev,imax).eq.0) THEN
c           ii=1
c           jj=jj+jjinc
c         endif
c       enddo
c     endif

C     TKE PBL HEIGHT (m)
                                                                                
      print*,'PBL HEIGHT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
c     j = 0
c     jpds = -1
c     jpds(3) = igdnum
c     jpds(5) = 221
c     jpds(6) = 1
c     jpds(7) = 0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,
c    +            iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=196
      jpdt(10)=1
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'pbl iret=',iret
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, maxpts
c         pbls(ii,jj) = grid(kk)
          pbls(ii,jj)=gfld%fld(kk)
          IF (mask(kk)) THEN
            mskpbl(ii,jj) = 1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 163
          END IF
        END DO
163     CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF

C     RI PBL HEIGHT (m)

c     j = 0
c     jpds = -1
c     jpds(3) = igdnum
c     jpds(5) = 7
c     jpds(6) = 220
c     jpds(7) = 0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,
c    +            iret)
c     print*,'ripbl iret=',iret
c     IF (iret.eq.0) THEN
c       ii = 1
c       jj = jj1
c       DO kk = 1, itot
c         pblris(ii,jj) = grid(kk)
c         IF (mask(kk)) THEN
c           mskpri(ii,jj) = 1
c         END IF
c         iprev = ii
c         ii = ii + 1
c         IF (mod(iprev,imax).eq.0) THEN
c           ii = 1
c           jj = jj + jjinc
c           IF (jj.gt.jmax.or.jj.lt.1) GO TO 164
c         END IF
c       END DO
c64     CONTINUE
c     ELSE
c       WRITE (6,1300) iv, jpds(5), l, iret
c     END IF

C     MIX LAYER HEIGHT (m)

c     j = 0
c     jpds = -1
c     jpds(3) = igdnum
c     jpds(5) = 67
c     jpds(6) = 1
c     jpds(7) = 0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,
c    +            iret)
      j = 0
      jids=-9999
      jpdtn=0
      jpdt=-9999
      jpdt(1)=19
      jpdt(2)=3
      jpdt(10)=1
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'mixht iret=',iret
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, maxpts
c         pblmxs(ii,jj) = grid(kk)
          pblmxs(ii,jj)=gfld%fld(kk)
          IF (mask(kk)) THEN
            mskpmx(ii,jj) = 1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 165
          END IF
        END DO
165     CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF

C     AEROSOL OPTICAL DEPTH - COLUMN

      j = 0
      jids=-9999
      jpdtn=48
      jpdt=-9999
      jpdt(1)=20
      jpdt(2)=102
      jpdt(3)=62000
      jpdt(20)=fhour
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'aod iret=',iret
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, maxpts
c         aod(ii,jj) = grid(kk)
          aod(ii,jj) = gfld%fld(kk)
c         print*,'ii,jj,aod(ii,jj)=',ii,jj,aod(ii,jj)
c         if(ifz.eq.1) aod(ii,jj)=aod(ii,jj)/100.
c         print*,'ii,jj,aod(ii,jj)=',ii,jj,aod(ii,jj)
c         IF (mask(kk)) THEN
            mskaod(ii,jj) = 1
c         END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 179
          END IF
        END DO
179     CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF

C     -== GET PRESSURE LEVEL VARIABLES Z, T, U, V & Q/RH ==-

      DO 100 iv = 1, 5
        j = 0
c       jpds = -1
c       jpds(3) = igdnum
c       jpds(5) = ivar(iv)
c       ETA OUTPUTS SPECIFIC HUMIDITY
c       IF ((modgen.eq.83.or.modgen.eq.84.or.modgen.eq.85.or.modgen.eq.
c    +              89.or.modgen.eq.110).and.iv.eq.5) jpds(5) = ivar(6)
c       if (iv.eq.5) jpds(5)=52
c       if (iv.eq.5) jpds(5)=51  ! specific humidity
c       jpds(6) = 100
c       jpds(13) = 1
c     j = 0
c     jids=-9999
c     jpdtn=0
c     jpdt=-9999
c     jpdt(1)=idisp(iv)
c     jpdt(2)=ivar(iv)
c     jpdt(10)=100
      call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
        DO 90 l = 1, numlev
c         jpds(7) = levs(l,levset)
          j=0
          jids=-9999
          jpdtn=0
          jpdt=-9999
          jpdt(1)=idisp(iv)
          jpdt(2)=ivar(iv)
          jpdt(10)=100
          jpdt(12) = levs(l,levset) * 100
          print*,'jpdt(12)=',jpdt(12)
          call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
          print*,'ivar(iv),idisp(iv),iret=',ivar(iv),idisp(iv),iret
c         CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +                grid,iret)
c
c Check to see if the winds are grid-relative or earth-relative
c
      if(btest(kgds(6),3) ) then
       iearth=0  !  grid relative
      else
       iearth=1  !  earth relative
      endif
C
C  1/31/05 : look for RH in Eta output if specific humidity is not present
C
      IF (iv.eq.5.and.iret.ne.0) THEN
c      jpds(5) = 52
c      jpds(6) = 100
c      jpds(7) = levs(l,levset)
c      jpds(13) = 1
       jpdt(2)=1
c      CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +            grid,iret)
       call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
       print*,'no q RH iret=',iret
      ENDIF
          IF (iret.eq.0) THEN
            ii = 1
            jj = jj1
            DO kk = 1, maxpts
              IF (iv.eq.1) THEN
c               z(ii,jj,l) = grid(kk)
                z(ii,jj,l) = gfld%fld(kk)
                IF (mask(kk)) THEN
                  maskz(ii,jj,l) = 1
                END IF
              ELSE IF (iv.eq.2) THEN
c               t(ii,jj,l) = grid(kk)
                t(ii,jj,l) = gfld%fld(kk)
c               if(ii.eq.189.and.jj.eq.80) then
c               print*,'ii,jj,l,t,mask=',
c    *            ii,jj,l,t(ii,jj,l),mask(kk)
c               endif
                IF (mask(kk)) THEN
                  maskt(ii,jj,l) = 1
c                 print*,'maskt is 1 for=',ii,jj,l
                END IF
              ELSE IF (iv.eq.3) THEN
c               u(ii,jj,l) = grid(kk)
                u(ii,jj,l) = gfld%fld(kk)
                IF (mask(kk)) THEN
                  masku(ii,jj,l) = 1
                END IF
              ELSE IF (iv.eq.4) THEN
c               v(ii,jj,l) = grid(kk)
                v(ii,jj,l) = gfld%fld(kk)
                IF (mask(kk)) THEN
                  maskv(ii,jj,l) = 1
                END IF
              ELSE IF (iv.eq.5) THEN
                IF (jpds(5).eq.52) THEN
                  IF (maskt(ii,jj,l).ge.1.) THEN
c                   dgrid(kk)=grid(kk)
                    dgrid(kk)=gfld%fld(kk)
c                   rh = max(0.01,grid(kk)*.01)
c                   rh = min(1.00,rh)
c                   vaps = w3fa09(t(ii,jj,l)) * 10.
                    rh = dmax1(D0001,dgrid(kk)*D001)
                    rh = dmin1(D1,rh)
                    vaps = w3fa09a(t(ii,jj,l)) * D10
                    vap = rh * vaps
c                   q(ii,jj,l) = .622 * vap / (pgd(ii,jj,l)-.378*vap)
                    q(ii,jj,l) = D622 * vap / (pgd(ii,jj,l)- D378*vap)
                  ELSE
                    q(ii,jj,l) = 0.0
                    mask(kk) = .false.
                  END IF
                ELSE
c                 q(ii,jj,l) = grid(kk)
                  q(ii,jj,l) = gfld%fld(kk)
                END IF
                IF (mask(kk)) THEN
                  maskq(ii,jj,l) = 1
                  alnq(ii,jj,l) = log(q(ii,jj,l))
                END IF
              END IF
              iprev = ii
              ii = ii + 1
              IF (mod(iprev,imax).eq.0) THEN
                ii = 1
                jj = jj + jjinc
                IF (jj.gt.jmax.or.jj.lt.1) GO TO 90
              END IF
            END DO
          ELSE
            WRITE (6,1300) iv, jpds(5), l, iret
          END IF
   90     CONTINUE
  100     CONTINUE

C  CLOUD BASE HEIGHT
                                                                                             
c     j=0
c     jpds=-1
c     jpds(5)=7
c     jpds(6)=2
c     jpds(7)=0
c     CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,
c    +              grid,iret)
c     IF (iret.eq.0) THEN
c       ii=1
c       jj=jj1
c       DO kk=1,itot
c         cdtz(ii,jj)=grid(kk)
c         IF (mask(kk)) THEN
c           mskcdtz(ii,jj)=1
c         ENDIF
c         iprev=ii
c         ii=ii+1
c         IF (mod(iprev,imax).eq.0) THEN
c           ii=1
c           jj=jj+jjinc
c           IF (jj.gt.jmax.or.jj.lt.1) GO TO 169
c          ENDIF
c        ENDDO
c69     CONTINUE
c       ELSE
c        WRITE(6,1300) iv,jpds(5),l,iret
c       ENDIF

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

      deallocate(grid)
      deallocate(grid2)
      deallocate(grid3)
      deallocate(grid4)

      deallocate(mask)

      call gf_free(gfld)

          RETURN
          END
