      SUBROUTINE getbak(lugb,lugi,numlev,iearth,
     *    ismart,ipcp,istat,fhour,src)
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
      real*8 d90
      real*8 D001,D0001,D00001
      real*8 D622
      real*8 D378
      real*8 D10
      real*8 D1
      real*8 vap, vaps
      real*8 rh
      real*8 dgrid(itot)

      DIMENSION RICDAT(5),IFDATE(8),IVDATE(8)
      INTEGER jpds(200), jgds(200), kpds(200), kgds(200)
c      INTEGER levs(40,6), ivar(6)
      INTEGER levs(40,10), ivar(6), idisp(6)
      integer fhour
      CHARACTER*8 modesc(110)
      character*10 src

C GRIB2
      type(gribfield) :: gfld
      integer ifile,j,jdisc,jpdtn,jgdtn,iret
      integer,dimension(300) :: jids,jpdt,jgdt
      real fldmin,fldmax,firstval,lastval
      logical :: unpack=.true.
c     character(255) :: cin
C GRIB2

      real,allocatable:: utrans(:,:),vtrans(:,:)
      real,allocatable:: u80(:,:),v80(:,:)
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
      if(src(:6).eq.'SRMEAN'.or.src(:5).eq.'NARRE') then
       jpdtn=2
      else if(src(:2).eq.'SR'.and.src(:6).ne.'SRMEAN') then
       jpdtn=1
      else
       jpdtn=0
      endif
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
      print*,'lugb,lugi=',lugb,lugi
C     
!for pdt, define catogory, parameter and level
!eg: tmp at 2 m
      print*,'fhour=',fhour
      j=0
      jpdt(1)=0      ! table 4.1
      jpdt(2)=0      ! table 4.2-0-0
      jpdt(10)=103   ! table 4.5
      jpdt(12)=2
      jpdt(20)=fhour
c     jpdt(20)=-9999
!-- set grid def template
      jgdtn=-1
!-- set product def array
      jgdt=-9999
        if(ismart.eq.1) then
          print*,'ismart'
          jpdt(10)=1
          jpdt(12)=-9999
        endif
        call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *         unpack,j,gfld,iret)
c       print*,'iret check=',iret
       if(iret.ne.0) then
          if(src(:5).eq.'NARRE') then
           print*,'setting for narre'
           jpdt(1)=2
           jpdt(2)=2
           jpdt(10)=103
           jpdt(12)=10
          endif
c         do i=1,20
c          print*,'i,jpdt(i)=',i,jpdt(i)
c         enddo
          call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *         unpack,j,gfld,iret)
          if(iret.ne.0) then
           print*,'iret=',iret
           print*,'No successful call to getgb2'
           stop
          endif
        endif
        ifdate(1)=gfld%idsect(6)
        ifdate(2)=gfld%idsect(7)
        ifdate(3)=gfld%idsect(8)
        ifdate(5)=gfld%idsect(9)
        ricdat(2)=gfld%ipdtmpl(9)
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

      allocate(z(imax,jmax,maxlev))
      allocate(t(imax,jmax,maxlev))
      allocate(u(imax,jmax,maxlev))
      allocate(v(imax,jmax,maxlev))
      allocate(q(imax,jmax,maxlev))
      allocate(alnq(imax,jmax,maxlev))

      allocate(pgd(imax,jmax,maxlev))
      allocate(maskz(imax,jmax,maxlev))
      allocate(maskt(imax,jmax,maxlev))
      allocate(masku(imax,jmax,maxlev))
      allocate(maskv(imax,jmax,maxlev))
      allocate(maskq(imax,jmax,maxlev))

      allocate(ps(imax,jmax))
      allocate(zs(imax,jmax))
      allocate(ts(imax,jmax))
      allocate(qs(imax,jmax))
      allocate(us(imax,jmax))
      allocate(vs(imax,jmax))
      allocate(pm(imax,jmax))

      allocate(cape(imax,jmax))
      allocate(cin(imax,jmax))
      allocate(bcape(imax,jmax))
      allocate(pli(imax,jmax))

      allocate(vis(imax,jmax))
      allocate(tmax(imax,jmax))
      allocate(tmin(imax,jmax))
      allocate(dpt(imax,jmax))
      allocate(tocc(imax,jmax))
      allocate(gust(imax,jmax))
      allocate(pbls(imax,jmax))

      allocate(pwo(imax,jmax))
      allocate(trop(imax,jmax))
      allocate(cdtz(imax,jmax))
      allocate(pblris(imax,jmax))
      allocate(pblmxs(imax,jmax))
      allocate(haines(imax,jmax))
      allocate(trans(imax,jmax))
      allocate(wnd80(imax,jmax))
      allocate(ceil(imax,jmax))

      allocate(mskps(imax,jmax))
      allocate(mskzs(imax,jmax))
      allocate(mskts(imax,jmax))
      allocate(mskqs(imax,jmax))
      allocate(mskus(imax,jmax))
      allocate(mskvs(imax,jmax))
      allocate(mskpm(imax,jmax))

      allocate(mskcp(imax,jmax))
      allocate(mskcn(imax,jmax))
      allocate(mskbc(imax,jmax))
      allocate(mskli(imax,jmax))

      allocate(mskvis(imax,jmax))
      allocate(msktmax(imax,jmax))
      allocate(msktmin(imax,jmax))
      allocate(mskdpt(imax,jmax))
      allocate(msktocc(imax,jmax))
      allocate(mskgust(imax,jmax))
      allocate(mskpbl(imax,jmax))

      allocate(mskpw(imax,jmax))
      allocate(msktrp(imax,jmax))
      allocate(mskcdtz(imax,jmax))
      allocate(mskpri(imax,jmax))
      allocate(mskpmx(imax,jmax))
      allocate(mskhaines(imax,jmax))
      allocate(msktrans(imax,jmax))
      allocate(mskwnd80(imax,jmax))
      allocate(mskceil(imax,jmax))

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
      vdate = (ivdate(1)*1000000+ivdate(2)*10000+ivdate(3)*100+
     1          ivdate(5))
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
      mskbc = 0
      mskvis = 0
      mskdpt = 0
      msktocc = 0
      mskgust = 0
      mskhaines = 0
      msktrans = 0
      mskwnd80 = 0
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

      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=0
      jpdt(10)=1
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'iret,maxpts=',iret,maxpts
      print*,'gfld%ibmap=',gfld%ibmap
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, itot
          ps(ii,jj) = gfld%fld(kk)*0.01
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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

      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=1
      jpdt(10)=101
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      if (iret.ne.0) then
       jpdt(2)=198
       call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      endif
      maxpts=gfld%ndpts
      print*,'pressure iret=',iret
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, itot
          pm(ii,jj) = gfld%fld(kk) * 0.01
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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

      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=5
      jpdt(10)=1
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'sfc ht iret=',iret
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, itot
          zs(ii,jj) = gfld%fld(kk)
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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

      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=0
      jpdt(2)=0
      jpdt(10)=103
      jpdt(12)=2
        if(ismart.eq.1) then
          jpdt(10)=1
          jpdt(12)=-9999
        endif

      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'ts iret=',iret
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, itot
          ts(ii,jj) = gfld%fld(kk)
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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

      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=1
      jpdt(2)=0
      jpdt(10)=103
      jpdt(12)=2
        if(ismart.eq.1) then
          jpdt(10)=1
          jpdt(12)=-9999
        endif

      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'Q iret=',iret
C
C  1/31/05 : look for RH in Eta output if specific humidity is not present
C
      IF (iret.ne.0) THEN
       jpdt(2)=1
       call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
       maxpts=gfld%ndpts
       print*,'RH iret=',iret
      ENDIF
c
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, itot
          IF (jpdt(2).eq.1) THEN
            IF (mskps(ii,jj).eq.1.and.mskts(ii,jj).eq.1) THEN
              rh = max(0.001,gfld%fld(kk)*.01)
              rh = min(1.00,rh)
              vaps = w3fa09(ts(ii,jj)) * 10.
              vap = rh * vaps
              qs(ii,jj) = .622 * vap / (ps(ii,jj)-.378*vap)
            ELSE
              qs(ii,jj) = 0.0
c             gfld%bmap(kk) = .false.
              gfld%ibmap = 255
            END IF
          ELSE
            qs(ii,jj) = gfld%fld(kk)
          END IF
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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

      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=2
      jpdt(2)=2
      jpdt(10)=103
      jpdt(12)=10
        if(ismart.eq.1) then
          jpdt(10)=1
          jpdt(12)=-9999
        endif

      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'u iret=',iret
c Check to see if the winds are grid-relative or earth-relative
c
      if(btest(gfld%igdtmpl(12),3) ) then
       iearth=0  !  grid relative
      else
       iearth=1  !  earth relative
      endif
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, itot
          us(ii,jj) = gfld%fld(kk)
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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
        jpdt(2) = 3
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'v iret=',iret
        IF (iret.eq.0) THEN
          ii = 1
          jj = jj1
          DO kk = 1, itot
            vs(ii,jj) = gfld%fld(kk)
            IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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
      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=7
      jpdt(2)=6
      jpdt(10)=1
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'cape iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,itot
          cape(ii,jj)=gfld%fld(kk)
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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
      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=7
      jpdt(2)=6
      jpdt(10)=108
      jpdt(12)=18000
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'best cape iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,itot
          bcape(ii,jj)=gfld%fld(kk)
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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

      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=7
      jpdt(2)=7
      jpdt(10)=1
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'cin iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,itot
          cin(ii,jj)=gfld%fld(kk)
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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

      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=7
      jpdt(2)=0
      jpdt(10)=108
      jpdt(12)=3000
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'pli iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,itot
          pli(ii,jj)=gfld%fld(kk)
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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
          
      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=0
      jpdt(2)=4
      jpdt(10)=1
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'tmax iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,itot
         tmax(ii,jj)=gfld%fld(kk)
         IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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
      
      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=0
      jpdt(2)=5
      jpdt(10)=1
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'tmin iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,itot
         tmin(ii,jj)=gfld%fld(kk)
         IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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
      
      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=19
      jpdt(2)=0
      jpdt(10)=1
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'vis iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,itot
         vis(ii,jj)=gfld%fld(kk)
         if(vis(ii,jj).gt.16090.0.and.vis(ii,jj).lt.bmiss) then
           vis(ii,jj)=16090.0
         endif
         IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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
   
      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=0
      jpdt(2)=6
      jpdt(10)=103
      jpdt(12)=2
        if(ismart.eq.1) then
          jpdt(10)=1
          jpdt(12)=-9999
        endif
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'dpt iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,itot
         dpt(ii,jj)=gfld%fld(kk)-tzero
         IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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
      
      j = 0
      jids=-9999
      jpdtnc=jpdtn
      jpdt=-9999
      jpdt(1)=6
      jpdt(2)=1
      jpdt(10)=200
        if(ismart.eq.1) then
          jpdt(10)=1
          jpdt(12)=-9999
        endif
      call getgb2(lugb,0,j,jdisc,jids,jpdtnc,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      if(iret.ne.0) then
        jpdtnc=8
        call getgb2(lugb,0,j,jdisc,jids,jpdtnc,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      endif
      if(iret.ne.0) then
        jpdt(10)=10
        call getgb2(lugb,0,j,jdisc,jids,jpdtnc,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      endif
      print*,'tocc iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,itot
         tocc(ii,jj)=gfld%fld(kk)
c        print*,'ii,jj,tocc(ii,jj)=',ii,jj,tocc(ii,jj)
         IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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
      
      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=5
      jpdt(10)=2
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'cldb iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,itot
          cdtz(ii,jj)=gfld%fld(kk)
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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

      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=2
      jpdt(2)=22
      jpdt(10)=1
c     jpdt(12)=10
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'gust iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,itot
         gust(ii,jj)=gfld%fld(kk)
         IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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

      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=1
      jpdt(2)=3
      jpdt(10)=200
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'pw iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,itot
         pwo(ii,jj)=gfld%fld(kk)
         IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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
                                                                                            
      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=0
      jpdt(10)=7
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'trop iret=',iret
      IF (iret.eq.0) THEN
        ii=1
        jj=jj1
        DO kk=1,itot
          trop(ii,jj)=gfld%fld(kk)/100.0
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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

C     TKE PBL HEIGHT (m)
                                                                                
c     print*,'PBL HEIGHT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=196
      jpdt(10)=1
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'pbl iret=',iret
      if (iret.ne.0) then
      jids=-9999
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=18
      jpdt(10)=1
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      print*,'second pbl iret=',iret
      endif
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, itot
          pbls(ii,jj)=gfld%fld(kk)
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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

      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=5
      jpdt(10)=220
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'ri pbl iret=',iret
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, itot
          pblris(ii,jj)=gfld%fld(kk)
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
            mskpri(ii,jj) = 1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 164
          END IF
        END DO
164     CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF

C     MIX LAYER HEIGHT (m)

      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=19
      jpdt(2)=3
      jpdt(10)=1
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'mixht iret=',iret
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, itot
          pblmxs(ii,jj)=gfld%fld(kk)
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
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

C     HAINES INDEX

      j=0 
      jids=-9999
      jpdt=-9999
      jdisc=2
      jpdt(1)=4
      jpdt(2)=2  
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'haines iret=',iret
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, itot
          haines(ii,jj) = gfld%fld(kk) 
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
            mskhaines(ii,jj) = 1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 166
          END IF
        END DO
166     CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF

C     TRANSPORT WIND - U and V COMPONENTS

      allocate(utrans(imax,jmax))
      allocate(vtrans(imax,jmax))
      j=0
      jids=-9999
      jpdt=-9999
      jdisc=0
      jpdt(1)=2
      jpdt(2)=2
      jpdt(10)=220
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'transport u iret=',iret
      if(iret.eq.0) then
       ii=1
       jj=jj1
       do kk=1,itot
        utrans(ii,jj)=gfld%fld(kk)
        iprev=ii
        ii=ii+1 
        IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 188
         END IF
        END DO
188     CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF 

      jpdt(2)=3
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'transport v iret=',iret
      if(iret.eq.0) then
       ii=1
       jj=jj1
       do kk=1,itot
        vtrans(ii,jj)=gfld%fld(kk)
        trans(ii,jj)=sqrt(utrans(ii,jj)**2+vtrans(ii,jj)**2)
        msktrans(ii,jj)=1 
        iprev=ii
        ii=ii+1
        IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 189
         END IF
        END DO
189     CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF

      deallocate(utrans)
      deallocate(vtrans)

C     80-M WIND - U and V COMPONENTS

      allocate(u80(imax,jmax))
      allocate(v80(imax,jmax))
      j=0
      jids=-9999
      jpdt=-9999
      jdisc=0
      jpdt(1)=2
      jpdt(2)=2
      jpdt(10)=103
      jpdt(12)=80
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'80-m u iret=',iret
      if(iret.eq.0) then
       ii=1
       jj=jj1
       do kk=1,itot
        u80(ii,jj)=gfld%fld(kk)
        iprev=ii
        ii=ii+1
        IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 198
         END IF
        END DO
198     CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF
      jpdt(2)=3
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      maxpts=gfld%ndpts
      print*,'80-m v iret=',iret
      if(iret.eq.0) then
       ii=1
       jj=jj1
       do kk=1,itot
        v80(ii,jj)=gfld%fld(kk)
        wnd80(ii,jj)=sqrt(u80(ii,jj)**2+v80(ii,jj)**2)
        mskwnd80(ii,jj)=1
        iprev=ii
        ii=ii+1
        IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 199
         END IF
        END DO
199     CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF

      deallocate(u80)
      deallocate(v80)

C     CEILING

      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=3
      jpdt(2)=5
      jpdt(10)=215
      call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
      print*,'ceiling iret=',iret
      maxpts=gfld%ndpts
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        DO kk = 1, itot
          ceil(ii,jj) = gfld%fld(kk)
          IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
            mskceil(ii,jj) = 1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 171
          END IF
        END DO
171     CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF


C     -== GET PRESSURE LEVEL VARIABLES Z, T, U, V & Q/RH ==-

      DO 100 iv = 1, 5
        j = 0
c       ETA OUTPUTS SPECIFIC HUMIDITY
c       IF ((modgen.eq.83.or.modgen.eq.84.or.modgen.eq.85.or.modgen.eq.
c    +              89.or.modgen.eq.110).and.iv.eq.5) jpds(5) = ivar(6)
        DO 90 l = 1, numlev
          j=0
          jdisc=0
          jids=-9999
          jpdt=-9999
          jpdt(1)=idisp(iv)
          jpdt(2)=ivar(iv)
          jpdt(10)=100
          jpdt(12) = levs(l,levset) * 100
          print*,'jpdt(12)=',jpdt(12)
          call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
          maxpts=gfld%ndpts
          print*,'ivar(iv),idisp(iv),iret=',ivar(iv),idisp(iv),iret
        
c
c Check to see if the winds are grid-relative or earth-relative
c
c     if(btest(kgds(6),3) ) then
c     if(btest(gfld%igdtmpl(12),3) ) then
c      iearth=0  !  grid relative
c     else
c      iearth=1  !  earth relative
c     endif
C
C  1/31/05 : look for RH in Eta output if specific humidity is not present
C
      IF (iv.eq.5.and.iret.ne.0) THEN
       jpdt(2)=1
       call getgb2(lugb,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *       unpack,j,gfld,iret)
       maxpts=gfld%ndpts
       print*,'no q RH iret=',iret
      ENDIF
          IF (iret.eq.0) THEN
            ii = 1
            jj = jj1
            print*,'maxpts=',maxpts
            DO kk = 1, itot
              IF (iv.eq.1) THEN
                z(ii,jj,l) = gfld%fld(kk)
                IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
                  maskz(ii,jj,l) = 1
                END IF
              ELSE IF (iv.eq.2) THEN
                t(ii,jj,l) = gfld%fld(kk)
                IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
                  maskt(ii,jj,l) = 1
                END IF
              ELSE IF (iv.eq.3) THEN
                u(ii,jj,l) = gfld%fld(kk)
                IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
                  masku(ii,jj,l) = 1
                END IF
              ELSE IF (iv.eq.4) THEN
                v(ii,jj,l) = gfld%fld(kk)
                IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
                  maskv(ii,jj,l) = 1
                END IF
              ELSE IF (iv.eq.5) THEN
                IF (jpdt(2).eq.1) THEN
                  IF (maskt(ii,jj,l).ge.1.) THEN
                    rh = max(0.01,gfld%fld(kk)*.01)
                    rh = min(1.00,rh)
                    vaps = w3fa09(t(ii,jj,l)) * 10.
                    vap = rh * vaps
                    q(ii,jj,l) = .622 * vap / (pgd(ii,jj,l)-.378*vap)
                  ELSE
                    q(ii,jj,l) = 0.0
                    gfld%bmap(kk) = .false.
                  END IF
                ELSE
                  q(ii,jj,l) = gfld%fld(kk)
                END IF
                IF (gfld%ibmap.eq.255.or.gfld%bmap(kk)) THEN
                  maskq(ii,jj,l) = 1
                  alnq(ii,jj,l) = alog(q(ii,jj,l))
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

c         deallocate(grid)
c         deallocate(grid2)
c         deallocate(grid3)
c         deallocate(grid4)

c         deallocate(mask)

          call gf_free(gfld)

          RETURN
          END
