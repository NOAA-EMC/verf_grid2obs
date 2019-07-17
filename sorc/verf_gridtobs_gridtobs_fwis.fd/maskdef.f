c     SUBROUTINE maskdef(imax,jmax,alat1,elon1,dxx,dyy,elonv,
c    +            alatan,latlong,lambert,polarstereo,namfcst,numfcst)
c     subroutine maskdef(namfcst,numfcst)
      subroutine maskdef(icyc)
C                .      .    .                                       .
C SUBPROGRAM:    GTGDEF      RETRIEVE grid definition parameters
C   PRGMMR: Geoff  DiMego    ORG: W/NP22     DATE: 97-12-29
C
C ABSTRACT: RETRIEVES grid definition parameters 
C
C PROGRAM HISTORY LOG:
C   98-01-04  Geoff DiMego   Brand new code
C
C USAGE:    CALL GTGDEF(IGRID,ISTAT,IMAX,JMAX,ALAT1,ELON1,
C    1    DXX,DYY,ELONV,ALATAN,LATLONG,LAMBERT,POLARSTEREO)
C
C   INPUT ARGUMENT LIST:
C     IGRID    - INTEGER NUMBER OF desired grid
C
C   OUTPUT ARGUMENT LIST:
C     ISTAT    - INTEGER =0  MEANS SUCCESSFUL COMPLETION
C                        = 1 MEANS GRID NOT currently supported
C     IMAX,JMAX- DIMENSIONS OF GRID
C     ALAT1,ELON1 LOCATION OF ORIGIN PT(1,1)
C     DXX,DYY  - MESH LENGTHS
C     ELONV    - VERTICAL LONGITUDE
C     ALATAN   - REFERENCE LATITUDE (FOR LAMBERT)
C     LATLONG,LAMBERT,POLARSTEREO  PROJECTION-TYPE SWITCHES
C            WHERE LATLONG,LAMBERT,POLARSTEREO are type LOGICAL
C
C REMARKS:

C ATTRIBUTES:
C   LANGUAGE: FORTRAN-77
C   MACHINE:  CRAY C-90
C$$$
c GRIB2
      use grib_mod
c GRIB2

            INCLUDE 'parm.inc'

      LOGICAL latlongf, lambertf, polarstereof
      logical*1 mask(ilim*jlim),masku(ilim*jlim)
      real grid(ilim*jlim),gridu(ilim*jlim)
c     character*24 namfcst(numfcst)
      integer inamfcst(7)
c     integer inamfcst(13)
      integer icount,icountr

C GRIB2
      type(gribfield) :: gfld
      integer ifile,j,jdisc,jpdtn,jgdtn,iret
      integer,dimension(200) :: jids,jpdt,jgdt
      real fldmin,fldmax,firstval,lastval
      logical :: unpack=.true.
C GRIB2

      common /firewx/ imaxf(7),iminf(7),jmaxf(7),jminf(7),alat1f(7),
     *                elon1f(7), dxxf(7), dyyf(7), elonvf(7),alatanf(7),
     *                latlongf(7),lambertf(7),polarstereof(7),
     *                fcsthrf(7),imask(ilim,jlim,7),imasku(ilim,jlim,7),
     *                icount,ifirst,icycfac

      INTEGER kgds(25),jpds(25),jgds(25),kpds(25)

      data inamfcst /0, 6, 12, 18, 24, 30, 36/ 
c     data inamfcst /0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36/

C     
C GRIB2
c Set GRIB2 field identification values to search for
      j=0              ! search from 0
      jdisc=0          ! for met field:0 hydro: 1, land: 2
c-- set id section
      jids=-9999
c-- set product def template, using template 4.0
      jpdtn=0
c-- set product def array
      jpdt=-9999
!-- set grid def template
      jgdtn=-1
!-- set product def array
      jgdt=-9999
C GRIB2
C
      istat = 0
      icount = 0
      icountr = 0
      ii = 0
      imask = 0
      inam = 0

C
C Add cycle factor
c
      icycfac=icyc
      if(icycfac.ge.6.and.icycfac.lt.12) icycfac=icycfac-6
      if(icycfac.ge.12.and.icycfac.lt.18) icycfac=icycfac-12
      if(icycfac.ge.18) icycfac=icycfac-18

      
C     USE W3FI71 WITH THE INPUT GRID NUMBER TO GET THE PDS
C     FROM WHICH TO EXTRACT THE GRID DEFINITION PARAMETERS
C     
c
       do 50 i=1,7
       icountr = icountr + 1
       lug=30+icountr
       lugi=lug+1
c      j = 0
       jf = ilim*jlim
c      jpds = -1
c      jpds(5) = 11
c      jpds(6) = 105
c      jpds(7) = 2
c      jpds(14) = inamfcst(i)
c      jpds(14) = i + icycfac
      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=0
      jpdt(2)=0
      jpdt(10)=103
      jpdt(12)=2
      jpdt(9)=inamfcst(i) + icycfac
c     jpdt(9)=i
       print*,'temp i,lug,inam=',i,lug,inamfcst(i)+icycfac
c      call getgb(lug,lugi,jf,j,jpds,jgds,kf,k,kpds,kgds,
c    *    mask,grid,iret)
       call getgb2(lug,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *         unpack,j,gfld,iret)
       if(iret.ne.0) then
         icountr = icountr - 1
         goto 50
       endif
       maxpts=gfld%ngrdpts
       print*,'gfld%ngrdpts,gfld%ipdtlen=',gfld%ngrdpts,gfld%ipdtlen
c      do k=1,gfld%ipdtlen
c       print*,'gfld%ipdtmpl(k)=',gfld%ipdtmpl(k)
c      enddo
c      print*,'maxpts=',maxpts
c      print*,'kpds=',kpds
c      print*,'kgds=',kgds
       print*,'temp iret=',iret
       do k=1,maxpts
         grid(k)=gfld%fld(k)
c        if(inamfcst(i).eq.30) print*,'k, grid(k)=',k,grid(k)
       enddo
c      grid=gfld%fld
c      if(iret.ne.0) then
c        icountr = icountr - 1
c        goto 50
c      endif
       j = 0
       jf = ilim*jlim
c      jpds = -1
c      jpds(5) = 33
c      jpds(6) = 105
c      jpds(7) = 10
c      jpds(14) = inamfcst(i)
c      jpds(14) = i
      j = 0
      jids=-9999
      jpdt=-9999
      jpdt(1)=2
      jpdt(2)=2
      jpdt(10)=103
      jpdt(12)=10
      jpdt(9)=inamfcst(i) + icycfac
c     jpdt(9)=i
       print*,'wind,i,lug,inam=',i,lug,inamfcst(i) + icycfac
c      call getgb(lug,lugi,jf,j,jpds,jgds,kf,k,kpds,kgds,
c    *    masku,gridu,iret)
        call getgb2(lug,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *         unpack,j,gfld,iret)
c      print*,'kpds=',kpds
c      print*,'kgds=',kgds
       print*,'wind iret=',iret
       maxpts=gfld%ngrdpts
c      print*,'maxpts=',maxpts
       do k=1,maxpts
         gridu(k)=gfld%fld(k)
c        print*,'gridu(k)=',gridu(k)
       enddo
c      if(iret.eq.0) then
       icount = icount + 1
       icountr = icountr + 1
c      inam = inam + 6
       print*,'icount=',icount
c      endif
C     
C     FILL IN GRIDEF COMMON BLOCK
C     W3FI71 RETURNS 2 EXTRA VALUES AT BEGINNING OF ARRAY KGDS
C     THE FOLLOWING DEFINED REGARDLESS OF GRID PROJECTION
C     
      if(icount.eq.1) ifirst=inamfcst(i) + icycfac
c     if(icount.eq.1) ifirst=i
c     imaxf(icount) = kgds(2)
c     jmaxf(icount) = kgds(3)
c     igrid = kpds(3)
      fcsthrf(icount) = inamfcst(i) + icycfac
c     fcsthrf(icount) = i
C     
C     USE KGDS(1+2) TO DETERMINE GRID PROJECTION
C     
C     KGDS(1+2) = 0 ----> LATITUDE/LONGITUDE
C     KGDS(1+2) = 1 ----> MERCATOR (NOT YET USED)
C     KGDS(1+2) = 3 ----> LAMBERT CONFORMAL
C     KGDS(1+2) = 5 ----> POLAR STEREOGRAPHIC
C     
      itmpl = gfld%igdtnum
      if(itmpl.eq.0) then
        imaxf(icount) = gfld%igdtmpl(8)
        jmaxf(icount) = gfld%igdtmpl(9)
        latlongf(icount) = .true.
        lambertf(icount) = .false.
        polarstereof(icount) = .false.
      else if (itmpl.eq.30) then
        imaxf(icount) = gfld%igdtmpl(8)
        jmaxf(icount) = gfld%igdtmpl(9)
        latlongf(icount) = .false.
        lambertf(icount) = .true.
        polarstereof(icount) = .false.
      else if (itmpl.eq.20) then
        imaxf(icount) = gfld%igdtmpl(8)
        jmaxf(icount) = gfld%igdtmpl(9)
        latlongf(icount) = .false.
        lambertf(icount) = .false.
        polarstereof(icount) = .true.
      ELSE
        iret = 99
c       WRITE (6,*) ' KGDS(1+2) = ', kgds(1+2)
        WRITE (6,*) ' GRID CAN NOT BE USED IN THIS CODE IRET= ', iret
        istat = 1
        RETURN
      END IF
C     
C     SET THE REST OF THE GRID PARAMETERS BASED ON PROJECTION TYPE
C     
C     Change has been made for LATLON definition --- Yuejian Zhu
C     Changed back with checking -- K. Brill
C     
      IF (latlongf(icount)) THEN
        alat1f(icount) = gfld%igdtmpl(12) * 0.000001
        elon1f(icount) = gfld%igdtmpl(13) * 0.000001
	IF ( elon1f(icount) .lt. 0.0 ) elon1f(icount) 
     *     = elon1f(icount) + 360
        elonvf(icount) = 0.0
        alatanf(icount) = 0.0
        dyyf(icount) = gfld%igdtmpl(17) * 0.000001
        dxxf(icount) = gfld%igdtmpl(18) * 0.000001
        if (igrid.eq.181.or.igrid.eq.182) alat1f(icount) = 
     *               alat1f(icount) - (
     +               jmaxf(icount)-1) * dyyf(icount)

      END IF
C     
      IF (lambertf(icount)) THEN
        alat1f(icount) = gfld%igdtmpl(10) * 0.000001
        elon1f(icount) = gfld%igdtmpl(11) * 0.000001
	IF ( elon1f(icount) .lt. 0.0 ) 
     *     elon1f(icount) = elon1f(icount) + 360
        elonvf(icount) = gfld%igdtmpl(14) * 0.000001
	IF ( elonvf(icount) .lt. 0.0 ) 
     *     elonvf(icount) = elonvf(icount) + 360
        alatanf(icount) = gfld%igdtmpl(19) * 0.000001
        dxxf(icount) = gfld%igdtmpl(15) * 0.000001
        dyyf(icount) = gfld%igdtmpl(16) * 0.000001
      END IF
C     
      IF (polarstereof(icount)) THEN
        alat1f(icount) = gfld%igdtmpl(10) * 0.000001
        elon1f(icount) = gfld%igdtmpl(11) * 0.000001
	IF ( elon1f(icount) .lt. 0.0 ) 
     *    elon1f(icount) = elon1f(icount) + 360
        elonvf(icount) = gfld%igdtmpl(14) * 0.000001
	IF ( elonvf(icount) .lt. 0.0 ) 
     *    elonvf(icount) = elonvf(icount) + 360
        alatanf(icount) = 0.0
        dxxf(icount) = gfld%igdtmpl(15) * 0.000001
        dyyf(icount) = gfld%igdtmpl(16) * 0.000001
      END IF
C     
      PRINT *, 'gridspecs ', lambertf(icount), alat1f(icount), 
     *           elon1f(icount), elonvf(icount), alatanf(icount)
     *           , dxxf(icount),
     +            dyyf(icount)
C     
      WRITE (6,*) ' GREETINGS FROM THE GRID-DEFINITION CODE! '
c     WRITE (6,*) ' THE GRID YOU HAVE CHOSEN IS NUMBER ', igrid
      IF (latlongf(icount)) THEN
        WRITE (6,*) ' A LAT/LON GRID WITH RES= ', dxxf(icount), ' BY ', 
     *              dyyf(icount), 
     +              ' DEG'
      ELSE IF (polarstereof(icount)) THEN
        WRITE (6,*) ' A POLAR STEREO GRID CENTERED AT ', elonvf(icount)
     *     , ' DEG E'
        WRITE (6,*) ' AND A HORIZONTAL RESOLUTION OF ', dxxf(icount)
     *      , ' KM'
      ELSE IF (lambertf(icount)) THEN
        WRITE (6,*) ' A LAMBERT CONFORMAL GRID CENTERED AT ', 
     *              elonvf(icount), 
     +              ' DEG E'
        WRITE (6,*) ' AND A HORIZONTAL RESOLUTION OF ', 
     *    dxxf(icount), ' KM'
      END IF
      WRITE (6,*) ' HORIZONTAL DIMENSIONS ARE ', imaxf(icount)
     *    , ' X', jmaxf(icount)

      ii = 1
      jj = 1
      do kk = 1,maxpts
       imask(ii,jj,icount) = grid(kk)
       imasku(ii,jj,icount) = gridu(kk) 
       iprev = ii
       ii = ii + 1
       if (mod(iprev,imaxf(icount)).eq.0) then
         ii = 1
         jj = jj + 1
         if (jj.gt.jmaxf(icount)) goto 40
       endif
      enddo
40    continue

50    continue
  
      call gf_free(gfld)

      RETURN
      END
