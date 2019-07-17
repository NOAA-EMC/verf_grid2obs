C-----------------------------------------------------------------------
C  EDITBUFR WILL RETAIN DATA LOCATED IN RETENTION GRID DEFINED
C  FROM STANDARD INPUT. THE FORMAT OF THE PARAMETER
C  LIST IS AS FOLLOWS:
C
C  FIRST WE READ THE GRID NUMBER THAT DEFINES THE RETENTION AREA
C
C  GRID#     - IS THE NUMBER OF THE RETENTION GRID
C
C  THEN WE READ THE DATE INFORMATION
C
C  YYMMDD    - IS THE (TWO-DIGIT) YEAR, MONTH, AND DAY TO BE KEPT
C              IF YYMMDD IS ZERO, THEN NO DATE CHECK IS PERFORMED
C              IF YYMMDD IS LESS THAN ZERO, THEN ITS MAGNITUDE IS
C              USED AS A LIMIT FOR THE MAGNITUDE OF THE OBSERVATION'S
C              DHR - WHICH IS THE TIME DIFFERENCE IN HUNDREDTH'S OF
C              AN HOUR BETWEEN THE OB TIME AND THE ANALYSIS TIME
C
C  THEN WE READ UP TO 30 DATATYPEs
C
C  DATATYPE  - IS THE OI REPORT TYPE (I.E. 120 IS A RAWINSONDE).
C              IT IS POSSIBLE TO SPECIFY A RANGE OF DATA TYPES
C              BY GIVING AN ABBREVIATED FORM. THAT IS, TO SPECIFY
C              TYPES 120 THROUGH 129, THE PARAMETER IS '12'. TYPES
C              100 THROUGH 199 WOULD BE SPECIFIED BY '1'. AND SO ON.
C
C
C  THE INPUT  PREPBUFR FILE IS ASSIGNED TO UNIT 20
C  THE OUTPUT PREPBUFR FILE IS ASSIGNED TO UNIT 50
C  THE OUTPUT PREPBUFR FILE FOR MAX VALUES IS ASSIGNED TO UNIT 70
C
C* Log:
C* K. Brill/EMC		11/98	Use INT for grid pt location & change
C*				check.  Use INT for TYP/TFAC.
C-----------------------------------------------------------------------
      PROGRAM editbufr
 
      PARAMETER (mxrt = 30,mxts = 4,mxtb = 1000000)
      PARAMETER (mxsta = 10000)
 
      dimension ozone1(mxsta,24),oz_max1(mxsta),oz_tmp1(24)
      dimension ozone8(mxsta,24),oz_max8(mxsta),oz_tmp8(24)
      DIMENSION xsta(mxsta), ysta(mxsta)
      REAL*8 bdate
      CHARACTER*8 subset, sslast
      DIMENSION rtyp(mxrt), tfac(mxrt)
      REAL*8 tab(mxts,mxtb) 
      DIMENSION ikep(mxtb)
      DIMENSION xob(mxtb), yob(mxtb), typ(mxtb), dhr(mxtb), oz(mxtb)
      DIMENSION xob1(mxtb), yob1(mxtb)
      DIMENSION nret(mxrt)
      LOGICAL keep, within(mxtb)
      LOGICAL latlong, lambert, polarstereo
      COMMON /gridef/ imax, jmax, kmax, alat1, elon1, dxx, dyy, elonv, 
     +            alatan, latlong, lambert, polarstereo
 
c     DATA lubfi /20/
c     DATA lubfp /21/
      DATA lundx /22/
      DATA lubfj /50/
      DATA lubfm /70/
c     CHARACTER*6 OZONE_AVE
      LOGICAL IFCOPY
 
      DATA bmiss /10E10/
      REAL*8 hdr(10), cat(255), obs(10,255), qms(10,255)
      CHARACTER*8 staid
      CHARACTER*80 headr, obstr, qmstr
      EQUIVALENCE (hdr(1),staid)
      headr='SID XOB YOB DHR ELV TYP T29 ITP           '
      obstr='TPHR QCIND COPOPM                         '
      qmstr='PQM QQM TQM ZQM WQM                       '

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      CALL DATELEN(10)

C     READ OZONE AVERAGING PERIOD (1HR or 8HR)
C     ----------------------------------------
c     READ (5,'(A6)') ozone_ave
 
C  READ THE RETENTION GRID NUMBER
C  ------------------------------
 
      READ (5,'(A80)')
      READ (5,*,END=10) iretgrid
      CALL gtgdef(iretgrid,istat)
      IF (istat.ne.0) THEN
        PRINT *, 'GTGDEF RETURNED NON-ZERO ISTAT=', istat
        STOP 'BAD RETENTION GRID'
      END IF
 
C     READ THE RETENTION DATE
C     -----------------------
 
      READ (5,'(A80)')
      READ (5,*,END=10) jdate
      IF (jdate.lt.0) THEN
        jwndo = -jdate
        jdate = 0
      ELSE
c       jwndo = 300
        jwndo = 1200
      END IF
      PRINT *, 'TIME-WINDOW=', jwndo, ' HUNDRETHS OF AN HOUR'
 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c         time=24
c         CALL raddate(jdate,time,bdate)
c         print*,'After call to raddate',jdate,bdate
c         jdate_dp1 = bdate
c         print*,jdate,jdate_dp1
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C     READ THE RETENTION OBTYPES
C     --------------------------
 
      READ (5,'(A80)')
 
      nobtyp = 0
 
      DO n = 1, mxrt
        READ (5,*,END=10) rtyp(n)
        nobtyp = n
      END DO
 
   10 IF (nobtyp.eq.0) THEN
        STOP 'EMPTY RETENTION LIST'
      ELSE
        PRINT *
        DO n = 1, nobtyp
          tfac(n) = 0
          IF (rtyp(n).lt.1000) tfac(n) = 1
          IF (rtyp(n).lt.100) tfac(n) = 10
          IF (rtyp(n).lt.10) tfac(n) = 100
          IF (rtyp(n).eq.0) tfac(n) = 1000
c         IF (tfac(n).eq.0) CALL abort('BAD OBTYPE PARAMETER INPUT')
          IF (tfac(n).eq.0) then
              print*,'BAD OBTYPE PARAMETER INPUT'
              CALL abort
          endif
          PRINT 1000, n, rtyp(n), iretgrid, jdate, jwndo
 1000     FORMAT ('RETENTION OBTYP ',I2,F6.0,3I12)
        END DO
      END IF
 
      do ifile=1,2  !For AQ we need 04z-04z, so we need the next day bufr

      if(ifile.eq.1) then
        lubfi=20
        lubfj=50
      elseif(ifile.eq.2) then
        lubfi=21
        lubfj=51
      endif

C     MAKE A TABLE OUT OF THE OBSERVATION TYPS, LATS, AND LONS
C     --------------------------------------------------------
 
      CALL datebf(lubfi,iy,im,id,ih,idate)
      if(ifile.eq.1) idate_day1=idate

      PRINT *, 'DATES OF INPUT BUFR FILES ', idate
      CALL ufbtab(lubfi,tab,mxts,mxtb,ntab,'XOB YOB TYP DHR')
      PRINT *, 'TABULATED ', ntab, ' REPORTS FROM INPUT BUFR FILE'

      
      if(ifile.eq.1) then
       ntab1=1
       ntab2=ntab
      elseif(ifile.eq.2) then
       ntab1=ntab2+1
       ntab2=ntab2+ntab
      endif
      print*,'ifile,ntab1,ntab2=== ',ifile,ntab1,ntab2

      DO n = 1, ntab
        if(ifile.eq.1) m=n
        if(ifile.eq.2) m=n+ntab1-1
        xob(n) = tab(1,n)
        yob(n) = tab(2,n)
        typ(n) = tab(3,n)
        dhr(n) = tab(4,n)
        xob1(m) = xob(n)
        yob1(m) = yob(n)
c       if(n.eq.10000) then
c       if(n.ge.1.and.n.le.500) then
c         print*,'xob,yob,typ,dhr=',xob(n),yob(n),typ(n),dhr(n)
c       endif
      END DO
    
 
C     MARK THE OBSERVATIONS LOCATED WITHIN BOTH GRID & TIME-WINDOW
C     ------------------------------------------------------------
 
      ikep = 0
      nwith = 0
      nwndo = 0
c     N = 1 
c     print *,N,ALAT1,ELON1,DXX,ELONV,ALATAN
c     print *,N,IMAX,JMAX,KMAX
      DO n = 1, ntab
 
C       Depending on the grid type,
C       CALCULATE Grid coordinates of obs lat,long 
C       ------------------------------------------
        IF (latlong) THEN
C         Latitiude - Longitude grid  - NOT global
          xi = (xob(n)-elon1) / dxx + 1.0
          yj = (yob(n)-alat1) / dyy + 1.0
        END IF
        IF (polarstereo) THEN
C         Polar Stereographic grid
C         W3FB06 wants grid spacing in meters
          dxm = dxx * 1000.
          CALL w3fb06(yob(n),xob(n),alat1,elon1,dxm,elonv,xi,yj)
        END IF
        IF (lambert) THEN
C         Lambert Conic Conformal grid
C         W3FB11 wants grid spacing in meters
          dxm = dxx * 1000.
          CALL w3fb11(yob(n),xob(n),alat1,elon1,dxm,elonv,alatan,xi,yj)
        END IF
        kxi = INT(xi)
        kyj = INT(yj)
c       IF(N.EQ.1.OR.MOD(N,500).EQ.0) THEN
c       print *,N,YOB(N),XOB(N),TYP(N),DHR(N)
c       print *,N,KXI,KYJ,IMAX,JMAX
c       ENDIF

C       CHECK IF OB IS WITHIN DOMAIN
C       ----------------------------
        IF ( kxi .lt. 1 .or. kxi .gt. imax .or.
     +	     kyj .lt. 1 .or. kyj .gt. jmax ) THEN
          within(n) = .false.
        ELSE
          within(n) = .true.
          nwith = nwith + 1
C         NOW CHECK IF OB IS ALSO WITHIN TIME WINDOW
C         ------------------------------------------
          IF (nint(abs(dhr(n)*100.)).gt.jwndo) THEN
            nwndo = nwndo + 1
            within(n) = .false.
          END IF
        END IF
      END DO
      PRINT *, 'GRID-AREA INCLUSIONS =', nwith
      PRINT *, 'TIME-WINDOW EXCLUSIONS =', nwndo
C     NOW CHECK OB TYPE
C     -----------------
      print*,'ntab,nobtyp=',ntab,nobtyp
      DO n = 1, ntab
        DO nr = 1, nobtyp
          keep = nint(rtyp(nr)) .eq. int(typ(n)/tfac(nr)) .or. 
     +                rtyp(nr) .eq. 0
c         if (keep.eqv..false.) then
c           print*,'rtyp(nr),typ(n),tfac(nr),n,nr=',
c    *        rtyp(nr),typ(n),tfac(nr),n,nr
c         endif
          IF (within(n).and.keep) ikep(n) = nr
        END DO
      END DO
 
C     OPEN THE OUTPUT FILE (POSITION MOD IF THE FILE EXISTS)
C     ------------------------------------------------------
C     OPEN THE INPUT FILE TO READ THROUGH
C     -----------------------------------
 
      print*,'---- 1 -----'
      CALL openbf(lubfi,'IN ',lubfi)
      print*,'---- 2 -----'
      CALL datebf(lubfj,iy,im,id,ih,idate)
      print*,'---- 3 -----'
c      print*,'===lubfj=== ',lubfj,idate
c--marIF (idate.ge.0) CALL openbf(lubfj,'APN',lubfi)
c--marIF (idate.lt.0) CALL openbf(lubfj,'OUT',lubfi)
      CALL openbf(lubfj,'OUT',lubfi)
      print*,'---- 4 -----'
      print*,'---- 5 -----'
       print*, '===lubfj=== ',ifile,lubfj

C     WRITE OUT ALL DATA MARKED FOR RETENTION
C     ---------------------------------------
 
      isub = 0
      nrep = 0
      nrej = 0
      nret = 0
 
      DO WHILE (ireadmg(lubfi,subset,idate).eq.0)
        nsub = nmsub(lubfi)
        IF (isub.eq.0) THEN
          sslast = subset
          PRINT *, 'SUBSET,NSUB,IDATE ', subset, nsub, idate
        END IF
        nrep = nrep + nsub
        adate = idate
c       print*,'nsub=',nsub
        DO n = 1, nsub
          isub = isub + 1
c         print*,'isub=',isub
          iobtyp = ikep(isub)
c         print*,'adate,dhr(isub)=',adate,dhr(isub)
c         print*,'Before call to raddate',idate,dhr(isub),bdate
          if(abs(dhr(isub)).gt.0.1E+20) goto 99


          CALL raddate(idate,dhr(isub),bdate)
c         print*,'After call to raddate',idate,dhr(isub),bdate
c         kdate = bdate * .01
          kdate = bdate
c         print*,idate,dhr(isub),kdate
c    +             ,xob(isub),yob(isub),jdate
c         print*,'iobtyp=',iobtyp

c--mar    IF (iobtyp.gt.0.and.(kdate.eq.jdate.or.jdate.eq.0)) THEN
c--mar    IF (iobtyp.gt.0.and.kdate.eq.jdate) THEN

          IF (iobtyp.gt.0) THEN
            IF (subset.ne.sslast) THEN
              PRINT *, 'SUBSET,ISUB,IDATE ', subset, isub, idate
              CALL closmg(lubfj)
              sslast = subset
            END IF
c           print*,'===idate_day1=== ',idate_day1
            CALL openmb(lubfj,subset,idate)
c--mar      CALL openmb(lubfj,subset,kdate)
            CALL copysb(lubfi,lubfj,iret)
            nret(iobtyp) = nret(iobtyp) + 1
          print*,idate,dhr(isub),kdate
     +             ,xob(isub),yob(isub),isub,nret(iobtyp),jdate
          ELSE
c           print*,'In ELSE block of IF iobtyp block'
c           print*,'lubfi=',lubfi
            CALL copysb(lubfi,0,iret)
c           print*,'After call to copysb, iret=',iret
            nrej = nrej + 1
99        END IF
        END DO
      END DO
 
      CALL closbf(lubfi)
      CALL closbf(lubfj)
 
      enddo   !end ifile loop

c     stop
 
c================================================================
C     MAKE A LIST OF STATIONS
C     -----------------------

      xsta=0.
      ysta=0.
      ista=1
      oz_max=0.

      do n=1,ntab2
       do nsta=1,mxsta
        if((xob1(n).ne.xsta(nsta)).or.(yob1(n).ne.ysta(nsta))) then
         iexist=0
        else
         iexist=1
         goto 500
        endif
       enddo
 500  continue
         if(iexist.eq.0)then
         xsta(ista)=xob1(n)
         ysta(ista)=yob1(n)
         ista=ista+1
         endif
      enddo

      num_sta=ista-1
      print*,'NUMBER OF STATIONS IN INPUT BUFR FILE   ',num_sta
c     stop

c     do ista=1,num_sta
c     print*,ista,xsta(ista),ysta(ista)
c     enddo
c================================================================

c===========================================================
C     OPEN INPUT FILE AGAIN TO READ ALL RECORDS AND MAKE MAX
C     ------------------------------------------------------

      do ifile=1,2

      if(ifile.eq.1) then
        lubfi=50
      elseif(ifile.eq.2) then
        lubfi=51
      endif


      CALL openbf(lubfi,'IN ',lubfi)
      DO WHILE (ireadmg(lubfi,subset,idate).eq.0)
       DO WHILE (ireadsb(lubfi).eq.0)
       CALL ufbint(lubfi,hdr,10,1,iret,headr)
       CALL ufbint(lubfi,obs,10,255,iret,obstr)

       xreal4=hdr(2)
       yreal4=hdr(3)
       do ista=1,num_sta
       if(xreal4.eq.xsta(ista).and.yreal4.eq.ysta(ista))then
        if(ifile.eq.1) ihr=hdr(4)+8
        if(ifile.eq.2) ihr=hdr(4)+32
         if (obs(1,1).eq.-1.and.ihr.ge.1.and.ihr.le.24) then
           ozone1(ista,ihr)=obs(3,1)
         elseif (obs(1,1).eq.-8.and.ihr.ge.8.and.ihr.le.31) then
           ihrm7=ihr-7
           ozone8(ista,ihrm7)=obs(3,1)
         endif
       endif
       enddo

       END DO
      END DO
      CALL closbf(lubfi)
      enddo   !end ifile loop

      oz_max1=bmiss
      oz_max8=bmiss

      do ista=1,num_sta
       count_bmiss1=0.
       count_bmiss8=0.
       iavecount=0
       ozave=0.
       do ihr=1,24
        oz_tmp1(ihr)=ozone1(ista,ihr)
        oz_tmp8(ihr)=ozone8(ista,ihr)
          if(oz_tmp1(ihr).eq.bmiss) count_bmiss1=count_bmiss1+1
          if(oz_tmp8(ihr).eq.bmiss) count_bmiss8=count_bmiss8+1
c       print*,ihr,oz_tmp1(ihr),oz_tmp8(ihr)
        if(oz_tmp1(ihr).ne.bmiss)then
         ozave=ozave+oz_tmp1(ihr)
         iavecount=iavecount+1
        endif
       enddo
      if(iavecount.ne.0) ozave=ozave/iavecount
      if(count_bmiss1.lt.3) oz_max1(ista)=ozave
      ii=maxloc(oz_tmp1,dim=1,mask=oz_tmp1.ne.bmiss)
      if(ii.ne.0) oz_max1(ista)=ozone1(ista,ii)
c--   if(oz_max1(ista).eq.0.) oz_max1(ista)=bmiss
      if(oz_max1(ista).eq.0.or.count_bmiss1.ge.3) 
     +    oz_max1(ista)=bmiss
c         print*,'max',oz_max1(ista),ii
      hr_max_edt=ii
      print*,'-mxsta1-',ista,xsta(ista),ysta(ista),oz_max1(ista),
     +       hr_max_edt,'count_bmiss1=',count_bmiss1
      ii=maxloc(oz_tmp8,dim=1,mask=oz_tmp8.ne.bmiss)
      if(ii.ne.0) oz_max8(ista)=ozone8(ista,ii)
c--   if(oz_max8(ista).eq.0.) oz_max8(ista)=bmiss
      if(oz_max8(ista).eq.0.or.count_bmiss8.ge.3) 
     +    oz_max8(ista)=bmiss
c         print*,'max',oz_max8(ista),ii
      hr_max_edt=ii+7
      print*,'-mxsta8-',ista,xsta(ista),ysta(ista),oz_max8(ista),
     +       hr_max_edt,'count_bmiss8=',count_bmiss8
      enddo


C     OPEN THE INPUT FILE AGAIN TO READ THROUGH
C     -----------------------------------------
C     OPEN OUTPUT FILE FOR MAX VALUES
C     -------------------------------

      lubfi=50
      CALL openbf(lubfi,'IN ',lubfi)
      CALL openbf(lubfm,'OUT',lundx)

      DO WHILE (ireadmg(lubfi,subset,idate).eq.0)
        nsub = nmsub(lubfi)

        DO WHILE (ireadsb(lubfi).eq.0)
        CALL ufbint(lubfi,hdr,10,1,iret,headr)
        CALL UFBINT(LUBFI,CAT, 1,255,iret,'CAT')
        CALL UFBINT(LUBFI,QMS,10,255,iret,QMSTR)
        CALL ufbint(lubfi,obs,10,255,iret,obstr)
 
        obs(3,1) = bmiss

        if(hdr(4).eq.0) then    !assign time 12Z to max data
         xreal4=hdr(2)
         yreal4=hdr(3)
         iexist=0
         do ista=1,num_sta
          if(xreal4.eq.xsta(ista).and.yreal4.eq.ysta(ista))then
            if (obs(1,1).eq.-1) then
             obs(3,1)=oz_max1(ista)
            elseif (obs(1,1).eq.-8) then
             obs(3,1)=oz_max8(ista)
            endif
            iexist=1
          else
            if(ista.eq.num_sta.and.iexist.eq.0)print*,'Station ',
     +             xreal4,yreal4,' is not in the list, ave =',
     +             obs(1,1)
          endif
         enddo
c        print*,'Station ',staid,xreal4,yreal4,obs(1,1),obs(3,1)
         if(obs(1,1).eq.-1)print*,'Station ',staid,xreal4,yreal4

C        WRITE OUTPUT
C        ------------

         CALL openmb(lubfm,subset,idate_day1)
         CALL ufbint(lubfm,hdr,10,1,iret,headr)
         CALL ufbint(lubfm,cat,1,1,iret,'CAT')
         CALL ufbint(lubfm,obs,10,1,iret,obstr)
         CALL ufbint(lubfm,qms,10,1,iret,qmstr)
         CALL writsb(lubfm)
        endif

        END DO
      END DO

      CALL closbf(lubfi)
      CALL closbf(lubfm)
c===========================================================
C     GENERATE REPORT AND EXIT
C     ------------------------
 
      PRINT *
      DO n = 1, nobtyp
        PRINT *, 'RETAINED REPORTS OF TYPE#', n, rtyp(n), '=', nret(n)
      END DO
      mret = nrep - nrej
      PRINT *
      PRINT *, 'TOTAL INCLUDED REPORTS=', mret
      PRINT *
      PRINT *, 'TOTAL EXCLUDED REPORTS=', nrej
      PRINT *
      PRINT *, 'EDITBUFR PROCESSED ', nrep, ' REPORTS '
      PRINT *
 
      STOP
      END
