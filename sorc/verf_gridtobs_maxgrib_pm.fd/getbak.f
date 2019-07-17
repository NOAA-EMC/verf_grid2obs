      SUBROUTINE getbak(lugb,lugi,numlev,vdate,fhr,istat)
C                .      .    .                                       .
      INCLUDE 'parm.inc'
      LOGICAL latlong, lambert, polarstereo
      COMMON /gridef/ imax, jmax, kmax, alat1, elon1, dxx, dyy, elonv, 
     +            alatan, latlong, lambert, polarstereo
C
C
      COMMON /surfce/ ps(ilim,jlim), mskps(ilim,jlim), zs(ilim,jlim), 
     +            mskzs(ilim,jlim), ts(ilim,jlim), mskts(ilim,jlim), 
     +            qs(ilim,jlim), mskqs(ilim,jlim), us(ilim,jlim), 
     +            mskus(ilim,jlim), vs(ilim,jlim), mskvs(ilim,jlim), 
     +            pm(ilim,jlim), mskpm(ilim,jlim), cape(ilim,jlim),
     +            mskcp(ilim,jlim), cin(ilim,jlim), mskcn(ilim,jlim),
     +            pli(ilim,jlim), mskli(ilim,jlim),
     +            oz1s(ilim,jlim), mskoz1(ilim,jlim),
     +            oz8s(ilim,jlim), mskoz8(ilim,jlim),
     +            oz1_max(ilim,jlim),oz8_max(ilim,jlim)
C
      COMMON /days/ day, cycle 
      DIMENSION grid(itot), grid1(itot)
      DIMENSION RICDAT(5),IFDATE(8),IVDATE(8)
      INTEGER jpds(400), jgds(400), kpds(400), kgds(400)
      CHARACTER*8 modesc(210)
      LOGICAL*1 mask(itot)
      CHARACTER*4 day
      CHARACTER*3 cycle
C
      DATA modesc /38 * 'xxxxxxxx', '80KM NGM', 37 * 'xxxxxxxx', 
     +            'T126 AVN', 'T126 MRF', 'xxxxxxxx', 'T62 MRF ', 2 * 
     +            'xxxxxxxx', '80KM ETA', '32KM ETA', '29KM ETA', 
     +            '60KM RUC', 2 * 'xxxxxxxx', '48KM ETA', 15 * 
     +            'xxxxxxxx', '40KM RUC', 4 * 'xxxxxxxx', '10KM ETA',
     +            99 * 'xxxxxxxx', '12KM CMQ'/
C
      CHARACTER cbuf(mbuf)
      LOGICAL*1 lb(jf)
      REAL f(jf)
      REAL*8 vdate
      REAL*8 bdate
      INTEGER jens(400), kens(400)

      istat = 0
C     
C     READ INDEX FILE TO GET GRID SPECS 
C     
      irgi = 1
      irgs = 1
      kmax = 0
      jr = 0
      kskip = 0
      CALL getgi(lugi,kskip,mbuf,cbuf,nlen,nnum,irgi)
      REWIND lugi
C     
C     NNUM IS THE NUMBER OF GRIB MESSAGES IN THE FILE. READ THE
C     FIRST PDS ONLY TO GET THE GRID NUMBER AND DATE
C     
C     
      DO k = 1, nnum
        jr = k - 1
        jpds = -1
        jgds = -1
        CALL getgb1s(cbuf,nlen,nnum,jr,jpds,jgds,jens,kr,kpds,kgds,kens,
     +              lskip,lgrib,irgs)
C       
C       USE  F I R S T  PDS TO GET START DATE & FORECAST LENGTH
C       
        IF (k.eq.1) THEN
          iys = kpds(8)
c----- For ozone, because it's averaging
          fhr = kpds(15)
c         print*,'==== fhr= ',kpds(14),fhr
          ifdate(2) = kpds(9)
          ifdate(3) = kpds(10)
          ifdate(5) = kpds(11)
          ricdat(2) = kpds(15)
          icent = kpds(21)
          if(ifdate(2).lt.100) then
            ifdate(1) = (icent-1)*100 + iys
          else
            ifdate(1) = icent * 100
          endif
        END IF
      END DO
c  10 CONTINUE
C     
C     FILL IN GRIDEF COMMON BLOCK
C     THE FOLLOWING DEFINED REGARDLESS OF GRID PROJECTION
C     
      imax = kgds(2)
      jmax = kgds(3)
C     
C     USE KGDS(1) TO DETERMINE GRID PROJECTION
C     
C     KGDS(1) = 0 ----> LATITUDE/LONGITUDE
C     KGDS(1) = 1 ----> MERCATOR (NOT YET UESD)
C     KGDS(1) = 3 ----> POLAR STEREOGRAPHIC
C     KGDS(1) = 5 ----> LAMBERT CONFORMAL
C     
      IF (kgds(1).eq.0) THEN
        latlong = .true.
        lambert = .false.
        polarstereo = .false.
      ELSE IF (kgds(1).eq.3) THEN
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
        alat1 = kgds(4) * 0.001
        elon1 = kgds(5) * 0.001
        IF (elon1.lt.0.0) elon1 = elon1 + 360.
        elonv = 0.0
        alatan = 0.0
        dxx = kgds(9) * 0.001
        dyy = kgds(10) * 0.001
C       IF we have a damn global grid with pt 1,1 at the NP, flip it
        IF (alat1.ge.89.999) THEN
          alat1 = alat1 - dyy * (jmax-1)
          alat1 = max(-90.,alat1)
          jj1 = jmax
          jjinc = -1
        END IF
      END IF
C     
      IF (lambert) THEN
        alat1 = kgds(4) * 0.001
        elon1 = kgds(5) * 0.001
        IF (elon1.lt.0.0) elon1 = elon1 + 360.
        elonv = kgds(7) * 0.001
        IF (elonv.lt.0.0) elonv = elonv + 360.
        alatan = kgds(12) * 0.001
        dxx = kgds(8) * 0.001
        dyy = kgds(9) * 0.001
      END IF
C     
      IF (polarstereo) THEN
        alat1 = kgds(4) * 0.001
        elon1 = kgds(5) * 0.001
        IF (elon1.lt.0.0) elon1 = elon1 + 360.
        elonv = kgds(7) * 0.001
        IF (elonv.lt.0.0) elonv = elonv + 360.
        alatan = 0.0
        dxx = kgds(8) * 0.001
        dyy = kgds(9) * 0.001
      END IF
C     
C     ADD IGF FORECAST HOURS TO START DATE TO GET VALID DATE
C     ADDDATE IS MODIFIED VERSION OF W3ADDATE
C     
      call w3movdat(RICDAT,IFDATE,IVDATE)
      vdate = (ivdate(1)*1000000+ivdate(2)*10000+ivdate(3)*100+
     1          ivdate(5))
      print *,ifdate,ricdat
      print *,ivdate
      print *,vdate
C     
C     GET GRID NUMBER FROM PDS
C     
      igdnum = kpds(3)
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
      WRITE (6,*) ' WHICH IS THE ', modesc(modgen)
C     
C     PROCESS THE GRIB FILE
C     
      numval = imax * jmax
C     initialize mask arrays
      mskoz1 = 0
      mskoz8 = 0
C     
 1200 FORMAT (' IMAX,JMAX,NUMLEV,KMAX,LEVSET ',5I4)
c     print*,'----------1----------'

C     -== GET SURFACE FIELDS ==-
      l = 0
      iv = 0

C     SURFACE 1HR AVE PM 2.5 (ACTUALLY OZONE AT FIRST SIGMA MID-LAYER)

c     to start each new file with its index, set J=-1
      j = -1
      jpds = -1
      jpds(3) = igdnum
      jpds(5) = 157
c     jpds(6) = 1
c     jpds(7) = 0
      jpds(6) = 107
      jpds(7) = 10000
         jpds(14) = fhr-1
         if(jpds(14).eq.-1) jpds(14)=0

c     print*,kpds(3),kpds(5),kpds(6),kpds(7),kpds(14),kpds(15),igdnum
c     print*,'----------2----------'
      CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,
     +            iret)
c     print*,kpds(3),kpds(5),kpds(6),kpds(7),kpds(14),kpds(15),igdnum

       print*,'iret= ', iret
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        print*,'itot= ', itot
        DO kk = 1, itot
          oz1s(ii,jj) = grid(kk)
          grid1(kk)=oz1_max(ii,jj)
c         print*,oz1s(ii,jj)
          IF (mask(kk)) THEN
            mskoz1(ii,jj) = 1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 155
          END IF
        END DO
  155   CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
 1300   FORMAT (' IV, IVAR, L, IRET:  ',4I5)
      END IF

C     SURFACE 1HR AVE PM 2.5 from Smoke
c     no need to convert units from log10 here; done later in prepfits

c     to start each new file with its index, set J=-1
      j = -1
      jpds = -1
      jpds(3) = igdnum
      jpds(5) = 164
      jpds(6) = 1
      jpds(7) = 0
         jpds(14) = fhr-1
         if(jpds(14).eq.-1) jpds(14)=0

c     print*,kpds(3),kpds(5),kpds(6),kpds(7),kpds(14),kpds(15),igdnum
c     print*,'----------2----------'
      CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,
     +            iret)
c     print*,kpds(3),kpds(5),kpds(6),kpds(7),kpds(14),kpds(15),igdnum

       print*,'iret= ', iret
      IF (iret.eq.0) THEN
        ii = 1
        jj = jj1
        print*,'itot= ', itot
        DO kk = 1, itot
          oz1s(ii,jj) = grid(kk)
          grid1(kk)=oz1_max(ii,jj)
c         print*,oz1s(ii,jj)
          IF (mask(kk)) THEN
            mskoz1(ii,jj) = 1
          END IF
          iprev = ii
          ii = ii + 1
          IF (mod(iprev,imax).eq.0) THEN
            ii = 1
            jj = jj + jjinc
            IF (jj.gt.jmax.or.jj.lt.1) GO TO 157
          END IF
        END DO
  157   CONTINUE
      ELSE
        WRITE (6,1300) iv, jpds(5), l, iret
      END IF
          print*,'------1',fhr,oz1s(1,1),lugb,lugi,grid(1),grid1(1)
          print*,kpds

          jdate = (ivdate(1)*1000000+ivdate(2)*10000+ivdate(3)*100+
     1          ivdate(5))
          if (cycle.eq.'12Z') then
            if (day.eq.'DAY1') time = -24
            if (day.eq.'DAY2') time = -48
          endif
          if (cycle.eq.'06Z') then
            if (day.eq.'DAY1') time = -30
            if (day.eq.'DAY2') time = -54
          endif
          print*, jdate
          CALL raddate(jdate,time,bdate)
          print*, bdate
          idate=bdate
          print*, idate
          mxyr=idate/1000000
          mxmo=mod((idate/10000),100)
          mxda=mod((idate/100),100)
          mxhr=mod(idate,100)
          if(mxyr.ge.100) then
             mxcen=mod((mxyr/100),100) + 1
             mxyr=mod(mxyr,100)
          endif
          print*,mxcen,mxyr,mxmo,mxda,mxhr,idate

          kpds(21) = mxcen
          kpds(8)  = mxyr
          kpds(9)  = mxmo
          kpds(10) = mxda
          kpds(11) = mxhr

          if (cycle.eq.'12Z') then
            if (day.eq.'DAY1') then
              kpds(15) = 24
              kpds(14) = 23
            elseif (day.eq.'DAY2') then
              kpds(15) = 48
              kpds(14) = 47
            endif  
          endif
          if (cycle.eq.'06Z') then
            if (day.eq.'DAY1') then
              kpds(15) = 30
              kpds(14) = 29
            elseif (day.eq.'DAY2') then
              kpds(15) = 54
              kpds(14) = 53
            endif  
          endif

       print*,'-------putgb-----1'
          print*,kpds
       call baopen(33,"outgribmax",iiii)
      CALL putgb(33,kf,kpds,kgds,mask,grid1,
     +            iret)
       print*,'iret= ', iret

c-C     SURFACE 8HR AVE OZONE (ACTUALLY OZONE AT FIRST SIGMA MID-LAYER)
c-
c-      j = 0
c-      jpds = -1
c-      jpds(3) = igdnum
c-      jpds(5) = 180
c-      jpds(6) = 107
c-      jpds(7) = 10000
c-         jpds(14) = fhr-8
c-         if(jpds(14).eq.-1) jpds(14)=0
c-
c-      CALL getgb(lugb,lugi,numval,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,
c-     +            iret)
c-
c-      IF (iret.eq.0) THEN
c-        ii = 1
c-        jj = jj1
c-        DO kk = 1, itot
c-          oz8s(ii,jj) = grid(kk)
c-          grid1(kk)=oz8_max(ii,jj)
c-          IF (mask(kk)) THEN
c-            mskoz8(ii,jj) = 1
c-          END IF
c-          iprev = ii
c-          ii = ii + 1
c-          IF (mod(iprev,imax).eq.0) THEN
c-            ii = 1
c-            jj = jj + jjinc
c-            IF (jj.gt.jmax.or.jj.lt.1) GO TO 156
c-          END IF
c-        END DO
c-  156   CONTINUE
c-      ELSE
c-        WRITE (6,1300) iv, jpds(5), l, iret
c-      END IF
c-          print*,'------8',fhr,oz8s(1,1),lugb,lugi,grid(1),grid1(1)
c-
c-          jdate = (ivdate(1)*1000000+ivdate(2)*10000+ivdate(3)*100+
c-     1          ivdate(5))
c-          if (cycle.eq.'12Z') then
c-            if (day.eq.'DAY1') time = -24
c-            if (day.eq.'DAY2') time = -48
c-          endif
c-          if (cycle.eq.'06Z') then
c-            if (day.eq.'DAY1') time = -30
c-            if (day.eq.'DAY2') time = -54
c-          endif
c-          print*, jdate
c-          CALL raddate(jdate,time,bdate)
c-          print*, bdate
c-          idate=bdate
c-          mxyr=idate/1000000
c-          mxmo=mod((idate/10000),100)
c-          mxda=mod((idate/100),100)
c-          mxhr=mod(idate,100)
c-          if(mxyr.ge.100) then
c-             mxcen=mod((mxyr/100),100) + 1
c-             mxyr=mod(mxyr,100)
c-          endif
c-          print*,macen,mxyr,mxmo,mxda,mxhr,idate
c-
c-          kpds(21) = mxcen
c-          kpds(8)  = mxyr
c-          kpds(9)  = mxmo
c-          kpds(10) = mxda
c-          kpds(11) = mxhr
c-
c-          if (cycle.eq.'12Z') then
c-            if (day.eq.'DAY1') then
c-              kpds(15) = 24
c-              kpds(14) = 16
c-            elseif (day.eq.'DAY2') then
c-              kpds(15) = 48
c-              kpds(14) = 40
c-            endif  
c-          endif
c-          if (cycle.eq.'06Z') then
c-            if (day.eq.'DAY1') then
c-              kpds(15) = 30
c-              kpds(14) = 22
c-            elseif (day.eq.'DAY2') then
c-              kpds(15) = 54
c-              kpds(14) = 46
c-            endif  
c-          endif
c-
c-       print*,'-------putgb-----8'
c-      CALL putgb(33,kf,kpds,kgds,mask,grid1,
c-     +            iret)
c-       print*,'iret= ', iret

          RETURN
          END
