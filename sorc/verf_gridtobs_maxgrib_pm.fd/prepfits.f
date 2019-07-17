c----------------------------------------------------------------------
C----------------------------------------------------------------------
      PROGRAM prepfits
C 
      INCLUDE 'parm.inc'
 
      REAL*8 bak(10,255)
      REAL*8 vdata,vdate,bdate

      CHARACTER*20 file1(mxb), fndx1(mxb), file8(mxb), fndx8(mxb)
      COMMON /vdates/ vdata, vdate(mxb), fhr(mxb), nofo(mxb)
 
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

      COMMON /days/ day, cycle

      DIMENSION oz1_all(24,ilim,jlim), oz1_tmp(24),
     +          oz8_all(24,ilim,jlim), oz8_tmp(24)

      CHARACTER*8 src(mxb)
      DIMENSION iges(mxb), irepv(mxb)
      CHARACTER*4 day
      CHARACTER*3 cycle
 
      DATA bmiss /10E10/
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C
C

      iread=11
      iwrite=6
 
      luges = 21
      lugbi = 23


      PRINT *
      PRINT *, '******  BEGINNING PREPFITS PROCESSING ******'
      PRINT *
 
      oz1_max=0.
      oz8_max=0.

C     READ A LIST OF INPUT BACKGROUND & INDEX FILES TO PROCESS
C     --------------------------------------------------------
 
      DO n = 1, mxb
        READ (5,'(A8,1X,A20,1X,A20)',END=10) src(n),file1(n)
        READ (5,'(A8,1X,A20,1X,A20)',END=10) src(n),fndx1(n)
      END DO
  
c=== for 12z run ====
c      if (file1(9)(5:6).eq.'01'.or.file1(9)(7:8).eq.'01') then
c        day = 'DAY1'
c      elseif (file1(9)(5:6).eq.'25'.or.file1(9)(7:8).eq.'25') then
c        day = 'DAY2'
c      else
c        print*,'=== CHECK INPUT LIST: prepfits.in* ====='
c      endif
c=== for both 06z and 12z runs ====
      if (file1(1)(7:8).eq.'17') then
c- 12z run
        cycle = '12Z'
        if (file1(9)(5:6).eq.'01'.or.file1(9)(7:8).eq.'01') then
          day = 'DAY1'
        elseif (file1(9)(5:6).eq.'25'.or.file1(9)(7:8).eq.'25') then
          day = 'DAY2'
        else
          print*,'=== CHECK INPUT LIST: prepfits.in* ====='
        endif
      endif
      if (file1(1)(7:8).eq.'23') then
c- 06z run
        cycle = '06Z'
        if (file1(3)(5:6).eq.'01'.or.file1(3)(7:8).eq.'01') then
          day = 'DAY1'
        elseif (file1(3)(5:6).eq.'25'.or.file1(3)(7:8).eq.'25') then
          day = 'DAY2'
        else
          print*,'=== CHECK INPUT LIST: prepfits.in* ====='
        endif
      endif
   
      print*,'== day  ',day

   10 nbak = n - 1
      PRINT *, 'NBAK=', nbak
 
C     MAKE SURE THE BACKGROUNDS ARE CHRONOLOGICAL AND TRANSFORM THEM
C     --------------------------------------------------------------
 
      do ifile=1,1

      DO ibak = 1, nbak
        CLOSE (luges)
        CLOSE (lugbi)
        if(ifile.eq.1) call baopen(21,file1(ibak),ishl)
        if(ifile.eq.1) call baopen(23,fndx1(ibak),ishi)
        if(ifile.eq.2) call baopen(21,file8(ibak),ishl)
        if(ifile.eq.2) call baopen(23,fndx8(ibak),ishi)
     
C       For each background, call GETBAK 
C       -------------------------------------------------
        IF (ishl.eq.0.and.ishi.eq.0) THEN
          CALL getbak(luges,lugbi,numlev,vdate(ibak),fhr(ibak),
     +                iges(ibak))

c       print*,'========ozs= ',oz1s(1,1),oz8s(1,1)
        do ii=1,ilim
        do jj=1,jlim
           if(ifile.eq.1) oz1_all(ibak,ii,jj)=oz1s(ii,jj)
           if(ifile.eq.2) oz8_all(ibak,ii,jj)=oz8s(ii,jj)
        enddo
        enddo

        PRINT*,'DONE WITH IBAK=',IBAK,FILE1(IBAK),IGES(IBAK)
        ELSE
          PRINT *
          PRINT *, 'MISSING: ', file1(ibak)
          PRINT *
          iges(ibak) = -1
        END IF
C       
C       
        nofo(ibak) = 1
      END DO

      enddo !end ifile loop
 
        do ii=1,ilim
        do jj=1,jlim
         do ihr=1,24
           oz1_tmp(ihr)=oz1_all(ihr,ii,jj)
           oz8_tmp(ihr)=oz8_all(ihr,ii,jj)
         enddo
           imaxhour=maxloc(oz1_tmp,dim=1,mask=oz1_tmp.ne.bmiss)
           oz1_max(ii,jj)=oz1_all(imaxhour,ii,jj)
c     print*,ii,jj,oz1_all(1,ii,jj),oz1_all(2,ii,jj),oz1_max(ii,jj),
c    +         imaxhour
c     print*,ii,jj,oz1_tmp(1),oz1_tmp(2),oz1_max(ii,jj),imaxhour
           imaxhour=maxloc(oz8_tmp,dim=1,mask=oz8_tmp.ne.bmiss)
           oz8_max(ii,jj)=oz8_all(imaxhour,ii,jj)
c     print*,ii,jj,oz8_all(1,ii,jj),oz8_all(2,ii,jj),oz8_max(ii,jj),
c    +         imaxhour
c     print*,ii,jj,oz8_tmp(1),oz8_tmp(2),oz8_max(ii,jj),imaxhour
        enddo
        enddo
 
C  Call getbak again to write the max grib file
C  --------------------------------------------
       ibak = 8 ! assign 12z to max value file
        CLOSE (luges)
        CLOSE (lugbi)
        call baopen(21,file1(ibak),ishl)
        call baopen(23,fndx1(ibak),ishi)

C       For each background, call GETBAK
C       -------------------------------------------------
        IF (ishl.eq.0.and.ishi.eq.0) THEN
          CALL getbak(luges,lugbi,numlev,vdate(ibak),fhr(ibak),
     +                iges(ibak))
        END IF

      PRINT *, '************************************'
      PRINT *, 'WROTE OUT FITS FOR ', irepo, ' REPORTS'
      PRINT *, '************************************'
      PRINT *, (irepv(i),I=1,nbak)
      PRINT *, '************************************'
      PRINT *, ' # of non-missing frcst reports for ONLYSF = ',
     +         kntgsf
      PRINT *, '************************************'
 
C     END OF PREPFITTING
C     ------------------
 
   20 CONTINUE
      STOP
      CALL abort
      END
