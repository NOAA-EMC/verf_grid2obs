C----------------------------------------------------------------------
C----------------------------------------------------------------------
      PROGRAM prepfits
C**********************************************************************
C	INPUTS:
C
C	UNIT		File description
C
C	 5	Records of model name  GRIB file / name  index file
C	11	LEVCAT namelist input file.  Example:
C                 &LEVCAT
C                 NUMLEV=39
C                 FIT=.T.,.T.,.F.,.F.,.T.,.T.,.T.,.F.,.F.,.F.
C                 &END
C	20	INPUT BUFR file created by editbufr
C	21	GRIB file for forecast
C	22	PREPFITS BUFR table file
C	23	GRIB index file for forecast
C	50	Output BUFR file with obs & forecast profiles
C***********************************************************************
C 
      use gridef
      use observ
      use backgv
      use guesfc
      use obstrs
      use debug
      use counts
      use vdates
      use guser
      use anl
      INCLUDE 'parm.inc'
 
      REAL*8 bdate
      real*8 hstid
      real*8 hi
      real*8 pob(255,mxr)
      integer fhours(mxb)

      character*80 file(mxb), fndx(mxb)
      character*80 obstrout,obstrout2,obstrout3,obstrout4
      character*8 sub(mxr),std(mxr)

      DIMENSION mate(mxr)
 
      character*8 cnf(2,255), sslast
      CHARACTER*10 src(mxb)
      character*10 source
      REAL*8 snf(2,255) 
      DIMENSION iges(mxb), irepv(mxb)
      EQUIVALENCE (cnf,snf)
      CHARACTER*8 staid
      EQUIVALENCE (hstid,staid)
      LOGICAL fit(0:9), valid, valix
      LOGICAL onlysf
c     REAL Q(200),P(200),T(200),PINT(201),Z(200),U(200),V(200)
c     REAL q1(100),p1(100),t1(100),z1(100),u1(100),v1(100),d1(100)
c     REAL q2(100),p2(100),t2(100),z2(100),u2(100),v2(100)
c     REAL a1(100),a2(100)
 
c     DATA headr /'SID XOB YOB DHR ELV TYP T29 ITP           '/
c     DATA obstr /'POB QOB TOB ZOB UOB VOB PMO CAPE CINH LI  '/
c     DATA qmstr /'PQM QQM TQM ZQM WQM                       '/
 
      DATA bmiss /10E10/
      DATA tzero /273.15/
      NAMELIST /levcat/ numlev, fit
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      target = 'KEOD'
C
c     NAMELIST /levcat/ numlev, fit
C
      headr='SID XOB YOB DHR ELV TYP T29 ITP           '
      obstr='POB QOB TOB ZOB UOB VOB PMO XDR YDR HRDR'
      obstr2='.REHOVI HOVI VTVI'
      obstr3='.DTHMXTM MXTM .DTHMITM MITM'
      obstr4='TDO'
c     obstr5='.DTHTOPC TOPC'
      obstr5='TOCC'
      obstr6='.DTMXGS MXGS'
c     obstr7='TP01 TP03 TP06 TP12 TP24'
      obstr7='PWO'
      obstr8='VSSO CLAM HOCB'
c     obstr8='CDBP'
c     obstr7='.DTHDOFS DOFS'
c     obstrout='POB QOB TOB ZOB UOB VOB PMO MXTM MITM TDO HOVI'
c     obstrout='POB QOB TOB ZOB UOB VOB PMO MXTM TDO'
c     obstrout='POB QOB TOB ZOB UOB VOB PMO MXTM'
c     obstrout='POB QOB TOB ZOB UOB VOB PMO PBL CAPE CINH LI TROP PWO'
      obstrout='POB QOB TOB ZOB UOB VOB PMO'
      obstrout2='MXTM MITM TDO HOVI TOCC MXGS THI TCH CDBZ CEIL'
      obstrout3=
     *  'CAPE CINH LI TROP PWO BCAPE SHR06 BCSVR SCAPE HAINES SCAPE2 
     *BCSVR2'
      print*,'obstrout3=',obstrout3
c  NEW PBL
      obstrout4='TKEPBL RIPBL TRANS VENT WIND80'
c  NEW PBL
      qmstr='PQM QQM TQM ZQM WQM                       '

 
      obsi=bmiss
      obsi2=bmiss
      obsi4=bmiss
      bak2=bmiss
      bak4=bmiss
      numlev = 6
C     categories to fit
C     0 - surface
      fit(0) = .true.
C     1 - mandatory level
      fit(1) = .true.
C     2 - sig level T & moisture
      fit(2) = .false.
C     3 - winds by pressure
      fit(3) = .false.
C     4 - winds by height
      fit(4) = .false.
C     5 - tropopause
      fit(5) = .false.
C     6 - ANY single level
      fit(6) = .false.
C     7 - auxiliary
      fit(7) = .true.
C     8 & 9 - reserved
      fit(8) = .false.
      fit(9) = .false.

      allocate(hdr(10))
      allocate(cat(255))
      allocate(qms(10,255))
      allocate(obs(10,255))
      allocate(obs7(10,255))
      allocate(obs2(10,255))
      allocate(obs3(10,255))
      allocate(obs4(10,255))
      allocate(obs5(10,255))
      allocate(obs6(10,255))
      allocate(obs8(3,5))
      allocate(obsi(10,255))
      allocate(obsi2(12,1))
      allocate(obsi4(5,2))
      allocate(bak(10,255))
      allocate(bak2(10,255))
      allocate(bak3(12,1))
      allocate(bak4(5,2))

      allocate(qtemp(255,mxr))

      allocate(xdft(3,255,mxr))
      allocate(xlev(mxr))

      allocate(fhr(mxb))
      allocate(vdate(mxb))
      allocate(nofo(mxb))

      allocate(xyz(2,mxr))

      qtemp=bmiss

      call getanl

      print*,'end of getanl'

      call datelen(10)

      iread=11
      iwrite=6
      rewind(iread)
      rewind(5)
      READ (iread,levcat,END=20)
      WRITE (iwrite,levcat)
      print*,'fit=',fit
      do i=1,9
       print*,'i,fit(i)=',i,fit(i)
      enddo
      print*,'fit(2),fit(7)=',fit(2),fit(7)

      ismart=0
      if(fit(8)) then
        ismart=1
        print*,'Smartinit run'
      endif

      ivirt=0
      if(fit(9)) then
        ivirt=1
        print*,'Virtual CAPE'
      else
        print*,'Sensible CAPE'
      endif
 
      lubfi = 20
      luges = 31
      lundx = 22
      lugbi = 61
      lubfo = 50

      kntgsf = 0
 
      PRINT *
      PRINT *, '******  BEGINNING PREPFITS PROCESSING ******'
      PRINT *
 
C     GET THE DATE FROM A PREPBUFR FILE AND A LIST OF LAT/LON'S
C     SUBROUTINE UFBXY3 RETURNS A LIST OF UNIQUE LON/LAT/LOCATIONS
C     ------------------------------------------------------------
 
      CALL datebf(lubfi,iy,im,id,ih,idate)
      PRINT '(" DATA VALID AT ",I4,3I2.2)', iy, im, id, ih
C     COMPUTE VALID DATE OF DATA ONCE AND FOR ALL !!!!
      dhr = 0.0
      CALL raddate(idate,dhr,vdata)
      PRINT *, 'VDATA=', vdata
C     
      IF (idate.ge.0) THEN
c       CALL ufbxy3(lubfi,xyz,mate,2,mxr,nrep)
        CALL ufbxy3(lubfi,mate,2,mxr) 
	PRINT *, 'Debug index = ', indux
        PRINT *, 'MAKING FITS FOR ', nrep, ' REPORTS'
c       PRINT *, 'MATE:', (mate(n),N=1,nrep)
      ELSE
        GO TO 30
      END IF
 
C     READ A LIST OF INPUT BACKGROUND & INDEX FILES TO PROCESS
C     --------------------------------------------------------
 
      DO n = 1, mxb+1
        READ (5,*,END=10) fhours(n), src(n), file(n)
        print*,'fhours, file=',fhours(n),file(n)
        source=trim(src(n))
        READ (5,*,END=10) fhours(n), src(n), fndx(n)
      END DO
  
      call errmsg('PREPFITS - TOO MANY BACKGROUND FIELDS')
      call errexit(1)
   10 nbak = n - 1
      PRINT *, 'NBAK=', nbak
c
cccccc 2/21/2006:  Read in the entire prepbufr file FIRST
c 2/22/2006:  Read in the prepbufr file's drift winds/time
c
C     OPEN THE INPUT FILE TO READ THROUGH
C     -----------------------------------
                                                                                          
      CALL openbf(lubfi,'IN ',lubfi)
      call maxout(100000)
      irep = 0
      irepo = 0
      irepv = 0
                                                                                          
C     READ THROUGH THE PREPBUFR MESSAGES
C     ----------------------------------
                                                                                          
      DO WHILE (ireadmg(lubfi,subset,idate).eq.0)
        nsub = nmsub(lubfi)
        PRINT *, 'SUBSET,NSUB,IDATE: ', subset, nsub, idate
                                                                                          
ccccC       READ A PREPBUFR REPORT
c 2/22/2006:  Read in the drift items here
C       ----------------------
                                                                                          
        sslast = subset
        DO WHILE (ireadpb(lubfi).eq.0)
          IF (subset.ne.sslast) THEN
            PRINT *, 'SUBSET CHANGED TO: ', subset,
     +                  ' BY READPB AT IREP', irep
            PRINT *, 'SUBSET,NSUB,IDATE: ', subset, nsub, idate
          END IF
          irep = irep + 1
21        continue
          IF (hdr(2).ne.xyz(1,irep).or.hdr(3).ne.xyz(2,irep)) THEN
            irep=irep-1
c           PRINT *, 'IREP: ob vs xyz', irep, hdr(2), xyz(1,irep),
c    +                  hdr(3), xyz(2,irep)
c           PRINT *, 'IREP-1 & IREP+1 xyz', xyz(1,irep-1),
c    +                  xyz(2,irep-1), xyz(1,irep+1), xyz(2,irep+1)
c           GO TO 20
            go to 21
          END IF
          hstid = hdr(1)
c
c Now, label the drift values by irep
c
      xlev(irep)=nlev
      sub(irep)=subset
      std(irep)=staid

      if(subset(:6).eq.'AIRSND') then
        if(obs(8,1).lt.0.) obs(8,1)=obs(8,1)+360.
      endif

      do i=1,nlev
        xdft(1,i,irep)=obs(8,i)
        xdft(2,i,irep)=obs(9,i)
        xdft(3,i,irep)=obs(10,i)
        obs(8,i)=bmiss
        obs(9,i)=bmiss
        obs(10,i)=bmiss
        pob(i,irep)=obs(1,i)
      enddo   !  i=1,nlev

      enddo   !  ireadpb
      enddo   !  ireadmg
 
C     MAKE SURE THE BACKGROUNDS ARE CHRONOLOGICAL AND TRANSFORM THEM
C     --------------------------------------------------------------
 
      DO ibak = 1, nbak
c       CLOSE (luges)
c       CLOSE (lugbi)
c       call baclose(luges,iret)
c       call baclose(lugbi,iret)
        luges=luges+1
        lugbi=lugbi+1
c       ishl = ishell('assign -a '//file(ibak)//' -s unblocked fort.21')
c       ishi = ishell('assign -a '//fndx(ibak)//' -s unblocked fort.23')
c       ishl = system('ln -s -f '//file(ibak)//'  fort.21')
c       ishi = system('ln -s -f '//fndx(ibak)//'  fort.23')
c       call system('ln -s -f '//file(ibak)//'  fort.21',ishl)
c       call system('ln -s -f '//fndx(ibak)//'  fort.23',ishi)
c       open(21,file=file(ibak),form='unformatted',iostat=ishl)
c       open(23,file=fndx(ibak),form='unformatted',iostat=ishi)
        call baopenr(luges,file(ibak),ishl)
        call baopenr(lugbi,fndx(ibak),ishi)
        print*,'luges,lugbi=',luges,lugbi
c       ishl=0
c       ishi=0
     
C       For each background, call GETBAK and then call GETPROF
C       To get a background profile at each ob location
C       -------------------------------------------------
        IF (ishl.eq.0.and.ishi.eq.0) THEN
          print*,'getbak for ',file(ibak)
          CALL getbak(luges,lugbi,numlev,
     +           iearth,ismart,ipcp,iges(ibak),fhours(ibak),source)
          print*,'between getbak and getprof'
          IF (iges(ibak).eq.0) CALL getprofupr(iearth,sub,pob,std)
        PRINT*,'DONE WITH IBAK=',IBAK,FILE(IBAK),IGES(IBAK)
        ELSE
          PRINT *
          PRINT *, 'MISSING: ', file(ibak)
          PRINT *
          iges(ibak) = -1
        END IF
C       
C       Here we should determine if time interpolation is desirable
C       BUT I'm WIRING it to NOT do ANY !!!!!
C       
        nofo(ibak) = 1
      END DO  
 
C     OPEN THE INPUT FILE TO READ THROUGH
C     -----------------------------------
      
      CALL closbf(lubfi)
      CALL openbf(lubfi,'IN ',lubfi)
      irep = 0
      irepo = 0
      irepv = 0

C     OPEN THE OUTPUT FILE (POSITION MOD IF THE FILE EXISTS)
C     ------------------------------------------------------
 
      CALL datebf(lubfo,iy,im,id,ih,idate)
      print*,'idate=',idate
      IF (idate.ge.0) CALL openbf(lubfo,'APN',lubfo)
      IF (idate.lt.0) CALL openbf(lubfo,'OUT',lundx)


C     READ THROUGH THE PREPBUFR MESSAGES
C     ----------------------------------
 
      DO WHILE (ireadmg(lubfi,subset,idate).eq.0)
        nsub = nmsub(lubfi)
        PRINT *, 'SUBSET,NSUB,IDATE: ', subset, nsub, idate
 
C       READ A PREPBUFR REPORT
C       ----------------------
  
        sslast = subset
        DO WHILE (ireadpb(lubfi).eq.0)
          IF (subset.ne.sslast) THEN
            PRINT *, 'SUBSET CHANGED TO: ', subset, 
     +                  ' BY READPB AT IREP', irep
            PRINT *, 'SUBSET,NSUB,IDATE: ', subset, nsub, idate
          END IF
          irep = irep + 1
22        continue
          IF (hdr(2).ne.xyz(1,irep).or.hdr(3).ne.xyz(2,irep)) THEN
            irep=irep-1
c           PRINT *, 'IREP: ob vs xyz', irep, hdr(2), xyz(1,irep), 
c    +                  hdr(3), xyz(2,irep)
c           PRINT *, 'IREP-1 & IREP+1 xyz', xyz(1,irep-1), 
c    +                  xyz(2,irep-1), xyz(1,irep+1), xyz(2,irep+1)
c           GO TO 20
            go to 22
          END IF
          hstid = hdr(1)
      if(subset(:6).eq.'AIRSND') then
        if(obs(8,1).lt.0.) obs(8,1)=obs(8,1)+360.
      endif
C-------------------------------------------------------------
C         READPB (ABOVE) REPLACES THE FOLLOWING BLOCK FOR READING DATA
C-------------------------------------------------------------
C         DO WHILE(IREADSB(LUBFI).EQ.0)
C-------------------------------------------------------------
C         CALL UFBINT(LUBFI,HDR,10,  1,IRET,HEADR)
C         CALL UFBINT(LUBFI,CAT, 1,255,NLEV,'CAT')
C         CALL UFBINT(LUBFI,QMS,10,255,NLEV,QMSTR)
C         CALL UFBINT(LUBFI,OBS,10,255,NLEV,OBSTR)
C------------------------------------------------------------
C         READPB RETURN A COMBINED MASS/WIND/REPORT IN THE OB ARRAYS
C------------------------------------------------------------
c         IF(CAT(1).GT.7) THEN
c         PRINT*,'IREP,HDR: ',IREP,STAID,(HDR(I),I=2,8)
c         PRINT*,'NLEV,CAT: ',NLEV,(CAT(N),N=1,NLEV)
c         PRINT*,'NLEV=1,OBS: ',(OBS(I,1),I=1,6)
c         PRINT*,'NLEV=1,QMS: ',(QMS(I,1),I=1,5)
c         ENDIF
c         DHR = HDR(4)
c         ADATE = IDATE
c         CALL RADDATE(ADATE,DHR,VDATA)
C         PRINT*,'IREP,VDATA: ',IREP,VDATA

c---- assign backgrond moisture to obs for ACARS ----
c         print*,'obs(2,1),qtemp',
c    *     obs(2,1),qtemp(1,irep),qtemp(2,irep),irep
          if(subset(:6).eq.'AIRSND') then
            do i=1,xlev(irep)
             obs(2,i)=qtemp(i,irep)
            enddo
c           print*,'obs(2,1),qtemp',
c    *        obs(2,1),qtemp(1,irep),qtemp(2,irep),irep
          endif
c         print*,'obs(2,1),qtemp',
c    *     obs(2,1),qtemp(1,irep),qtemp(2,irep),irep,subset
c         print*,'obs(2,2),qtemp',
c    *     obs(2,2),qtemp(1,irep),qtemp(2,irep),irep

      if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.'AIRSND'
     *    .or.subset(:6).eq.'AIRCAR'.or.subset(:6).eq.'AIRCFT') 
     *    then
       do i=1,nlev
         obs(8,i)=bmiss
         obs(9,i)=bmiss
         obs(10,i)=bmiss
c        obs(11,i)=bmiss
       enddo
      endif

      if(subset(:6).eq.'GPSIPW') then
       obs=bmiss
       obsi=bmiss
c      obsi(9,2)=obs7(1,nlev7)
c      obs(13,1)=obs7(1,nlev7)
       if(nlev7.gt.0) obsi2(5,2)=obs7(1,nlev7)
      endif

      if(subset(:6).eq.'ADPSFC'.or.subset(:6).eq.'MSONET') then
      obsi=bmiss
      do i=1,nlev8
       if(obs8(3,i).ne.bmiss.and.obs8(3,i).ne.0.0) then
        obsi(9,2)=obs8(3,i)
        goto 312
       endif
      enddo
312   continue
c
c Heat index calculation
c
      i1=0
      i2=0
      if(obs(3,1).ne.bmiss.and.obs(2,1).ne.bmiss.and.obs(1,1).ne.bmiss
     *  .and.obs(3,1).gt.26.0) then
        call htindex(obs(3,1),obs(1,1),obs(2,1),hi)
        i1=1
      endif
      if(obs(3,1).ne.bmiss.and.obs(5,1).ne.bmiss.and.obs(6,1).ne.bmiss
     *  .and.obs(3,1).lt.5.0) then
        call windchill(obs(3,1),obs(5,1),obs(6,1),wc)
        i2=1
      endif
      if(nlev3.gt.0) obsi(1,2)=obs3(2,nlev3)-tzero   ! TMAX
      if(nlev3.gt.0) obsi(2,2)=obs3(4,nlev3)-tzero   ! TMIN
      if(nlev4.gt.0) obsi(3,2)=obs4(1,nlev4)        ! DPT
      if(nlev2.gt.0) obsi(4,2)=obs2(2,nlev2)        ! VIS
c
c  The following code determines the ceiling in the observation file.  
c  We are looking at several levels of the atmosphere
c  Using the VSSO (obs8(1,i)), where 7=low, 8=middle, and 9=high
c  And CLAM (obs8(2,i)), where the percentage of cloud is given thru a
c  look-up table.  The value of HOCB (obs8(3,i)) is the ceiling
c  only if CLAM is > than 50%.  If not we go up a level.  PS 28 Jan 2015
c  
c     print*,'ceiling staid=',staid
      do i=1,5
c      print*,'obs8(1,i),obs8(2,i),obs8(3,i)=',
c    *    obs8(1,i),obs8(2,i),obs8(3,i)
c      if((obs8(2,i).gt.4.and.obs8(2,i).lt.9).and.
c    *   obs8(2,i).ne.bmiss.and.obs8(3,i).ne.bmiss.
c    *   and.(nint(obs8(1,i)).eq.7.or.nint(obs8(1,i)).eq.8.
c    *   .or.nint(obs(8,i)).eq.9)) then
       if(((obs8(2,i).gt.4.and.obs8(2,i).lt.9).or.
     *   nint(obs8(2,i)).eq.12).and.
     *   obs8(3,i).ne.bmiss)
     *    then
        obsi(10,2)=obs8(3,i)
        goto 15
       else
        obsi(10,2)=bmiss
       endif       
      enddo
15    continue
c     if(i.ne.6) then
c     print*,'ceiling staid,i=',staid,i
c     do j=1,i
c     print*,'j,obs8(1,j),obs8(2,j),obs8(3,j)=',
c    *    j,obs8(1,j),obs8(2,j),obs8(3,j)
c     enddo
c     print*,'obsi(9,2)=',obsi(9,2)  
c     endif
      if(obsi(4,2).gt.16090.0.and.obsi(4,2).lt.bmiss) then
        obsi(4,2)=16090.0
      endif
      if(nlev5.gt.0) obsi(5,2)=obs5(1,nlev5)
      if(nlev6.gt.0) obsi(6,2)=obs6(2,nlev6)
      if(i1.eq.1) then
        obsi(7,2)=hi
      else
        obsi(7,2)=bmiss 
      endif
      if(i2.eq.1) then
        obsi(8,2)=wc
      else
        obsi(8,2)=bmiss
      endif
      endif
      obs2=bmiss
      obs3=bmiss
      obs4=bmiss
      obs5=bmiss
      obs6=bmiss
      obs7=bmiss
      obs8=bmiss

      ipw=0
c     if(subset(:6).eq.'ADPUPA') then
      if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.'AIRSND'
     *          .or.subset(:6).eq.'PROFLR'.or.
     *          subset(:6).eq.'VADWND') then
       call obscapepw(ivirt,cape,cin,pli,pw,bcape)
       obsi=bmiss
c      obsi(9,2)=pw
c      do i=1,nlev
c       obs(13,i)=bmiss
c      enddo
c      obs(13,1)=pw
      obsi2(1,1)=cape
      obsi2(2,1)=cin
      obsi2(3,1)=pli
       obsi2(5,1)=pw
       obsi2(6,1)=bcape
c      bak(9,1)=cpi(irep,1)
c      bak(10,1)=cni(irep,1)
c      bak(11,1)=pxi(irep,1)

       if(subset(:6).eq.'AIRSND')then
       obs(4,:)=bmiss
       endif

       call obspbl(hpbl,trans,vent,w80,idate)
       obsi4(3,2)=trans
       obsi4(4,2)=vent
       if(w80.eq.0) w80=bmiss
       obsi4(5,2)=w80
       call calshear(1,shr)
c      obsi2(7,1)=shr
c      obsi2(8,1)=cape*shr
       bcsvrob=cape*shr
       if(obsi2(7,1).eq.bmiss) obsi2(8,1)=bmiss
       if(obsi2(1,1).eq.0.0) then
c        obsi2(1,1)=bmiss
         obsi2(7,1)=bmiss
         obsi2(8,1)=bmiss
         obsi2(9,1)=bmiss
       else
         obsi2(7,1)=shr
         obsi2(8,1)=cape*shr
         if(obsi2(7,1).eq.bmiss) obsi2(8,1)=bmiss
         obsi2(9,1)=obsi2(1,1)
         obsi2(11,1)=obsi2(1,1)
         obsi2(12,1)=obsi2(8,1)
       endif

c      bak(8,1)=pbli(irep,1)
c      bak3(1,2)=pbli(irep,1)
c      obsi2(1,2)=hpbl
       obsi4(1,2)=hpbl
       obsi4(2,2)=hpbl
c      obsi4(3,2)=hpbl
c      
       do i=1,nlev
        if(cat(i).eq.5.0) then
          obsi2(5,2)=obs(1,i)
          goto 123
        endif
       enddo
c      do i=2,nlev
c        obs(12,i)=bmiss
c      enddo
123    continue

       call calhaines(hindex,subset)
       obsi2(10,1)=hindex
c
c
      endif    !  if ADPUPA 

C         FILTER THE DATA TO FIT BY CATEGORY
C         ----------------------------------
 
          mlev = 0
          DO l = 1, nlev
            kat = nint(cat(l))
            IF (kat.gt.9.or.kat.lt.0) kat = 9
c           if(staid.eq.'71111') print*,'71111 fit=',kat,fit(kat)
            IF (fit(kat)) THEN
              mlev = mlev + 1
              cat(mlev) = cat(l)
              IF (kat.eq.6.and.(subset(:6).eq.'ADPSFC'.or.subset(:6).eq.
     +         'SFCSHP'.or.subset(:6).eq.'MSONET')) cat(mlev) = 0
              DO i = 1, 7
                obs(i,mlev) = obs(i,l)
              enddo
              do i = 1,5
                qms(i,mlev) = qms(i,l)
c               if(qms(i,mlev).eq.bmiss) qms(i,mlev)=9.0
              END DO
            END IF
          END DO
c          if(staid.eq.'71111') print*,'mlev=',mlev
           IF ( mlev .eq. 1 .and. subset (:6) .eq. 'ADPSFC' .or.
     +	   subset(:6).eq.'MSONET'.or.subset (:6) .eq. 'SFCSHP'.or.
     *     subset(:6).eq.'GPSIPW' ) THEN
	      mlev = 2
	      cat (mlev) = cat (1)
              DO i = 1, 10
                obs(i,mlev) = obs(i,1)
              enddo
              do i = 1,5
                qms(i,mlev) = qms(i,1)
              END DO
	      onlysf = .true.
	   ELSE
	      onlysf = .false.
	   END IF
	   
C----------------------------------------------------------------------
C         IF(MOD(IREP,350).EQ.0 .OR. MLEV.EQ.0) THEN
c         IF (mlev.eq.0) THEN
c           PRINT *, 'IREP,NLEV,MLEV: ', irep, nlev, mlev
c           PRINT *, 'IRET,HDR: ', iret, staid, (hdr(i),I=2,8)
c           PRINT *, 'NLEV,CAT: ', nlev, (cat(n),N=1,nlev)
c           PRINT *, 'NLEV=1,OBS: ', (obs(i,1),I=1,7)
c           PRINT *, 'NLEV=1,QMS: ', (qms(i,1),I=1,5)
c         END IF
 
          nlev = mlev
          print*,'are we here?'

C         MAKE SURE THERE IS DATA TO WRITE OUT
          IF (nlev.gt.0) THEN
            irepo = irepo + 1
C           MAKE SURE A MESSAGE OF THE PROPER TYPE AND DATE IS OPEN FOR OUTPUT
C           USING THE CURRENT SUBSET NAME (SSLAST) IN CASE IREADPB CHANGED IT
C           ------------------------------------------------------------------

            call datelen(10)
            print*,'openmb'
            CALL openmb(lubfo,sslast,idate)

C           WRITE OUT THE REGISTRATION DATA FOR EACH REPORT
C           -----------------------------------------------
 
            print*,'ufbint'
            CALL ufbint(lubfo,hdr,10,1,iret,headr)
 
C           WRITE OUT THE REGISTRATION DATA FOR EACH LEVEL
C           ----------------------------------------------
 
            CALL ufbint(lubfo,cat,1,nlev,iret,'CAT')
            CALL ufbint(lubfo,obs,10,nlev,iret,'PRC')
            CALL ufbint(lubfo,qms,10,nlev,iret,qmstr)

C           WRITE THE OBSERVED VALUES FOR EACH LEVEL
C           ----------------------------------------
 
            DO l = 1, nlev
              cnf(1,l) = 'PRP'
              snf(2,l) = bmiss
            END DO
                                                                                             
c          if(subset(:6).eq.'MSONET'.or.subset(:6).eq.'ADPSFC'
c    *         .or.subset(:6).eq.'GPSIPW') then
c           print*,'subset,staid=',subset,staid
c           do i=1,6
c           print*,'i,obsi2(i,2)=',i,obsi2(i,2)
c           enddo
c           endif
            
c           if(subset(:6).eq.'ADPUPA') then
c            print*,'subset,staid=',subset,staid
c            do j=1,7
c            do i=1,nlev
c             print*,'i,j,obs(j,i)=',i,j,obs(j,i),qms(j,i)
c             print*,'i,obs(1,i),obs(2,i),obs(3,i),obs(4,i),obs(5,i)=',
c    *         i,obs(1,i),obs(2,i),obs(3,i),obs(4,i),obs(5,i)
c            enddo
c            enddo
c           endif

c           if(subset(:6).eq.'ADPUPA') then
c           if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.'AIRSND')
c    &          then
c            print*,'staid,subset=',staid,subset,fhour
c            do i=1,12
c            obsi2(i,1)=obsi2(i,2)
c            print*,'i,obsi2,bak3=',i,obsi2(i,1),bak3(i,1)
c            enddo
c           endif
c           if(subset(:6).eq.'ADPUPA') then
            if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.'AIRSND') then
             do i=1,5
             obsi4(i,1)=obsi4(i,2)
             enddo
            endif
c 

            CALL ufbint(lubfo,snf,2,nlev,iret,'SRC FHR')
            CALL ufbint(lubfo,obs,10,nlev,iret,obstrout)
            if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq. 
     *         'GPSIPW'.or.subset(:6).eq.'AIRSND') then
            call ufbint(lubfo,obsi2,12,1,iret,obstrout3)
            endif
            if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.'AIRSND')
     &          then
            call ufbint(lubfo,obsi4,5,1,iret,obstrout4)
            endif
            if(subset(:6).eq.'MSONET'.or.subset(:6).eq.'ADPSFC'.
     *        or.subset(:6).eq.'GPSIPW') then
             call ufbint(lubfo,obsi,10,nlev,iret,obstrout2)
            endif

C           WRITE OUT THE INTERPOLATED BACKGROUND(S) INTO THE FIT FILE
C           ----------------------------------------------------------
            valix = .false.
            DO n = 1, nbak
c             if(staid.eq.'71111') print*,'n,nbak=',n,nbak
              IF (n.lt.nbak.and.nofo(n).ne.1) THEN
                icdate = idint(vdate(n))
                CALL raddate(icdate,abs(fhr(n)-fhr(n+1)),bdate)
                valid = vdate(n) .le. vdata .and. vdate(n+1) .ge. vdata
     +                      .and. iges(n) .eq. 0 .and. iges(n+1) .eq. 0
     +                      .and. vdate(n+1) .eq. bdate .and. .not. 
     +                      valix
c               if(staid.eq.'71111') then
c                print*,'vdate(n)=',vdate(n)
c                print*,'vdata=',vdata
c                print*,'vdate(n+1)=',vdate(n+1)
c                print*,'iges(n)=',iges(n)
c                print*,'iges(n+1)=',iges(n+1)
c                print*,'bdate=',bdate
c                print*,'valix=',valix
c                print*,'valid=',valid
c               endif
                valix = valid
              ELSE IF (n.eq.nbak.or.nofo(n).eq.1) THEN
                valid = vdate(n) .eq. vdata .and. iges(n) .eq. 0
c               if (staid.eq.'71111') then
c                print*,'vdate(n),vdata,iges(n)=',vdate(n),vdata,iges(n)
c                print*,'valid=',valid
c               endif
              END IF
              IF (valid) THEN
                irepv(n) = irepv(n) + 1
                CALL getfct(irep,n,fhour,onlysf,subset,source)

c     if(subset(:6).eq.'ADPUPA') then
      if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.'AIRSND')
     &          then
        call calshear(2,shr6)
        bak3(7,1)=shr6
        bak3(8,1)=shr6*bak3(1,1)
        bak3(9,1)=bak3(1,1)
        if(bak3(1,1).eq.0.0.and.obsi2(1,1).eq.0.0) then
          bak3(11,1)=bmiss
          bak3(12,1)=bmiss
        else
c         if(fhour.eq.0.0) print*,'staid,bak3(1,1),obsi2(1,1)=',
c    *     staid,bak3(1,1),obsi2(1,1)
          bak3(11,1)=bak3(1,1)
          bak3(12,1)=bak3(8,1)
        endif
        if(bak3(7,1).eq.bmiss) bak3(8,1)=bmiss
      endif
c
c Ventilation rate:  transport wind times the PBL height
c
      bak4(4,2)=bak4(3,2)*bak4(2,2)

      if(subset(:6).eq.'ADPSFC'.or.subset(:6).eq.'MSONET') then
      i1=0
      i2=0
      if(bak(3,1).ne.bmiss.and.bak(2,1).ne.bmiss.and.bak(1,1).ne.bmiss.
     *  and.bak(3,1).gt.26.0) then
        call htindex(bak(3,1),bak(1,1),bak(2,1),hi)
        i1=1
      endif
      if(bak(3,1).ne.bmiss.and.bak(5,1).ne.bmiss.and.bak(6,1).ne.bmiss
     *  .and.bak(3,1).lt.5.0) then
        call windchill(bak(3,1),bak(5,1),bak(6,1),wc)
        i2=1
       endif
       if(i1.eq.1) then
         bak2(7,2)=hi
       else 
         bak2(7,2)=bmiss
       endif
       if(i2.eq.1) then
         bak2(8,2)=wc
       else
         bak2(8,2)=bmiss
       endif
       endif
                if (subset(:6).ne.'ADPUPA') then
                  do iii=8,10
                    bak(iii,1)=bmiss
                  enddo
                endif
                DO l = 1, nlev
                  cnf(1,l) = src(n)
                  snf(2,l) = fhours(n)
                END DO

c       print 1555, HPBL,PBLI(IREP,1),jpbl,(p1(kk)/100,kk=1,jpbl)
c       if(hpbl.gt.0.)
c    +    write(100,1556) hdr(3), hdr(2), hpbl, staid
c       write(101)(p1(kk)/100,kk=1,100),
c    +            (z1(kk),kk=1,100),
c    +            (t1(kk),kk=1,100),
c    +            (d1(kk),kk=1,100),
c    +            (q1(kk),kk=1,100),
c    +            (u1(kk),kk=1,100),
c    +            (v1(kk),kk=1,100)
1555    format('OB PBL=',f5.0,' MDL PBL=',f5.0,i4,15f6.0)
1556    format(f6.2,1x,f7.2,2x,f5.0,3x,a8)

                CALL ufbint(lubfo,snf,2,nlev,iret,'SRC FHR')
      
c          if(subset(:6).eq.'MSONET'.or.subset(:6).eq.'ADPSFC'
c    *         .or.subset(:6).eq.'GPSIPW') then
c           print*,'subset,staid=',subset,staid
c           do i=1,6
c           print*,'i,bak3(i,2)=',i,bak3(i,2)
c           enddo
c           endif

c           if(subset(:6).eq.'ADPSFC') then
c            print*,'subset,staid,nlev=',subset,staid,nlev
c            do i=7,7
c            do i=1,nlev
c             bak(12,i)=bmiss
c             print*,'i,j,bak(j,1)=',i,j,bak(j,1)
c             print*,'i,bak(1,i),bak(2,i),bak(3,i),bak(4,i),bak(5,i)=',
c    *         i,bak(1,i),bak(2,i),bak(3,i),bak(4,i),bak(5,i)
c             print*,'i,bak(3,i),obs(3,i)=',i,bak(3,i),obs(3,i)
c             if(obs(1,i).eq.500.0.and.n.eq.1) then
c            print*,'subset,staid,nlev=',subset,staid,nlev
c              print*,'obs(3,i),obs(5,i),bak(3,i),bak(5,i)=',
c    *          obs(3,i),obs(5,i),bak(3,i),bak(5,i)
c            if(qms(1,1).ne.bmiss) print*,'i,obs(i,1),bak(i,1)=',
c    *         i,obs(i,1),bak(i,1),qms(i,1)
c            if(qms(1,2).ne.bmiss) print*,'i,obs(i,2),bak(i,2)=',
c    *         i,obs(i,2),bak(i,2),qms(i,2)
c             endif
c            enddo
c            enddo
c           endif

c           if(subset(:6).eq.'ADPUPA') then
c           if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.'AIRSND')
c    &         then
c            print*,'staid,subset,fhour=',staid,subset,fhours(n)
c            do i=1,12
c            bak3(i,1)=bak3(i,2)
c            print*,'i,obsi2,bak3=',i,obsi2(i,1),bak3(i,1)
c            enddo
c           endif

            if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.'AIRSND')
     &        then
             do i=1,5
             bak4(i,1)=bak4(i,2)
             enddo
c      if(obsi4(1,1).ne.bmiss) then
c        print 1557, hdr(3),hdr(2),hdr(5),obsi4(1,2),
c    *   bak4(1,1),bak4(2,1),bak4(3,1)
c        write(102,1557) hdr(3),hdr(2),hdr(5),obsi4(1,2),
c    *   bak4(1,1),bak4(2,1),bak4(3,1),nlev
c      endif
c1557    format(2f8.2,' elv=',f8.2,' obspbl=',g12.5,' TKEPBL=',g12.5,
c    *  ' RIPBL=',g12.5,' MIXHT=',g12.5,' nlev=',i3)
            endif
c            
             CALL ufbint(lubfo,bak,10,nlev,iret,obstrout)
            if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.
     *         'GPSIPW'.or.subset(:6).eq.'AIRSND') then
            call ufbint(lubfo,bak3,12,1,iret,obstrout3)
            endif
            if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.'AIRSND')
     &          then
             call ufbint(lubfo,bak4,5,1,iret,obstrout4)
            endif
             if(subset(:6).eq.'ADPSFC'.or.subset(:6).eq.'MSONET'.
     *         or.subset(:6).eq.'GPSIPW') then
                do i=1,9
                 if(obsi(i,2).eq.bmiss) bak2(i,2)=bmiss
                enddo
                call ufbint(lubfo,bak2,10,nlev,iret,obstrout2)
              endif
              ELSE
                PRINT *, 'VALID IS FALSE!!  IREP,N,IGES(N)', irep, n, 
     +                      iges(n)
                PRINT *, 'VALID IS FALSE!!  VDATE(N),VDATA', vdate(n), 
     +                      vdata
              END IF
            END DO
 
C           BEFORE END OF READ LOOPS - WRITE THE PREPFITS SUBSET
C           ----------------------------------------------------

            CALL writsb(lubfo)
          END IF
          sslast = subset
 
C       END OF READ LOOPS
        END DO
      END DO
 
C     CLOSE THE BUFR FILES
C     --------------------
 
      CALL closbf(lubfi)
      CALL closbf(lubfo)
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
   30 call errmsg('PREPFITS - BAD OR MISSING PREPBUFR FILE')
      call errexit(0)
      END

