c----------------------------------------------------------------------
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
      INCLUDE 'parm.inc'
 
      REAL*8 hdr(10), cat(255), obs(10,255), qms(10,255)
      REAL*8 bak(10,255)
      REAL*8 vdata,vdate,bdate

      COMMON /observ/ hdr, cat, obs, qms, nlev
      COMMON /obstrs/ headr, obstr, qmstr, subset, idate, nsub
      COMMON /backgv/ bak, nbak
      COMMON /counts/ kntgsf

      CHARACTER*8 target
      integer fhours(mxb)
      COMMON /debug/ target, indux
      
      CHARACTER*80 headr, obstr, qmstr, file(mxb), fndx(mxb)
      DIMENSION mate(mxr)
      COMMON /guesfc/ psi(mxr,mxb), zsi(mxr,mxb), tsi(mxr,mxb),
     +            usi(mxr,mxb), vsi(mxr,mxb), qsi(mxr,mxb),
     +            pmi(mxr,mxb), cpi(mxr,mxb), cni(mxr,mxb),
     +            pxi(mxr,mxb), oz1i(mxr,mxb), oz8i(mxr,mxb)
      real*8 psi,zsi,tsi,usi,vsi,qsi,pmi,cpi,cni,pxi,oz1i,oz8i
      REAL*8 xyz(2,mxr)
      COMMON /guser/ xyz, nrep, ibak
      COMMON /vdates/ vdata, vdate(mxb), fhr(mxb), nofo(mxb)
 
      CHARACTER*8 subset, cnf(2,255), sslast
      CHARACTER*8 src(mxb)
      REAL*8 snf(2,255) 
      DIMENSION iges(mxb), irepv(mxb)
      EQUIVALENCE (cnf,snf)
      CHARACTER*8 staid
c     EQUIVALENCE (hstid,staid)
      EQUIVALENCE (hdr(1),staid)
      LOGICAL fit(0:9), valid, valix
      LOGICAL onlysf
      REAL Q(100),P(100),T(100),PINT(101)
      REAL q1(100),p1(100),t1(100)
c     CHARACTER*6 OZONE_AVE
 
c     DATA headr /'SID XOB YOB DHR ELV TYP T29 ITP           '/
c     DATA obstr /'POB QOB TOB ZOB UOB VOB PMO CAPE CINH LI  '/
c     DATA qmstr /'PQM QQM TQM ZQM WQM                       '/
 
      DATA bmiss /10E10/
      NAMELIST /levcat/ numlev, fit
 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      target = 'KEOD'
C
c     NAMELIST /levcat/ numlev, fit
C
      headr='SID XOB YOB DHR ELV TYP T29 ITP           '
c     obstr='POB QOB TOB ZOB UOB VOB PMO CAPE CINH LI  '
      obstr='TPHR QCIND COPO                           '
      qmstr='PQM QQM TQM ZQM WQM                       '

      IRET1=ISETPRM('MXLCC', 12)
      IRET2=ISETPRM('MAXMEM', 75000000)
      IRET3=ISETPRM('MXCDV', 50000)
      IRET4=ISETPRM('MAXJL', 128000)
      IRET5=ISETPRM('MAXSS', 150000)
      IRET6=ISETPRM('MXMSGL', 2500000)
      IRET7=ISETPRM('MXRST', 500)

      call openbf (1, "FIRST", 2)

      numlev = 19
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
      fit(6) = .true.
C     7 - auxiliary
      fit(7) = .true.
C     8 & 9 - reserved
      fit(8) = .false.
      fit(9) = .false.

      call datelen(10)

      iread=11
      iwrite=6
      READ (iread,levcat,END=20)
      WRITE (iwrite,levcat)
 
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
        CALL ufbxy3(lubfi,xyz,mate,2,mxr,nrep)
	PRINT *, 'Debug index = ', indux
        PRINT *, 'MAKING FITS FOR ', nrep, ' REPORTS'
c       PRINT *, 'MATE:', (mate(n),N=1,nrep)
      ELSE
        GO TO 30
      END IF
 
C     READ OZONE AVERAGING PERIOD
C     --------------------------------------------------------
c     READ (5,'(A6)') ozone_ave

C     READ A LIST OF INPUT BACKGROUND & INDEX FILES TO PROCESS
C     --------------------------------------------------------
 
      print*,'mxb=',mxb
      DO n = 1, mxb+1
c       READ (5,'(A8,1X,A50)',END=10) src(n), file(n)
c       READ (5,'(A8,1X,A50)',END=10) src(n), fndx(n)
        READ (5,*,END=10) fhours(n), src(n), file(n)
        READ (5,*,END=10) fhours(n), src(n), fndx(n)
        print*,'fhours,src(n),file(n)=',fhours(n),src(n),file(n)
      END DO
  
      CALL errmsg("PREPFITS - TOO MANY BACKGROUND FIELDS")
      CALL ERREXIT(20)
   10 nbak = n - 1
      PRINT *, 'NBAK=', nbak
 
C     MAKE SURE THE BACKGROUNDS ARE CHRONOLOGICAL AND TRANSFORM THEM
C     --------------------------------------------------------------
 
      DO ibak = 1, nbak
c       CLOSE (luges)
c       CLOSE (lugbi)
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
        call baopen(luges,file(ibak),ishl)
        call baopen(lugbi,fndx(ibak),ishi)
c       ishl=0
c       ishi=0
     
C       For each background, call GETBAK and then call GETPROF
C       To get a background profile at each ob location
C       -------------------------------------------------
        IF (ishl.eq.0.and.ishi.eq.0) THEN
          CALL getbak(luges,lugbi,numlev,vdate(ibak),fhours(ibak),
     +                src(ibak),iges(ibak))
          IF (iges(ibak).eq.0) CALL getprof
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

      CALL openbf(lubfi,'IN ',lubfi)
      irep = 0
      irepo = 0
      irepv = 0

C     OPEN THE OUTPUT FILE (POSITION MOD IF THE FILE EXISTS)
C     ------------------------------------------------------
 
      CALL datebf(lubfo,iy,im,id,ih,idate)
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
          IF (hdr(2).ne.xyz(1,irep).or.hdr(3).ne.xyz(2,irep)) THEN
            PRINT *, 'IREP: ob vs xyz', irep, hdr(2), xyz(1,irep), 
     +                  hdr(3), xyz(2,irep)
            PRINT *, 'IREP-1 & IREP+1 xyz', xyz(1,irep-1), 
     +                  xyz(2,irep-1), xyz(1,irep+1), xyz(2,irep+1)
            GO TO 20
          END IF
          hstid = hdr(1)
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
c         PRINT*,'NLEV=1,OBS: ',(OBS(I,1),I=1,10)
c         PRINT*,'NLEV=1,QMS: ',(QMS(I,1),I=1,5)
c         ENDIF
c         DHR = HDR(4)
c         ADATE = IDATE
c         CALL RADDATE(ADATE,DHR,VDATA)
c         PRINT*,'IREP,VDATA: ',IREP,VDATA

      ILEV=0
      do k=1,100
       p(k)=0.
       p1(k)=0.
       q(k)=0.
       q1(k)=0.
       t(k)=0.
       t1(k)=0.
      enddo
      slindx=0.
      cape=0.
      cin=0.
      peql=0.
      plcl=0.
      ifcape=0.
      do iii=8,10
        bak(iii,1)=bmiss
        obs(iii,1)=bmiss
      enddo
      if(subset(:6).eq.'ADPUPA') then
      DO K=1,100
      if(obs(2,k).ne.BMISS.and.
     *    obs(3,k).ne.BMISS) then
        ilev=ilev+1
        P1(ilev)=OBS(1,k)*100.
        Q1(ilev)=OBS(2,k)/1000000.
        T1(ilev)=OBS(3,k)+273.15
      endif
      ENDDO
      do k=1,ilev
       p(k)=p1(ilev-k+1)
       q(k)=q1(ilev-k+1)
       t(k)=t1(ilev-k+1)
      enddo
      if(ilev.gt.1) then
c     print*,'Entering CALCAPE'
      CALL CALCAPE(T,Q,P,pint,ILEV,1,1,ILEV,
     *   CAPE,CIN,PLCL,PEQL,PLI)
c     print*,'irep=',irep
c     print*,'cpi(irep,1),cni(irep,1),pxi(irep,1)=',
c    * cpi(irep,1),cni(irep,1),pxi(irep,1)
c     if(cpi(irep,1).ne.100000000000..and.
c    *   cni(irep,1).ne.100000000000..and.
c    *   pxi(irep,1).ne.100000000000.) then
      if(cpi(irep,1).ne.BMISS.and.
     *   cni(irep,1).ne.BMISS.and.
     *   pxi(irep,1).ne.BMISS) then
      print*,'OB CAPE,CIN,SLINDX=',CAPE,CIN,PLI
      print*,'IREP=',IREP
c     print*,'HDR(2),HDR(3)=',HDR(2),HDR(3)
      print*,'MDL CAPE,CIN,SLINDX=',CPI(IREP,1),CNI(IREP,1),
     *    PXI(IREP,1)
      obs(8,1)=cape
      obs(9,1)=cin
      obs(10,1)=pli
      bak(8,1)=cpi(irep,1)
      bak(9,1)=cni(irep,1)
      bak(10,1)=pxi(irep,1)
      endif
      endif
      endif


C         FILTER THE DATA TO FIT BY CATEGORY
C         ----------------------------------
 
          mlev = 0
          DO l = 1, nlev
            kat = nint(cat(l))
            IF (kat.gt.9.or.kat.lt.0) kat = 9
            IF (fit(kat)) THEN
              mlev = mlev + 1
              cat(mlev) = cat(l)
              IF (kat.eq.6.and.(subset(:6).eq.'ADPSFC'.or.subset(:6).eq.
c    +                    'SFCSHP'.or.subset(:6).eq.'AIRNOW')) 
     +                    'SFCSHP')) 
     +                     cat(mlev) = 0
              DO i = 1, 10
                obs(i,mlev) = obs(i,l)
                qms(i,mlev) = qms(i,l)
              END DO
            END IF
          END DO
	  IF ( mlev .eq. 1 .and. subset (:6) .eq. 'ADPSFC' .or.
c    +   subset (:6) .eq. 'SFCSHP'.or.subset(:6).eq.'AIRNOW' ) THEN
     +   subset (:6) .eq. 'SFCSHP') THEN
	      mlev = 2
	      cat (mlev) = cat (1)
              DO i = 1, 10
                obs(i,mlev) = obs(i,1)
                qms(i,mlev) = qms(i,1)
              END DO
	      onlysf = .true.
	   ELSE
	      onlysf = .false.
	   END IF
	   
C----------------------------------------------------------------------
C         IF(MOD(IREP,350).EQ.0 .OR. MLEV.EQ.0) THEN
          IF (mlev.eq.0) THEN
            PRINT *, 'IREP,NLEV,MLEV: ', irep, nlev, mlev
            PRINT *, 'IRET,HDR: ', iret, staid, (hdr(i),I=2,8)
            PRINT *, 'NLEV,CAT: ', nlev, (cat(n),N=1,nlev)
            PRINT *, 'NLEV=1,OBS: ', (obs(i,1),I=1,7)
            PRINT *, 'NLEV=1,QMS: ', (qms(i,1),I=1,5)
          END IF
 
          nlev = mlev

C         MAKE SURE THERE IS DATA TO WRITE OUT
          IF (nlev.gt.0) THEN
            irepo = irepo + 1
C           MAKE SURE A MESSAGE OF THE PROPER TYPE AND DATE IS OPEN FOR OUTPUT
C           USING THE CURRENT SUBSET NAME (SSLAST) IN CASE IREADPB CHANGED IT
C           ------------------------------------------------------------------

            call datelen(10)
            CALL openmb(lubfo,sslast,idate)

C           WRITE OUT THE REGISTRATION DATA FOR EACH REPORT
C           -----------------------------------------------
 
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
  
            if(subset(:6).ne.'ADPUPA') then
               obs(8,1)=bmiss
               obs(9,1)=bmiss
               obs(10,1)=bmiss
            endif

c==         if(subset(:6).eq.'AIRNOW'.and.ozone_ave.eq.'8HR_AV') then
c==            obs(:,1)=obs(:,2)
c==         endif

         PRINT*,'==NLEV=1,OBS:',(OBS(I,1),I=1,3),staid,(HDR(I),I=2,8)

        dhr=real(hdr(4))

        if(dhr.ge.0) then
          CALL raddate(idate,dhr,vdata1)
          minutes=dhr*60.
          hfrac=dhr
        else
          CALL raddate(idate,-1.0,vdata1)
          minutes=(dhr+1.)*60.
          hfrac=dhr+1
        endif

        print*,'vdata1=',vdata1

        ihr = mod(int(vdata1),100)
        obshr = real(ihr)+hfrac

       icent = mod(idate/100000000,100)
       iy = mod(int(vdata1)/1000000,100)
       im = mod(int(vdata1)/10000,100)
       id = mod(int(vdata1)/100,100)
       ihr = mod(int(vdata1),100)

       if(hdr(2).lt.0..and.obs(3,1).ne.bmiss) then
         write(103,1557) hdr(3), hdr(2)+360., obs(3,1)*1e9,staid,hdr(5),
     +   hdr(4),idate,int(vdata1),int(minutes),obshr,
     +   icent,iy,im,id,ihr,int(minutes)
       else
         write(103,1557) hdr(3), hdr(2), obs(3,1)*1e9,staid,hdr(5),
     +   hdr(4),idate,int(vdata1),int(minutes),obshr,
     +   icent,iy,im,id,ihr,int(minutes)
       endif
1557    format(f6.2,1x,f7.2,2x,f5.0,3x,a8,2x,f5.0,2x,f6.2,4x,2i12,i2.2
     +         ,2x,f5.2,2x,2i2.2,"/",i2.2,"/",i2.2,1x,i2.2,":",i2.2)
            CALL ufbint(lubfo,snf,2,nlev,iret,'SRC FHR')
            CALL ufbint(lubfo,obs,3,nlev,iret,obstr)

C           WRITE OUT THE INTERPOLATED BACKGROUND(S) INTO THE FIT FILE
C           ----------------------------------------------------------
            valix = .false.
            DO n = 1, nbak
              IF (n.lt.nbak.and.nofo(n).ne.1) THEN
                icdate = idint(vdate(n))
                CALL raddate(icdate,abs(fhr(n)-fhr(n+1)),bdate)
                valid = vdate(n) .le. vdata .and. vdate(n+1) .ge. vdata
     +                      .and. iges(n) .eq. 0 .and. iges(n+1) .eq. 0
     +                      .and. vdate(n+1) .eq. bdate .and. .not. 
     +                      valix
                valix = valid
              ELSE IF (n.eq.nbak.or.nofo(n).eq.1) THEN
                valid = vdate(n) .eq. vdata .and. iges(n) .eq. 0
              END IF
              IF (valid) THEN
                irepv(n) = irepv(n) + 1
	      if (subset(:6).eq.'AIRNOW') onlysf = .true.
                CALL getfct(irep,n,fhour,onlysf,subset)
                if (subset(:6).ne.'ADPUPA') then
                  do iii=8,10
                    bak(iii,1)=bmiss
                  enddo
                endif

                DO l = 1, nlev
                  cnf(1,l) = src(n)
                  snf(2,l) = fhours(n)
                END DO
                CALL ufbint(lubfo,snf,2,nlev,iret,'SRC FHR')
c               do j=1,50
c               if(bak(3,j).ne.bmiss) then
c               print*,'j,bak(3,j)=',j,bak(3,j)
c               endif
c               enddo
            if(subset(:6).eq.'AIRNOW')then
              if(obs(1,1).eq.-1) bak(3,1)=bak(3,1)
              if(obs(1,1).eq.-8) bak(3,1)=bak(3,2)
              print*,'AIRNOW==',obs(1,1),bak(3,1),bak(3,2)

              if(obs(1,1).eq.-1.and.obs(3,1).ne.bmiss)
     +          write(101,1556) hdr(3), hdr(2), obs(3,1)*1e9, staid
              if(obs(1,1).eq.-8.and.obs(3,1).ne.bmiss)
     +          write(108,1556) hdr(3), hdr(2), obs(3,1)*1e9, staid
1556    format(f6.2,1x,f7.2,2x,f5.0,3x,a8)

            endif

                CALL ufbint(lubfo,bak,10,nlev,iret,obstr)
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
   30 CALL errmsg("PREPFITS - BAD OR MISSING/INCOMPLETE AIRNOW DATA")
      CALL ERREXIT(0)
      END
