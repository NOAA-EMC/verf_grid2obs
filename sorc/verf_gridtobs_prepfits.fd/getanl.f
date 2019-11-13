C----------------------------------------------------------------------
C GETANL interpolates analysis if obs missing (moisture for ACARS)
C----------------------------------------------------------------------
      SUBROUTINE getanl
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
      use surfce
      use grid3d
      use gridef
      use guess
      use weights
      use vrtfac
      use anl
      INCLUDE 'parm.inc'
 
      REAL*8 bdate
      real*8 hstid
      real*8 hi
      real*8 pob(255,mxr)
c     real*8 xdft(3,255,32000)
c     real*8 pob(255,32000)
c     integer xlev(32000)
c     real*8 qtemp(255,32000)

      character*80 file(mxb), fndx(mxb)
      character*80 obstrout,obstrout2,obstrout3,obstrout4
      character*8 sub(mxr),std(mxr)
c     character*8 subset,sub(32000)

c     DIMENSION mate(mxr)
      DIMENSION mate(32000)
 
      character*8 cnf(2,255), sslast
      CHARACTER*10 src(mxb)
      character*10 source,sourcer
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
      obstr8='HOCB'
c     obstr8='CDBP'
c     obstr7='.DTHDOFS DOFS'
c     obstrout='POB QOB TOB ZOB UOB VOB PMO MXTM MITM TDO HOVI'
c     obstrout='POB QOB TOB ZOB UOB VOB PMO MXTM TDO'
c     obstrout='POB QOB TOB ZOB UOB VOB PMO MXTM'
c     obstrout='POB QOB TOB ZOB UOB VOB PMO PBL CAPE CINH LI TROP PWO'
      obstrout='POB QOB TOB ZOB UOB VOB PMO'
      obstrout2='MXTM MITM TDO HOVI TOCC MXGS THI TCH CDBZ'
      obstrout3='CAPE CINH LI TROP PWO BCAPE'
c  NEW PBL
      obstrout4='TKEPBL RIPBL MIXHT'
c  NEW PBL
      qmstr='PQM QQM TQM ZQM WQM                       '

      print*,'====== BEGIN BACKGROUND (ANALYSIS) PROCESSING ======='
 
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

      print*,'start of getanl'

      call datelen(10)

      iread=11
      iwrite=6
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
 
      lubfi = 20
      luges = 21
c     lundx = 22
      lugbi = 23
c     lubfo = 51

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
 
c--- for analysis
      nbak=1
      file(1)='ANL.tm00'
      fndx(1)='ANLi.tm00'
c     file(1)='AWIP12.tm00'
c     fndx(1)='AWIP12i.tm00'

c==      DO n = 1, mxb+1
c==        READ (5,'(A8,1X,A50)',END=10) src(n), file(n)
           READ (5,*) fhour, sourcer
           rewind(5)
           source=trim(sourcer)
c          print*,'file=',file(n)
           print*,'file=',source
c==        READ (5,'(A8,1X,A50)',END=10) src(n), fndx(n)
c==      END DO
  
c==      call errmsg('PREPFITS - TOO MANY BACKGROUND FIELDS')
c==      call errexit(1)
c==   10 nbak = n - 1
      PRINT *, 'NBAK=', nbak
c
cccccc 2/21/2006:  Read in the entire prepbufr file FIRST
c 2/22/2006:  Read in the prepbufr file's drift winds/time
c
C     OPEN THE INPUT FILE TO READ THROUGH
C     -----------------------------------
                                                                                          
      CALL openbf(lubfi,'IN ',lubfi)
c     call maxout(30000)
c     call maxout(2500000)
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
c         print*,'subset,sslast=',subset,sslast
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
c       print*,'After ireadpb(anl) ',staid,'nlev= ',nlev
c     print*,'subset,staid(anl)=',subset,staid,irep

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
        CLOSE (luges)
        CLOSE (lugbi)
c       ishl = ishell('assign -a '//file(ibak)//' -s unblocked fort.21')
c       ishi = ishell('assign -a '//fndx(ibak)//' -s unblocked fort.23')
c       ishl = system('ln -s -f '//file(ibak)//'  fort.21')
c       ishi = system('ln -s -f '//fndx(ibak)//'  fort.23')
c       call system('ln -s -f '//file(ibak)//'  fort.21',ishl)
c       call system('ln -s -f '//fndx(ibak)//'  fort.23',ishi)
c       open(21,file=file(ibak),form='unformatted',iostat=ishl)
c       open(23,file=fndx(ibak),form='unformatted',iostat=ishi)
        call baopenr(21,file(ibak),ishl)
        call baopenr(23,fndx(ibak),ishi)
c       ishl=0
c       ishi=0
     
C       For each background, call GETBAK and then call GETPROF
C       To get a background profile at each ob location
C       -------------------------------------------------
        IF (ishl.eq.0.and.ishi.eq.0) THEN
         print*,'  ibak= ', ibak
          CALL getbak(luges,lugbi,numlev,
     +             iearth,ismart,ipcp,iges(ibak),fhour,source)
          print*,'after getbak','  iges(ibak)=',iges(ibak),'  ibak= ', ibak
c         print*,'=======',iearth,pob
          IF (iges(ibak).eq.0) CALL getprofupr(iearth,sub,pob,std)
          print*,'+++pbl',pbli(1,1),pbli(2,1),pbli(3,1)
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

      if(subset(:6).eq.'AIRSND')
     *  obs(2,:)=0.

      if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.'AIRSND')
     *    then
       do i=1,nlev
         obs(8,i)=bmiss
         obs(9,i)=bmiss
         obs(10,i)=bmiss
c        obs(11,i)=bmiss
       enddo
      endif

      if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.'AIRSND')
     &          then
c      print*,'before obspbl'

       if(subset(:6).eq.'AIRSND')then
       obs(4,:)=bmiss
       endif

c      call obspbl(hpbl,idate)
c      bak(8,1)=pbli(irep,1)
c      bak3(1,2)=pbli(irep,1)
c      obsi2(1,2)=hpbl
       obsi4(1,2)=hpbl
       obsi4(2,2)=hpbl
       obsi4(3,2)=hpbl
c      
       do i=1,nlev
        if(cat(i).eq.5.0) then
          obsi2(4,2)=obs(1,i)
          goto 123
        endif
       enddo
c      do i=2,nlev
c        obs(12,i)=bmiss
c      enddo
123    continue

c
c
      endif    !  if ADPUPA 

C         FILTER THE DATA TO FIT BY CATEGORY
C         ----------------------------------
c-- No filter for analysis fields
          onlysf = .false.
          mlev = nlev
 
C         MAKE SURE THERE IS DATA TO WRITE OUT
          IF (nlev.gt.0) THEN
            irepo = irepo + 1
C           MAKE SURE A MESSAGE OF THE PROPER TYPE AND DATE IS OPEN FOR OUTPUT
C           USING THE CURRENT SUBSET NAME (SSLAST) IN CASE IREADPB CHANGED IT
C           ------------------------------------------------------------------

            
c           if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.'AIRSND')
c    &          then
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
c            do i=1,6
c            obsi2(i,1)=obsi2(i,2)
c            enddo
c           endif
            if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.'AIRSND')
     &         then
             do i=1,3
             obsi4(i,1)=obsi4(i,2)
             enddo
            endif
c 

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
                CALL getfct(irep,n,fhour,onlysf,subset,source)
c-- For backgroung moisture instead of obs moisture---
            DO l = 1, nlev
              qtemp(l,irep) = bak(2,l)
c             print*,'==anl,irep=',irep,qtemp(l,irep),l
            END DO

              ELSE
                PRINT *, 'VALID IS FALSE!!  IREP,N,IGES(N)', irep, n, 
     +                      iges(n)
                PRINT *, 'VALID IS FALSE!!  VDATE(N),VDATA', vdate(n), 
     +                      vdata
              END IF
            END DO
 
C           BEFORE END OF READ LOOPS - WRITE THE PREPFITS SUBSET
C           ----------------------------------------------------

          END IF
          sslast = subset
 
C       END OF READ LOOPS
        END DO
      END DO
 
C     CLOSE THE BUFR FILES
C     --------------------
 
      CALL closbf(lubfi)
c     CALL closbf(lubfo)
      PRINT *, '************************************'
      PRINT *, 'WROTE OUT FITS FOR ', irepo, ' REPORTS'
      PRINT *, '************************************'
      PRINT *, (irepv(i),I=1,nbak)
      PRINT *, '************************************'
      PRINT *, ' # of non-missing frcst reports for ONLYSF = ',
     +         kntgsf
      PRINT *, '************************************'
      print*,'++++bak(2,1)',bak(2,1)
 
C     END OF PREPFITTING
C     ------------------
 
      print*,'====== END BACKGROUND (ANALYSIS) PROCESSING ======='
   20 CONTINUE
 
       IF (ishl.eq.0.and.ishi.eq.0) THEN
        deallocate(pi)
        deallocate(zi)
        deallocate(ti)
        deallocate(ui)
        deallocate(vi)
        deallocate(qi)
        deallocate(ai)

        deallocate(psi)
        deallocate(zsi)
        deallocate(tsi)
        deallocate(usi)
        deallocate(vsi)
        deallocate(qsi)
        deallocate(bcpi)
        deallocate(pxi)
        deallocate(tmxi)
        deallocate(dpti)
        deallocate(visbi)
        deallocate(tmni)
        deallocate(tocci)
        deallocate(gusti)
        deallocate(pbli)
        deallocate(pwi)
        deallocate(tropi)
        deallocate(cdti)
        deallocate(pblrii)
        deallocate(pblmxi)
        deallocate(pmi)
        deallocate(cpi)
        deallocate(cni)
        deallocate(hni)
        deallocate(transi)
        deallocate(wnd80i)
        deallocate(ceili)

        deallocate(kxid)
        deallocate(kyjd)
        deallocate(wtswd)
        deallocate(wtsed)
        deallocate(wtnwd)
        deallocate(wtned)

        deallocate(kxi)
        deallocate(kyj)
        deallocate(wtsw)
        deallocate(wtse)
        deallocate(wtnw)
        deallocate(wtne)

        deallocate(rm1)
        deallocate(rm2)
        deallocate(vrterp)
c       call baclose(luges,iret)
c       call baclose(lugbi,iret)
       endif
      RETURN
   30 call errmsg('PREPFITS - BAD OR MISSING PREPBUFR DATE')
      call errexit(1)
      END

