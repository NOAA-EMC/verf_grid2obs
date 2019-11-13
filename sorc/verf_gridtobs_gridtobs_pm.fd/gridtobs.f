C********************************************************************
C  GRIDTOBS   -  CREATE VERIFICATION STATS BETWEEN OBSERVATIONS AND
C                BACKGROUNDS.  THE BACKGROUND & OBS VALUES ARE READ
C                FROM A PREPBTIM (GLOBAL) OR PREPFITS (MESO) BUFR FILE
C
      PROGRAM gridtobs

      INCLUDE 'parm.inc'

      DIMENSION sumdata(mxfcst,mxvrbl,maxlvl,mxarea,maxobs), 
     +            sumgrid(mxfcst,mxvrbl,maxlvl,mxarea,maxobs), 
     +            sumprod(mxfcst,mxvrbl,maxlvl,mxarea,maxobs), 
     +            ssqdata(mxfcst,mxvrbl,maxlvl,mxarea,maxobs), 
     +            ssqgrid(mxfcst,mxvrbl,maxlvl,mxarea,maxobs), 
     +            count(mxfcst,mxvrbl,maxlvl,mxarea,maxobs,mxstat)

      DIMENSION nchrmodel(maxmod), nchrfcst(mxfcst), nchrvfdate(mxdate),
     +            nchrvfyobs(maxobs), nchrarea(mxarea), 
     +            nchrstat(mxstat), nchrvarbl(mxvrbl), 
     +            nchrlevel(maxlvl)
      CHARACTER*24 namodel(maxmod), namfcst(mxfcst), namvfdate(mxdate),
     +            namvfyobs(maxobs), namarea(mxarea), namstat(mxstat), 
     +            namvarbl(mxvrbl), namlevel(maxlvl)

      COMMON /names/ namodel, namfcst, namvfdate, namvfyobs, namarea, 
     +            namstat, namvarbl, namlevel
      COMMON /nchrs/ nchrmodel, nchrfcst, nchrvfdate, nchrvfyobs, 
     +            nchrarea, nchrstat, nchrvarbl, nchrlevel
      LOGICAL	  vtflg, nmbflg
      COMMON /cnvrsns/ vtflg, nmbflg (maxmod), concon (maxmod),
     +		       cenlon (maxmod)
      CHARACTER*3 regions (30)
      COMMON /grdef/ mode(mxarea), imax(mxarea), imin(mxarea), 
     +            jmax(mxarea), jmin(mxarea), alat1(mxarea), 
     +            elon1(mxarea), dxx(mxarea), dyy(mxarea), 
     +            elonv(mxarea), alatan(mxarea), latlong(mxarea), 
     +            lambert(mxarea), polarstereo(mxarea), numreg(mxarea),
     +            ig104(147,110), regions

      CHARACTER*3 namversion
      CHARACTER*132 vdbhdr132, input, substr (3)

      CHARACTER*1 blank, equal
C
      CHARACTER*50 headr, obstr, qmstr
      CHARACTER*8 subset, stnid
      CHARACTER*3 field
      
C
      INTEGER	lstart (maxobs), lstop (maxobs)
      REAL*8 hdr(10), obs(12,255,mxb), qms(8,255), obs1(12,255)
      real*8 dhr,probs
      EQUIVALENCE (hdr(1),stnid)
C
C...   STRING FOR HEADER PARAMETERS
C
      DATA headr /'SID XOB YOB DHR ELV TYP T29 ITP'/
C
C...   STRING FOR THE OB, GUESS, ANALYSIS ....
C
c     DATA obstr /'SRC FHR POB QOB TOB ZOB UOB VOB PMO CAPE CINH LI'/
      DATA obstr /'SRC FHR TPHR QCIND COPOPM                       '/
C
C...   STRING FOR THE QUALITY MARKS ....
C
      DATA qmstr /'CAT PRC PQM QQM TQM ZQM WQM    '/
      DATA stnid /'        '/
C----------------------------------------------------
      DATA blank /' '/
      DATA equal /'='/
      DATA namversion /'V01'/
C
      DATA bmiss /10E10/
      DATA rmiss /99999./

C--AIRNOW ozone verification using thresholds.
c==   PARAMETER (NUMTHR=6,IHR=10)
      parameter (numthr=7)
c==   REAL*8 THRESH(NUMTHR)
      dimension thresh(numthr)
c     CHARACTER THR(6)*4
c==   REAL FG(NUMTHR,IHR),HG(NUMTHR,IHR),OG(NUMTHR,IHR),TG(IHR)
c==   REAL*8  FOT(NUMTHR,IHR),HOT(NUMTHR,IHR),OOT(NUMTHR,IHR)
      dimension fot(mxfcst,mxvrbl,maxlvl,mxarea,maxobs,numthr),
     +          hot(mxfcst,mxvrbl,maxlvl,mxarea,maxobs,numthr),
     +          oot(mxfcst,mxvrbl,maxlvl,mxarea,maxobs,numthr),
     +          fg(mxfcst,mxvrbl,maxlvl,mxarea,maxobs,numthr),
     +          hg(mxfcst,mxvrbl,maxlvl,mxarea,maxobs,numthr),
     +          og(mxfcst,mxvrbl,maxlvl,mxarea,maxobs,numthr),
     +          tg(mxfcst,mxvrbl,maxlvl,mxarea,maxobs)
c      DATA THRESH /50., 65., 85., 105., 125., 150./
       data thresh /5., 10., 12., 15., 20., 25., 35./
C--
      call datelen(10)
   10 CONTINUE

C     
C     READ VERIFICATION DATABASE VERSION AND INPUT UNIT NUMBER
C
C     Also read in logical flag (T or F) to convert virtual temperature
C     observed into actual temperature.
C     
c     READ (5,'(A)',END=160) input
c     CALL ST_CLST ( input, ' ', ' ', 3, substr, num, ier )
c     namversion = substr (1)
c     CALL ST_NUMB ( substr (2), lunin, ier )
c     vtflg = ( substr (3) .eq. 't' .or. substr (3) .eq. 'T' )
      READ (5,*,END=160) namversion,lunin
      PRINT '(A3,I5)', namversion, lunin
      CALL datebf(lunin,iy,im,id,ih,idate)
      print*,'iy,im,id,ih=',iy,im,id,ih
      CALL openbf(lunin,'IN ',lunin)
      PRINT *, ' DATE OF INPUT BUFR FILE', idate, ' UNIT=', lunin
C     
C     READ REST OF THIS CONTROL-FILE GROUP
C     TO GET NUMBER OF THINGS TO BE VERIFIED
C     
      CALL readcntl(numodel,numfcst,numvfdate,numvfyobs,numarea,numstat,
     +            numvarbl,numlevel,numvector)

      IF (numvfdate.eq.1.and.nchrvfdate(1).eq.2) THEN
        WRITE (namvfdate(1)(1:10),'(I10)') idate
        nchrvfdate(1) = 10
        PRINT '(" NEW NAMVFDATE =",A24)', namvfdate(1)
        PRINT '(" CHARACTER COUNT =",I5)', nchrvfdate(1)
      END IF

C     
C     LET'S DO THE LOOPS
C     
C     OUTERMOST LOOP OVER VERIFYING DATE
C     
c     print*,'numvfdate=',numvfdate
      DO 150 ivfdate = 1, numvfdate
          FG      = 0.
          HG      = 0.
          OG      = 0.
          TG      = 0.
C       
C       OUTER LOOP OVER FORECAST MODEL TO BE VERIFIED
C       
c       print*,'numodel=',numodel
        DO 140 imodel = 1, numodel
C         
C         ZERO THE SUMS
C         
          count = 0.0
          sumdata = 0.0
          ssqdata = 0.0
          sumgrid = 0.0
          ssqgrid = 0.0
          sumprod = 0.0
C         
C         HERE WE GET A BUFR FILE AND START LOOP OVER EACH OB
C         
          DO WHILE (ireadmg(lunin,subset,jdate).eq.0)
            nsub = nmsub(lunin)
            PRINT *, 'SUBSET,NSUB,JDATE', subset, nsub, jdate
            ichk = 0
            print*,'numvfyobs=',numvfyobs
            DO 20 iob = 1, numvfyobs
	      lstart (iob) = 1
	      lstop (iob) = 0
              IF (namvfyobs(iob)(:6).eq.'ANYSFC') THEN
                IF (subset(:6).eq.'ADPSFC'.or.subset(:6).eq.'SFCSHP'.or.
     +                      subset(:6).eq.'ADPUPA'.or.subset(:6).eq.
     +                      'PROFLR'.or.subset(:6).eq.'ANOWPM') THEN
			    ichk = 1
			    IF (subset(:6).eq.'ADPSFC'.or.
     +                          subset(:6).eq.'ANOWPM'.or.
     +			        subset(:6).eq.'SFCSHP') lstop (iob) = 1
		END IF
              ELSE IF (namvfyobs(iob)(:6).eq.'ANYAIR') THEN
                IF (subset(:6).eq.'AIRCAR'.or.subset(:6).eq.'AIRCFT') 
     +                      ichk = 1
	      ELSE IF (namvfyobs(iob)(:6).eq.'ONLYSF') THEN
                IF (subset(:6).eq.'ADPSFC'.or.subset(:6).eq.'SFCSHP'.
     +              or.subset(:6).eq.'ANOWPM')
     +		THEN
	            ichk = 1
		    lstart (iob) = 2
		END IF
              ELSE IF (subset(:6).eq.namvfyobs(iob)(:6)) THEN
                ichk = 1
		IF (subset(:6).eq.'ADPSFC'.or.
     +		    subset(:6).eq.'SFCSHP'.or.subset(:6).eq.'ANOWPM') 
     +              lstop (iob) = 1
              END IF
   20       CONTINUE
            DO WHILE (ireadsb(lunin).eq.0.and.ichk.eq.1)
              CALL ufbint(lunin,hdr,10,1,nlev,headr)
              PRINT *, 'NLEV,HDR', nlev, hdr
              dhr = hdr(4)
              adate = jdate
c             CALL raddate(adate,dhr,vdata)
              CALL raddate(idate,dhr,vdata)
              print*,'======idate,vdata',idate,vdata
C             
C             INNER LOOP FOR ACTUAL OB-TYPE TO BE VERIFIED
C             
              DO 80 iob = 1, numvfyobs
                IF (isitob(subset,hdr(6),iob).eq.0) THEN
C                 
C                 INNER LOOP OVER VERIFYING AREA
C                 
                  DO 70 iar = 1, numarea
                 print*,'===iar1'
                    IF (inarea(imodel,stnid,hdr(2),hdr(3),iar,rm1,rm2)
     +			.eq.0) THEN

                 print*,'===iar2','lunin=',lunin,'nlev=',nlev,'nevn=',
     +                  nevn,'obstr=',obstr
                      CALL ufbin3(lunin,obs,12,255,mxb,nlev,nevn,obstr)
               do kk=1,2
               print*,'===',kk,stnid,(obs(ii,1,kk),ii=2,5)
               enddo

c              do kk=1,mxb
c                     if (iar.eq.1) then
c                       k=1
c                       do j=1,50
c                       if (obs(4,j,k).ne.bmiss) then
c                       print*,'j,k,obs(4,j,k)=',j,k,obs(4,j,k)
c                       endif
c                       enddo
c                     endif
                      CALL ufbint(lunin,qms,8,255,nlev,qmstr)
C                     
C                     INNER LOOP OVER OBSERVATION LEVELS
C                     
		      IF ( lstop (iob) .eq. 0 ) THEN
			nstop = nlev
		      ELSE
			nstop = lstop (iob)
		      END IF
                      DO 60 nlv = lstart (iob), nstop
                        probs = obs(3,nlv,nevn)
                        kat = nint(qms(1,nlv))
C                       INNER LOOP OVER REQUESTED LEVELS
                        DO 50 ilv = 1, numlevel
c                         print*,'subset,probs,kat,ilv=',
c   *                       subset,probs,kat,ilv
                        iilayer=inlayer(subset,probs,kat,ilv)
                        print*,'iilayer=',iilayer
                        iilayer=0
                        IF (iilayer.eq.0) THEN
c                         IF (inlayer(subset,probs,kat,ilv).eq.0) THEN
C                           
C                           INNER LOOP OVER VARIABLE 
C                           
                            DO 40 ivr = 1, numvarbl
                            print*,namvarbl(ivr)
                            if(namvarbl(ivr).eq.'OZON/1'.or.
     +                         namvarbl(ivr).eq.'OZON/8') then
                                if(namvarbl(ivr)(5:6).eq.'/1')
     +                             ozaverage=-1
                                if(namvarbl(ivr)(5:6).eq.'/8')
     +                             ozaverage=-8
                                namvarbl(ivr)='COPO'
                                nchrvarbl(ivr)=4
                            endif
                            print*,ozaverage,namvarbl(ivr)

                            if(namvarbl(ivr).eq.'PM25/1'.or.
     +                         namvarbl(ivr).eq.'PM25MX'.or.
     +                         namvarbl(ivr).eq.'PM25AV') then
                               ozaverage=-1
                               if(namvarbl(ivr)(5:6).eq.'/1')
     +                           field="/1" 
                               if(namvarbl(ivr)(5:6).eq.'MX')
     +                           field="MX" 
                               if(namvarbl(ivr)(5:6).eq.'AV')
     +                           field="AV" 
                               namvarbl(ivr)='COPOPM'
                               nchrvarbl(ivr)=6
                            endif

C                             
C                             INNER LOOP OVER NUMBER OF FORECAST HOUR (EVENT)
C                             
c      icnt=0
                              DO 30 ifh = 1, numfcst
c      print*,'igotdata=',igotdata(obsval,forcst,obs,qms,nevn,
c    +                                      nlv,imodel,ifh,ivr,ilv,
c    +                                      iob,rm1,rm2,ozaverage)
                            print*,'---',ozaverage,namvarbl(ivr)
                                IF (igotdata(obsval,forcst,obs,qms,nevn,
     +                                      nlv,imodel,ifh,ivr,ilv,
     +					    iob,rm1,rm2,ozaverage)
     +                                      .eq.0) THEN
c      icnt=icnt+1
c      print*,'icnt=',icnt
       print*,'obsval=',obsval,'forcst=',forcst,stnid

                                 if(subset(:6).eq.'ANOWPM') then
                                  obsval=obsval*1.e9
                                  forcst=forcst*1.e9
                                 endif
        print*,'WE ARE HERE!!!', numstat,namstat

                     DO 200 ist=1,numstat

                     IF(namstat(ist).eq.'FHO') THEN

c==  +                                  OG(LTHR,ifh)=OG(LTHR,ifh)+1.
c==  +                                  FG(LTHR,ifh)=FG(LTHR,ifh)+1.
c==  +                                  HG(LTHR,ifh)=HG(LTHR,ifh)+1.
c==                                   TG(ifh)=TG(ifh)+1.
                                   DO LTHR=1,NUMTHR
                                      IF(obsval.GT.THRESH(LTHR))
     +                               og(ifh,ivr,ilv,iar,iob,lthr)=
     +                               og(ifh,ivr,ilv,iar,iob,lthr)+1.
                                      IF(forcst.GT.THRESH(LTHR))
     +                                  fg(ifh,ivr,ilv,iar,iob,lthr)=
     +                                  fg(ifh,ivr,ilv,iar,iob,lthr)+1.
                                      IF(obsval.GT.THRESH(LTHR).AND.
     +                                   forcst.GT.THRESH(LTHR))
     +                                  hg(ifh,ivr,ilv,iar,iob,lthr)=
     +                                  hg(ifh,ivr,ilv,iar,iob,lthr)+1.
                                   ENDDO
                                        tg(ifh,ivr,ilv,iar,iob)=
     +                                  tg(ifh,ivr,ilv,iar,iob)+1.
   
                                  count(ifh,ivr,ilv,iar,iob,ist) = 
     +                                    count(ifh,ivr,ilv,iar,iob,ist)
     +                                    + 1.0

                     ELSEIF(namstat(ist).eq.'SL1L2') THEN

                                  sumd = sumdata(ifh,ivr,ilv,iar,iob) *
     +                                    count(ifh,ivr,ilv,iar,iob,ist)
     +                                    + obsval
                                  ssqd = ssqdata(ifh,ivr,ilv,iar,iob) *
     +                                    count(ifh,ivr,ilv,iar,iob,ist)
     +                                    + obsval * obsval
                                  sumg = sumgrid(ifh,ivr,ilv,iar,iob) *
     +                                    count(ifh,ivr,ilv,iar,iob,ist)
     +                                    + forcst
                                  ssqg = ssqgrid(ifh,ivr,ilv,iar,iob) *
     +                                    count(ifh,ivr,ilv,iar,iob,ist)
     +                                    + forcst * forcst
                                  prod = forcst * obsval
                                  sump = sumprod(ifh,ivr,ilv,iar,iob) *
     +                                    count(ifh,ivr,ilv,iar,iob,ist)
     +                                    + prod
                                  count(ifh,ivr,ilv,iar,iob,ist) = 
     +                                    count(ifh,ivr,ilv,iar,iob,ist)
     +                                    + 1.0
      print*,'count=',ifh,ivr,ilv,iar,iob,ist,
     +        count(ifh,ivr,ilv,iar,iob,ist),obsval,forcst,stnid 
                                  sumdata(ifh,ivr,ilv,iar,iob) = sumd /
     +                                    count(ifh,ivr,ilv,iar,iob,ist)
                                  ssqdata(ifh,ivr,ilv,iar,iob) = ssqd /
     +                                    count(ifh,ivr,ilv,iar,iob,ist)
                                  sumgrid(ifh,ivr,ilv,iar,iob) = sumg /
     +                                    count(ifh,ivr,ilv,iar,iob,ist)
                                  ssqgrid(ifh,ivr,ilv,iar,iob) = ssqg /
     +                                    count(ifh,ivr,ilv,iar,iob,ist)
                                  sumprod(ifh,ivr,ilv,iar,iob) = sump /
     +                                    count(ifh,ivr,ilv,iar,iob,ist)

                     END IF  !stat type
  200     CONTINUE
                                END IF !igotodata
C                             
C                             END INNER LOOP OVER NUMBER OF FORECAST HOUR (EVENT)
   30                         CONTINUE

                            if(namvarbl(ivr).eq.'COPO'.and.
     +                                     ozaverage.eq.-1)then
c                                          namvarbl(ivr)='COPO-1'
                                           namvarbl(ivr)='OZON/1'
                                           nchrvarbl(ivr)=6
                            endif
                            if(namvarbl(ivr).eq.'COPO'.and.
     +                                     ozaverage.eq.-8)then
c                                          namvarbl(ivr)='COPO-8'
                                           namvarbl(ivr)='OZON/8'
                                           nchrvarbl(ivr)=6
                            endif
                            if(namvarbl(ivr).eq.'COPOPM'.and.
     +                                     ozaverage.eq.-1)then
                                           if(field.eq.'/1')
     +                                       namvarbl(ivr)='PM25/1'
                                           if(field.eq.'MX')
     +                                       namvarbl(ivr)='PM25MX'
                                           if(field.eq.'AV')
     +                                       namvarbl(ivr)='PM25AV'
                                           nchrvarbl(ivr)=6
                            endif

C                           END INNER LOOP OVER VARIABLE
   40                       CONTINUE
                          END IF
C                       END INNER LOOP OVER REQUESTED LEVELS
   50                   CONTINUE
C                     END INNER LOOP OVER OBSERVATION LEVELS
   60                 CONTINUE
                    END IF
C                 END INNER LOOP OVER VERIFYING AREA
   70             CONTINUE
                END IF
C             END LOOP FOR ACTUAL OB-TYPE
   80         CONTINUE
C           END DO WHILE LOOP OVER BUFR REPORT
            END DO
C         END DO WHILE LOOP OVER BUFR MESSAGE
          END DO
C         
C         NOW IS THE TIME TO WRITE OUT THE STAT RECORDS
C         
c         ist = 1
         DO 300 ist=1,numstat
          DO 130 iob = 1, numvfyobs
            DO 120 iar = 1, numarea
              DO 110 ifh = 1, numfcst
                numv = numvarbl
                IF (numvector.gt.0) numv = numvarbl - 1
                DO 100 ivr = 1, numv
                  DO 90 ilv = 1, numlevel
c                   rtest = count(ifh,ivr,ilv,iar,iob,ist)
c                   print *, 'COUNT(IFH,IVR,ILV,IAR,IOB,IST)=',RTEST
                    IF (count(ifh,ivr,ilv,iar,iob,ist).gt.0.0) THEN
                      iend = 3
                      vdbhdr132(1:iend) = namversion(:3)
                      iend = iend + 1
                      vdbhdr132(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrmodel(imodel)
                      vdbhdr132(istrt:iend) = 
     +                            namodel(imodel)(:nchrmodel(imodel))
                      iend = iend + 1
                      vdbhdr132(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrfcst(ifh)
                      vdbhdr132(istrt:iend) = 
     +                            namfcst(ifh)(:nchrfcst(ifh))
                      iend = iend + 1
                      vdbhdr132(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrvfdate(ivfdate)
                      vdbhdr132(istrt:iend) = 
     +                            namvfdate(ivfdate)(:
     +                            nchrvfdate(ivfdate))
                      iend = iend + 1
                      vdbhdr132(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrvfyobs(iob)
                      vdbhdr132(istrt:iend) = 
     +                            namvfyobs(iob)(:nchrvfyobs(iob))
                      iend = iend + 1
                      vdbhdr132(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrarea(iar)
                      vdbhdr132(istrt:iend) = 
     +                            namarea(iar)(:nchrarea(iar))
                      iend = iend + 1
                      vdbhdr132(iend:iend) = blank
                      istrt = iend + 1
                      ivstrt = istrt
                      iend = iend + nchrstat(ist)
                      vdbhdr132(istrt:iend) = 
     +                            namstat(ist)(:nchrstat(ist))
                      iend = iend + 1
c-------------------------------------------------------------
                        if(namstat(ist).eq.'FHO')then
                          vdbhdr132(iend:iend) = '>'
                          ibreak1=iend
                          ibreak2=iend + 1
                          iend = iend + 1
                        endif
c-------------------------------------------------------------
                      vdbhdr132(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrvarbl(ivr)
                      vdbhdr132(istrt:iend) = 
     +                            namvarbl(ivr)(:nchrvarbl(ivr))


                      print*,'@@@@@@@@@@@@',namvarbl(ivr),TPHR
c                       if(namvarbl(ivr).eq.'COPO') then
c                         if(namstat(ist).eq.'FHO')
c    +                       vdbhdr132(istrt:iend) = 'OZON'
c                         istrt = iend + 1
c                         iend = iend + 2
c                          if(ozaverage.eq.-1.) then
c                            vdbhdr132(istrt:iend) = '-1'
c                            if(namstat(ist).eq.'FHO')
c    +                          vdbhdr132(istrt:iend) = '/1'
c                          elseif(ozaverage.eq.-8.) then
c                            vdbhdr132(istrt:iend) = '-8'
c                            if(namstat(ist).eq.'FHO')
c    +                          vdbhdr132(istrt:iend) = '/8'
c                          endif
c                       endif

                      iend = iend + 1
                      vdbhdr132(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrlevel(ilv)
                      vdbhdr132(istrt:iend) = 
     +                            namlevel(ilv)(:nchrlevel(ilv))
                      iend = iend + 1
                      vdbhdr132(iend:iend) = blank
                      iend = iend + 1
                      vdbhdr132(iend:iend) = equal
                      iend = iend + 1
                      vdbhdr132(iend:iend) = blank

                     IF(namstat(ist).eq.'FHO') THEN

c------ ozone FHO
c==                       print*,'FG=',FG
c==                       print*,'OG=',OG
c==                       print*,'HG=',HG
c==                       print*,'TG=',TG
c==                       IF(TG(ifh).gt.0)THEN
                          IF(TG(ifh,ivr,ilv,iar,iob).gt.0)THEN
                           DO LTHR=1,NUMTHR
c==                FOT(LTHR,ifh)=FG(LTHR,ifh)/TG(ifh)
c==                HOT(LTHR,ifh)=HG(LTHR,ifh)/TG(ifh)
c==                OOT(LTHR,ifh)=OG(LTHR,ifh)/TG(ifh)
                   fot(ifh,ivr,ilv,iar,iob,lthr)=
     +             fg(ifh,ivr,ilv,iar,iob,lthr)/tg(ifh,ivr,ilv,iar,iob)
                   hot(ifh,ivr,ilv,iar,iob,lthr)=
     +             hg(ifh,ivr,ilv,iar,iob,lthr)/tg(ifh,ivr,ilv,iar,iob)
                   oot(ifh,ivr,ilv,iar,iob,lthr)=
     +             og(ifh,ivr,ilv,iar,iob,lthr)/tg(ifh,ivr,ilv,iar,iob)

                        PRINT'(A,F7.0,7E18.9)',VDBHDR132(:IEND)
c                  WRITE (50,1201) vdbhdr132(:26),
c    +               vdbhdr132(35:38),THRESH(LTHR),
c    +                             vdbhdr132(53:60),
c    +                       TG(ifh),FOT(LTHR,ifh)
c    +                      ,HOT(LTHR,ifh),OOT(LTHR,ifh)
c1201        FORMAT (A,1x,'ANOWPM',1x,A,1x,' FHO>',f4.0,1x,
c    +               'OZON/1',1x,A,f7.0,1x,3f8.5)
                        WRITE (50,1202) vdbhdr132(:ibreak1), 
     +                       THRESH(LTHR),
     +                       vdbhdr132(ibreak2:iend), 
     +                       count(ifh,ivr,ilv,iar,iob,ist), 
c==  +                       FOT(LTHR,ifh),HOT(LTHR,ifh),OOT(LTHR,ifh)
     +                       fot(ifh,ivr,ilv,iar,iob,lthr),
     +                       hot(ifh,ivr,ilv,iar,iob,lthr),
     +                       oot(ifh,ivr,ilv,iar,iob,lthr)
                           ENDDO
                         ENDIF
 1202                   FORMAT (A,f0.0,A,F7.0,3f8.5)
c------

                     ELSEIF(namstat(ist).eq.'SL1L2') THEN

c                     print *, 'IVR=',IVR,' and NUMVECTOR=',NUMVECTOR
                      IF (ivr.ne.numvector) THEN

                        PRINT'(A,F7.0,5E18.9)',VDBHDR132(:IEND)
                        WRITE (50,1000) vdbhdr132(:iend), 
     +                              count(ifh,ivr,ilv,iar,iob,ist), 
     +                              sumgrid(ifh,ivr,ilv,iar,iob), 
     +                              sumdata(ifh,ivr,ilv,iar,iob), 
     +                              sumprod(ifh,ivr,ilv,iar,iob), 
     +                              ssqgrid(ifh,ivr,ilv,iar,iob), 
     +                              ssqdata(ifh,ivr,ilv,iar,iob)
 1000                   FORMAT (A,F7.0,5E18.9)
                      END IF
C                     
                      IF (ivr.eq.numvector) THEN
                        sump = sumprod(ifh,ivr,ilv,iar,iob) + 
     +                              sumprod(ifh,ivr+1,ilv,iar,iob)
                        sumg = ssqgrid(ifh,ivr,ilv,iar,iob) + 
     +                              ssqgrid(ifh,ivr+1,ilv,iar,iob)
                        sumd = ssqdata(ifh,ivr,ilv,iar,iob) + 
     +                              ssqdata(ifh,ivr+1,ilv,iar,iob)
                        vdbhdr132(ivstrt:ivstrt) = namversion(1:)
C                       PRINT'(A,F7.0,7E18.9)',VDBHDR132(:IEND),
                        WRITE (50,1100) vdbhdr132(:iend), 
     +                              count(ifh,ivr,ilv,iar,iob,ist), 
     +                              sumgrid(ifh,ivr,ilv,iar,iob), 
     +                              sumgrid(ifh,ivr+1,ilv,iar,iob), 
     +                              sumdata(ifh,ivr,ilv,iar,iob), 
     +                              sumdata(ifh,ivr+1,ilv,iar,iob), 
     +                              sump, sumg, sumd
 1100                   FORMAT (A,F7.0,7E18.9)
                      END IF !vector
                    END IF   !stat type
                   END IF    !count
   90             CONTINUE
  100           CONTINUE
  110         CONTINUE
  120       CONTINUE
  130     CONTINUE
C         END LOOP OVER STATISTIC TYPE
  300     CONTINUE
C         END OUTER LOOP OVER FORECAST MODEL
  140     CONTINUE
C         END OUTERMOST LOOP OVER VERIFYING DATE
  150     CONTINUE
          CALL closbf(lunin)
C         GO BACK AND DO IT ALL AGAIN   WHOOPEE!
          GO TO 10
C23456789012345678901234567890123456789012345678901234567890123456789012
  160     CONTINUE
          STOP
          END
