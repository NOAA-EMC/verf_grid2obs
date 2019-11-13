C********************************************************************
C  GRIDTOBS   -  CREATE VERIFICATION STATS BETWEEN OBSERVATIONS AND
C                BACKGROUNDS.  THE BACKGROUND & OBS VALUES ARE READ
C                FROM A PREPBTIM (GLOBAL) OR PREPFITS (MESO) BUFR FILE
C
C  Modification history:
C   Add PBS score, by Binbin Zhou, SAIC. Jan 18, 2008
C
C
      PROGRAM gridtobs

      INCLUDE 'parm.inc'

      DIMENSION sumdata(mxfcst,mxvrbl,maxlvl,mxarea,maxobs), 
     +            sumgrid(mxfcst,mxvrbl,maxlvl,mxarea,maxobs), 
     +            sumprod(mxfcst,mxvrbl,maxlvl,mxarea,maxobs), 
     +            ssqdata(mxfcst,mxvrbl,maxlvl,mxarea,maxobs), 
     +            ssqgrid(mxfcst,mxvrbl,maxlvl,mxarea,maxobs), 
     +            sumensvar(mxfcst,mxvrbl,maxlvl,mxarea,maxobs), 
     +            sumfracle3(mxfcst,mxvrbl,maxlvl,mxarea,maxobs),            
     +            count(mxfcst,mxvrbl,maxlvl,mxarea,maxobs)

      DIMENSION nchrmodel(maxmod),nchrfcst(mxfcst),nchrvfdate(mxdate),
     +            nchrvfyobs(maxobs), nchrarea(mxarea), 
     +            nchrstat(mxstat), nchrvarbl(mxvrbl), 
     +            nchrlevel(maxlvl)
      CHARACTER*24 namodel(maxmod), namfcst(mxfcst), namvfdate(mxdate),
     +            namvfyobs(maxobs), namarea(mxarea), namstat(mxstat), 
     +            namvarbl(mxvrbl), namlevel(maxlvl)

      character*(80) msg
      character*(2) alunin,alunin2
      character*(10) aobs,aobs2

      DIMENSION lunin(maxmod),obsvaln(maxmod),forcstn(maxmod)
      DIMENSION count1(mxfcst,mxvrbl,maxlvl,mxarea,maxobs),
     +          rhnt(mxfcst,mxvrbl,maxlvl,mxarea,maxobs,maxmod),
     +          dif(maxmod)
      DIMENSION count2(mxfcst,mxvrbl,maxlvl,mxarea,maxobs),
     +          rhet(mxfcst,mxvrbl,maxlvl,mxarea,maxobs,maxmodp1),
     +          frcst_rng(maxmod)

      COMMON /names/ namodel, namfcst, namvfdate, namvfyobs, namarea, 
     +            namstat, namvarbl, namlevel
      COMMON /nchrs/ nchrmodel, nchrfcst, nchrvfdate, nchrvfyobs, 
     +            nchrarea, nchrstat, nchrvarbl, nchrlevel
      LOGICAL	  vtflg, nmbflg, mask
      COMMON /cnvrsns/ vtflg, nmbflg (maxmod), concon (maxmod),
     +		       cenlon (maxmod)
      CHARACTER*3 regions (29)
      COMMON /grdef/ mode(mxarea), imax(mxarea), imin(mxarea), 
     +            jmax(mxarea), jmin(mxarea), alat1(mxarea), 
     +            elon1(mxarea), dxx(mxarea), dyy(mxarea), 
     +            elonv(mxarea), alatan(mxarea), latlong(mxarea), 
     +            lambert(mxarea), polarstereo(mxarea), numreg(mxarea),
     +            ig104(147,110), regions

      !add PBS score 
      DIMENSION fcstsave(maxmod)
      integer pbsmrk(mxvrbl)
      CHARACTER*24 op(mxvrbl),pbsthr(mxvrbl,10)
      integer  nchrpbsthr(mxvrbl,10)
      real     rpbsthr(mxvrbl,10)
      real, allocatable, dimension(:,:,:,:,:,:,:) :: pfield
      CHARACTER*420 vdbhdr420
      COMMON /pbs/ pbsmrk,op,pbsthr,nchrpbsthr,rpbsthr 

      CHARACTER*3 namversion
      CHARACTER*132 vdbhdr132, input, substr (4)
      CHARACTER*24 ens_name,ens_name_tmp

      CHARACTER*1 blank, equal
C
      CHARACTER*50 headr, obstr, obstr2,obstr3,obstr4,qmstr
      CHARACTER*8 subset, stnid
C
      INTEGER	lstart (maxobs), lstop (maxobs)
     
      REAL*8 hdr(10), obs(9,255,mxb), qms(7,255),
     *       obs2(10,255,mxb),obs4(6,2,mxb),obs5(3,2,mxb)

      real*8 dhr,probs
      EQUIVALENCE (hdr(1),stnid)

      REAL*8 oobs(9,255,mxb,maxmod), qqms(7,255,maxmod)
      REAL*8 oobs2(10,255,mxb,maxmod)
      real*8 oobs4(6,2,mxb,maxmod),oobs5(3,2,mxb,maxmod)

      REAL*8 hdr_prev_mem(10),obs_prev_mem(9,255,mxb)
      CHARACTER*50 headr_prev_mem
      CHARACTER*8 subset_prev_mem
      real timef,elapsed

      elapsed=timef()
C
C...   STRING FOR HEADER PARAMETERS
C
      DATA headr /'SID XOB YOB DHR ELV TYP T29 ITP'/
C
C...   STRING FOR THE OB, GUESS, ANALYSIS ....
C
      DATA obstr /'SRC FHR POB QOB TOB ZOB UOB VOB PMO'/
      data obstr2 /'MXTM MITM TDO HOVI TOCC MXGS THI TCH CDBZ'/
      data obstr3 /'CAPE CINH LI TROP PWO'/
      data obstr4 /'TKEPBL RIPBL MIXHT'/
C
C...   STRING FOR THE QUALITY MARKS ....
C
      DATA qmstr /'CAT PRC PQM QQM TQM ZQM WQM    '/
      DATA stnid /'        '/
C----------------------------------------------------
      DATA blank /' '/
      DATA equal /'='/
c     DATA namversion /'V01'/
C
      DATA bmiss /10E10/
      DATA rmiss /99999./
      call datelen(10)
      n1=0 
   10 CONTINUE
C     
C     READ VERIFICATION DATABASE VERSION AND INPUT UNIT NUMBER
C
C     Also read in logical flag (T or F) to convert virtual temperature
C     observed into actual temperature.
C     
c     READ (5,'(A)',END=160) input
c     print*,'input=',input
c     CALL ST_CLST ( input, ' ', ' ', 3, substr, num, ier )
c     namversion = substr (1)
c     CALL ST_NUMB ( substr (2), lunin(1), ier )
c     vtflg = ( substr (4) .eq. 't' .or. substr (4) .eq. 'T' )
C      PRINT '(A3,I5)', namversion, lunin(1)

      READ (5,'(A3,I5,2x,A12)',END=160) namversion,lunin(1),ens_name
      lng=len(trim(ens_name))
      print*,'lng=',lng
      ens_name_tmp=ens_name(1:lng)
      ens_name=ens_name_tmp 
c     read (5,*,end=160) namversion,lunin
c143   format (A3,I5,2x,a11)
      PRINT '(A3,I5,2x,a12)', namversion, lunin(1),ens_name
c     PRINT *, namversion, lun

c     CALL ST_RMBL ( substr (3), ens_name, lng, ier )
      nchrensemble = lng
       PRINT *, 'ENSEMBLE NAME:    ', ens_name
C     
C     READ REST OF THIS CONTROL-FILE GROUP
C     TO GET NUMBER OF THINGS TO BE VERIFIED
C     

      CALL readcntl(numodel,numfcst,numvfdate,numvfyobs,numarea,numstat,
     +            numvarbl,numlevel,numvector)

      print*, 'num model= ',numodel
      print*, 'num fcsts= ',numfcst
      print*, 'num vfydates= ',numvfdate
      print*, 'num vfyobs= ',numvfyobs
      print*, 'num area= ',numarea
      print*, 'num stats= ',numstat
      print*, 'num varbls= ',numvarbl
      print*, 'num levels= ',numlevel
      print*, 'num vectors= ',numvector
      print*, 'num  pbsmrk=', pbsmrk

      
      allocate(pfield(numfcst,numvarbl,numlevel,numarea,numvfyobs,
     +                          10,24))                              !maximum of thresholds is set to 10 to save space


    
C     
C     LET'S DO THE LOOPS
C     
C       
C     LOOP OVER ENSEMBLE MEMBER TO BE VERIFIED
C      
      DO 150 imodel = 1, numodel
      if(imodel.gt.1) lunin(imodel)=lunin(imodel-1)+1

      CALL datebf(lunin(imodel),iy,im,id,ih,idate)

      print*,'iy,im,id,ih=',iy,im,id,ih,' idate=',idate

      CALL openbf(lunin(imodel),'IN ',lunin(imodel))
       PRINT *, ' DATE OF INPUT BUFR FILE', idate, 
     +         ' UNIT=', lunin(imodel)
  150     CONTINUE


      IF (numvfdate.eq.1.and.nchrvfdate(1).eq.2) THEN
        WRITE (namvfdate(1)(1:10),'(I10)') idate
        nchrvfdate(1) = 10
c        PRINT '(" NEW NAMVFDATE =",A24)', namvfdate(1)
c        PRINT '(" CHARACTER COUNT =",I5)', nchrvfdate(1)
      END IF
C         
C       OUTER LOOP OVER VERIFYING DATE
C     
        DO 140 ivfdate = 1, numvfdate
C         
C         ZERO THE SUMS
C         
          count = 0.0
          sumdata = 0.0
          ssqdata = 0.0
          sumgrid = 0.0
          ssqgrid = 0.0
          sumprod = 0.0

	  sumensvar=0.0
	  sumfracle3=0.0
          rhnt   = 0.0
          rhet   = 0.0
          pfield = 0.

          count1  = 0.0
          count2  = 0.0
          dif     = 9999.

          print*,'ivfdate=',ivfdate
C         
C         HERE WE GET A BUFR FILE AND START LOOP OVER EACH OB
C         

          do while (ireadns(lunin(1),subset,jdate,iret).eq.0)
              subset_prev_mem     = subset
              jdate_prev_mem      = jdate

           DO 151 imodel = 2, numodel
            if(ireadns(lunin(imodel),subset,jdate,iret).ne.0) then
               write(alunin,'(i2)') lunin(imodel)
               write(alunin2,'(i2)') lunin(imodel-1)
               msg='NO BUFR RECORDS ON UNIT' // alunin
               call errmsg(msg)
               call errexit(1)
            endif
c------------- check if all ens members have the same hdr,nlev,headra--
               if(subset.ne.subset_prev_mem)          goto 501
               if(jdate.ne.jdate_prev_mem)            goto 501
               goto 601
 501           continue
                 write(alunin,'(i2)') lunin(imodel)
                 write(alunin2,'(i2)') lunin(imodel-1)
                 msg='SUBSET OR JDATE ON UNIT ' // alunin //
     *                   ' IS DIFFERENT FROM UNIT ' // alunin2
                 call errmsg(msg)
                 call errexit(1)
 601           continue
c----------------------------------------------------------------------
              subset_prev_mem     = subset
              jdate_prev_mem      = jdate
c           nsub = nmsub(lunin(imodel))
c           PRINT *, 'SUBSET,NSUB,JDATE ', lunin(imodel),subset, nsub, jdate
  151      CONTINUE

c           nsub = nmsub(lunin(1))
c           PRINT *, 'SUBSET,JDATE ', subset, jdate
            ichk = 0
            DO 20 iob = 1, numvfyobs
	      lstart (iob) = 1
	      lstop (iob) = 0
              IF (namvfyobs(iob)(:6).eq.'ANYSFC') THEN
                IF (subset(:6).eq.'ADPSFC'.or.subset(:6).eq.'SFCSHP'.or.
     +                      subset(:6).eq.'ADPUPA'.or.subset(:6).eq.
     +                      'PROFLR') THEN
			    ichk = 1
			    IF (subset(:6).eq.'ADPSFC'.or.
     +			        subset(:6).eq.'SFCSHP') lstop (iob) = 1
		END IF
              ELSE IF (namvfyobs(iob)(:6).eq.'ANYAIR') THEN
                IF (subset(:6).eq.'AIRCAR'.or.subset(:6).eq.'AIRCFT') 
     +                      ichk = 1
	      ELSE IF (namvfyobs(iob)(:6).eq.'ONLYSF') THEN
                IF (subset(:6).eq.'ADPSFC'.or.subset(:6).eq.'SFCSHP')
     +		THEN
	            ichk = 1
		    lstart (iob) = 2
		END IF
              ELSE IF (subset(:6).eq.namvfyobs(iob)(:6)) THEN
                ichk = 1
		IF (subset(:6).eq.'ADPSFC'.or.
     +		    subset(:6).eq.'SFCSHP') lstop (iob) = 1
              END IF
   20       CONTINUE

              CALL ufbint(lunin(1),hdr,10,1,nlev,headr)
              hdr_prev_mem     = hdr
              nlev_prev_mem    = nlev
              headr_prev_mem   = headr

              DO 152 imodel = 2, numodel

               CALL ufbint(lunin(imodel),hdr,10,1,nlev,headr)

c------------- check if all ens members have the same hdr,nlev,headr --

               do i=1,10
                if (hdr(i).ne.hdr_prev_mem(i)) then
                 goto 502
                endif
               enddo
               if (nlev.ne.nlev_prev_mem) then
                goto 502
               endif
               if (headr.ne.headr_prev_mem) then
                goto 502
               endif
               goto 602
 502           continue
                 write(alunin,'(i2)') lunin(imodel)
                 write(alunin2,'(i2)') lunin(imodel-1)
                 msg='HDR OR NLEV OR HEADR ON UNIT' // alunin
     *                 // 'IS DIFFERENT FROM UNIT' // alunin2
                 call errmsg(msg)
                 call errexit(1)
 602           continue
c-----------------------------------------------------------------------
              hdr_prev_mem     = hdr
              nlev_prev_mem    = nlev
              headr_prev_mem   = headr
               
  152        CONTINUE
c             PRINT*, 'NLEV,HDR ',nlev,hdr
              dhr = hdr(4)
              adate = jdate
C             
C             INNER LOOP FOR ACTUAL OB-TYPE TO BE VERIFIED
C             
              DO 80 iob = 1, numvfyobs
                IF (isitob(subset,hdr(6),iob).eq.0) THEN
C                 
C                 INNER LOOP OVER VERIFYING AREA
C                 
                  DO 70 iar = 1, numarea
                    IF (inarea(1,stnid,hdr(2),hdr(3),iar,rm1,rm2)
     +			.eq.0) THEN

                    DO 153 imodel = 1, numodel

                      CALL ufbin3(lunin(imodel),obs,9,255,mxb,
     +                            nlev,nevn,obstr)
           if((subset(:6).eq.'ADPSFC'.or.subset(:6).eq.'SFCSHP').
     *        and.obs(9,1,nevn).eq.bmiss)
     *       then
             obs(9,2,nevn)=bmiss
           endif
                      oobs(:,:,:,imodel)=obs(:,:,:)
                  
c                     call ufbin3(lunin(imodel),obs2,10,255,mxb,
c    *                            lnev,nevn,obstr2)
c                     oobs2(:,:,:,imodel)=obs2(:,:,:)
c               if(subset(:6).eq.'ADPUPA'.or.subset(:6).eq.
c    *             'GPSIPW') then
c                     obs4=bmiss
c                     nlev2=nlev
c                     nlev=1
c                     call ufbin3(lunin(imodel),obs4,6,2,mxb,2,
c    *                  nevn,obstr3)
c                     oobs4(:,:,:,imodel)=obs4(:,:,:)
c                     print*,'nlev=',nlev
c                     print*,'after special ufbin3,subset=',subset
c                     nlev=nlev2
c                endif

              if(subset(:6).eq.'ADPUPA') then
                      obs5=bmiss
                      nlev2=nlev
                      nlev=1
c                     print*,'before ufbin3, lunin=',lunin(imodel)
                      call ufbin3(lunin(imodel),obs5,3,2,mxb,nlev,nevn,
     *                  obstr4)
c                     print*,'after ufbin3, lunin=',lunin(imodel)
                      nlev=nlev2
                      oobs5(:,:,:,imodel)=obs5(:,:,:)
              endif

                      CALL ufbint(lunin(imodel),qms,7,255,nlev,qmstr)
                      qqms(:,:,imodel)=qms(:,:)

  153               CONTINUE
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
                         IF (inlayer(subset,probs,kat,ilv).eq.0) THEN
C                           
C                           INNER LOOP OVER VARIABLE 
C                           
                            DO 40 ivr = 1, numvarbl
C                             
C                             INNER LOOP OVER NUMBER OF FORECAST HOUR (EVENT)
C                             
c      icnt=0
                              DO 30 ifh = 1, numfcst

c                               write(*,*) 'ivr=',ivr, 'ifh=',ifh

                                obs(:,:,:)=oobs(:,:,:,1)
                                obs2(:,:,:)=oobs2(:,:,:,1)
                                obs4(:,:,:)=oobs4(:,:,:,1)
                                obs5(:,:,:)=oobs5(:,:,:,1)
                                qms(:,:)=qqms(:,:,1)

                                IF (igotdata(obsvaln(1),forcstn(1),
     +                                       obs,obs2,obs4,obs5,qms,
     +                                   nevn,nlv,1,ifh,ivr,ilv,
     +                                       iob,rm1,rm2,subset,stnid)
     +                                       .eq.0) THEN

                                    DO 154 imodel = 2, numodel
            write(alunin,'(i2)') lunin(imodel)
            write(alunin2,'(i2)') lunin(imodel-1)
                                     obs(:,:,:)=oobs(:,:,:,imodel)
                                     obs2(:,:,:)=oobs2(:,:,:,imodel)
                                     obs4(:,:,:)=oobs4(:,:,:,imodel)
                                     obs5(:,:,:)=oobs5(:,:,:,imodel)
                                     qms(:,:)=qqms(:,:,imodel)

                                     if (igotdata(obsvaln(imodel),
     +                                              forcstn(imodel),
     +                                       obs,obs2,obs4,obs5,qms,
     +                                       nevn, nlv,imodel,ifh,ivr
     +                                             ,ilv, iob,rm1,rm2
     *                                             ,subset,stnid)
     +                                              .ne.0) then
C            goto 30
            write(alunin,'(i2)') lunin(imodel)
            msg='IGOTDATA IS -1 FOR UNIT' // alunin // ' - skipping ...'
                                         call errmsg(msg)
            goto 30
                                         call errexit(1)
                                      endif

c---------------------- check if obs is the same in all files -------

                                        if(obsvaln(imodel).ne.
     +                                     obsvaln(imodel-1))then
                              write(alunin,'(i2)') lunin(imodel)
                              write(alunin2,'(i2)') lunin(imodel-1)
                              write(aobs,'(f8.2)') obsvaln(imodel)
                              write(aobs2,'(f8.2)') obsvaln(imodel-1)
                                           msg=
     *                                     'OBS IS DIFFERENT:' //
     *                                     aobs2 //
     *                                     'ON UNIT=' // alunin2
     +                                     // aobs //
     +                                     'ON UNIT=' // alunin //
     +                                     ' STNID= ' // stnid
                                           call errmsg(msg)
                                           call errexit(1)
                                        endif
c--------------------------------------------------------------------
  154                                CONTINUE

Chuang: adding ensemble variance and fraction to ESL1L2

                     fcstsave = forcstn    !Binbin Add for PBS

c                     write(*,*) 'fcstsave=',(fcstsave(im),im=1,numodel)

                     DO 200 ist=1,numstat

                     IF(namstat(ist).eq.'ESL1L2') THEN
                                     obsval=0.
                                     forcst=0.
                    

                                     DO 155 imodel = 1, numodel
                                       obsval=obsval+obsvaln(imodel)
                                       forcst=forcst+forcstn(imodel)
  155                                CONTINUE

                                  obsval=obsval/numodel
                                  forcst=forcst/numodel
				  fcstva=0.0
				  fcstdiffsq=0.0
				  do imodel = 1, numodel
				    fcstdiffsq=fcstdiffsq+(forcstn(imodel)-
     +                                       forcst)**2.0
                                  end do  
				  fcstva=(fcstdiffsq/float(numodel-1))**0.5
				  fracle3=0.0   
				  do imodel = 1, numodel
				    fcstdiff= abs(forcstn(imodel)-forcst)
				    if(fcstdiff .ge. (3.0*fcstva))fracle3=
     +				        fracle3+1.0
                                  end do 
				  fracle3=fracle3/float(numodel)
c				  print*,'Debug:fcstva fracle3=',fcstva,fracle3 				       
c      icnt=icnt+1
                                  sumd = sumdata(ifh,ivr,ilv,iar,iob) *
     +                                        count(ifh,ivr,ilv,iar,iob)
     +                                        + obsval
                                  ssqd = ssqdata(ifh,ivr,ilv,iar,iob) *
     +                                        count(ifh,ivr,ilv,iar,iob)
     +                                        + obsval * obsval
                                  sumg = sumgrid(ifh,ivr,ilv,iar,iob) *
     +                                        count(ifh,ivr,ilv,iar,iob)
     +                                        + forcst
                                  ssqg = ssqgrid(ifh,ivr,ilv,iar,iob) *
     +                                        count(ifh,ivr,ilv,iar,iob)
     +                                        + forcst * forcst
                                  prod = forcst * obsval
                                  sump = sumprod(ifh,ivr,ilv,iar,iob) *
     +                                        count(ifh,ivr,ilv,iar,iob)
     +                                        + prod
                                  sume=sumensvar(ifh,ivr,ilv,iar,iob)*
     +                                        count(ifh,ivr,ilv,iar,iob)
     +                                        +fcstva**2.0
                                  sumf=sumfracle3(ifh,ivr,ilv,iar,iob)*
     +                                        count(ifh,ivr,ilv,iar,iob)
     +                                        +fracle3                               
                                  count(ifh,ivr,ilv,iar,iob) = 
     +                                        count(ifh,ivr,ilv,iar,iob)
     +                                        + 1.0
c				  print*,'Debug:sume sumf=',sume,sumf 
                                  sumdata(ifh,ivr,ilv,iar,iob) = sumd /
     +                                        count(ifh,ivr,ilv,iar,iob)
                                  ssqdata(ifh,ivr,ilv,iar,iob) = ssqd /
     +                                        count(ifh,ivr,ilv,iar,iob)
                                  sumgrid(ifh,ivr,ilv,iar,iob) = sumg /
     +                                        count(ifh,ivr,ilv,iar,iob)
                                  ssqgrid(ifh,ivr,ilv,iar,iob) = ssqg /
     +                                        count(ifh,ivr,ilv,iar,iob)
                                  sumprod(ifh,ivr,ilv,iar,iob) = sump /
     +                                        count(ifh,ivr,ilv,iar,iob)
                                  sumensvar(ifh,ivr,ilv,iar,iob) = 
     +				          sume/count(ifh,ivr,ilv,iar,iob)
                                  sumfracle3(ifh,ivr,ilv,iar,iob) = sumf
     +				              /count(ifh,ivr,ilv,iar,iob)				  
c 				  print*,'Debug:sumensvar,sumfracle3=',
c    +				    sumensvar(ifh,ivr,ilv,iar,iob)
c    +                              ,sumfracle3(ifh,ivr,ilv,iar,iob)
                     ELSEIF(namstat(ist).eq.'RHNT') THEN
                                  DO 156 imodel = 1, numodel
                                    dif(imodel)=abs(obsvaln(imodel)
     +                                              -forcstn(imodel))
  156                             CONTINUE
                    ibest_model=minloc(dif, dim=1, mask = dif.ne.9999.)
                          rhnt(ifh,ivr,ilv,iar,iob,ibest_model)=
     +                     rhnt(ifh,ivr,ilv,iar,iob,ibest_model)
     +                     + 1.0
                                  count1(ifh,ivr,ilv,iar,iob) = 
     +                               count1(ifh,ivr,ilv,iar,iob)
     +                               + 1.0
                     ELSEIF(namstat(ist).eq.'RHET') THEN
                       frcst_rng=9999.

                       DO 158 imodel = 1, numodel
                        ii=minloc(forcstn,dim=1,mask=forcstn.ne.9999.)
                        frcst_rng(imodel)=forcstn(ii)
                        forcstn(ii)=9999.
  158                  CONTINUE

                       DO 162 imodel = 1, numodel
                        if(obsvaln(1).lt.frcst_rng(imodel)) then
                          rhet(ifh,ivr,ilv,iar,iob,imodel)=
     +                    rhet(ifh,ivr,ilv,iar,iob,imodel)
     +                    + 1.0
                          goto 163
                       endif
c     print*,'===',imodel,obsvaln(1),frcst_rng(imodel)
  162                 CONTINUE
  163                 CONTINUE
                    if(obsvaln(1).ge.frcst_rng(numodel)) then
                          rhet(ifh,ivr,ilv,iar,iob,numodel+1)=
     +                     rhet(ifh,ivr,ilv,iar,iob,numodel+1)
     +                     + 1.0
                    endif
c           print*,'rhet===',rhet(ifh,ivr,ilv,iar,iob,imodel)
                                  count2(ifh,ivr,ilv,iar,iob) = 
     +                               count2(ifh,ivr,ilv,iar,iob)
     +                               + 1.0


       ELSE IF(namstat(ist).eq.'PBS' .and. pbsmrk(ivr).gt.0) THEN


         !Loop for all thresholds for the irv-th variable
         do 170 ith=1,pbsmrk(ivr)


          pfield(ifh,ivr,ilv,iar,iob,ith,1)= 
     +                  pfield(ifh,ivr,ilv,iar,iob,ith,1)+1.0

          !(1) Computer prob for different thresholds at each station
          num=0
          do imodel = 1, numodel
           if(op(ivr).eq.'>') then 
            if(fcstsave(imodel).gt.rpbsthr(ivr,ith)) num=num+1
           else if(op(ivr).eq.'<') then
            if(fcstsave(imodel).lt.rpbsthr(ivr,ith)) num=num+1                  
           else
            msg='Wrong setting in PBS threshold!'
            call errmsg(msg)
            call errexit(1)
           end if
          end do


          prob = (1.0*num)/numodel 

          
          !(2) Then to see&count the prob belongs to which probability range pfield(2~24)

          if(prob.eq.1.0) then

           pfield(ifh,ivr,ilv,iar,iob,ith,12)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,12)+1.    
      
          else if(prob.gt.0.0 .and. prob.le.0.1) then 
            pfield(ifh,ivr,ilv,iar,iob,ith,2)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,2)+1.
          else if(prob.gt.0.1 .and. prob.le.0.2) then
            pfield(ifh,ivr,ilv,iar,iob,ith,3)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,3)+1.
          else if(prob.gt.0.2 .and. prob.le.0.3) then
            pfield(ifh,ivr,ilv,iar,iob,ith,4)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,4)+1.
          else if(prob.gt.0.3 .and. prob.le.0.4) then
            pfield(ifh,ivr,ilv,iar,iob,ith,5)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,5)+1.
          else if(prob.gt.0.4 .and. prob.le.0.5) then
            pfield(ifh,ivr,ilv,iar,iob,ith,6)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,6)+1.
          else if(prob.gt.0.5 .and. prob.le.0.6) then
            pfield(ifh,ivr,ilv,iar,iob,ith,7)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,7)+1.
          else if(prob.gt.0.6 .and. prob.le.0.7) then
            pfield(ifh,ivr,ilv,iar,iob,ith,8)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,8)+1.
          else if(prob.gt.0.7 .and. prob.le.0.8) then
            pfield(ifh,ivr,ilv,iar,iob,ith,9)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,9)+1.
          else if(prob.gt.0.8 .and. prob.le.0.9) then
            pfield(ifh,ivr,ilv,iar,iob,ith,10)= 
     +         pfield(ifh,ivr,ilv,iar,iob,ith,10)+1.
          else if(prob.gt.0.9 .and. prob.le.1.0) then
            pfield(ifh,ivr,ilv,iar,iob,ith,11)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,11)+1.
          end if

         if(op(ivr).eq.'>') then
           if(prob.eq.1.0 .and. obsvaln(1).gt.rpbsthr(ivr,ith)) then     
             pfield(ifh,ivr,ilv,iar,iob,ith,24)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,24)+1. 

           else if(prob.eq.0.0.and.obsvaln(1).gt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,14)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,14)+1.
           else if(prob.gt.0.0.and.prob.le.0.1.and.
     +             obsvaln(1).gt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,14)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,14)+1.  
           else if(prob.gt.0.1.and.prob.le.0.2.and.
     +             obsvaln(1).gt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,15)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,15)+1.
           else if(prob.gt.0.2.and.prob.le.0.3.and.
     +             obsvaln(1).gt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,16)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,16)+1.
           else if(prob.gt.0.3.and.prob.le.0.4.and.
     +             obsvaln(1).gt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,17)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,17)+1.
           else if(prob.gt.0.4.and.prob.le.0.5.and.
     +             obsvaln(1).gt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,18)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,18)+1.
           else if(prob.gt.0.5.and.prob.le.0.6.and.
     +             obsvaln(1).gt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,19)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,19)+1.
           else if(prob.gt.0.6.and.prob.le.0.7.and.
     +             obsvaln(1).gt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,20)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,20)+1.
           else if(prob.gt.0.7.and.prob.le.0.8.and.
     +             obsvaln(1).gt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,21)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,21)+1.
           else if(prob.gt.0.8.and.prob.le.0.9.and.
     +             obsvaln(1).gt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,22)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,22)+1.
           else if(prob.gt.0.9.and.prob.le.1.0.and.
     +             obsvaln(1).gt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,23)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,23)+1.
           end if
         else 
           if(prob.eq.1.0 .and. obsvaln(1).lt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,24)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,24)+1.
           else if(prob.eq.0.0.and.obsvaln(1).lt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,14)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,14)+1.
           else if(prob.gt.0.0.and.prob.le.0.1.and.
     +             obsvaln(1).lt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,14)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,14)+1.
           else if(prob.gt.0.1.and.prob.le.0.2.and.
     +             obsvaln(1).lt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,15)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,15)+1.
           else if(prob.gt.0.2.and.prob.le.0.3.and.
     +             obsvaln(1).lt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,16)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,16)+1.
           else if(prob.gt.0.3.and.prob.le.0.4.and.
     +             obsvaln(1).lt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,17)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,17)+1.
           else if(prob.gt.0.4.and.prob.le.0.5.and.
     +             obsvaln(1).lt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,18)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,18)+1.
           else if(prob.gt.0.5.and.prob.le.0.6.and.
     +             obsvaln(1).lt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,19)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,19)+1.
           else if(prob.gt.0.6.and.prob.le.0.7.and.
     +             obsvaln(1).lt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,20)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,20)+1.
           else if(prob.gt.0.7.and.prob.le.0.8.and.
     +             obsvaln(1).lt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,21)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,21)+1.
           else if(prob.gt.0.8.and.prob.le.0.9.and.
     +             obsvaln(1).lt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,22)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,22)+1.
           else if(prob.gt.0.9.and.prob.le.1.0.and.
     +             obsvaln(1).lt.rpbsthr(ivr,ith)) then
             pfield(ifh,ivr,ilv,iar,iob,ith,23)=
     +         pfield(ifh,ivr,ilv,iar,iob,ith,23)+1.
           end if
         end if
  
c        write(*,*) 'fcstsave=',(fcstsave(imodel),imodel=1,numodel),
c     +   ' obsvaln=',obsvaln(1), op(ivr), rpbsthr(ivr,ith), 
c     +   ' prob=',prob
c        write(*,'(24f7.4)') 
c     +     (pfield(ifh,ivr,ilv,iar,iob,ith,i),i=1,24) 

 170    continue

       END IF  !stat type
  200     CONTINUE
                                END IF !igotodata
C                             
C                             END INNER LOOP OVER NUMBER OF FORECAST HOUR (EVENT)
   30                         CONTINUE
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
c           END DO
C         END DO WHILE LOOP OVER BUFR MESSAGE
          END DO

         write(*,*) "Read Bufr and data accumulation done!" 

          do ifh = 1, numfcst
           do ivr = 1, numvarbl
            do ilv = 1, numlevel
             do iar = 1, numarea
               do iob = 1, numvfyobs
                 do ith = 1, pbsmrk(ivr)
                  do ip = 2,24
                    if (pfield(ifh,ivr,ilv,iar,iob,ith,1).gt.0) then
                     pfield(ifh,ivr,ilv,iar,iob,ith,ip) = 
     +               pfield(ifh,ivr,ilv,iar,iob,ith,ip)/
     +               pfield(ifh,ivr,ilv,iar,iob,ith,1)
                    end if
                  end do
                end do
               end do
              end do
             end do
            end do
          end do

 
C         
C         NOW IS THE TIME TO WRITE OUT THE STAT RECORDS
C         
c         ist = 1
         DO 300 ist=1,numstat
     
         if (namstat(ist).ne.'PBS' ) then

          DO 130 iob = 1, numvfyobs
            DO 120 iar = 1, numarea
              DO 110 ifh = 1, numfcst
                numv = numvarbl
                IF (numvector.gt.0) numv = numvarbl - 1
                DO 100 ivr = 1, numv
                  DO 90 ilv = 1, numlevel
                    rtest = count(ifh,ivr,ilv,iar,iob)
C                    print *, 'COUNT(IFH,IVR,ILV,IAR,IOB)=',RTEST
                    IF (count(ifh,ivr,ilv,iar,iob).gt.0.0) THEN
                      iend = 3
                      vdbhdr132(1:iend) = namversion(:3)
                      iend = iend + 1
                      vdbhdr132(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrensemble
                      vdbhdr132(istrt:iend) = 
     +                            ens_name(:nchrensemble)
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
                      ivend = iend
                      vdbhdr132(istrt:iend) = 
     +                            namstat(ist)(:nchrstat(ist))

c ESL1L2 needs to have ensemble member numbers info attached 
                      IF(namstat(ist).eq.'ESL1L2') THEN
                       iend = iend + 1
                       vdbhdr132(iend:iend) = '/'
		       istrt = iend + 1
		       if(numodel .lt. 10)then
		        iend = iend + 1
		        write(vdbhdr132(istrt:iend),'(i1)')numodel
		       else if(numodel .lt. 100)then
		        iend = iend + 2 
		        write(vdbhdr132(istrt:iend),'(i2)')numodel
		       else
		        iend = iend + 3 
		        write(vdbhdr132(istrt:iend),'(i3)')numodel
		       end if 
		      END IF

                      iend = iend + 1
                      vdbhdr132(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrvarbl(ivr)
                      vdbhdr132(istrt:iend) = 
     +                            namvarbl(ivr)(:nchrvarbl(ivr))
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

                 IF(namstat(ist).eq.'ESL1L2') THEN

                      print *, 'IVR=',IVR,' and NUMVECTOR=',NUMVECTOR
                      IF (ivr.ne.numvector) THEN
                        PRINT'(A,F7.0,5E18.9)',VDBHDR132(:IEND) 
                        WRITE (50,1000) vdbhdr132(:iend), 
     +                              count(ifh,ivr,ilv,iar,iob), 
     +                              sumgrid(ifh,ivr,ilv,iar,iob), 
     +                              sumdata(ifh,ivr,ilv,iar,iob), 
     +                              sumprod(ifh,ivr,ilv,iar,iob), 
     +                              ssqgrid(ifh,ivr,ilv,iar,iob), 
     +                              ssqdata(ifh,ivr,ilv,iar,iob),
     +                              sumensvar(ifh,ivr,ilv,iar,iob), 
     +                              sumfracle3(ifh,ivr,ilv,iar,iob)      
C 1000                   FORMAT (A,F7.0,7E18.9)
 1000                   FORMAT (A,F7.0,7(1X,1PE13.6))
                      END IF
C                     
                      IF (ivr.eq.numvector) THEN
                        sump = sumprod(ifh,ivr,ilv,iar,iob) + 
     +                              sumprod(ifh,ivr+1,ilv,iar,iob)
                        sumg = ssqgrid(ifh,ivr,ilv,iar,iob) + 
     +                              ssqgrid(ifh,ivr+1,ilv,iar,iob)
                        sumd = ssqdata(ifh,ivr,ilv,iar,iob) + 
     +                              ssqdata(ifh,ivr+1,ilv,iar,iob)
                        vdbhdr132(ivstrt+1:ivstrt+1) = namversion(1:)
c                       PRINT'(A,F7.0,7E18.9)',VDBHDR132(:IEND)
                        WRITE (50,1100) vdbhdr132(:iend), 
     +                              count(ifh,ivr,ilv,iar,iob), 
     +                              sumgrid(ifh,ivr,ilv,iar,iob), 
     +                              sumgrid(ifh,ivr+1,ilv,iar,iob), 
     +                              sumdata(ifh,ivr,ilv,iar,iob), 
     +                              sumdata(ifh,ivr+1,ilv,iar,iob), 
     +                              sump, sumg, sumd,
     +                              sumensvar(ifh,ivr,ilv,iar,iob),
     +                              sumensvar(ifh,ivr+1,ilv,iar,iob),
     +                              sumfracle3(ifh,ivr,ilv,iar,iob),
     +                              sumfracle3(ifh,ivr+1,ilv,iar,iob)
C 1100                   FORMAT (A,F7.0,11E18.9)
 1100                   FORMAT (A,F7.0,11(1X,1PE13.6))
                      END IF

                     ELSEIF(namstat(ist).eq.'RHNT') THEN
                      vdbhdr132(ivstrt:ivend) = 
     +                            namstat(ist)(:nchrstat(ist))
                        DO 157 imodel = 1, numodel
                        rhnt(ifh,ivr,ilv,iar,iob,imodel) =
     +                       rhnt(ifh,ivr,ilv,iar,iob,imodel)/
     +                       count1(ifh,ivr,ilv,iar,iob) 
  157                             CONTINUE
                        WRITE (50,1200) vdbhdr132(:iend), 
     +                   count1(ifh,ivr,ilv,iar,iob), 
     +    (rhnt(ifh,ivr,ilv,iar,iob,imodel),imodel=1,numodel-1) 
C 1200                   FORMAT (A,F7.0,<numodel-1>E18.9)
 1200                   FORMAT (A,F7.0,<numodel-1>(1X,F7.5))
                     ELSEIF(namstat(ist).eq.'RHET') THEN
                      vdbhdr132(ivstrt:ivend) = 
     +                            namstat(ist)(:nchrstat(ist))
                        DO 167 imodel = 1, numodel+1
                        rhet(ifh,ivr,ilv,iar,iob,imodel) =
     +                       rhet(ifh,ivr,ilv,iar,iob,imodel)/
     +                       count2(ifh,ivr,ilv,iar,iob) 
  167                             CONTINUE
                        WRITE (50,1300) vdbhdr132(:iend), 
     +                   count2(ifh,ivr,ilv,iar,iob), 
     +    (rhet(ifh,ivr,ilv,iar,iob,imodel),imodel=1,numodel) 
C 1300                   FORMAT (A,F7.0,<numodel>E18.9)
 1300                   FORMAT (A,F7.0,<numodel>(1X,F7.5))


                     END IF  !stat type
                    END IF  !count
   90             CONTINUE
  100           CONTINUE
  110         CONTINUE
  120       CONTINUE
  130     CONTINUE

       else    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   
c   PBS case:
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          DO 2001 iob = 1, numvfyobs
            DO 2002 iar = 1, numarea
              DO 2003 ifh = 1, numfcst
                numv = numvarbl
                IF (numvector.gt.0) numv = numvarbl - 1
                DO 2004 ivr = 1, numv
 
                 DO 2005 ith = 1, pbsmrk(ivr)      !only for those variables with pbsmrk(ivr) > 0
                
                  DO  2006 ilv = 1, numlevel
                    IF (pfield(ifh,ivr,ilv,iar,iob,ith,1).gt.0.0) THEN
                      iend = 3
                      vdbhdr420(1:iend) = namversion(:3)
                      iend = iend + 1
                      vdbhdr420(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrensemble
                      vdbhdr420(istrt:iend) =
     +                            ens_name(:nchrensemble)
                      iend = iend + 1
                      vdbhdr420(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrfcst(ifh)
                      vdbhdr420(istrt:iend) =
     +                            namfcst(ifh)(:nchrfcst(ifh))
                      iend = iend + 1
                      vdbhdr420(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrvfdate(ivfdate)
                      vdbhdr420(istrt:iend) =
     +                            namvfdate(ivfdate)(:
     +                            nchrvfdate(ivfdate))

                      iend = iend + 1
                      vdbhdr420(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrvfyobs(iob)
                      vdbhdr420(istrt:iend) =
     +                            namvfyobs(iob)(:nchrvfyobs(iob))
                      iend = iend + 1
                      vdbhdr420(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrarea(iar)
                      vdbhdr420(istrt:iend) =
     +                            namarea(iar)(:nchrarea(iar))
                      iend = iend + 1
                      vdbhdr420(iend:iend) = blank
                      istrt = iend + 1
                      ivstrt = istrt
                      iend = iend + nchrstat(ist)
                      ivend = iend
                      vdbhdr420(istrt:iend) =
     +                            namstat(ist)(:nchrstat(ist))

                       istrt = iend + 1
                       iend = iend + 5
                       vdbhdr420(istrt:iend) = '_ENS:'
                      
                       istrt = iend + 1
                       if(numodel .lt. 10)then
                        iend = iend + 1
                        write(vdbhdr420(istrt:iend),'(i1)')numodel
                       else if(numodel .lt. 100)then
                        iend = iend + 2
                        write(vdbhdr420(istrt:iend),'(i2)')numodel
                       else
                        iend = iend + 3
                        write(vdbhdr420(istrt:iend),'(i3)')numodel
                       end if

                      istrt = iend + 1
                      iend = iend + 1
                      vdbhdr420(istrt:iend) = '/'

                      istrt = iend + 1
                      iend = iend + 1
                      vdbhdr420(istrt:iend) = op(ivr)
                      
                      istrt = iend + 1
                      iend = iend + nchrpbsthr(ivr,ith)
                      print*,'pbsthr(ivr,ith)=',ivr,ith,pbsthr(ivr,ith)
                      vdbhdr420(istrt:iend) = 
     +                    pbsthr(ivr,ith)(:nchrpbsthr(ivr,ith))
                      

                      iend = iend + 1
                      vdbhdr420(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrvarbl(ivr)
                      vdbhdr420(istrt:iend) =
     +                            namvarbl(ivr)(:nchrvarbl(ivr))
                      iend = iend + 1
                      vdbhdr420(iend:iend) = blank
                      istrt = iend + 1
                      iend = iend + nchrlevel(ilv)
                      vdbhdr420(istrt:iend) =
     +                            namlevel(ilv)(:nchrlevel(ilv))
                      iend = iend + 1
                      vdbhdr420(iend:iend) = blank
                      iend = iend + 1
                      vdbhdr420(iend:iend) = equal
                      iend = iend + 1
                      vdbhdr420(iend:iend) = blank


                      WRITE (50,1555) vdbhdr420(:iend),
     +                  int(pfield(ifh,ivr,ilv,iar,iob,ith,1)),
     +                   (pfield(ifh,ivr,ilv,iar,iob,ith,ip),ip=2,24)
1555                  FORMAT (A,I10, 23(1X,F7.5))

                   END IF
2006            continue
2005            continue
2004            continue
2003            continue
2002            continue
2001            continue

       end if  ! end of PBS

C         END LOOP OVER STATISTIC TYPE
  300     CONTINUE

       

C         END OUTER LOOP OVER VERIFYING DATE
  140     CONTINUE

         do imodel=1,numodel
          CALL closbf(lunin(imodel))
         enddo

          deallocate(pfield)

C         GO BACK AND DO IT ALL AGAIN   WHOOPEE!
          GO TO 10
C23456789012345678901234567890123456789012345678901234567890123456789012
  160     CONTINUE
          elapsed=timef()
          PRINT*, 'TOTAL TIME (min)= ',elapsed/6.0e4
          STOP
          END
