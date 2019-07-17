      SUBROUTINE readcntl(numodel,numfcst,numvfdate,numvfyobs,numarea,
     +            numstat,numvarbl,numlevel,numvector)
C
      INCLUDE 'parm.inc'
C
      DIMENSION nchrmodel(maxmod), nchrfcst(mxfcst), nchrvfdate(mxdate),
     +            nchrvfyobs(maxobs), nchrarea(mxarea), 
     +            nchrstat(mxstat), nchrvarbl(mxvrbl), 
     +            nchrlevel(maxlvl)
      CHARACTER*24 namodel(maxmod), namfcst(mxfcst), 
     +            namvfdate(mxdate), namvfyobs(maxobs), 
     +            namarea(mxarea), namstat(mxstat), 
     +            namvarbl(mxvrbl), namlevel(maxlvl)
      CHARACTER*80 input
      character*80 msg
      character*10 a,b
      CHARACTER*10 substr (20)
C
      COMMON /names/ namodel, namfcst, namvfdate, namvfyobs, namarea, 
     +            namstat, namvarbl, namlevel
      COMMON /nchrs/ nchrmodel, nchrfcst, nchrvfdate, nchrvfyobs, 
     +            nchrarea, nchrstat, nchrvarbl, nchrlevel
      LOGICAL	  vtflg
      COMMON /cnvrsns/ vtflg, nmbgrd (maxmod), concon (maxmod),
     +		       cenlon (maxmod)


      integer pbsmrk(mxvrbl)
      CHARACTER*24 op(mxvrbl),pbsthr(mxvrbl,10)
      integer  nchrpbsthr(mxvrbl,10)
      real     rpbsthr(mxvrbl,10)
      COMMON /pbs/ pbsmrk,op,pbsthr,nchrpbsthr,rpbsthr

      CHARACTER*1 blank
      DATA blank /' '/

      CHARACTER*80  fhovar(mxthr+2)
      CHARACTER*80  fhovar1,tempstr
      INTEGER       idelim(mxthr+1)

      integer numthr(mxvrbl)
      real thresh(mxvrbl,mxthr)
      character*1 arrow


       pbsmrk = 0
       pbsthr = ''
       nchrpbsthr = 0
       rpbsthr = 0.0
C
C   READ NUMBER OF MODELS TO BE VERIFIED, first model name, and
C   optional wind rotation grid # for the model.
C
C   Note:  no adherence to format is necessary for this input.
C
c       READ (5,'(A)') input
c       CALL ST_CLST ( input, ' ', ' ', 3, substr, num, ier )
cCALL ST_NUMB ( substr (1), numodel, ier )
cCALL ST_RMBL ( substr (2), namodel (1), lng, ier )
cnchrmodel (1) = lng
cIF ( substr (3) .ne. ' ' ) THEN
c   CALL ST_NUMB ( substr (3), nmbgrd (1), ier )
cELSE
c   nmbgrd (1) = -1
c       END IF

      read(5,'(i5,2x,a24)') numodel,namodel(1)
      nchrmodel(1) = LEN(TRIM(namodel(1)))
      print*,'numodel=',numodel
      nmbgrd(1) = -1
C
C     NUMBER OF VERIFYING MODELS IS LIMITED TO MAXMOD
      IF (numodel.gt.maxmod) THEN
        write(a,'(i2)') maxmod
        write(b,'(i2)') numodel
        msg="  NUMBER OF VERIFYING MODELS EXCEEDS LIMIT OF " // a // b
        call errmsg(msg)
        call errexit(1)
      END IF
C     
C     READ NUMODEL MODEL PNEMONICS AND GET CHARACTER COUNT FOR EACH
C     Also read in optional wind rotation flag (false if missing).
C     
C   Note:  no adherence to format is necessary for this input.
C
c     DO 30 n = 1, numodel
c       IF ( n .ne. 1 ) THEN
c    READ (5,'(A)') input
c     	    CALL ST_CLST ( input, ' ', ' ', 2, substr, num, ier )
c    CALL ST_RMBL ( substr (1), namodel (n), lng, ier )
c    nchrmodel (n) = lng
c    IF ( substr (2) .ne. ' ' ) THEN
c       CALL ST_NUMB ( substr (2), nmbgrd (n), ier )
c    ELSE
c       nmbgrd (n) = -1
c           END IF
cEND IF
cCALL SETMODEL(N,NAMODEL(N),NCHRMODEL(N))
cPRINT *, ' GRD # for wind rotation = ', nmbgrd (n)
c  30 CONTINUE

      IF (numodel.gt.1) THEN
         DO i = 2,numodel
            READ(5,'(a24)') namodel(i)
            nmbgrd(i) = -1
         ENDDO
      ENDIF

      DO 30 n = 1, numodel
        CALL SETMODEL(N,NAMODEL(N),NCHRMODEL(N))
        PRINT *, ' GRD # for wind rotation = ', nmbgrd (n)
   30 CONTINUE

        write (*,*) "read model done!"
C     
C     READ NUMBER OF FORECAST HOURS TO BE VERIFIED AND first hour.
C     
cCALL ST_READ  ( namfcst, nchrfcst, numfcst, iret )
cIF ( iret .ne. 0 ) THEN
c        call errmsg('EOF encountered reading forecast hours.')
c        call errexit(1)
cEND IF

      READ(5,*) numfcst,namfcst(1)
      nchrfcst(1) = LEN(TRIM(namfcst(1)))
C
C*      NUMBER OF FORECAST HOURS IS LIMITED TO MXFCST
        IF (numfcst.gt.mxfcst) THEN
          write(a,'(i2)') maxmod
          write(b,'(i2)') numodel
          msg="  NUMBER OF FORECAST HOURS EXCEEDS LIMIT OF" // a // b
          call errmsg(msg)
          call errexit(1)
        END IF

C     READ THE REST OF THE FCST HRS IF NECESSARY
C
      IF ( numfcst.gt.1) THEN
         DO i = 2,numfcst
C           READ(5,'(a4)') namfcst(i)
            READ(5,*) namfcst(i)
            print*,'i,namfcst(i)=',i,namfcst(i)
            nchrfcst(i) = LEN(TRIM(namfcst(i)))
         ENDDO
      ENDIF

C     
C       SET NUMFCST FORECAST HOUR PNEMONICS
C     
        DO 60 n = 1, numfcst
          CALL setfcst(n,namfcst(n),nchrfcst(n))
   60   CONTINUE

       write (*,*) "read fcst hour done!"

C     
C     READ NUMBER OF VERIFICATION DATES
C     
cCALL ST_READ  ( namvfdate, nchrvfdate, numvfdate, iret )
cIF ( iret .ne. 0 ) THEN
c        call errmsg('End of file encountered reading verifying dates')
c        call errexit(1)
cEND IF

      READ(5,*) numvfdate,namvfdate(1)
      nchrvfdate(1) = LEN(TRIM(namvfdate(1)))
C
C*      NUMBER OF VERIFYING DATES IS LIMITED TO MXDATE
        IF (numvfdate.gt.mxdate) THEN
          write(a,'(i2)') mxdate
          write(b,'(i2)') numvfdate
          msg='NUMBER OF VERIFYING DATES EXCEEDS LIMIT OF' // a // b
          call errmsg(msg)
          call errexit(1)
        END IF

C     READ THE REST OF THE VERIF DATES IF NECESSARY
C
      IF ( numvfdate.gt.1) THEN
         DO i = 2,numvfdate
            READ(5,*) namvfdate(i)
            nchrvfdate(i) = LEN(TRIM(namvfdate(i)))
         ENDDO
      ENDIF

        write (*,*) "read verification dates  done!"

C     
C       READ NUMVFDATE VERIFICATION PNEMONICS AND GET CHARACTER COUNT
C     
C     READ NUMBER OF VERIFYING OB TYPES
C     
cCALL ST_READ  ( namvfyobs, nchrvfyobs, numvfyobs, iret )
cIF ( iret .ne. 0 ) THEN
c           msg=' End of file encountered reading verifying ob types.'
c           call errmsg(msg)
c           call errexit(1)
cEND IF

       READ(5,*) numvfyobs,namvfyobs(1)
      nchrvfyobs(1) = LEN(TRIM(namvfyobs(1)))
C
C*      NUMBER OF VERIFYING OBS IS LIMITED TO MAXOBS
        IF (numvfyobs.gt.maxobs) THEN
          write(a,'(i2)') maxobs
          write(b,'(i2)') numvfyobs
          msg='NUMBER OF VERIFYING OBS EXCEEDS LIMIT OF ' // a // b
          call errmsg(msg)
          call errexit(1)
        END IF

      IF ( numvfyobs.gt.1 ) THEN
         DO i = 2,numvfyobs
            READ(5,*) namvfyobs(i)
            nchrvfyobs(i) = LEN(TRIM(namvfyobs(i)))
c           print*,'length of ob typ=',nchrvfyobs(i)
c           print*,'i,namvfyobs(i)=',i,namvfyobs(i)
         ENDDO
      ENDIF
C     
C       SET NUMVFYOBS VERIFYING OBS PNEMONICS
C     
        DO 110 n = 1, numvfyobs
          CALL setobtyp(n,namvfyobs(n),nchrvfyobs(n))
  110   CONTINUE

        write (*,*) "read obstype done!"

C     
C     READ NUMBER OF VERIFICATION AREAS
C     
cCALL ST_READ  ( namarea, nchrarea, numarea, iret )
cIF ( iret .ne. 0 ) THEN
c        call errmsg('EOF encountered reading verification areas')
c        call errexit(1)
cEND IF

      read(5,'(i5,2x,a24)') numarea,namarea(1)
      print*,'namarea(1)=',namarea(1)
      nchrarea(1) = LEN(TRIM(namarea(1)))
C
C*      NUMBER OF VERIFYING AREAS IS LIMITED TO MXAREA
        IF (numarea.gt.mxarea) THEN
          write(a,'(i2)') mxarea
          write(b,'(i2)') numarea
          msg='NUMBER OF VERIFYING AREAS EXCEEDS LIMIT OF ' // a // b
          call errmsg(msg)
          call errexit(1)
        END IF

      IF ( numarea.gt.1 ) THEN
         DO i = 2,numarea
c           READ(5,*) namarea(i)
            read(5,'(2x,a24)') namarea(i)
            print*,'i,namarea(i)=',i,namarea(i)
            nchrarea(i) = LEN(TRIM(namarea(i)))
c           print*,'length of area string=',nchrarea(i)
         ENDDO
      ENDIF
C     
C       SET NUMAREA DOMAIN PNEMONICS
C     
        DO 140 n = 1, numarea
          CALL setarea(n,namarea(n),nchrarea(n))
  140   CONTINUE

        write (*,*) "read area done!"

C     
C     READ NUMBER OF STATISTICAL SCORES
C     
cCALL ST_READ  ( namstat, nchrstat, numstat, iret )
cIF ( iret .ne. 0 ) THEN
c        call errmsg('End of file encountered reading statistic types')
c        call errexit(1)
cEND IF

      READ(5,*) numstat,namstat(1)
      print *,'numstat,namstat(1)=',numstat,namstat(1)
      nchrstat(1) = LEN(TRIM(namstat(1)))
C
C       NUMBER OF STATISTICS IS LIMITED TO MXSTAT
        IF (numstat.gt.mxstat) THEN
          write(a,'(i2)') mxstat
          write(b,'(i2)') numstat
          msg='NUMBER OF STATISTICS EXCEEDS LIMIT OF ' // a // b
          call errmsg(msg)
          call errexit(1)
        END IF

      IF (numstat.gt.1) THEN
         DO i = 2,numstat
            READ(5,*) namstat(i)
            print*,'i,namstat(i)=',i,namstat(i)
            nchrstat(i) = LEN(TRIM(namstat(i)))
c           print*,'length of stat string=',nchrstat(i)
         ENDDO
      ENDIF

       write (*,*) "read statics-type done!"

C     
C     READ VARIABLES TO BE VERIFIED
C     

c	CALL ST_READ  ( namvarbl, nchrvarbl, numvarbl, iret )
c	IF ( iret .ne. 0 ) THEN
c	    WRITE (6,*)
c     +		' End of file encountered reading parameters.'
c	    STOP
c	END IF

c       READ (5, '(A)') input
c       CALL ST_CLST( input, ' ', ' ', 10, substr, num, ier)
c       CALL ST_NUMB( substr(1), numvarbl, ier)
c       CALL ST_RMBL( substr(2), namvarbl(1), nchrvarbl(1), ier)

      READ(5,*) numvarbl, fhovar1
      fhovar(1) = fhovar1
      print*,'FIRST VARIABLE',numvarbl,fhovar1
C
C       NUMBER OF VARIABLES IS LIMITED TO MXVRBL

        IF (numvarbl.gt.mxvrbl) THEN
          write(a,'(i2)') mxvrbl
          write(b,'(i2)') numvarbl
          msg='NUMBER OF VARIABLES EXCEEDS LIMIT OF ' // a // b
          call errmsg(msg)
          call errexit(1)
        END IF

      DO i = 2, numvarbl
         READ(5,'(2x,a80)') fhovar1
         fhovar(i)=fhovar1
      ENDDO
      DO i = 1, numvarbl
         print*,'VARIABLES= ',i,fhovar(i)
         tempstr = fhovar(i)
         olen = LEN(fhovar(i))
         nlen = LEN(TRIM(fhovar(i)))
         print*,'olen,nlen= ',olen,nlen
         print*,'tempstr=',tempstr
         nindex = INDEX(tempstr(1:nlen),' ')
         print*, 'nindex= ',nindex
         IF (nindex.eq.0) THEN
            namvarbl(i)=fhovar(i)
            nchrvarbl(i)=LEN(TRIM(namvarbl(i)))
c           print*,'NO THRESHOLDS FOR VARIABLE ',namvarbl(i)
            if(i.eq.1) then
             num=2
            else
             num=1
            endif
            numthr(i) = 0
            thresh(i,:) = 0
            rpbsthr(i,:) = 0 
         ELSE
            namvarbl(i)=tempstr(1:nindex-1)
            nchrvarbl(i)=LEN(TRIM(namvarbl(i)))
            print*,'namvarbl(i)=',namvarbl(i)
            ndelim = 0
C     FIND THE LOCATIONS OF ALL THE SPACES
            DO j = 1,nlen
               IF (tempstr(j:j) == ' ') THEN
                  ndelim = ndelim + 1
                  idelim(ndelim) = j
                  print*,'ndelim= ',ndelim
                  print*,'idelim= ',idelim(ndelim)
               ENDIF
            ENDDO
            num=ndelim+1
C     USE THE LOCATIONS OF THE SPACES TO EXTRACT THE INFO
C     WE NEED ABOUT THE NUMBER OF THRESHOLDS AND WHAT THE
C     THRESHOLD VALUES ARE
            READ(tempstr(idelim(1)+1:idelim(2)-1),*) arrow
            op(i)=arrow
            print*,'arrow=',op(i)
            read(tempstr(idelim(2)+1:idelim(3)-1),*) numthr(i)
            pbsmrk(i)=numthr(i)
            print*,'i,pbsmrk(i)=',i,pbsmrk(i)
            print*,'THERE ARE ',numthr(i),
     +             ' THRESHOLDS FOR VARIABLE ',namvarbl(i)
            DO j = 1,numthr(i)-1
c           do j = 1,3
               READ(tempstr(idelim(j+2)+1:idelim(j+3)-1),*) thresh(i,j)
               pbsthr(i,j)=tempstr(idelim(j+2)+1:idelim(j+3)-1)
c              print*, 'THRESHOLD VALUES = ',thresh(i,j)
            ENDDO
            ibeg = numthr(i)+1
c           ibeg = 4
            print*,'ibeg,nlen = ',ibeg,nlen
            READ(tempstr(idelim(ibeg+1)+1:nlen),*) thresh(i,numthr(i))
            pbsthr(i,numthr(i))=tempstr(idelim(ibeg+1)+1:nlen)
            DO j = 1,numthr(i)
c           do j = 1,4
               print*, 'THRESHOLD VALUES = ',thresh(i,j)
               rpbsthr(i,j)=thresh(i,j)
               print*,'pbsthr=',i,j,pbsthr(i,j)
               nchrpbsthr(i,j)=len(trim(pbsthr(i,j)))
            ENDDO
         ENDIF
      ENDDO


        write(*,*) numvarbl,namvarbl(1),num

         if (num.gt.3) then        !has PBS threholds
c          CALL ST_RMBL( substr(3), op(1), n1, ier)
c          pbsmrk(1)=num-3
c          
c          do nx = 4, num
c           CALL ST_RMBL(substr(nx),pbsthr(1,nx-3),
c    +               nchrpbsthr(1,nx-3),ier)
c           rpbsthr(1,nx-3) = ChartoReal(pbsthr(1,nx-3))
c          end do

         write(*,*) namvarbl(1),op(1),(rpbsthr(1,nx),nx=1,pbsmrk(1))

          end if

c       write (*,*) "read 1st variable done!"


       IF(numvarbl.GE.2) THEN
         DO n=2,numvarbl
c          READ (5, '(A)') input
c          CALL ST_CLST( input, ' ', ' ', 10, substr, num, ier)
c          CALL ST_RMBL( substr(1), namvarbl(n), nchrvarbl(n), ier)

           write(*,*) namvarbl(n),num

           if (num.gt.3) then
c            CALL ST_RMBL( substr(2), op(n), n1, ier)
c            pbsmrk(n)=num-2

c            do nx = 3, num
c             CALL ST_RMBL(substr(nx),pbsthr(n,nx-2),
c    +               nchrpbsthr(n,nx-2),ier)
c             rpbsthr(n,nx-2) = ChartoReal(pbsthr(n,nx-2))
c            end do

           write(*,*) namvarbl(n),op(n),(rpbsthr(n,nx),nx=1,pbsmrk(n))

           end if

          END DO
        END IF

c       write (*,*) "read 2nd variabels done!"


C
C       NUMBER OF VARIABLES IS LIMITED TO MXVRBL
c       IF (numvarbl.gt.mxvrbl) THEN
c         write(a,'(i2)') mxvrbl
c         write(b,'(i2)') numvarbl
c         msg='NUMBER OF VARIABLES EXCEEDS LIMIT OF ' // a // b
c         call errmsg(msg)
c         call errexit(1)
c       END IF
C     
C       FIND VECTOR VARIABLE (must be last in list)
C     
        numvector = 0
	
	DO ivr = 1, numvarbl
          IF (nchrvarbl(ivr).eq.4.and.namvarbl(ivr).eq.'VWND') THEN
            numvector = ivr
	    IF ( numvector .ne. numvarbl ) THEN
              call errmsg('VWND must be last in parm list')
              call errexit(1)
	    END IF
            numvarbl = numvarbl + 1
            namvarbl(ivr+1) = namvarbl(ivr) (1:1)
            nchrvarbl(ivr+1) = 1
          END IF
	END DO
C     
C     READ NUMBER OF LEVEL DESCRIPTIONS
C     
cCALL ST_READ  ( namlevel, nchrlevel, numlevel, iret )
cIF ( iret .ne. 0 ) THEN
c        call errmsg('End of file encountered reading verifying levels')
c        call errexit(1)
cEND IF

        READ(5,*) numlevel,namlevel(1)
        nchrlevel(1) = LEN(TRIM(namlevel(1)))
C
C*      NUMBER OF VERIFYING LEVELS IS LIMITED TO MAXLVL
        IF (numlevel.gt.maxlvl) THEN
          write(a,'(i2)') maxlvl
          write(b,'(i2)') numlevel
          msg='NUMBER OF VERIFYING LEVELS EXCEEDS LIMIT OF ' // a // b
          call errmsg(msg)
          call errexit(1)
        END IF

        IF (numlevel.gt.1) THEN
           DO i = 2, numlevel
              READ(5,*) namlevel(i)
              nchrlevel(i) = LEN(TRIM(namlevel(i)))
           ENDDO
        ENDIF
C     
C       READ NUMLEVEL LEVEL PNEMONICS AND GET CHARACTER COUNT
C     
        DO 220 n = 1, numlevel
          CALL setlevel(n,namlevel(n),nchrlevel(n))
  220   CONTINUE

        write (*,*) "read level done!"

        write (*,*) "readcntl  done!"

C*
      RETURN
      END
