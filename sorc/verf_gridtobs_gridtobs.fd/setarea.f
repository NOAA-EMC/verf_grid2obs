      SUBROUTINE setarea(iar,namarea,nchr)
C************************************************************************
      use grdefs
      INCLUDE 'parm.inc'
      CHARACTER*24 namarea, nam24
      CHARACTER*1 nam1(24)
      EQUIVALENCE (nam24,nam1(1))
c     CHARACTER*3 regnam(mxarea), regions(30)
      CHARACTER*3 regnam(mxarea)
c     COMMON /grdef/ mode(mxarea), imax(mxarea), imin(mxarea), 
c    +            jmax(mxarea), jmin(mxarea), alat1(mxarea), 
c    +            elon1(mxarea), dxx(mxarea), dyy(mxarea), 
c    +            elonv(mxarea), alatan(mxarea), latlong(mxarea), 
c    +            lambert(mxarea), polarstereo(mxarea), numreg(mxarea),
c    +            ig104(147,110), regions
c     LOGICAL latlong, lambert, polarstereo
      CHARACTER*8 dumstn
      CHARACTER*8 stnlist
      character*3 agrid
      COMMON /stndef/ nstns (mxarea), stnlist (mxarea,maxj)
C-------------------------------------------------------------------------
      nam24 = namarea
      mode(iar) = 0
      iqstn = INDEX ( nam24, '.STNS' )
      IF (nam1(1).eq.'G' .and. iqstn .eq. 0 ) THEN
        READ (nam24(2:4),'(a3)') agrid
        if(agrid.eq.'AFR') then
         print*,'agrid=',agrid
         mode(iar)=2
         return
        endif
        read(agrid,'(i3)') igrid
        print*,'igrid=',igrid
        CALL gtgdef(igrid,istat,imax(iar),jmax(iar),alat1(iar),
     +              elon1(iar),dxx(iar),dyy(iar),elonv(iar),alatan(iar),
     +              latlong(iar),lambert(iar),polarstereo(iar))
        imin(iar) = 1
        jmin(iar) = 1
        IF (nam24(5:5).eq.'/') THEN
          regnam(iar) = nam24(6:8)
          IF (ig104(75,55).eq.0) THEN
            READ (20,'(20I4)') ig104
C           PRINT*,'IG104: ',IG104
            READ (21,'(4X,A3)') regions
C         PRINT*,'REGIONS: ',REGIONS
          END IF
          DO nr = 1, 30
            IF (regnam(iar).eq.regions(nr)) THEN
              numreg(iar) = nr
              mode(iar) = 2
              PRINT *, regnam(iar), ' is subset # ', numreg(iar)
              RETURN
            END IF
          END DO
          mode(iar) = 1
        ELSE
          mode(iar) = 1
        END IF
      ELSE IF ( iqstn .ne. 0 ) THEN
	nstns (iar) = 0
	OPEN ( UNIT=22, FILE=nam24, STATUS='OLD',
     +         IOSTAT=ios )
	IF ( ios .ne. 0 ) WRITE (6,*)
     +	  'Cannot open file ', nam24
	DO WHILE ( ios .eq. 0 )
	  READ ( 22, '(A)', IOSTAT=ios ) dumstn
	  IF ( ios .eq. 0 ) THEN
	    nstns (iar) = nstns (iar) + 1
	    stnlist ( iar, nstns(iar) ) = dumstn
          END IF
	END DO
	mode (iar) = 3
	namarea = nam24 (1:iqstn-1)
	nchr = iqstn - 1
	CLOSE ( UNIT = 22 )
      ELSE
        mode(iar) = 0
        PRINT *, ' NON-CONFORMING AREA PARAMETER NEEDS TO BE Gxxx'
	PRINT *, ' OR xxx.STNS, where xxx is an upper case name.'
      END IF
      RETURN
      END
