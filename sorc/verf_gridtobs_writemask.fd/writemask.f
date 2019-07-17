      program writemask

c GRIB2
      use grib_mod
C GRIB2

      character*12 fname
      character*13 fnamei
      character*11 maskfile

      parameter (ilim=1000,jlim=1000)

      logical*1 mask(ilim*jlim)
      dimension grid(ilim*jlim)
      real imask(ilim*jlim)
      real imask2(ilim*jlim)
      integer jpds(200),jgds(200),kpds(200),kgds(200)
c GRIB2
      type(gribfield) :: gfld
      integer ifile,j,jdisc,jpdtn,jgdtn,iret
      integer,dimension(200) :: jids,jpdt,jgdt
      real fldmin,fldmax,firstval,lastval
      logical :: unpack=.true.
C
    
      read(5,'(a12)') fname
      read(5,'(a13)') fnamei
      read(5,'(a11)') maskfile

      print*,'fname=',fname
      print*,'fnamei=',fnamei
      print*,'maskfile=',maskfile

      call baopen(30,fname,iret)
      call baopen(31,fnamei,iret)
      call baopen(32,maskfile,iret)

C GRIB2
c Set GRIB2 field identification values to search for
      j=0              ! search from 0
      jdisc=0          ! for met field:0 hydro: 1, land: 2
c-- set id section
      jids=-9999
c-- set product def template, using template 4.0
      jpdtn=0
c-- set product def array
      jpdt=-9999
C 
    
      j=0
      j=0
      jpdt(1)=0      ! table 4.1
      jpdt(2)=0      ! table 4.2-0-0
      jpdt(10)=103   ! table 4.5
      jpdt(12)=2
c-- set grid def template
      jgdtn=-1
c-- set product def array
      jgdt=-9999
c     call getgb(30,31,jf,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,iret)
      call getgb2(30,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *         unpack,j,gfld,iret)
      if(iret.ne.0) then
        print*,'temp iret=',iret
c       print*,'kpds=',kpds
c       print*,'kgds=',kgds
        stop
      endif

      maxpts=gfld%ndpts
      do i=1,maxpts
c      if(mask(i)) then
       if(gfld%ibmap.eq.255.or.gfld%bmap(i)) THEN
        gfld%fld(i)=1
       else
        gfld%fld(i)=0
       endif
      print*,'i,gfld%fld(i),gfld%bmap(i)=',i,gfld%fld(i),gfld%bmap(i)
c     if(imask(i).eq.0.and.grid(i).gt.200) print*,'grid(i),mask(i)=',
c    *   grid(i),mask(i)
      enddo

c     call putgb(32,kf,kpds,kgds,mask,imask,iret)
      call putgb2(32,gfld,iret)
      
c     print*,'putgb iret=',iret

      j=0
      j = 0
      jids=-9999
c     jpdtn=0
      jpdt=-9999
      jpdt(1)=2
      jpdt(2)=2
      jpdt(10)=103
      jpdt(12)=10
c     call getgb(30,31,jf,j,jpds,jgds,kf,k,kpds,kgds,mask,grid,iret)
      call getgb2(30,0,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     *         unpack,j,gfld,iret)
c     print*,'kpds(14)=',kpds(14)
      if(iret.ne.0) then
        print*,'wind iret=',iret
c       print*,'kpds=',kpds
c       print*,'kgds=',kgds
        stop
      endif

      maxpts=gfld%ndpts
      do i=1,maxpts
c      if(mask(i)) then
       if(gfld%ibmap.eq.255.or.gfld%bmap(i)) THEN
        gfld%fld(i)=1
       else
        gfld%fld(i)=0
       endif
c     print*,'i,mask(i),imask(i)=',i,mask(i),imask(i)
      print*,'i,gfld%fld(i),gfld%bmap(i)=',i,gfld%fld(i),gfld%bmap(i)
c     if(imask(i).eq.0.and.grid(i).gt.200) print*,'grid(i),mask(i)=',
c    *   grid(i),mask(i)
      enddo

c     call putgb(32,kf,kpds,kgds,mask,imask,iret)
      call putgb2(32,gfld,iret)
c     print*,'putgb iret=',iret

      stop
      end
