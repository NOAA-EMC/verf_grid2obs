      subroutine obspbl(hpbl,transwnd,vent,w80,idate)

      use observ

      INCLUDE 'parm.inc'

      REAL Q(100),P(100),T(100),PINT(101),Z(100),U(100),V(100)
      REAL q1(100),p1(100),t1(100),z1(100),u1(100),v1(100),d1(100)
      REAL q2(100),p2(100),t2(100),z2(100),u2(100),v2(100)
      REAL a1(100),a2(100)

      real*8 hstid
      CHARACTER*8 staid
      EQUIVALENCE (hstid,staid)

      DATA bmiss /10E10/

      parameter (G=9.8E0,CP=1004.6E0,CAPA=0.28589641E0,ROG=287.04/G)

      hstid = hdr(1)

      print*,'--------- obspbl ----------'

      do k=1,100
      p1(k)=obs(1,k)
      t1(k)=obs(3,k)
      q1(k)=obs(2,k)
      if(obs(2,k).eq.bmiss) q1(k)=bmiss
      z1(k)=obs(4,k)
      u1(k)=obs(5,k)
      v1(k)=obs(6,k)
c     print*,'k,u1(k),v1(k)=',k,u1(k),v1(k),z1(k)
      enddo
                                                                                
c--- Order all observations by pressure from bottom to top
      p2=bmiss
      t2=bmiss
      q2=bmiss
      z2=bmiss
      u2=bmiss
      v2=bmiss
      allmiss=1.
        nump=0
      do k=1,100
                                                                                
         do ii=1,100
          if(p1(ii).ne.bmiss) then
            allmiss=0.
            goto 500
          else
            allmiss=1.
          endif
         enddo
500       continue
                                                                                
      if(allmiss.ne.1)then
        imx=maxloc(p1, dim=1, mask = p1.ne.bmiss)
        p2(k)=p1(imx)
        p1(imx)=bmiss
        t2(k)=t1(imx)
        q2(k)=q1(imx)
        z2(k)=z1(imx)
        u2(k)=u1(imx)
        v2(k)=v1(imx)
        nump=k
      endif
      enddo

c     print*,'nump= ',nump
c     nump1=nump+1
c     do k=1,nump1
c     print*,'t2,k=',k,p2(k),t2(k),q2(k),z2(k),u2(k),v2(k)
c     enddo
                                                                                
c--- Make interpolation (temp to wind levels; wind to temp levels)
      do ivar=1,5
                                                                                
        if(ivar.eq.1) a2=t2
        if(ivar.eq.2) a2=q2
        if(ivar.eq.3) a2=u2
        if(ivar.eq.4) a2=v2
        if(ivar.eq.5) a2=z2
                                                                                
      a1=bmiss
      do k=1,100
        if(k.gt.nump) goto 501
         if(a2(k).eq.bmiss) then
          km=k
          kp=k
188       km=km-1
           if(km.eq.0) goto 501
           if(a2(km).eq.bmiss) goto 188
189       kp=kp+1
           if(kp.gt.nump) goto 501
           if(a2(kp).eq.bmiss) goto 189
            coef=(alog(p2(k)/p2(km)))/(alog(p2(kp)/p2(km)))
            a1(k)=a2(km)+coef*(a2(kp)-a2(km))
         else
           a1(k)=a2(k)
         endif
501   continue
      enddo
                                                                                
        if(ivar.eq.1) t1=a1
        if(ivar.eq.2) q1=a1
        if(ivar.eq.3) u1=a1
        if(ivar.eq.4) v1=a1
        if(ivar.eq.5) z1=a1
                                                                                
      enddo
                                                                                
      p1=p2

c-- If height is missing derive from hydrostatic

      if(z1(1).eq.bmiss) z1(1)=hdr(5)
c-    if(z1(1).eq.bmiss) z1(1)=0.

      do k=2,nump
      if(z1(k).eq.bmiss) then
       if(t1(k).ne.bmiss.and.t1(k-1).ne.bmiss.and.
     &    p1(k).ne.bmiss.and.p1(k-1).ne.bmiss) then
c          if(q1(k).eq.bmiss) q1(k)=0.
c          if(q1(k-1).eq.bmiss) q1(k-1)=0.
           ttop=(t1(k)+273.15)*(1.+0.608*q1(k)/1000000.)
           tbot=(t1(k-1)+273.15)*(1.+0.608*q1(k-1)/1000000.)
           ptop=p1(k)*100.
           pbot=p1(k-1)*100.
           z1(k)=z1(k-1)-ROG*(ttop+tbot)*(ptop-pbot)/(ptop+pbot)
       endif
      endif
      enddo

c     do k=1,nump1
c     print*,'t1,k=',k,p1(k),t1(k),q1(k),z1(k),u1(k),v1(k)
c     enddo
                                                                                
      e0=611
      xll=2500000
      rv=461
                                                                                
      iilev=0
      DO K=1,100
      if(p1(k).ne.BMISS.and.
     *    q1(k).ne.BMISS.and.
     *    t1(k).ne.BMISS.and.
     *    z1(k).ne.BMISS.and.
     *    u1(k).ne.BMISS.and.
     *    v1(k).ne.BMISS) then
        iilev=iilev+1
        P1(iilev)=p1(k)*100.
        Q1(iilev)=q1(k)/1000000.
        T1(iilev)=t1(k)+273.15
        Z1(iilev)=z1(k)
        U1(iilev)=u1(k)
        V1(iilev)=v1(k)
        e=q1(iilev)*p1(iilev)/0.622
        d1(iilev)=1/(1/273.16-(rv/xll)*alog(e/611.))
c       print*,'k,p1,t1=',k,p1(k),z1(k),t1(k)
c        if(iilev.gt.1.and.p1(iilev).gt.p1(iilev-1))print*,
c    +   '!!!PBL:lev=',iilev-1,'p=',p1(iilev-1),
c    +   'lev=',iilev,'p=',p1(iilev)
      endif
      ENDDO
      print*,'iilev=',iilev
                                                                                
        hpbl=bmiss
c       if(iilev.gt.2) then
        if(iilev.ge.10.and.p1(10).gt.50000.and.p1(iilev).le.50000.) then
c       print*,'before calpbl (go)',iilev,' p=',(p1(ii),ii=1,iilev)
        CALL CALPBL(T1,Q1,P1,Z1,U1,V1,IILEV,HPBL,JPBL)
c---- Substract station's elevation from HPBL
        hpbl=hpbl-hdr(5)
c--QC
        if(hpbl.gt.5000.and.hpbl.ne.bmiss) hpbl=5000.
c--
c       obs(8,1)=hpbl
c       else
c       print*,'before calpbl(not go)',iilev,' p=',(p1(ii),ii=1,iilev)
        endif

c
c  Now that we have our PBL height, we can now find the transport wind
c  and ventilation rate
c
c  For the transport winds, just average up the winds in the boundary
c  layer.
c

        sumwndu=0.0
        sumwndv=0.0
        icountu=0
        icountv=0

        do i=1,iilev
         
         if(z1(i).le.hpbl.and.u1(i).ne.bmiss.
     *     and.v1(i).ne.bmiss) then
          print*,'i,z1(i),u1(i),v1(i),hpbl=',
     *     i,z1(i),u1(i),v1(i),hpbl
          sumwndu=u1(i)+sumwndu
          sumwndv=v1(i)+sumwndv
          print*,'i,sumwndu=',i,sumwndu
          print*,'i,sumwndv=',i,sumwndv
          icountu=icountu+1
          icountv=icountv+1
         else
          if(icountu.gt.0) then
            utrans=sumwndu/icountu
          else
            utrans=0
          endif
          if(icountv.gt.0) then
            vtrans=sumwndv/icountv
          else
            vtrans=0
          endif
          transwnd=sqrt(utrans**2+vtrans**2)
          goto 346
          endif

         enddo
346      continue
         print*,'sumwndu,utrans,transwnd=',sumwndu,utrans,transwnd
         print*,'sumwndv,vtrans,transwnd=',sumwndv,vtrans,transwnd
c
c Ventilation rate:  the transport wind times the pbl height
c
         vent=transwnd*hpbl
c
c Calculate 80-m agl wind
c
         do i=1,iilev
c         print*,'i,z1(i),z1(i)-hdr(5),u1(i),v1(i)=',
c    *     i,z1(i),z1(i)-hdr(5),u1(i),v1(i)
          if((z1(i)-hdr(5)).gt.80.0) then
           hgt1=z1(i)-hdr(5)
           hgt2=z1(i-1)-hdr(5)
           uwnd1=u1(i)
           uwnd2=u1(i-1)
           vwnd1=v1(i)
           vwnd2=v1(i-1)

           coef=(80.0-hgt2)/(hgt1-hgt2)

           u80=uwnd2+coef*(uwnd1-uwnd2)
           v80=vwnd2+coef*(vwnd1-vwnd2)

c          print*,'hgt1,hgt2=',hgt1,hgt2
c          print*,'uwnd1,uwnd2,u80=',uwnd1,uwnd2,u80
c          print*,'vwnd1,vwnd2,u80=',vwnd1,vwnd2,v80

           goto 555
          endif

         enddo

555     continue
      
        w80=sqrt(u80*u80+v80*v80) 

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

        ihr = mod(int(vdata1),100)
        obshr = real(ihr)+hfrac

       icent = mod(idate/100000000,100)
       iy = mod(int(vdata1)/1000000,100)
       im = mod(int(vdata1)/10000,100)
       id = mod(int(vdata1)/100,100)
       ihr = mod(int(vdata1),100)

        if(hpbl.gt.0.and.hpbl.ne.bmiss)then
        if(hdr(2).lt.0.) then
        write(100,1556) hdr(3), hdr(2)+360., hpbl, staid, hdr(5)
        write(103,1557) hdr(3), hdr(2)+360., hpbl, staid, hdr(5),
     +       hdr(4),idate,int(vdata1),int(minutes),obshr,
     +    icent,iy,im,id,ihr,int(minutes)
        else
        write(100,1556) hdr(3), hdr(2), hpbl, staid, hdr(5)
c       write(103,*) int(vdata1)
        write(103,1557) hdr(3), hdr(2), hpbl, staid, hdr(5),
     +    hdr(4),idate,int(vdata1),int(minutes),obshr,
     +    icent,iy,im,id,ihr,int(minutes)
        endif
        write(101)(p1(kk)/100,kk=1,100),
     +            (z1(kk),kk=1,100),
     +            (t1(kk),kk=1,100),
     +            (d1(kk),kk=1,100),
     +            (q1(kk),kk=1,100),
     +            (u1(kk),kk=1,100),
     +            (v1(kk),kk=1,100)

1556    format(f6.2,1x,f7.2,2x,f5.0,3x,a8,2x,f5.0)
1557    format(f6.2,1x,f7.2,2x,f5.0,3x,a8,2x,f5.0,2x,f6.2,4x,2i12,i2.2
     +         ,2x,f5.2,2x,2i2.2,"/",i2.2,"/",i2.2,1x,i2.2,":",i2.2)
        endif

        return
        end
