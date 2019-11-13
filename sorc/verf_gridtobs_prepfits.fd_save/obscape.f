      subroutine obscapepw(ivirt,cape,cin,pli,pw,bcape)

c
c   Purpose:  Set up the observed profile of temperature, pressure
c   and moisture for the processing of the convective indices CAPE, CIN,
c   and LI.  Returns the observed values of CAPE, CIN, and LI and puts then
c   in the obs array. 
c
c   Author:  Perry Shafran - 19 May 2006
c

      use observ
      use backgv
      use obstrs

      real p(200),p1(200),q(200),q1(200),t(200),t1(200),z(200),z1(200)
      real pint(201)

      character*80 obstrout,obstrout2

      real*8 pbnd(6),qbnd(6),tbnd(6),pbint(7),psum(200),nsum(200)
      integer lvlbnd(6)

      DATA bmiss /10E10/

      ILEV=0
      do k=1,200
       p(k)=bmiss
       p1(k)=bmiss
       q(k)=bmiss
       q1(k)=bmiss
       t(k)=bmiss
       t1(k)=bmiss
       z(k)=bmiss
       z1(k)=bmiss
      enddo
      slindx=bmiss
      cape=bmiss
      cin=bmiss
      peql=bmiss
      plcl=bmiss
      ifcape=bmiss
c     do iii=9,11
c       bak(iii,2)=bmiss
c       obs(iii,2)=bmiss
c       bak(iii,1)=bmiss
c       obs(iii,1)=bmiss
c     enddo
      DO K=1,200
      if(obs(2,k).ne.bmiss.and.
     *    obs(3,k).ne.bmiss) then
        ilev=ilev+1
        P1(ilev)=OBS(1,k)*100.
        Q1(ilev)=OBS(2,k)/1000000.
        T1(ilev)=OBS(3,k)+273.15
        Z1(ilev)=OBS(4,k)
      endif
      ENDDO
      if(ilev.eq.0) return
      do k=1,ilev
       p(k)=p1(ilev-k+1)
       q(k)=q1(ilev-k+1)
       t(k)=t1(ilev-k+1)
       z(k)=z1(ilev-k+1)
       print*,'k,p(k),t(k),q(k)=',k,p(k),t(k),q(k)
      enddo
      itype=1
      p1d=bmiss
      t1d=bmiss
      q1d=bmiss
      p(24)=1012
      t(24)=26.2
      q(24)=19231
      p(23)=1007
      t(23)=27.4
      q(23)=20282
      p(22)=1000
      t(22)=26.8
      q(22)=19347
      p(21)=999
      t(21)=26.6
      q(21)=19133
      p(20)=965
      t(20)=25
      q(20)=18420
      p(19)=925
      t(19)=23
      q(19)=15676
      p(18)=850
      t(18)=18.6
      q(18)=12014
      p(17)=811
      t(17)=16.2
      q(17)=10265
      p(16)=775
      t(16)=14.2
      q(16)=10253
      p(15)=706
      t(15)=11.2
      q(15)=6799
      p(14)=700
      t(14)=10.6
      q(14)=6570
      p(13)=582
      t(13)=-0.3
      q(13)=5449
      p(12)=552
      t(12)=-1.9
      q(12)=3519
      p(11)=547
      t(11)=-2.3
      q(11)=4183
      p(10)=537
      t(10)=-3.1
      q(10)=3037
      p(9)=504
      t(9)=-6.1
      q(9)=4301
      p(8)=500
      t(8)=-6.5
      q(8)=4138
      p(7)=493
      t(7)=-7.5
      q(7)=3556
      p(6)=483
      t(6)=-6.7
      q(6)=2735
      p(5)=462
      t(5)=-8.9
      q(5)=2012
      p(4)=404
      t(4)=-16.3
      q(4)=1921
      p(3)=400
      t(3)=-16.1
      q(3)=2079
      p(2)=343
      t(2)=-23.7
      q(2)=1361
      p(1)=300
      t(1)=-30.7
      q(1)=734

      do i=1,24
       p(i)=p(i)*100
       q(i)=q(i)/1000000.
       t(i)=t(i)+273.15
       print*,'i,p(i),q(i),t(i)=',i,p(i),q(i),t(i)
      enddo
      
      ivirt=1
      itype=1
      ilev=24
 
      CALL CALCAPE(ivirt,itype,T,Q,P,p1d,t1d,q1d,pint,ILEV,1,1,ILEV,
     *   CAPE,CIN,PLCL,PEQL,PLI)
      print*,'test cape=',cape
c     obs(9,1)=cape
c     obs(10,1)=cin
c     obs(11,1)=pli
      obsi2(2,2)=cape
      obsi2(3,2)=cin
      obsi2(4,2)=pli

c
c Calculate the 0-6 km wind shear AND the product of CAPE and the wind
c shear - the product is known as the Craven-Brooks Severe Weather
c Parameter - 14 May 2013

c     do k=1,200
c      print*,'obs(4,k)=',obs(4,k)
c     enddo 

      umean5=ust5/count5
      vmean5=vst5/count5

      ushr6=umean5-umean1
      vshr6=vmean5-vmean1
c
c  Now calculate the Best Cape - 11 July 2011
c
c  In order to calculate Best Cape the arrays P1D, T1D, and Q1D need to be obtained - borrowing
c  code from the Post
c
c  From BNDLYR
c

       nbnd=6
       dpbnd=30.e2
       pbnd=0.0
       tbnd=0.0
       qbnd=0.0
       lvlbnd=0
       psum=0.0
       nsum=0.0
       print*,'ilev=',ilev
       pbint(1)=p(ilev)
       do lbnd=2,nbnd+1
        pbint(lbnd)=pbint(lbnd-1)-dpbnd
       enddo

c      do lbnd=1,nbnd+1
c        print*,'lbnd,pbint(lbnd)=',lbnd,pbint(lbnd)
c      enddo

      tbnd=0.0
      qbnd=0.0

      do lbnd=1,nbnd
       do l=1,ilev-1
        pm=(p(l)+p(l+1))*0.5
        if((pbint(lbnd).ge.pm).and.(pbint(lbnd+1).le.pm)) then
c         print*,'l,p(l),p(l+1)=',l,p(l),p(l+1)
          dp=p(l+1)-p(l)
c         print*,'dp=',dp
          nsum(lbnd)=nsum(lbnd)+1
          psum(lbnd)=psum(lbnd)+dp
          lvlbnd(lbnd)=lvlbnd(lbnd)+l
c         print*,'l,lbnd,t(l),tbnd=',l,lbnd,t(l),tbnd(lbnd)
          tbnd(lbnd)=tbnd(lbnd)+t(l)*dp
c         print*,'after: l,lbnd,t(l),tbnd=',l,lbnd,t(l),tbnd(lbnd)
          qbnd(lbnd)=qbnd(lbnd)+q(l)*dp
c         print*,'after: l,lbnd,q(l),qbnd=',l,lbnd,q(l),qbnd(lbnd)
         endif
        enddo
       enddo

      do lbnd=1,nbnd
        if(psum(lbnd).ne.0) then
c         print*,'in here'
c         print*,'psum(lbnd)=',psum(lbnd)
          rpsum=1./psum(lbnd)
c         print*,'rpsum=',rpsum
          lvlbnd(lbnd)=lvlbnd(lbnd)/nsum(lbnd)
          pbnd(lbnd)=(pbint(lbnd)+pbint(lbnd+1))*0.5
          tbnd(lbnd)=tbnd(lbnd)*rpsum
          qbnd(lbnd)=qbnd(lbnd)*rpsum
         endif
       enddo

       do lbnd=1,nbnd
        if(psum(lbnd).eq.0) then
         l=200
         pmin=9999999.
         pbnd(lbnd)=(pbint(lbnd)+pbint(lbnd+1))*0.5
c
         do ll=1,ilev-1
          pm=(p(ll)+p(ll+1))*0.5
          delp=abs(pm-pbnd(lbnd))
          if(delp.lt.pmin) then
           pmin=delp
           l=ll
          endif
         enddo

         dp=p(l+1)-p(l)
         lvlbnd(lbnd)=l
         tbnd(lbnd)=t(l)
         qbnd(lbnd)=q(l)
        endif

       enddo


c      do lbnd=1,nbnd
c       print*,'lbnd,pbnd,tbnd,qbnd,lvlbnd(lbnd)=',
c    *    lbnd,pbnd(lbnd),tbnd(lbnd),qbnd(lbnd),lvlbnd(lbnd)
c      enddo

C  Having acquired the boundary level values, calculate the theta-E 

       thte=0.0
       egrid1=-99999.
       egrid2=-99999.
       do lbnd=1,nbnd

        eps=18.015/28.964
        oneps=1.-eps
        
        evp=pbnd(lbnd)*qbnd(lbnd)/(eps+oneps*qbnd(lbnd))
        rmx=eps*evp/(pbnd(lbnd)-evp)
        ckapa=0.2845*(1.-0.28*rmx)
        rkapa=1./ckapa
        arg=evp*0.01
        arg=amax1(1.e-12,arg)
        tt1=tbnd(lbnd)
        denom=3.5*alog(tt1)-alog(evp*0.01)-4.805
        tlcl=2840./denom+55.
        plcl=p1d*(tlcl/tbnd(lbnd))**rkapa
        fac=(1000.e2/pbnd(lbnd))**ckapa
        eterm=(3.376/tlcl-0.00254)*(rmx*1.e3*(1.0+0.81*rmx))
        thte=tbnd(lbnd)*fac*exp(eterm)
c       print*,'egrid2,thte=',egrid2,thte
        if(thte.gt.egrid2) then
         egrid2=thte
         lb2=lvlbnd(lbnd)
         p1d=pbnd(lbnd)
         q1d=qbnd(lbnd)
         t1d=tbnd(lbnd)
        endif
       enddo

c      print*,'p1d,t1d,q1d,lb2=',p1d,t1d,q1d,lb2

       itype=2
      CALL CALCAPE(ivirt,itype,T,Q,P,p1d,t1d,q1d,pint,ILEV,1,1,ILEV,
     *   BCAPE,CIN,PLCL,PEQL,PLI)

       ilev=0
       do k=1,200
       if(obs(2,k).ne.bmiss.and.
     *    obs(3,k).ne.bmiss) then
       ilev=ilev+1
       p1(ilev)=obs(1,k)*100
       q1(ilev)=obs(2,k)/1000000.
       endif
       enddo
       do k=1,ilev
        p(k)=p1(ilev-k+1)
        q(k)=q1(ilev-k+1)
       enddo
       htm=1
       call calpw(pw,p,q,ilev)
       

c     obsi(9,2)=pw
      return
      end
