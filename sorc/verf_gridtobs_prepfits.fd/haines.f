       subroutine calhaines(haines,subset)

       use observ

       real psfc, t950,t850,t700,t500
       real q850,q700,rh850,rh700,haines

       character*8 subset

       DATA bmiss /10E10/

c      print*,'obs(2,1)=',obs(7,1)

       print*,'in haines, subset=',subset

       psfc=bmiss
       tsfc=bmiss
       t850=bmiss
       q850=bmiss
       t950=bmiss
       t700=bmiss
       q700=bmiss
       t500=bmiss

c      if(obs(1,1).ne.bmiss) psfc=obs(1,1)
c      if(obs(3,1).ne.bmiss) tsfc=obs(3,1)
c      print*,'psfc=',psfc
c
c  Find the surface pressure
c

      do i=1,nlev

c     print*,'obs(4,i),hdr(5)=',obs(4,i),hdr(5)
       if(obs(4,i).eq.hdr(5)) then
        psfc=obs(1,i)
        tsfc=obs(3,i)
        print*,'psfc=',psfc
        goto 444
       endif

      enddo
444   continue
      if(psfc.eq.bmiss) then
       psfc=obs(1,1)
       print*,'psfc=',psfc
      endif
c
c If no height at all,use the hydrostatic (as in the obspbl routine)
c to calculate it
c
c      if(psfc.eq.bmiss)
c       obs(4,1)=hdr(5)
c       do i=2,nlev
c        if(obs(3,i).ne.bmiss.and.obs(3,i-1).ne.bmiss.and.
c    *      obs(1,i).ne.bmiss.and.obs(1,i-1).ne.bmiss) then
c         ttop=(obs(3,i)+273.15)*(1.+0.608*obs(2,i)/1000000.)
c         tbot=(obs(3,i-1)+273.15)*(1+0.608*obs(3,i-1)/1000000.)
c         ptop=obs(1,i)*100
c         pbot=obs(1,i-1)*100
c         obs(4,i)=obs(4,i-1)-ROG*(ttop+tbot)*(ptop-pbot)/(ptop+pbot)
c         if(obs(4,i).eq.hdr(5) tsfc=obs(3,i)
c       enddo
       
c
c  Many times the 950-mb temperature is unavailable..interpolate to
c  calculate the 950-mb temp if the surface pressure is > 950
c
      do i=1,nlev

       if(obs(1,i).lt.950.0) then
        coef=(log(950.0/obs(1,i)))/(log(obs(1,i+1)/obs(1,i)))
        t950=obs(3,i)+coef*(obs(3,i+1)-obs(3,i))
        t950=t950+273.15
        if(t950.gt.500.) t950=bmiss
        goto 333
       endif
      enddo

333   continue

       do i=1,nlev
c       print*,'obs(1,i),obs(3,i),obs(4,i),elv=',
c    *    obs(1,i),obs(3,i),obs(4,i),hdr(5)
        if(obs(1,i).eq.850.0) then
         if(obs(3,i).ne.bmiss) t850=obs(3,i)+273.15
         if(obs(2,i).ne.bmiss) q850=obs(2,i)/1.e6
        endif
c       if(obs(1,i).eq.950.0) then
c        if(obs(3,i).ne.bmiss) t950=obs(3,i)+273.15
c       endif
        if(obs(1,i).eq.700.0) then
         if(obs(3,i).ne.bmiss) t700=obs(3,i)+273.15
         if(obs(2,i).ne.bmiss) q700=obs(2,i)/1.e6
        endif
        if(obs(1,i).eq.500.0)  then
         if(obs(3,i).ne.bmiss) t500=obs(3,i)+273.15
        endif
       enddo

       print*,'t500, t700, t850, t950=',
     *  t500, t700, t850, t950 
       print*,'q850, q700=',q850,q700

       if(q850.ne.bmiss) then
        ratio=q850/(1.0-q850)
        vap=850.0*ratio/(0.622+ratio)
        vaps=w3fa09(t850)*10.
        rh850=(vap/vaps)*100.
       else
        rh850=bmiss
       endif

       if(q700.ne.bmiss) then
        ratio=q700/(1.0-q700)
        vap=700.0*ratio/(0.622+ratio)
        vaps=w3fa09(t700)*10.
        rh700=(vap/vaps)*100.
       else
        rh700=bmiss
       endif

       print*,'rh850, rh700=',rh850,rh700
    
c      if(psfc.ne.bmiss.and.t950.ne.bmiss.and.t850.ne.bmiss.
c    *     and.t700.ne.bmiss.and.rh850.ne.bmiss.and.rh700.ne.
c    *     bmiss) then

c       print*,'psfc, psfc.gt.950.=',psfc, psfc.gt.950.
        if(psfc.gt.950.0.and.t950.ne.bmiss.and.t850
     *     .ne.bmiss.and.rh850.ne.bmiss) then
         print*,'psfc.gt.950'
         hat=t950-t850
         tmois=t850-273.15
         rhmois=rh850
         st1=8
         st2=3
         mt1=10
         mt2=5
        elseif(psfc.gt.850.0.and.t850.ne.bmiss.and.t700
     *    .ne.bmiss.and.rh850.ne.bmiss) then
         print*,'psfc.gt.850'
         hat=t850-t700
         tmois=t850-273.15
         rhmois=rh850
         st1=11
         st2=5
         mt1=13
         mt2=5
        elseif(t700.ne.bmiss.and.t500.ne.bmiss.and
     *   .rh700.ne.bmiss) then
         print*,'else'
         hat=t700-t500
         tmois=t700-273.15
         rhmois=rh700
         st1=22
         st2=17
         mt1=21
         mt2=14
        else
         haines=bmiss
         goto 123
        endif

        term=log10(rhmois)/7.5+(tmois/(tmois+237.3))
        dpmois=(term*237.3)/(1.0-term)

        hainesm=tmois-dpmois
        slopet=1/(st1-st2)
        intt=1.5-(st2-0.5)*slopet
        hainest=(slopet*hat)+intt
        slopem=1/(mt1-mt2)
        intm=1.5-(mt2-0.5)*slopem
        hainesm=(slopem*dpmois)+intm
        haines=hainest+hainesm
c      else
c       haines=bmiss 
c      endif  

123    print*,'haines=',haines
       
       return
       end
