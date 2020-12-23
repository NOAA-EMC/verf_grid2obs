       subroutine calshear(iobmod,shr6)

c
c Purpose: Given a profile of height, and u and v wind, calculate the
c 0-6 km wind shear
c
c Author:  Perry Shafran
c Date: 17 May 2013
c
       use observ
       use backgv

       real z1(100),u1(100),v1(100),p1(100),zmid(100+1)

       real z2(100),u2(100),v2(100),p2(100)

       real a1(100),a2(100)
        
       integer count1, count5

       DATA bmiss /10E10/

       p1=bmiss
       z1=bmiss
       u1=bmiss
       v1=bmiss

       if(iobmod.eq.1) then
c         print*,'obs'
          do i=1,100
           p1(i)=obs(1,i)
           z1(i)=obs(4,i)
           u1(i)=obs(5,i)
           v1(i)=obs(6,i)
c          print*,'i,z1(i)=',i,z1(i),u1(i),v1(i)
          enddo
        elseif(iobmod.eq.2) then
c         print*,'model'
          do i=1,100
           p1(i)=bak(1,i)
           z1(i)=bak(4,i)
           u1(i)=bak(5,i)
           v1(i)=bak(6,i)
c          print*,'i,z1(i)=',i,z1(i),u1(i),v1(i)
          enddo
        endif
        
        if(iobmod.eq.1) then

        z2=bmiss
        u2=bmiss
        v2=bmiss
        p2=bmiss

        do k=1,100
         imx=maxloc(p1, dim=1, mask = p1.ne.bmiss)
         p2(k)=p1(imx)
         p1(imx)=bmiss
         z2(k)=z1(imx)
         u2(k)=u1(imx)
         v2(k)=v1(imx)
         nump=k
        enddo

        do ivar=1,5

           if(ivar.eq.1) a2=z2
           if(ivar.eq.2) a2=u2
           if(ivar.eq.3) a2=v2

           a1=bmiss

           do k=1,100
            if(k.ge.nump) goto 501
            if(a2(k).eq.bmiss) then
             km=k
             kp=k
188          km=km-1
             if(km.eq.0) goto 501
             if(a2(km).eq.bmiss) goto 188
189          kp=kp+1
             if(kp.eq.nump) goto 501
             if(a2(kp).eq.bmiss) goto 189
             coef=(alog(p2(k)/p2(km)))/(alog(p2(kp)/p2(km)))
             a1(k)=a2(km)+coef*(a2(kp)-a2(km))
            else
             a1(k)=a2(k)
            endif
501         continue
           enddo

           if(ivar.eq.1) z1=a1
           if(ivar.eq.2) u1=a1
           if(ivar.eq.3) v1=a1

        enddo
        
        endif

c       print*,'ordered'
c       do i=1,100
c         print*,'i,z1(i)=',i,z1(i),u1(i),v1(i)
c       enddo

        ust5=0.0
        vst5=0.0
        ust1=0.0
        vst1=0.0
        count1=0
        count5=0

        do i=1,100

         if(z1(i).eq.bmiss.or.u1(i).eq.bmiss.or.
     *      v1(i).eq.bmiss) goto 320
         if(z1(i).eq.0.and.u1(i).eq.0.and.v1(i).eq.0) 
     *      goto 320
      
         zmid(i)=(z1(i)+z1(i+1))*0.5
         dzabv=zmid(i)-z1(1)

         if(dzabv.lt.6000.and.dzabv.ge.5500) then
          ust5=ust5+u1(i)
          vst5=vst5+v1(i)
          count5=count5+1
c         print*,'ust5,vst5,count5=',ust5,vst5,count5
         endif

         if(dzabv.lt.500) then
          ust1=ust1+u1(i)
          vst1=vst1+v1(i)
          count1=count1+1
         endif
320     continue
        enddo

        if(count5.eq.0) then
         do i=100,1,-1
          zmid(i)=(z1(i)+z1(i+1))*0.5
          dzabv=zmid(i)-z1(1)

          if(dzabv.lt.7000.and.dzabv.ge.6000) then
           ust5=ust5+u1(i)
           vst5=vst5+v1(i)
           count5=1
c          print*,'ust5,vst5,count5=',ust5,vst5,count5
           goto 321
          endif
         enddo
        endif
321    continue

        if(count5.gt.0.and.count1.gt.0) then
          umean5=ust5/count5
          vmean5=vst5/count5
          umean1=ust1/count1
          vmean1=vst1/count1

          ushr6=umean5-umean1
          vshr6=vmean5-vmean1

          shr6=sqrt(ushr6*ushr6+vshr6*vshr6)

         else
          shr6=bmiss
         endif

c        print*,'ushr6,vshr6=',ushr6,vshr6,shr6,iobmod

         return
         end
