      module surfce

      save

      real,allocatable :: ps(:,:),zs(:,:),ts(:,:),qs(:,:),us(:,:)
     *   ,vs(:,:),pm(:,:)
      real,allocatable :: cape(:,:),cin(:,:),bcape(:,:),pli(:,:)
      real,allocatable :: vis(:,:),tmax(:,:),tmin(:,:),dpt(:,:),
     *   tocc(:,:),gust(:,:)
      real,allocatable :: pbls(:,:),pwo(:,:),trop(:,:),cdtz(:,:),
     *   pblris(:,:),pblmxs(:,:),haines(:,:),trans(:,:),wnd80(:,:)
      real,allocatable :: ceil(:,:)

      integer,allocatable :: mskps(:,:),mskzs(:,:),mskts(:,:),
     *    mskqs(:,:),mskus(:,:),mskvs(:,:),mskpm(:,:)
      integer,allocatable :: mskcp(:,:),mskcn(:,:),mskbc(:,:),mskli(:,:)
      integer,allocatable :: mskvis(:,:),msktmax(:,:),msktmin(:,:),
     *    mskdpt(:,:),msktocc(:,:),mskgust(:,:)
      integer,allocatable :: mskpbl(:,:),mskpw(:,:),msktrp(:,:),
     *   mskcdtz(:,:),mskpri(:,:),mskpmx(:,:),mskhaines(:,:)
      integer,allocatable :: msktrans(:,:),mskwnd80(:,:),mskceil(:,:)

      end module surfce
