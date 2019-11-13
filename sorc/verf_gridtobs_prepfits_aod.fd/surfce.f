      module surfce

      save

      real(8),allocatable :: ps(:,:),zs(:,:),ts(:,:),qs(:,:),us(:,:)
     *   ,vs(:,:),pm(:,:)
      real,allocatable :: cape(:,:),cin(:,:),bcape(:,:),pli(:,:)
      real,allocatable :: vis(:,:),tmax(:,:),tmin(:,:),dpt(:,:),
     *   tocc(:,:),gust(:,:)
      real,allocatable :: pbls(:,:),pwo(:,:),trop(:,:),cdtz(:,:),
     *   pblris(:,:),pblmxs(:,:)

      real,allocatable :: aod(:,:)

      integer,allocatable :: mskps(:,:),mskzs(:,:),mskts(:,:),
     *    mskqs(:,:),mskus(:,:),mskvs(:,:),mskpm(:,:)
      integer,allocatable :: mskcp(:,:),mskcn(:,:),mskbc(:,:),mskli(:,:)
      integer,allocatable :: mskvis(:,:),msktmax(:,:),msktmin(:,:),
     *    mskdpt(:,:),msktocc(:,:),mskgust(:,:)
      integer,allocatable :: mskpbl(:,:),mskpw(:,:),msktrp(:,:),
     *   mskcdtz(:,:),mskpri(:,:),mskpmx(:,:)

      integer,allocatable :: mskaod(:,:)

      end module surfce
