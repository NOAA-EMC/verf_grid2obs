      module grdefs

      save
   
      real,allocatable :: alat1(:), elon1(:), dxx(:), dyy(:)
      real,allocatable :: elonv(:), alatan(:)
      logical, allocatable :: latlong(:), lambert(:), polarstereo(:)
      character(3),allocatable :: regions(:)
      integer,allocatable :: mode(:), imax(:), imin(:), jmax(:), jmin(:)
      integer,allocatable :: ig104(:,:),numreg(:)

      end module grdefs
    
