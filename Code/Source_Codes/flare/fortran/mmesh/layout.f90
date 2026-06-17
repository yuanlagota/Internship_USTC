module flare_mmesh_layout
  use iso_fortran_env
  implicit none


  integer, parameter :: &
     LOWER_SIDE = 1, &
     UPPER_SIDE = 2


  integer, parameter :: &
     PLATE    = -3, &
     VACUUM   = -2, &
     CORE     = -1, &
     PERIODIC = 1, &
     UPDOWN   = 2, &
     MAPPING  = 3


  type, private :: bound
     integer :: radial(2), poloidal(2), toroidal(2)
  end type bound


  ! number of layers per (toroidal) block, number of zones
  integer :: layers, nz
  ! map zone index to block index and layer index
  integer, allocatable :: block_index(:), layer_index(:)


  ! connectivity of flux surface contours at X-points
  integer, allocatable :: connectX(:)
  ! connectivity of zones
  type(bound), allocatable :: connectZ(:)


   contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init_layout()
  use flare_mmesh_parameters

  integer :: iblock, ilayer, iz


  select case(layout)
  ! generic layout
  case(LAYOUT_GENERIC)
     layers = 1
     allocate (connectX, source=[0]) ! irrelevant here, will be deallocated later
     allocate (connectZ(0:layers-1))
     connectZ(0) = bound([CORE, VACUUM], [PERIODIC, PERIODIC], [MAPPING, MAPPING])


  ! (lower) single null
  case(LAYOUT_LSN)
     layers = 3
     allocate (connectX, source=[1])
     allocate (connectZ(0:layers-1))
     connectZ(0) = bound([CORE, MAPPING],   [PERIODIC, PERIODIC], [MAPPING, MAPPING])
     connectZ(1) = bound([MAPPING, VACUUM], [PLATE,    PLATE],    [MAPPING, MAPPING])
     connectZ(2) = bound([VACUUM, MAPPING], [PLATE,    PLATE],    [MAPPING, MAPPING])


  ! disconnected double null
  case(LAYOUT_DDN)
     layers = 6
     allocate (connectX, source=[-2, -2])
     allocate (connectZ(0:layers-1))
     connectZ(0) =     bound([CORE, MAPPING],    [PERIODIC, PERIODIC], [MAPPING, MAPPING])
     connectZ(1) =     bound([MAPPING, MAPPING], [PLATE,    PLATE],    [MAPPING, MAPPING])
     do iz=2,3
        connectZ(iz) = bound([MAPPING, VACUUM],  [PLATE,    PLATE],    [MAPPING, MAPPING])
     enddo
     do iz=4,5
        connectZ(iz) = bound([VACUUM, MAPPING],  [PLATE,    PLATE],    [MAPPING, MAPPING])
     enddo


  ! connected double null
  case(LAYOUT_CDN)
     layers = 5
     allocate (connectX, source=[2, 1])
     allocate (connectZ(0:layers-1))
     connectZ(0) =     bound([CORE, MAPPING],    [PERIODIC, PERIODIC], [MAPPING, MAPPING])
     do iz=1,2
        connectZ(iz) = bound([MAPPING, VACUUM],  [PLATE,    PLATE],    [MAPPING, MAPPING])
     enddo
     do iz=3,4
        connectZ(iz) = bound([VACUUM, MAPPING],  [PLATE,    PLATE],    [MAPPING, MAPPING])
     enddo

  end select


  ! set up block and layer indices for zones
  nz = blocks * layers
  allocate (block_index(0:nz-1), layer_index(0:nz-1))

  iz = 0
  do iblock=0,blocks-1
     do ilayer=0,layers-1
        block_index(iz) = iblock
        layer_index(iz) = ilayer
        iz = iz + 1
     enddo
  enddo

  end subroutine init_layout
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function toroidal_mapping_type(iblock, iside) result(itype)
  use flare_mmesh_parameters
  integer, intent(in) :: iblock, iside
  integer             :: itype


  ! default mapping type
  itype = MAPPING
  if (.not. stellarator_symmetry) return


  ! 1st mapping surface
  if (iside == LOWER_SIDE  .and.  iblock == 0) then
     if (updown_symmetry == -1) then
        itype = 2
     elseif (updown_symmetry == 0) then
        itype = -MAPPING
     endif

  ! last mapping surface
  elseif (iside == 2  .and.  iblock == blocks-1) then
     if (updown_symmetry == 1) then
        itype = 2
     elseif (updown_symmetry == 0) then
        itype = -MAPPING
     endif

  endif

  end function toroidal_mapping_type
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_layout()


  if (allocated(connectX)) deallocate (connectX)
  if (allocated(connectZ)) deallocate (connectZ)

  end subroutine free_layout
  !-----------------------------------------------------------------------------

end module flare_mmesh_layout
