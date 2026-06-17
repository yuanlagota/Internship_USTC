module base_mesh
  use iso_fortran_env
  use moose_geometry, only: interp_curve
  use moose_contours
  use flare_rpath2d
  use flare_mmesh_parameters, only: toroidal_discretization
  use mesh_interface
  use elements
  use mfs_mesh
  implicit none
  private


  character(len=*), public, parameter :: &
     TOPO_LSN    = "lsn", &
     TOPO_DDN    = "ddn", &
     TOPO_CDN    = "cdn", &
     TOPO_DDNSFP = "ddnsf+"


  integer, parameter :: &
     LEFT   =  1, &
     CENTER =  0, &
     RIGHT  = -1

  type t_layer
     ! number of (poloidal) elements in layer
     integer   :: nz

     ! element indices
     integer, dimension(:), allocatable :: iz

     ! base element index (i0), poloidal layer (ipl) and side (ipl_side)
     integer   :: i0

     ! radial and poloidal resolution in layer
     integer   :: nr = UNDEFINED, np = UNDEFINED

     ! toroidal discretization
     type(toroidal_discretization) :: T


     contains
     procedure :: initialize
     procedure :: setup_resolution
     procedure :: map_poloidal_resolution
  end type t_layer
  type(t_layer), dimension(:), allocatable :: L


  ! magnetic axis
  real(real64)  :: Pmag(2)

  ! X-point(s), separatrix(ces) and radial paths
  type(xcontour), dimension(:),   allocatable :: S
  type(rpath2d_curve),      dimension(:,:), allocatable, target :: R
  integer,            dimension(:),   allocatable :: connectX
  type(interp_curve), target :: S0, S0L, S0R
  integer       :: nX, nrpath, layers

  type(t_mfs_mesh), dimension(:), allocatable :: M
  type(t_mfs_mesh), dimension(:), allocatable, target :: Mtmp

  logical :: bugfix_ddn

  public :: setup_topology
  public :: setup_geometry
  public :: setup_interfaces
  public :: generate_equi2d_base_mesh
  public :: init_equi2d_base_mesh_generator

  contains
!=======================================================================



!=======================================================================
  subroutine initialize(this, nz, iz, i0)
  class(t_layer)      :: this
  integer, intent(in) :: nz, iz(nz), i0


  if (nz <= 0) then
     write (6, *) 'error in t_layer%initialize: number of elements must be positive!'
     stop
  endif


  this%nz = nz
  allocate (this%iz(nz))
  this%iz = iz

  ! set base element index
  this%i0 = i0

  end subroutine initialize
!=======================================================================



!=======================================================================
  subroutine setup_resolution(this, il)
  use flare_mmesh_parameters, only: nr, np, npR, npL
  class(t_layer)      :: this
  integer, intent(in) :: il

  integer :: i, i1, iz, ipl, ipl0


  ! 1. set radial resolution in this layer..............................
  this%nr = nr(il)
  !.....................................................................


  ! 2. set resolution in elements.......................................
  ! case A: innermost domain
  if (il == 0) then
     if (this%nz == 1) then
        iz = this%iz(1);  Z(iz)%nr = nr(il);  Z(iz)%np = np(il)
        Z(iz)%ipl      = 0
        Z(iz)%ipl_side = CENTER
     elseif (this%nz == 2) then
        iz = this%iz(1);  Z(iz)%nr = nr(il);  Z(iz)%np = npL(il)
        Z(iz)%ipl      = 0
        Z(iz)%ipl_side = LEFT

        iz = this%iz(2);  Z(iz)%nr = nr(il);  Z(iz)%np = npR(il)
        Z(iz)%ipl      = 0
        Z(iz)%ipl_side = RIGHT
     else
        write (6, 9000);  write(6, 9001);  stop
     endif

  ! case B: outer layers -> poloidal resolution is already defined in at least one element
  else
     ! 1. set radial resolution throughout layer
     do i=1,this%nz
        iz = this%iz(i);  Z(iz)%nr = nr(il)
     enddo

     ! 2. set poloidal resolution on lower/right side of layer
     ! 2.1 find index of first element with defined poloidal resolution
     do i1=1,this%nz
        iz = this%iz(i1)
        if (Z(iz)%np /= UNDEFINED) exit
     enddo
     if (i1 > this%nz) then
        write (6, 9000);  write(6, 9002);  stop
     endif
     ipl0 = Z(iz)%ipl; if (Z(iz)%ipl_side == LEFT) ipl0 = ipl0 + 1
     ! 2.2 go backwards and set up poloidal layer and corresponding resolution
     do i=i1-1,1,-1
        iz             = this%iz(i)
        ipl            = ipl0 + i1-i
        Z(iz)%np       = npR(ipl)
        Z(iz)%ipl      = ipl
        Z(iz)%ipl_side = RIGHT
     enddo

     ! 3. set poloidal resolution on upper/left side of layer
     ! 3.1 find index of first element with defined poloidal resolution
     do i1=this%nz,1,-1
        iz = this%iz(i1)
        if (Z(iz)%np /= UNDEFINED) exit
     enddo
     if (i1 < 1) then
        write (6, 9000);  write(6, 9002);  stop
     endif
     ipl0 = Z(iz)%ipl; if (Z(iz)%ipl_side == RIGHT) ipl0 = ipl0 + 1
     ! 3.2 go backwards and set up poloidal layer and corresponding resolution
     do i=i1+1,this%nz
        iz             = this%iz(i)
        ipl            = ipl0 + i-i1
        Z(iz)%np       = npL(ipl)
        Z(iz)%ipl      = ipl
        Z(iz)%ipl_side = LEFT
     enddo
  endif
  !.....................................................................


  ! 3. set poloidal resolution in layer.................................
  this%np = 0
  do i=1,this%nz
     iz      = this%iz(i)
     this%np = this%np + Z(iz)%np
  enddo
  !.....................................................................


 9000 format('error in t_layer%setup_resolution:')
 9001 format('innermost domain with ', i0, ' > 2 elements not supported!')
 9002 format('poloidal resolution undefined in all elements!')
  end subroutine setup_resolution
!=======================================================================



!=======================================================================
  subroutine map_poloidal_resolution(this)
  class(t_layer)      :: this

  integer :: i, iside, iz, iz_map


  do i=1,this%nz
     iz = this%iz(i)

     do iside=-1,1,2
        iz_map = Z(iz)%map_r(iside)
        ! map to another element?
        if (iz_map < 0) cycle

        ! map poloidal resolution
        if (Z(iz_map)%np == UNDEFINED) then
           Z(iz_map)%np       = Z(iz)%np
           Z(iz_map)%ipl      = Z(iz)%ipl
           Z(iz_map)%ipl_side = Z(iz)%ipl_side

        ! poloidal resolution in mapped element already defined?
        elseif (Z(iz_map)%np /= Z(iz)%np) then
           write (6, 9000) iz, iz_map, Z(iz)%np, Z(iz_map)%np;  stop
        endif
     enddo
  enddo

 9000 format('error in t_layer%map_poloidal_resolution:'//, 'element ', i0, ' maps to ', i0, &
             ', but poloidal resolution is ', i0, ' vs. ', i0, '!')
  end subroutine map_poloidal_resolution
!=======================================================================



!=======================================================================
  subroutine setup_topology(nx_out, nrpath_out)
  use flare_mmesh_parameters, only: topology => layout
  integer, intent(out), optional :: nx_out, nrpath_out


  bugfix_ddn = .false.
  ! 1. initialize topology
  layers = -1
  select case(topology)
  ! lower single null (LSN)
  case(TOPO_LSN, "lower_single_null")
     topology = TOPO_LSN
     nrpath = 3
     call initialize_elements(6)
     call initialize_interfaces(nrpath) ! radial interfaces
     nX = 1;  allocate(connectX(nX))
     connectX(1) = 1

  ! disconnected double null (DDN)
  case(TOPO_DDN, "disconnected_double_null")
     topology = TOPO_DDN
     bugfix_ddn = .true.
     nrpath = 10
     call initialize_elements(16)
     call initialize_interfaces(nrpath) ! radial interfaces
     nX = 2;  allocate(connectX(nX))
     connectX(1) = -2
     connectX(2) = -2

  ! connected double null (CDN)
  case(TOPO_CDN, "connected_double_null")
     topology = TOPO_CDN
     nrpath = 10
     call initialize_elements(12)
     call initialize_interfaces(nrpath) ! radial interfaces
     nX = 2;  allocate(connectX(nX))
     connectX(1) = 2
     connectX(2) = 1

  ! snowflake + (in disconnected double null)
  case(TOPO_DDNSFP, "disconnected_double_null_snowflake_plus")
     topology = TOPO_DDNSFP
     nrpath = 14
     call initialize_elements(22)
     call initialize_interfaces(nrpath) ! radial interfaces
     nX = 3;  allocate(connectX(nX))
     connectX(1) = -3
     connectX(2) = -1
     connectX(3) = -3

  case default
     write (6, 9000) trim(topology)
     stop
  end select
  call initialize_poloidal_interfaces(4*nX)



  ! 2. setup element topology
  select case(topology)
  ! lower single null (LSN)
  case(TOPO_LSN)
     ! innermost domain
     call Z(1)%setup_boundary(LOWER, POLOIDAL, PERIODIC, 1) ! periodic poloidal boundaries at interface R1
     call Z(1)%setup_boundary(UPPER, POLOIDAL, PERIODIC, 1) ! periodic poloidal boundaries
     call Z(1)%setup_boundary(LOWER, RADIAL,   CORE)     ! core boundary
     call Z(1)%setup_mapping (UPPER, RADIAL,   Z(2), 1)  ! connect to main SOL at interface I1

     ! main SOL
     call Z(2)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(2)%setup_mapping (LOWER, POLOIDAL, Z(3), 2)     ! connect to right divertor leg at interface R2
     call Z(2)%setup_mapping (UPPER, POLOIDAL, Z(4))  ! connect to left divertor leg

     ! right divertor leg (SOL)
     call Z(3)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(3)%setup_boundary(LOWER, POLOIDAL, DIVERTOR) ! divertor target
     call Z(3)%setup_mapping (LOWER, RADIAL,   Z(5), 2)  ! connect to right PFR at interface I2

     ! left divertor leg (SOL)
     call Z(4)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(4)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! divertor target
     call Z(4)%setup_mapping (LOWER, RADIAL,   Z(6), 3)  ! connect to right PFR at interface I3

     ! right divertor leg (PFR)
     call Z(5)%setup_boundary(LOWER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(5)%setup_boundary(LOWER, POLOIDAL, DIVERTOR) ! divertor target
     call Z(5)%setup_mapping (UPPER, POLOIDAL, Z(6), 4)  ! connect to left PFR at interface R4

     ! left divertor leg (PFR)
     call Z(6)%setup_boundary(LOWER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(6)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! divertor target


  ! DDN
  case(TOPO_DDN)
     ! innermost domain
     call Z(1)%setup_boundary (LOWER, RADIAL,   CORE)     ! core boundary
     call Z(1)%setup_mapping  (UPPER, RADIAL,   Z(3), 1)  ! connect to main SOL at interface 1
     call Z(2)%setup_boundary (LOWER, RADIAL,   CORE)     ! core boundary
     call Z(2)%setup_mapping  (UPPER, RADIAL,   Z(4), 2)  ! connect to main SOL at interface 2
     call Z(1)%setup_mapping  (UPPER, POLOIDAL, Z(2))     ! connect left and right segments
     call Z(2)%setup_mapping  (UPPER, POLOIDAL, Z(1), 1)  ! connect left and right segments at interface R1

     ! primary SOL
     call Z(3)%setup_mapping  (UPPER, RADIAL,   Z(5), 8)  ! connect to right secondary SOL at PARTIAL interface I6 (this should only be used to find the poloidal side for the generating radial path!)
     call Z(4)%setup_mapping  (UPPER, RADIAL,   Z(6), 7)  ! connect to left secondary SOL at interface I5 (SAME NOTE AS ABOVE)
     call Z(3)%setup_mapping  (UPPER, POLOIDAL, Z(4))     ! connect left and right segments
     call Z(3)%setup_mapping  (LOWER, POLOIDAL, Z(7), 2)  ! connect to right divertor leg at interface R2
     call Z(7)%setup_boundary (LOWER, POLOIDAL, DIVERTOR) ! right divertor target
     call Z(4)%setup_mapping  (UPPER, POLOIDAL, Z(8))     ! connect to left divertor leg
     call Z(8)%setup_boundary (UPPER, POLOIDAL, DIVERTOR) ! left divertor target
     call Z(7)%setup_mapping  (UPPER, RADIAL,   Z(9), 5)  ! VIRTUAL interface I5
     call Z(8)%setup_mapping  (UPPER, RADIAL,   Z(12),6)  ! VIRTUAL interface I6

     ! secondary SOL
     ! right branch
     call Z(5)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(5)%setup_mapping  (LOWER, POLOIDAL, Z(9))     ! connect to right divertor leg (right branch)
     call Z(9)%setup_boundary (LOWER, POLOIDAL, DIVERTOR) ! right divertor target
     call Z(5)%setup_mapping  (UPPER, POLOIDAL, Z(10), 6)    ! connect to left divertor leg  (right branch) at interface R6
     call Z(10)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! left divertor target
     call Z(9)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(10)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     ! left branch
     call Z(6)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(6)%setup_mapping  (LOWER, POLOIDAL, Z(11), 7)    ! connect to right divertor leg (left branch) at interface R7
     call Z(11)%setup_boundary(LOWER, POLOIDAL, DIVERTOR) ! right divertor target
     call Z(6)%setup_mapping  (UPPER, POLOIDAL, Z(12))    ! connect to left divertor leg  (left branch)
     call Z(12)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! left divertor target
     call Z(11)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(12)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain

     ! primary PFR
     call Z(7)%setup_mapping  (LOWER, RADIAL,   Z(13), 3)  ! connect to right primary PFR at interface I3
     call Z(13)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(8)%setup_mapping  (LOWER, RADIAL,   Z(14), 4)  ! connect to left primary PFR at interface I4
     call Z(14)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(13)%setup_boundary(LOWER, POLOIDAL, DIVERTOR)  ! right divertor target
     call Z(13)%setup_mapping (UPPER, POLOIDAL, Z(14), 4)  ! connect to left primary PFR at interface R4
     call Z(14)%setup_boundary(UPPER, POLOIDAL, DIVERTOR)  ! left divertor target

     ! secondary PFR
     call Z(10)%setup_mapping (LOWER, RADIAL,   Z(15),10)  ! connect to right secondary PFR at interface I10
     call Z(15)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(11)%setup_mapping (LOWER, RADIAL,   Z(16), 9)  ! connect to left secondary PFR at interface I9
     call Z(16)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(16)%setup_mapping (UPPER, POLOIDAL, Z(15), 8)  ! connect to left secondary PFR at interface R8
     call Z(16)%setup_boundary(LOWER, POLOIDAL, DIVERTOR)  ! right divertor target
     call Z(15)%setup_boundary(UPPER, POLOIDAL, DIVERTOR)  ! left divertor target


  ! connected double null (CDN)
  case(TOPO_CDN)
     ! innermost domain
     call Z(1)%setup_boundary (LOWER, RADIAL,   CORE)     ! core boundary
     call Z(1)%setup_mapping  (UPPER, RADIAL,   Z(3), 1)  ! connect to main SOL at interface 1
     call Z(2)%setup_boundary (LOWER, RADIAL,   CORE)     ! core boundary
     call Z(2)%setup_mapping  (UPPER, RADIAL,   Z(6), 2)  ! connect to main SOL at interface 2
     call Z(1)%setup_mapping  (UPPER, POLOIDAL, Z(2))     ! connect left and right segments
     call Z(2)%setup_mapping  (UPPER, POLOIDAL, Z(1), 1)  ! connect left and right segments at interface R1

     ! rhs SOL
     call Z(3)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(3)%setup_mapping  (LOWER, POLOIDAL, Z(4), 2)  ! connect to right divertor leg (right branch)
     call Z(4)%setup_boundary (LOWER, POLOIDAL, DIVERTOR) ! right lower divertor target
     call Z(3)%setup_mapping  (UPPER, POLOIDAL, Z(5))     ! connect to upper divertor leg  (right branch) at interface R6
     call Z(5)%setup_boundary (UPPER, POLOIDAL, DIVERTOR) ! upper divertor target
     call Z(4)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(5)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain

     ! lhs SOL
     call Z(6)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(6)%setup_mapping  (LOWER, POLOIDAL, Z(7))     ! connect to right divertor leg (left branch) at interface R7
     call Z(7)%setup_boundary (LOWER, POLOIDAL, DIVERTOR) ! right divertor target
     call Z(6)%setup_mapping  (UPPER, POLOIDAL, Z(8), 3)  ! connect to left divertor leg  (left branch)
     call Z(8)%setup_boundary (UPPER, POLOIDAL, DIVERTOR) ! left divertor target
     call Z(7)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(8)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain

     ! lower PFR
     call Z(4)%setup_mapping  (LOWER, RADIAL,   Z(9), 3)  ! connect to right lower PFR at interface I3
     call Z(9)%setup_boundary (LOWER, RADIAL,   VACUUM)   !
     call Z(8)%setup_mapping  (LOWER, RADIAL,   Z(10), 4) ! connect to left lower PFR at interface I4
     call Z(10)%setup_boundary(LOWER, RADIAL,   VACUUM)   !
     call Z(9)%setup_boundary (LOWER, POLOIDAL, DIVERTOR) ! right divertor target
     call Z(9)%setup_mapping  (UPPER, POLOIDAL, Z(10), 4) ! connect to left lower PFR at interface R4
     call Z(10)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! left divertor target

     ! upper PFR
     call Z(5)%setup_mapping  (LOWER, RADIAL,   Z(12), 6) ! connect to right upper PFR at interface I10
     call Z(12)%setup_boundary(LOWER, RADIAL,   VACUUM)   !
     call Z(7)%setup_mapping  (LOWER, RADIAL,   Z(11), 5) ! connect to left upper PFR at interface I9
     call Z(11)%setup_boundary(LOWER, RADIAL,   VACUUM)   !
     call Z(11)%setup_mapping (UPPER, POLOIDAL, Z(12), 8) ! connect to left upper PFR at interface R8
     call Z(11)%setup_boundary(LOWER, POLOIDAL, DIVERTOR) ! right divertor target
     call Z(12)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! left divertor target


  ! snowflake + (in disconnected double null)
  case(TOPO_DDNSFP)
     ! innermost domain
     call Z(1)%setup_boundary (LOWER, RADIAL,   CORE)     ! core boundary
     call Z(1)%setup_mapping  (UPPER, RADIAL,   Z(3), 1)  ! connect to main SOL at interface 1
     call Z(2)%setup_boundary (LOWER, RADIAL,   CORE)     ! core boundary
     call Z(2)%setup_mapping  (UPPER, RADIAL,   Z(4), 2)  ! connect to main SOL at interface 2
     call Z(1)%setup_mapping  (UPPER, POLOIDAL, Z(2))     ! connect left and right segments
     call Z(2)%setup_mapping  (UPPER, POLOIDAL, Z(1), 1)  ! connect left and right segments at interface R1

     ! primary SOL
     call Z(3)%setup_mapping  (UPPER, RADIAL,   Z(5), 12)  ! connect to right secondary SOL at PARTIAL interface I6 (this should only be used to find the poloidal side for the generating radial path!)
     call Z(4)%setup_mapping  (UPPER, RADIAL,   Z(6), 11)  ! connect to left secondary SOL at interface I5 (SAME NOTE AS ABOVE)
     call Z(3)%setup_mapping  (UPPER, POLOIDAL, Z(4))     ! connect left and right segments
     call Z(3)%setup_mapping  (LOWER, POLOIDAL, Z(7), 2)  ! connect to right divertor leg at interface R2
     call Z(7)%setup_boundary (LOWER, POLOIDAL, DIVERTOR) ! right divertor target
     call Z(4)%setup_mapping  (UPPER, POLOIDAL, Z(8))     ! connect to left divertor leg
     call Z(8)%setup_boundary (UPPER, POLOIDAL, DIVERTOR) ! left divertor target
     call Z(7)%setup_mapping  (UPPER, RADIAL,   Z(9), 5)  ! VIRTUAL interface I5
     call Z(8)%setup_mapping  (UPPER, RADIAL,   Z(12),6)  ! VIRTUAL interface I6

     ! secondary SOL
     ! right branch
     call Z(5)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(5)%setup_mapping  (LOWER, POLOIDAL, Z(9))     ! connect to right divertor leg (right branch)
     call Z(9)%setup_boundary (LOWER, POLOIDAL, DIVERTOR) ! right divertor target
     call Z(5)%setup_mapping  (UPPER, POLOIDAL, Z(10),10) ! connect to left divertor leg  (right branch) at interface R10
     call Z(10)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! left divertor target
     call Z(9)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(10)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     ! left branch
     call Z(6)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(6)%setup_mapping  (LOWER, POLOIDAL, Z(11),11) ! connect to right divertor leg (left branch) at interface R11
     call Z(11)%setup_boundary(LOWER, POLOIDAL, DIVERTOR) ! right divertor target
     call Z(6)%setup_mapping  (UPPER, POLOIDAL, Z(12))    ! connect to left divertor leg  (left branch)
     call Z(12)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! left divertor target
     call Z(11)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(12)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain

     ! primary PFR
     call Z(7)%setup_mapping  (LOWER, RADIAL,   Z(13), 3)  ! connect to right primary PFR at interface I3
     call Z(8)%setup_mapping  (LOWER, RADIAL,   Z(14), 4)  ! connect to left primary PFR at interface I4
     call Z(13)%setup_boundary(LOWER, POLOIDAL, DIVERTOR)  ! right divertor target
     call Z(13)%setup_mapping (UPPER, POLOIDAL, Z(14), 4)  ! connect to left primary PFR at interface R4
     call Z(14)%setup_boundary(UPPER, POLOIDAL, DIVERTOR)  ! left divertor target
     ! right snowflake PFR
     call Z(13)%setup_mapping (LOWER, RADIAL,   Z(17), 9)  ! connect to right primary PFR at interface I9
     call Z(17)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(18)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(17)%setup_boundary(LOWER, POLOIDAL, DIVERTOR)
     call Z(17)%setup_mapping (UPPER, POLOIDAL, Z(18), 6)  ! connect to left primary PFR at interface R6
     call Z(18)%setup_boundary(UPPER, POLOIDAL, DIVERTOR)
     ! left snowflake PFR
     call Z(14)%setup_mapping (LOWER, RADIAL,   Z(19), 7)  ! connect to right primary PFR at interface I7
     call Z(19)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(20)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(19)%setup_boundary(UPPER, POLOIDAL, DIVERTOR)
     call Z(19)%setup_mapping (LOWER, POLOIDAL, Z(20), 7)  ! connect to left primary PFR at interface R7
     call Z(20)%setup_boundary(LOWER, POLOIDAL, DIVERTOR)
     ! center snowflake PFR
     call Z(18)%setup_mapping (UPPER, RADIAL,   Z(21), 10) ! connect to right primary PFR at interface I10
     call Z(20)%setup_mapping (UPPER, RADIAL,   Z(22), 8)  ! connect to right primary PFR at interface I10
     call Z(21)%setup_boundary(UPPER, RADIAL,   VACUUM)    !
     call Z(22)%setup_boundary(UPPER, RADIAL,   VACUUM)    !
     call Z(21)%setup_boundary(UPPER, POLOIDAL, DIVERTOR)
     call Z(21)%setup_mapping (LOWER, POLOIDAL, Z(22), 8)  ! connect to left primary PFR at interface R6
     call Z(22)%setup_boundary(LOWER, POLOIDAL, DIVERTOR)


     ! secondary PFR
     call Z(10)%setup_mapping (LOWER, RADIAL,   Z(15),14)  ! connect to right secondary PFR at interface I10
     call Z(15)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(11)%setup_mapping (LOWER, RADIAL,   Z(16),13)  ! connect to left secondary PFR at interface I9
     call Z(16)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(16)%setup_mapping (UPPER, POLOIDAL, Z(15),12)  ! connect to left secondary PFR at interface R12
     call Z(16)%setup_boundary(LOWER, POLOIDAL, DIVERTOR)  ! right divertor target
     call Z(15)%setup_boundary(UPPER, POLOIDAL, DIVERTOR)  ! left divertor target

  end select
  call undefined_element_boundary_check(.false.)


  if (present(nx_out)) nx_out = nX
  if (present(nrpath_out)) nrpath_out = nrpath

 9000 format('error: invalid topology ', a, '!')
  end subroutine setup_topology
!=======================================================================



!=======================================================================
! Set up geometry of computational domain:
! Magnetic axis, X-points, separatrix(ces) and radial paths from X-points
! SHOULD THIS BE MOVED TO MODULE fieldline_grid?
!=======================================================================
  subroutine setup_geometry()
  use moose_utils
  use moose_math
  use moose_geometry, only: polygon2d, polygon2d_hypersurf, hypersurf2d
  use moose_interp_curve, only: arclength_parametrization
  use flare_model, only: boundary2d, equi2d
  use flare_fluxsurf2d, only: separatrix2d
  use flare_mmesh_parameters, only: d_SOL, d_PFR, auto_adjust_mesh
  use flare_mmesh_inner_boundary

  type(rpath2d_trace) :: Rtmp
  character(len=2) :: Sstr
  logical          :: reverse
  real(real64)     :: dx, tmp(3), theta_cut, psiN, x(2), x1(2), tau
  integer          :: ix, ierr, iSOL, iPFR, ipi, jx, k, orientation, o(4)


  ! 1. set up guiding surface for construction of base mesh(s) ---------
  call initialize_guiding_contour(0)
  ! ... done later for each block


  ! 2. CRITICAL POINTS -------------------------------------------------
  ! 2.a set up magnetic axis (Pmag) ------------------------------------
  Pmag = equi2d%r0
  write (6, 2000) Pmag
 2000 format(8x,'Magnetic axis at: ',2f10.4)


  ! 2.b check X-points -------------------------------------------------
  tmp(3) = 0.d0
  do ix=1,nX
     tmp(1:2) = equi2d%xpoint(ix)
     write (6, 2001) ix, tmp(1:2)
     write (6, 2002) equi2d%poloidal_angle(tmp)/pi*180.d0
  enddo
  write (6, *)
 2001 format(8x,i0,'. X-point at: ',2f10.4)
 2002 format(11x,'-> poloidal angle [deg]: ',f10.4)
 9001 format('error: ',i0,'. X-point is not defined!')


  ! 3. inner boundaries for plasma transport domain (EMC3) -------------
  ! ... already done


  ! 4.1 generate separatrix(ces) ---------------------------------------
  write (6, 4100)
  allocate (S(nX), R(nX,4))
  reverse = equi2d%Bp_sign == 1
  do ix=1,nX
     write (6, 4101) ix
     ! find cut-off poloidal angle for guiding X-point
     jx = connectX(ix)
     if (jx < 0  .and.  abs(jx).ne.ix  .and.  connectX(abs(jx)) == jx) then
        write (6, 4102) abs(jx), ix
        Rtmp = rpath2d_traceX(abs(jx), -1, RPATH2D_PSIN, 1.d0)
        R(abs(jx), DESCENT_CORE) = rpath2d_curve(Rtmp)
        poloidal_interface((abs(jx)-1)*4+1)%C => R(abs(jx), DESCENT_CORE)
        S(ix) = separatrix2d(ix, boundary=guiding_contour, reverse=reverse, Xadd=Rtmp%element(0))
     else
        S(ix) = separatrix2d(ix, boundary=guiding_contour, reverse=reverse)
     endif
     !write (Sstr, 2003) ix; call S(ix)%save(trim(Sstr), split=.true.)
  enddo
 2003 format('S',i0)
  if (connectX(1) == 1) then
     ! core segments of separatrix
     S0 = S(1)%branch(XCONTOUR_BRANCH(equi2d%Bp_sign,-1))%interp()

  elseif (connectX(1) > 1) then
     ! left and right core segment of separatrix
     S0R = S(1)%branch(XCONTOUR_BRANCH(equi2d%Bp_sign,-1))%interp()
     S0L = S(1)%branch(XCONTOUR_BRANCH(-equi2d%Bp_sign,-1))%interp()
  else
     ! setup left and right branch of core separatrix if second X-point is used for guidance
     S0R = S(1)%branch(XCONTOUR_BRANCH(equi2d%Bp_sign,-1))%interp()
     S0L = S(1)%branch(XCONTOUR_BRANCH(-equi2d%Bp_sign,-1))%interp()
  endif
  write (6, *)


  ! 4.2 generate radial paths from X-points ----------------------------
     ! and set up interfaces between elements
     !call Iface(1)%set_curve(S0)
  write (6, 3000)
  iSOL = 0
  iPFR = 0
  ipi  = 0
  x1   = equi2d%xpoint(1)
  do ix=1,nX
     jx = connectX(ix)

     o(ASCENT_LEFT ) = ASCENT_LEFT
     o(ASCENT_RIGHT) = ASCENT_RIGHT
     o(DESCENT_CORE) = DESCENT_CORE
     o(DESCENT_PFR ) = DESCENT_PFR
     ! orientation for secondary X-points
     orientation = 0
     if (ix > 1) then
        ! X-point is in private flux region of main X-point
        if (equi2d%psiN_X(ix) < equi2d%psiN_X(1)  .and.  jx /= 1) then

           x = equi2d%xpoint(ix)
           ! set direction of PFR and core
           if (x(1) < x1(1)) then
              o(DESCENT_PFR ) = ASCENT_LEFT
              o(DESCENT_CORE) = ASCENT_RIGHT
           else
              o(DESCENT_PFR ) = ASCENT_RIGHT
              o(DESCENT_CORE) = ASCENT_LEFT
           endif

           ! ...
           o(ASCENT_LEFT ) = DESCENT_CORE
           o(ASCENT_RIGHT) = DESCENT_PFR
        endif
     endif


     ! "core-interface"
     if (ix == 1  .or.  jx > 0) then
        write (6, 3010) ix
        R(ix, DESCENT_CORE) = rpath2d_curveX(ix, -1, RPATH2D_PSIN, psiN_1)
        !poloidal_interface(ipi+1)%C => R(ix, DESCENT_CORE)

        ! re-map to arclength
        allocate (poloidal_interface(ipi+1)%C, source=arclength_parametrization(R(ix, DESCENT_CORE)))
     else
        ! secondary X-point in radial connection to another X-point
!        if (abs(jx) .ne. ix) then
!           write (6, 3032) ix, abs(jx)
!           R(ix, DESCENT_CORE) = ...
!           call R(ix, DESCENT_CORE)%setup_linear(Xp(abs(jx))%X, Xp(ix)%X-Xp(abs(jx))%X)
!        endif
        if (R(ix, DESCENT_CORE)%interp_type == 0) then
           write (6, *) "ERROR: rpath segment ", ix, DESCENT_CORE, " is not set up!"
           stop
        endif
     endif


     ! scrape-off layer
     if (jx == ix) then
        ! single SOL, left and right branch on same flux surface
        iSOL = iSOL + 1
        write (6, 3020) ipi+2, ipi+3, ix, d_SOL(iSOL) * 1.d2
        R(ix, ASCENT_LEFT) = rpath2d_curveX(ix, -1, RPATH2D_ARCLENGTH, d_SOL(iSOL))
        poloidal_interface(ipi+3)%C => R(ix, ASCENT_LEFT)

        psiN = equi2d%psiN(R(ix, ASCENT_LEFT)%xb())
        R(ix, ASCENT_RIGHT) = rpath2d_curveX(ix, 1, RPATH2D_PSIN, psiN)
        !poloidal_interface(ipi+2)%C => R(ix, ASCENT_RIGHT)
        allocate (poloidal_interface(ipi+2)%C, source=arclength_parametrization(R(ix, ASCENT_RIGHT)))

     elseif (jx > ix) then
        ! connected X-points, left and right branch on individual flux surfaces
        iSOL = iSOL + 1
        write (6, 3021) ipi+2, ix, jx, d_SOL(iSOL) * 1.d2
        R(ix, ASCENT_LEFT) = rpath2d_curveX(ix, -1, RPATH2D_ARCLENGTH, d_SOL(iSOL))
        poloidal_interface(ipi+2)%C => R(ix, ASCENT_LEFT)
        psiN = equi2d%psiN(R(ix, ASCENT_LEFT)%xb())
        R(jx, ASCENT_RIGHT) = rpath2d_curveX(jx, 1, RPATH2D_PSIN, psiN)

        iSOL = iSOL + 1
        write (6, 3022) ipi+3, ix, jx, d_SOL(iSOL) * 1.d2
        R(ix, ASCENT_RIGHT) = rpath2d_curveX(ix, 1, RPATH2D_ARCLENGTH, d_SOL(iSOL))
        poloidal_interface(ipi+3)%C => R(ix, ASCENT_RIGHT)
        psiN = equi2d%psiN(R(ix, ASCENT_RIGHT)%xb())
        R(jx, ASCENT_LEFT) = rpath2d_curveX(jx, -1, RPATH2D_PSIN, psiN)

     elseif (jx == -ix) then
        ! outer SOL with left and right branch on individual flux surfaces
        iSOL = iSOL + 1
        write (6, 3023) ipi+2, ix, d_SOL(iSOL) * 1.d2
        k = 1;   if (bugfix_ddn) k = -1
        R(ix, ASCENT_LEFT) = rpath2d_curveX(ix, k, RPATH2D_ARCLENGTH, d_SOL(iSOL))
        poloidal_interface(ipi+2)%C => R(ix, ASCENT_LEFT)

        iSOL = iSOL + 1
        write (6, 3024) ipi+3, ix, d_SOL(iSOL) * 1.d2
        k = -1;   if (bugfix_ddn) k = 1
        R(ix, ASCENT_RIGHT) = rpath2d_curveX(ix, k, RPATH2D_ARCLENGTH, d_SOL(iSOL))
        poloidal_interface(ipi+3)%C => R(ix, ASCENT_RIGHT)

     elseif (jx < 0  .and.  connectX(abs(jx)) == jx) then
        ! this SOL's boundary is another separatrix
        psiN = equi2d%psiN_X(abs(jx))
        write (6, 3025) ipi+2, ipi+3, ix, abs(jx), psiN
        R(ix, ASCENT_LEFT) = rpath2d_curveX(ix, 1, RPATH2D_PSIN, psiN)
        !poloidal_interface(ipi+2)%C => R(ix, ASCENT_LEFT)
        allocate (poloidal_interface(ipi+2)%C, source=arclength_parametrization(R(ix, ASCENT_LEFT)))
        R(ix, ASCENT_RIGHT) = rpath2d_curveX(ix, -1, RPATH2D_PSIN, psiN)
        poloidal_interface(ipi+3)%C => R(ix, ASCENT_RIGHT)
     endif
     poloidal_interface(ipi+2)%inode(-1) = ix
     poloidal_interface(ipi+3)%inode(-1) = ix


     ! additional PFRs from X-point in PFR of the primary X-point
     if (jx < 0  .and.  connectX(abs(jx)) /= jx) then
        iPFR = iPFR + 1
        write (6, 3033) ipi+2, ix, d_PFR(iPFR) * 1.d2
        R(ix, ASCENT_LEFT) = rpath2d_curveX(ix, 1, RPATH2D_ARCLENGTH, d_PFR(iPFR))
        !call R(ix, ASCENT_LEFT)%flip()
        poloidal_interface(ipi+2)%C => R(ix, ASCENT_LEFT)

        iPFR = iPFR + 1
        write (6, 3034) ipi+3, ix, d_PFR(iPFR) * 1.d2
        R(ix, ASCENT_RIGHT) = rpath2d_curveX(ix, -1, RPATH2D_ARCLENGTH, d_PFR(iPFR))
        !call R(ix, ASCENT_RIGHT)%flip()
        poloidal_interface(ipi+3)%C => R(ix, ASCENT_RIGHT)

        poloidal_interface(ipi+2)%inode(-1) = UNDEFINED
        poloidal_interface(ipi+3)%inode(-1) = UNDEFINED
     endif
     ! private flux region
     iPFR = iPFR + 1
     do jx=1,nX
        if (jx == ix) cycle
        if (connectX(jx) == -ix) then
           if (connectX(ix) == -ix) cycle
           exit
        endif
     enddo
     ! regular PFR
     if (jx > nX) then
        write (6, 3030) ipi+4, ix, d_PFR(iPFR) * 1.d2
        R(ix, DESCENT_PFR) = rpath2d_curveX(ix, 1, RPATH2D_ARCLENGTH, -d_PFR(iPFR))
        !if (o(DESCENT_PFR) == DESCENT_PFR) call R(ix, DESCENT_PFR)%flip()
        poloidal_interface(ipi+4)%C => R(ix, DESCENT_PFR)

     ! connect to another X-points
     else
        write (6, 3032) ix, jx
        !call R(ix, DESCENT_PFR)%setup_linear(Xp(jx)%X, Xp(ix)%X-Xp(jx)%X)
        !call poloidal_interface(ipi+4)%set_curve(R(ix, DESCENT_PFR)%t_curve)
     endif
     poloidal_interface(ipi+4)%inode(1) = ix


     ! plot paths
!     do k=1,4
!        call R(ix, k)%plot(filename='rpath_'//trim(str(ix))//'_'//trim(str(k))//'.plt')
!     enddo

     ipi = ipi + 4
  enddo
  write (6, *)

  do ipi=1,poloidal_interfaces
     !call poloidal_interface(ipi)%C%plot(filename='R'//trim(str(ipi))//'.plt')
     call poloidal_interface(ipi)%C%savetxt('R'//trim(str(ipi))//'.dat')

     ! check intersection with guiding surface
     if (auto_adjust_mesh) then
     if (boundary2d%intersect(polygon2d(transpose(poloidal_interface(ipi)%C%u)), posi=tau)) then
        write (6, 9420) ipi

        dx = poloidal_interface(ipi)%C%b
        if (poloidal_interface(ipi)%inode(-1) > 0) then
           write (6, 9421) dx * tau
        elseif (poloidal_interface(ipi)%inode( 1) > 0) then
           write (6, 9421) dx * (1.d0 - tau)
        endif
        stop
     endif
     endif
  enddo


 4100 format(3x,'- Setting up block-structured decomposition')
 4101 format(8x,'generating separatrix for X-point ', i0)
 4102 format(8x,'connecting X-point ', i0, ' along radial path to separatrix ', i0)
 3000 format(3x,'- Setting up radial paths for block-structured decomposition')
 3010 format(8x,'generating core segment for X-point ', i0)
 3020 format(8x,i0,', ',i0,': generating SOL segment for X-point ', i0, ' (length = ', f0.3, ' cm)')
 3021 format(8x,i0,': generating left SOL segment for X-points ', i0, ', ', i0, ' (length = ', f0.3, ' cm)')
 3022 format(8x,i0,': generating right SOL segment for X-points ', i0, ', ', i0, ' (length = ', f0.3, ' cm)')
 3023 format(8x,i0,': generating left SOL segment for X-point ', i0, ' (length = ', f0.3, ' cm)')
 3024 format(8x,i0,': generating right SOL segment for X-point ', i0, ' (length = ', f0.3, ' cm)')
 3025 format(8x,i0,',',i0,': generating SOL segments for X-point ', i0, ' up to separatrix from X-point ', &
                i0, ' at psiN = ', f0.3)
 3030 format(8x,i0,': generating PFR segment for X-point ', i0, ' (length = ', f0.3, ' cm)')
 3032 format(8x,': generating PFR segment from X-point ', i0, ' to X-point ', i0)
 3033 format(8x,i0,': generating left PFR segment for X-point ', i0, ' (length = ', f0.3, ' cm)')
 3034 format(8x,i0,': generating right PFR segment for X-point ', i0, ' (length = ', f0.3, ' cm)')
 9420 format('ERROR: reference path for radial discretization (R',i0,'.plt) crosses boundary!')
 9421 format('d < ', f0.3, ' is required!')
  end subroutine setup_geometry
!=======================================================================



!=======================================================================
  subroutine setup_interfaces(ni)
  use moose_utils
  use flare_model, only: equi2d
  use mesh_interface
  integer, intent(out), optional :: ni


  integer :: ix, ix1, jx, iri, iconnect, k


  write (6, 1000)
  iri = 0
  do ix=1,nX
     jx = connectX(ix)

     ! connect back to same X-point
     if (jx == ix) then
        if (ix .ne. 1) then
           write (6, 9000) ix
           stop
        endif

        iri = iri + 1
        radial_interface(iri)%C => S0
        call radial_interface(iri)%setup(ix, ix)
        write (6, 1001) iri, ix

     ! all branches connect to divertor targets
     elseif (jx == -ix  .or.  (jx < 0  .and.  connectX(abs(jx)) /= jx)) then
        iconnect  = STRIKE_POINT
        ! are these "upstream" branches?
        do ix1=1,ix-1
           if (abs(connectX(ix1)) == ix) then
              iconnect = -ix1
              exit
           endif
        enddo
        iri = iri + 1
        if (iconnect == STRIKE_POINT) then
           k = XCONTOUR_BRANCH(equi2d%Bp_sign, 1)
           allocate (radial_interface(iri)%C, source=S(ix)%branch(k)%interp())
        endif
        call radial_interface(iri)%setup(ix, iconnect)
        write (6, 1002) iri, ix
        iri = iri + 1
        if (iconnect == STRIKE_POINT) then
           k = XCONTOUR_BRANCH(-equi2d%Bp_sign, 1)
           allocate (radial_interface(iri)%C, source=S(ix)%branch(k)%interp())
        endif
        call radial_interface(iri)%setup(iconnect, ix)
        write (6, 1002) iri, ix

     ! connect to other X-point OR
     ! main separatrix decomposition is guided by secondary X-point
     elseif (jx > ix  .or.  jx < 0) then
        if (ix .ne. 1) then
           if (jx > ix) then
              write (6, 9001) ix
           else
              write (6, 9002) ix, abs(jx)
           endif
           stop
        endif

        ! right core interface
        iri = iri + 1
        radial_interface(iri)%C => S0R
        call radial_interface(iri)%setup(ix, jx)
        write (6, 1003) iri, ix, jx

        ! left core interface
        iri = iri + 1
        radial_interface(iri)%C => S0L
        call radial_interface(iri)%setup(jx, ix)
        write (6, 1003) iri, ix, jx

     ! nothing to be done here anymore
     else

     endif

     ! divertor branches
     iri = iri + 1
     k = XCONTOUR_BRANCH(-equi2d%Bp_sign, 1)
     allocate (radial_interface(iri)%C, source=S(ix)%branch(k)%interp())
     call radial_interface(iri)%setup(STRIKE_POINT, ix)
     write (6, 1004) iri, ix
     iri = iri + 1
     k = XCONTOUR_BRANCH(equi2d%Bp_sign, 1)
     allocate (radial_interface(iri)%C, source=S(ix)%branch(k)%interp())
     call radial_interface(iri)%setup(ix, STRIKE_POINT)
     write (6, 1004) iri, ix
     ! add divertor branches for outer separatrix
     if (jx < -ix) then
        iri = iri + 1
        call radial_interface(iri)%setup(STRIKE_POINT, -ix)
        write (6, 1004) iri, abs(jx)
        iri = iri + 1
        call radial_interface(iri)%setup(-ix, STRIKE_POINT)
        write (6, 1004) iri, abs(jx)
     endif


  enddo
  print *


  !write (6, *) 'radial interfaces:'
  do iri=1,radial_interfaces
     !write (6, *) iri, radial_interface(iri)%inode(-1), radial_interface(iri)%inode(1)
     if (associated(radial_interface(iri)%C)) then
        call radial_interface(iri)%C%plot(filename='I'//trim(str(iri))//'.plt')
        call radial_interface(iri)%C%savetxt('I'//trim(str(iri))//'.dat')
     endif
  enddo
  if (present(ni)) ni = radial_interfaces

 1000 format(3x,'- Setting up radial interfaces:')
 1001 format(8x,i0,': X-point ',i0,' is connected back to itself')
 1002 format(8x,i0,': main branch of separatrix ',i0,' is connected to divertor target')
 1003 format(8x,i0,': main branch of separatrix ',i0,' is connected to/guided by X-point ',i0)
 1004 format(8x,i0,': divertor branch for separatrix ',i0)
 9000 format('error: seconday X-point ', i0, ' connects back to itself!')
 9001 format('error: seconday X-point ', i0, ' does not connect back to primary one!')
 9002 format('error: seconday X-point ', i0, ' used as guiding point for separatrix ', i0, '!')
  end subroutine setup_interfaces
!=======================================================================



!=======================================================================
  subroutine initialize_guiding_contour(iblock)
  use moose_error
  use moose_utils,    only: split
  use moose_math,     only: pi
  use moose_r3grid, only: length_scale
  use moose_geometry, only: polygon2d, hypersurf2d
  use flare_boundary, only: firstwall_rzslice
  use flare_model
  use flare_mmesh_parameters, only: T, filename => guiding_contour, &
                                       filename_L => guiding_contour_L, &
                                       filename_R => guiding_contour_R
  integer, intent(in) :: iblock

  character(len=256) :: gc, units


  gc = filename(iblock);   if (gc == "") gc = filename(-1)
  ! guiding contour from model boundary
  if (gc == "") then
     print 1001
     guiding_contour = hypersurf2d(firstwall_rzslice(T(iblock)%phi_base / 180.d0 * pi))
     guiding_contour_L => guiding_contour
     guiding_contour_R => guiding_contour


  ! user defined guiding contour
  else
     call split(gc, gc, units, set=':', default='m')
     print 1002, trim(gc)
     guiding_contour = hypersurf2d(polygon2d(gc, scale=length_scale(units)))
     guiding_contour_L => guiding_contour
     guiding_contour_R => guiding_contour

     gc = filename_L(iblock);   if (gc == "") gc = filename_L(-1)
     if (gc /= "") then
        print 1012, "L", trim(gc)
        call split(gc, gc, units, set=':', default='m')
        allocate (guiding_contour_L, source = hypersurf2d(polygon2d(gc, scale=length_scale(units))))
     endif

     gc = filename_R(iblock);   if (gc == "") gc = filename_R(-1)
     if (gc /= "") then
        print 1012, "R", trim(gc)
        call split(gc, gc, units, set=':', default='m')
        allocate (guiding_contour_R, source = hypersurf2d(polygon2d(gc, scale=length_scale(units))))
     endif
  endif
 1001 format(3x,"- Mesh generation is guided by first model boundary"/)
 1002 format(3x,"- Mesh generation is guided by ",a/)
 1012 format(3x,"- Mesh generation (",a,") is guided by ",a/)


  end subroutine initialize_guiding_contour
!=======================================================================



!=======================================================================
  subroutine initialize_layers(iblock)
  use moose_math, only: pi
  use flare_mmesh_parameters
  integer, intent(in) :: iblock

  integer, parameter :: COUNT_RUN = 1, SETUP_RUN = 2
  integer, dimension(:), allocatable :: markz, izl
  character(len=1) :: cside(-1:1) = (/'R','C','L'/)
  real(real64) :: phi
  integer      :: i, il, iz, iz0, iz_map, idir, irun, nzl(-1:1)


  ! initialize block
  do il=0,max_layers-1
     nr(il) = tblock(iblock)%nr(il)
     np(il) = tblock(iblock)%np(il)
     npL(il) = tblock(iblock)%npL(il)
     npR(il) = tblock(iblock)%npR(il)
  enddo
  phi = T(iblock)%phi_base / 180.d0 * pi


  write (6, 1000)
  allocate (markz(nelement), izl(-nelement:nelement))

  do irun=COUNT_RUN,SETUP_RUN
  markz = 0
  il    = 0
  if (irun == SETUP_RUN) allocate(L(0:layers-1))
  layer_loop: do
     ! exit if no more elements are unmarked
     if (sum(markz) == nelement) exit

     ! start new layer
     il  = il + 1
     izl = 0;  nzl = 0

     ! find first unmarked element
     do iz=1,nelement
        if (markz(iz) == 0) exit
     enddo

     ! set base element in layer
     iz0 = iz;  markz(iz) = 1;  izl(0) = iz0

     ! scan in both poloidal directions
     dir_loop: do idir=-1,1,2
        ! start poloidal scan at base element
        iz = iz0
        poloidal_scan: do
           iz_map = Z(iz)%map_p(idir)
           ! poloidal scan in both directions finished when returning to base element
           if (iz_map == PERIODIC) exit dir_loop
           if (iz_map == iz0) exit dir_loop
           ! poloidal scan in this direction finished at divertor targets
           if (iz_map == DIVERTOR) exit

           iz = iz_map;  markz(iz) = 1
           nzl(idir) = nzl(idir) + 1;  izl(idir*nzl(idir)) = iz
        enddo poloidal_scan
     enddo dir_loop

     if (irun == SETUP_RUN) then
        ! now set up element indices in this layer
        call L(il-1)%initialize(nzl(-1)+1+nzl(1), izl(-nzl(-1):nzl(1)), nzl(-1)+1)

        ! set up resolution in each element
        call L(il-1)%setup_resolution(il-1)

        ! map poloidal resolution in radial direction
        call L(il-1)%map_poloidal_resolution()
     endif
  enddo layer_loop
  layers = il
  enddo


  ! initialize toroidal discretization
  do il=0,layers-1;  L(il)%T = T(iblock);  enddo
  do iz=1,nelement;  Z(iz)%T = T(iblock);  enddo


  ! initialize mesh for each domain element
  allocate (Mtmp(nelement))
  allocate (M(0:layers-1))
  do il=0,layers-1
     write (6, 1001) il
     do i=1,L(il)%nz
        iz = L(il)%iz(i)
        call Mtmp(iz)%initialize(Z(iz)%nr, Z(iz)%np, phi)
        write (6, 1002) iz, Z(iz)%np, Z(iz)%nr, Z(iz)%ipl, cside(Z(iz)%ipl_side)
     enddo
     call M(il)%initialize(L(il)%nr, L(il)%np, phi)
  enddo


  ! cleanup
  deallocate (markz, izl)

 1000 format(3x,'- Setting up radial layers:')
 1001 format(8x,'Layer ',i0)
 1002 format(8x,i3,': ',i5,' x ',i3,5x,'(',i0,a1,')')
  end subroutine initialize_layers
!=======================================================================



!=======================================================================
  subroutine generate_equi2d_base_mesh(iblock)
  use moose_quantiles
  use flare_mmesh_parameters
  integer, intent(in) :: iblock

  logical, parameter :: Debug = .false.

  integer, dimension(:), allocatable :: markz
  class(qfunc), allocatable :: Sp, SpL, SpR, Sr
  integer         :: i, il, il0, iz, iz0, iz_map, iside


  ! initialize block
  !call initialize_guiding_contour(iblock)
  call initialize_layers(iblock)


  ! generate core-interface
  il = 0
  if (connectX(1) == 1) then
     ! single element
     allocate (Sp, source=generate_quantile_function(poloidal_spacing(il)))
     call Mtmp(1)%setup_boundary_nodes(UPPER, RADIAL, S0, Sp)
  else
     ! left and right elements
     allocate (SpR, source=generate_quantile_function(poloidal_spacing_R(il)))
     call Mtmp(1)%setup_boundary_nodes(UPPER, RADIAL, S0R, SpR)

     allocate (SpL, source=generate_quantile_function(poloidal_spacing_L(il)))
     call Mtmp(2)%setup_boundary_nodes(UPPER, RADIAL, S0L, SpL)
  endif


  ! main loop: generate mesh for each layer
  write (6, *)
  allocate (markz(nelement));  markz = 0
  do il=0,layers-1
     ! base element index
     iz0 = L(il)%iz(L(il)%i0)

     ! set up radial spacings for this layer
     allocate (Sr, source=generate_quantile_function(radial_spacing(il)))

     ! generate mesh for this layer
     call generate_layer(il, iz0, iblock, Sr)

     ! map radial interface to next element
     do i=1,L(il)%nz
        iz = L(il)%iz(i)
        markz(iz) = 1 ! mark elements in this layer

        do iside=-1,1,2
           iz_map = Z(iz)%map_r(iside)
           ! map to another element?
           if (iz_map < 0) cycle

           ! map poloidal resolution
           if (markz(iz_map) == 0) then
              call Mtmp(iz)%connect_to(Mtmp(iz_map), RADIAL, iside)
              if (Debug) then
                 call Mtmp(iz)%plot_boundary(RADIAL, iside, 'MAP', iz)
                 call Mtmp(iz_map)%plot_boundary(RADIAL, -iside, 'MAPPED', iz_map)
              endif

              ! status of radial interface
              !write (6, *) 'radial interface in element ', iz, ' side ', iside, ' is ', Z(iz)%rad_bound(iside)
           endif
        enddo
     enddo
     call Sr%free();   deallocate (Sr)
  enddo


  ! write output files
  il0 = iblock * layers
  do il=0,layers-1
     iz = il0 + il
     call M(il)%save("base", iz)
  enddo


  ! cleanup
  deallocate (M, Mtmp, markz, L)

  end subroutine generate_equi2d_base_mesh
!=======================================================================



!=======================================================================
  subroutine generate_layer(il, iz0, iblock, Sr)
  use moose_quantiles
  use moose_analysis, only: interp, loadtxt_interp, pchip, INTERP_PCHIP
  use moose_utils, only: str
  use flare_mmesh_parameters, only: poloidal_spacing, poloidal_spacing_L, poloidal_spacing_R, &
                            upstream_adjust_L, upstream_adjust_R
  integer,         intent(in) :: il, iz0, iblock
  class(qfunc),    intent(in) :: Sr

  logical, parameter :: Debug = .false.

  class(qfunc), allocatable :: Sp
  type(interp)   :: U(-1:1)
  integer :: i, idir, irside, iri, ir1, ipside, ipi, ip0, iz, iz_map, npz(-1:1)


  write (6, 1000) il, iz0
  !if (Debug) write (6, 1001) Mtmp(iz0)%fixed_coord_value
  write (6, *) 'reference discretization at:'


  ! select upper or lower radial boundary for reference nodes
  if (Mtmp(iz0)%ir0 == 0) then
     irside = LOWER
  elseif (Mtmp(iz0)%ir0 > 0) then
     irside = UPPER
  else
     write (6, 9000) il, iz0
     write (6, 9001)
     stop
  endif
  iri = Z(iz0)%rad_bound(irside)
  write (6, *) 'radial boundary ', irside, ' which is interface ', iri
  if (iri == UNDEFINED) then
     write (6, 9000) il, iz0
     write (6, 9002)
     stop
  endif


  ! select upper or lower poloidal boundary -> start from an X-point
  do ipside=-1,1,2
     ipi = Z(iz0)%pol_bound(ipside)
     if (radial_interface(iri)%inode(ipside) > 0  .and.  ipi /= UNDEFINED) exit
  enddo
  if (ipi == UNDEFINED) then
     write (6, 9000) il, iz0
     write (6, 9003)
     stop
  endif


  ! initialize radial discretization in layer
  ir1 = 0
  if (Z(iz0)%map_r(-1) == CORE) ir1 = 1
  call Mtmp(iz0)%setup_boundary_nodes(ipside, POLOIDAL, poloidal_interface(ipi)%C, Sr, ir1, debug=Debug)
  U(LEFT)  = load_upstream_adjust(upstream_adjust_L(il))
  U(RIGHT) = load_upstream_adjust(upstream_adjust_R(il))
  ! diagnostic output
  ip0 = 0;   if (ipside == 1) ip0 = Mtmp(iz0)%np
  open  (newunit=iz, file="TMP2_INIT_LAYER_"//str(il))
  do i=ir1,Mtmp(iz0)%nr
     write (iz, *) Mtmp(iz0)%mesh(i,ip0,:)
  enddo
  close (iz)


  ! generate mesh in base element
  ! no Sp needed in base element -> only for divertor legs
  select case(Z(iz0)%ipl_side)
  case(LEFT)
     allocate (Sp, source=generate_quantile_function(poloidal_spacing_L(Z(iz0)%ipl)))
  case(CENTER)
     allocate (Sp, source=generate_quantile_function(poloidal_spacing(Z(iz0)%ipl)))
  case(RIGHT)
     allocate (Sp, source=generate_quantile_function(poloidal_spacing_R(Z(iz0)%ipl)))
  end select
  call Z(iz0)%generate_mesh(Mtmp(iz0), irside, ipside, iblock, Sr, Sp, U(ipside), debug=Debug)
  call upstream_adjust(iz0)
  !if (Debug) call Mtmp(iz0)%plot_mesh('Mtmp'//trim(str(iz0))//'.plt')


  ! scan through poloidal elements in this layer
  npz    = 0
  npz(0) = 1
  idir_loop: do idir=1,-1,-2
     iz = iz0
     poloidal_scan: do
        iz_map = Z(iz)%map_p(idir)
        ! 1. poloidal boundary of layer?
        ! 1.1 periodic boundaries: connect element back to itself
        if (iz_map == PERIODIC) then
           call Mtmp(iz)%connect_to(Mtmp(iz), POLOIDAL, LOWER_TO_UPPER)
           call Mtmp(iz)%save('Mtmp', iz)
           exit idir_loop
        endif
        ! 1.2 back to initial/base element
        if (iz_map == iz0) then
           call Mtmp(iz)%connect_to(Mtmp(iz0), POLOIDAL, LOWER_TO_UPPER)
           call Mtmp(iz)%save('Mtmp', iz)
           exit idir_loop
        endif
        ! 1.3 divertor targets
        if (iz_map == DIVERTOR) then
           call Mtmp(iz)%save('Mtmp', iz)
           exit
        endif

        ! 2. connect mesh to next element
        write (6, *) 'connect element ', iz, ' to ', iz_map
        call Mtmp(iz)%connect_to(Mtmp(iz_map), POLOIDAL, idir)
        call Mtmp(iz)%save('Mtmp', iz)
        iz        = iz_map

        ! 3. generate mesh in next element
        call Sp%free();   deallocate (Sp)
        select case(Z(iz)%ipl_side)
        case(LEFT)
           allocate (Sp, source=generate_quantile_function(poloidal_spacing_L(Z(iz)%ipl)))
        case(CENTER)
           allocate (Sp, source=generate_quantile_function(poloidal_spacing(Z(iz)%ipl)))
        case(RIGHT)
           allocate (Sp, source=generate_quantile_function(poloidal_spacing_R(Z(iz)%ipl)))
        end select
        call Z(iz)%generate_mesh(Mtmp(iz), irside, ipside, iblock, Sr, Sp, U(Z(iz)%ipl_side))
        call upstream_adjust(iz)
        !if (Debug) call Mtmp(iz)%plot_mesh('Mtmp'//trim(str(iz))//'.plt')
        npz(idir) = npz(idir) + 1
     enddo poloidal_scan
  enddo idir_loop
  write (6, *)


  ! merge elements
  ip0 = 0
  do i=1,L(il)%nz
     iz  = L(il)%iz(i)
     call M(il)%copy(0, ip0, Mtmp(iz))
     ip0 = ip0 + Z(iz)%np
  enddo

 1000 format('Generate layer ', i0, ' from base element ', i0)
 1001 format(3x,'- reference discretization at: ',f0.3,' deg')
 9000 format('error in generate_layer for il, iz0 = ', i0, ', ', i0)
 9001 format('undefined reference nodes on radial boundary!')
 9002 format('undefined radial interface!')
 9003 format('undefined poloidal interface!')
  contains
  !.....................................................................
  subroutine upstream_adjust(iz)
  use moose_table
  integer, intent(in) :: iz

  character(len=256) :: filename
  type(table) :: D
  type(interp) :: U, V
  integer       :: iz_map, idir


  do idir=-1,1,2
     iz_map = Z(iz)%map_p(idir)

     ! periodic boundary or already at divertor mesh?
     if (iz_map == PERIODIC  .or.  iz_map == DIVERTOR) cycle

     ! next mesh is for divertor leg
     if (Z(iz_map)%map_p(idir) /= DIVERTOR) cycle

     select case(Z(iz_map)%ipl_side)
     case(LEFT)
        filename = upstream_adjust_L(il)
     case(RIGHT)
        filename = upstream_adjust_R(il)
     case default
        write (6, *) "ERROR: THIS SHOULD NOT HAPPEN"
        stop
     end select
     if (filename == '') cycle
     if (Debug) write (6, *) 'UPSTREAM ADJUST:', Z(iz_map)%ipl, iz, iz_map, &
        Z(iz_map)%map_p(idir), Z(iz_map)%ipl_side, trim(filename)

     D = table(filename, 3)
     U = pchip(D%values(:,1), D%values(:,2) / 1.d2)
     V = pchip(D%values(:,1), D%values(:,3) / 1.d2)
     call Mtmp(iz)%upstream_adjust(Z(iz_map)%ipl_side, U, V, Sr)
  enddo

  end subroutine upstream_adjust
  !.....................................................................
  function load_upstream_adjust(filename) result(U)
  character(len=*), intent(in) :: filename
  type(interp) :: U


  ! NOTE: this is only required for t_mfs_mesh%upstream_adjust_divertor_leg in order
  ! to find the radial index at which an upstream adjustment needs to be done in
  ! the quasi-orthogonal part of the divertor leg
  if (filename == '') then
     U = pchip([0.d0, 1.d0*Z(iz0)%nr], [0.d0, 0.d0])
  else
     U = loadtxt_interp(filename, INTERP_PCHIP)
  endif

  end function load_upstream_adjust
  !.....................................................................
  end subroutine generate_layer
!=======================================================================



!=======================================================================
  subroutine init_equi2d_base_mesh_generator()
  use flare_model, only: assert_equi2d
  use flare_mmesh_parameters


  ! set up geometry of computational domain
  call assert_equi2d("init_base_mesh_generator")
  call setup_topology()
  call setup_geometry()
  call setup_interfaces()

  end subroutine init_equi2d_base_mesh_generator
!=======================================================================

end module base_mesh
