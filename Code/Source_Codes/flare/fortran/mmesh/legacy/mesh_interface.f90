module mesh_interface
!module zone_interface
  use iso_fortran_env
  use moose_geometry, only: interp_curve
  implicit none
  private


  integer, parameter, public :: &
     STRIKE_POINT = 0, &
     UNDEFINED    = -100


  type, public :: t_mesh_interface
     integer       :: id

     ! geometric definition of interface
     class(interp_curve), pointer :: C => null()
     integer       :: inode(-1:1) = UNDEFINED ! lower and upper end node type:  > 0 X-point
                                  !                                 = 0 strike point
                                  !                                 < guiding point (other X-point)

     ! adjacent zones:
     ! Z(-1) = zone on lower side of interface, i.e. interface is upper zone boundary
     ! Z( 1) = zone on upper side of interface, i.e. interface is lower zone boundary
     integer :: Z(-1:1) = -1


     ! discretization of interface
     real(real64), dimension(:,:), allocatable :: x
     ! number of nodes
     integer :: n

     contains
     procedure :: setup
     procedure :: setup_discretization
     procedure :: geometry_undefined
  end type t_mesh_interface
  !@public :: setup_interfaces

  ! radial interfaces between zones
  type(t_mesh_interface), dimension(:), allocatable, public :: radial_interface
  integer, public :: radial_interfaces
  ! poloidal interfaces between zones
  type(t_mesh_interface), dimension(:), allocatable, public :: poloidal_interface
  integer, public :: poloidal_interfaces


  ! number of X-points in computation domain, and their connectivity
  integer, dimension(:), allocatable :: connectX
  integer :: nX




  public :: initialize_interfaces
  public :: initialize_poloidal_interfaces

  contains
!=======================================================================



!=======================================================================
  subroutine setup(this, lower_boundary, upper_boundary)
  use flare_model, only: equi2d
  class(t_mesh_interface) :: this
  integer, intent(in)     :: lower_boundary, upper_boundary

  real(real64) :: X(2)
  integer      :: ix, iside, inode(-1:1), ierr


  inode(-1) = lower_boundary
  inode( 1) = upper_boundary

  do iside=-1,1,2
     ix = inode(iside)
     this%inode(iside) = ix

     ! check input
     if (abs(ix) > equi2d%X%nelements()) then
        write (6, 9001) lower_boundary, upper_boundary
        stop
     endif
  enddo
   

  ! 2. generate interface geometry
  ! call F%generate_branch
  ! for stability: connect 2 X-points by joining curves!


 9001 format('error in t_mesh_interface%setup: invalid X-point IDs = ', 2i0)
 9002 format('error in t_mesh_interface%setup: X-point ', i0, 'is not defined!')
  end subroutine setup
!=======================================================================



!=======================================================================
  subroutine setup_discretization(this, n, S)
  use moose_quantiles
  class(t_mesh_interface)   :: this
  integer,       intent(in) :: n
  class(qfunc)              :: S

  integer      :: i


  this%n = n
  do i=0,n
     this%x(i,1:2) = this%C%eval(i, n, S)
  enddo

  end subroutine setup_discretization
!=======================================================================



!=======================================================================
  function geometry_undefined(this)
  class(t_mesh_interface) :: this
  logical                 :: geometry_undefined


  geometry_undefined = .false.
  if (.not.associated(this%C)) geometry_undefined = .true.

  end function geometry_undefined
!=======================================================================



!=======================================================================
  subroutine initialize_interfaces(n)
  integer, intent(in) :: n


  allocate (radial_interface(n))
  radial_interfaces = n

  end subroutine initialize_interfaces
!=======================================================================



!=======================================================================
  subroutine initialize_poloidal_interfaces(n)
  integer, intent(in) :: n

  integer :: i


  allocate (poloidal_interface(n))
  poloidal_interfaces = n
  do i=1,n
     poloidal_interface(i)%id = i
  enddo


  end subroutine initialize_poloidal_interfaces
!=======================================================================



!=======================================================================
  subroutine setup_zones

  integer :: i, iz

  ! set up zone numbers for interfaces
  iz = 1
  do i=1,radial_interfaces
  enddo

  end subroutine setup_zones
!=======================================================================

end module mesh_interface
