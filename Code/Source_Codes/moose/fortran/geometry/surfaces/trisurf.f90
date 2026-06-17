!===============================================================================
! Linear approximation of torosurf (each quadrilateral surface patch is split
! into 2 triangles)
!===============================================================================
module moose_trisurf
  use iso_fortran_env
  use moose_surface, only: surface_coords
  implicit none
  private


  ! triangular approximation of torosurf
  type, public :: trisurf
     real(real64), allocatable :: phi(:), x(:,:,:)
     real(real64), allocatable :: tedge(:,:,:,:), pedge(:,:,:,:), diag(:,:,:)
     real(real64), allocatable :: nvect(:,:,:,:), area(:,:,:)
     logical, allocatable :: masked(:,:,:)
     integer :: nu, nv

     contains
     procedure :: broadcast, free
     procedure :: rzslice, winding_number, ray_intersect, interp, normal_vector
     procedure :: tmesh => export_tmesh
  end type trisurf


  interface trisurf
     procedure :: init
  end interface trisurf



  ! surface coordinate constructor for triangular representation
  interface surface_coords
     procedure :: construct_surface_coords
  end interface surface_coords



  public :: &
     kindex

  contains
  !-----------------------------------------------------------------------------


! auxiliary procedures:
  !-----------------------------------------------------------------------------
  function construct_surface_coords(i, j, k, u, v) result(S)
  !
  ! generate standard surface coordinates from triangular representation
  !
  integer,      intent(in) :: i, j, k
  real(real64), intent(in) :: u, v
  type(surface_coords)     :: S


  if (k == 0) then
     S = surface_coords(i, j, u, v)
  else
     S = surface_coords(i, j, 1.d0-u, 1.d0-v)
  endif

  end function construct_surface_coords
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function kindex(u, v)
  !
  ! return k-index (sub-triangle number) for surface coordinates (u, v)
  !
  real(real64), intent(in) :: u, v
  integer                  :: kindex


  kindex = 0
  if (u + v > 1.d0) kindex = 1

  end function kindex
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(T) result(this)
  !
  ! initialize trisurf from torosurf *T*
  !
  !  i, j+1  O----O   i+1, j+1
  !          |\k=1|
  !          | \  |
  !          |  \ |
  !          |k=0\|
  !  i, j    O----O   i+1, j
  !
  use moose_math,     only: cross_product
  use moose_torosurf, only: torosurf
  type(torosurf), intent(in) :: T
  type(trisurf)              :: this

  integer, parameter :: eps10 = 1.d-10

  real(real64) :: a, cos_phi, sin_phi
  integer :: i, j, k


  ! set up nodes
  this%nu = T%nu
  this%nv = T%nv
  allocate (this%phi(0:this%nu), source = T%phi)
  allocate (this%x(3, 0:this%nv, 0:this%nu))
  do i=0,this%nu
     cos_phi = cos(T%phi(i))
     sin_phi = sin(T%phi(i))
     do j=0,this%nv
        this%x(1, j, i) = cos_phi * T%rz(1, j, i)
        this%x(2, j, i) = sin_phi * T%rz(1, j, i)
        this%x(3, j, i) = T%rz(2, j, i)
     enddo
  enddo


  ! set up edges
  allocate (this%tedge(3, 0:1, 0:this%nv, 0:this%nu-1))
  allocate (this%pedge(3, 0:1, 0:this%nv-1, 0:this%nu))
  allocate (this%diag(3, 0:this%nv-1, 0:this%nu-1))
  do i=0,this%nu-1
  do j=0,this%nv-1
     ! set up toroidal edges
     this%tedge(:, 0, j, i) = this%x(:, j, i+1) - this%x(:, j, i)
     this%tedge(:, 1, j, i) = this%x(:, j+1, i) - this%x(:, j+1, i+1)

     ! set up poloidal edges
     this%pedge(:, 0, j, i) = this%x(:, j+1, i) - this%x(:, j, i)
     this%pedge(:, 1, j, i) = this%x(:, j, i+1) - this%x(:, j+1, i+1)

     ! set up diagonals
     this%diag(:, j, i) = this%x(:, j+1, i) - this%x(:, j, i+1)
  enddo
  enddo


  ! set up normal vectors and surface area
  allocate (this%nvect(3, 0:1, 0:this%nv-1, 0:this%nu-1), source = 0.d0)
  allocate (this%area(0:1, 0:this%nv-1, 0:this%nu-1), source = 0.d0)
  allocate (this%masked(0:1, 0:this%nv-1, 0:this%nu-1), source = .true.)
  do i=0,this%nu-1
  do j=0,this%nv-1
  do k=0,1
     this%nvect(:, k, j, i) = cross_product(this%tedge(:, k, j, i), this%pedge(:, k, j, i))
     a = norm2(this%nvect(:, k, j, i))
     if (a > eps10 * norm2(this%tedge(:, k, j, i)) * norm2(this%pedge(:, k, j, i))) then
        this%nvect(:, k, j, i) = this%nvect(:, k, j, i) / a
        this%area(k, j, i) = abs(a) / 2
        this%masked(k, j, i) = .false.
     endif
  enddo
  enddo
  enddo

  end function init
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(trisurf), intent(inout) :: this


  call proc(0)%broadcast(this%nu)
  call proc(0)%broadcast(this%nv)
  call proc(0)%broadcast_allocatable(this%phi)
  call proc(0)%broadcast_allocatable(this%x)
  call proc(0)%broadcast_allocatable(this%tedge)
  call proc(0)%broadcast_allocatable(this%pedge)
  call proc(0)%broadcast_allocatable(this%diag)
  call proc(0)%broadcast_allocatable(this%nvect)
  call proc(0)%broadcast_allocatable(this%area)
  call proc(0)%broadcast_allocatable(this%masked)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(trisurf), intent(inout) :: this


  deallocate (this%phi, this%x, this%tedge, this%pedge, this%diag, this%nvect, this%area, this%masked)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function rzslice(this, phi) result(x)
  !
  ! generate slice at phi [rad]
  !
  use moose_algorithms, only: binary_search_L
  class(trisurf), intent(in) :: this
  real(real64),   intent(in) :: phi
  real(real64)               :: x(2, 0:2*this%nv)

  real(real64) :: phi_mod, t, nvec(2)
  integer :: i, j


  ! toroidal bounds check
  if (phi < this%phi(0)  .or.  phi > this%phi(this%nu)) then
     x = 0.d0
     return
  endif


  ! find (i, t) for phi
  if (phi == this%phi(0)) then
     i = 0
  else
     i = binary_search_L(this%phi, phi) - 1
  endif
  t = (phi - this%phi(i)) / (this%phi(i+1) - this%phi(i))


  ! set normal vector for r-z plane at phi
  nvec = [-sin(phi), cos(phi)]


  do j=0,this%nv-1
     ! 1. toroidal edges
     x(:, 2*j) = slice_edge(this%x(:, j, i), this%tedge(:, 0, j, i))

     ! 2. diagonals
     x(:, 2*j+1) = slice_edge(this%x(:, j, i+1), this%diag(:, j, i))
  enddo
  x(:, 2*this%nv) = slice_edge(this%x(:, j, i), -this%tedge(:, 1, this%nv-1, i))

  contains
  !.............................................................................
  pure function slice_edge(p, v) result(rz)
  real(real64), intent(in) :: p(3), v(3)
  real(real64)             :: rz(2)

  real(real64) :: den, num, x(3)


  den = dot_product(v(1:2), nvec)
  num = dot_product(p(1:2), nvec)
  if (den == 0.d0) then
     if (num == 0.d0) then
        rz = [sqrt(p(1)**2 + p(2)**2), p(3)]
     else
        rz = 0.d0
     endif
  else
     x = p - num / den * v
     rz = [sqrt(x(1)**2 + x(2)**2), x(3)]
  endif

  end function slice_edge
  !.............................................................................
  end function rzslice
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function winding_number(this, p) result(wn)
  !
  ! compute winding number for p = (r, z, phi)
  !
  use moose_polygon2d, only: backend => winding_number
  class(trisurf), intent(in) :: this
  real(real64),   intent(in) :: p(3)
  integer                    :: wn


  wn = backend(this%rzslice(p(3)), p(1:2), .false.)

  end function winding_number
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function ray_intersect(this, p0, d, p0_on_surface, px, t, S)
  !
  ! check for intersection of ray: p0 + t * d with trisurf
  !
  ! if p0_on_surface, S is expected as input so that the intersection check with
  ! this particular triangle will be skipped
  !
  use moose_tmesh, only: ray_intersects_triangle
  class(trisurf),        intent(in   ) :: this
  real(real64),          intent(in   ) :: p0(3), d(3)
  logical,               intent(in   ) :: p0_on_surface
  real(real64),          intent(  out) :: px(3), t
  type(surface_coords),  intent(inout) :: S
  logical                              :: ray_intersect

  real(real64) :: tijk, u, v
  integer :: i, j, k, nskip(3)


  ! avoid spurious intersection if ray starts on surface
  nskip = -1
  if (p0_on_surface) nskip = [S%uindex, S%vindex, kindex(S%u, S%v)]


  ray_intersect = .false.
  t = huge(1.d0)
  do i=0,this%nu-1
  do j=0,this%nv-1
  do k=0,1
     if (all([i, j, k] == nskip)  .or.  this%masked(k, j, i)) cycle

     if (ray_intersects_triangle(p0, d, this%x(:,j+k,i+k), this%tedge(:,k,j,i), this%pedge(:,k,j,i), tijk, u, v)) then
        if (tijk < t) then
           ray_intersect = .true.
           t = tijk
           S = surface_coords(i, j, k, u, v)
           px = p0 + t * d
        endif
     endif
  enddo
  enddo
  enddo

  end function ray_intersect
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function interp(this, S) result(x)
  !
  ! evaluate Cartesian coordinates for *S*
  !
  class(trisurf),        intent(in) :: this
  class(surface_coords), intent(in) :: S
  real(real64)                      :: x(3)

  real(real64) :: u, v
  integer :: i, j, k


  i = S%uindex
  j = S%vindex
  k = kindex(S%u, S%v)


  if (k == 0) then
     x = this%x(:,j,i)
     u = S%u
     v = S%v
  else
     x = this%x(:,j+1,i+1)
     u = 1.d0 - S%u
     v = 1.d0 - S%v
  endif
  x = x + u * this%tedge(:,k,j,i) + v * this%pedge(:,k,j,i)

  end function interp
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function normal_vector(this, S)
  !
  ! return normal vector for surface point *S*
  !
  class(trisurf),        intent(in) :: this
  class(surface_coords), intent(in) :: S
  real(real64)                      :: normal_vector(3)

  integer :: k


  k = kindex(S%u, S%v)
  normal_vector = this%nvect(:, k, S%vindex, S%uindex)

  end function normal_vector
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function export_tmesh(this) result(T)
  !
  ! construct tmesh from trisurf
  !
  use moose_grids, only: tmesh
  class(trisurf), intent(in) :: this
  type(tmesh)                :: T

  integer, allocatable :: inode(:,:)
  integer :: i, j, k


  T = tmesh((this%nu + 1) * (this%nv + 1), this%nu * this%nv * 2, 3)
  allocate (inode(0:this%nv, 0:this%nu))


  ! define nodes
  k = 0
  do j=0,this%nv
  do i=0,this%nu
     inode(j, i) = k
     T%x(:, k) = this%x(:, j, i)
     k = k + 1
  enddo
  enddo


  ! define cells
  k = 0
  do j=0,this%nv-1
  do i=0,this%nu-1
     ! lower triangle
     if (.not.this%masked(0, j, i)) then
        T%triangles(:, k) = [inode(j, i), inode(j, i+1), inode(j+1, i)]
        k = k + 1
     endif

     ! upper triangle
     if (.not.this%masked(1, j, i)) then
        T%triangles(:, k) = [inode(j+1, i+1), inode(j+1, i), inode(j, i+1)]
        k = k + 1
     endif
  enddo
  enddo
  deallocate (inode)

  end function export_tmesh
  !-----------------------------------------------------------------------------

end module moose_trisurf
