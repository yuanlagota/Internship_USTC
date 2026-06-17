module moose_hypermesh3d
  use iso_fortran_env
  use moose_grids, only: block_structured, tpzmesh3d, r3grid
  implicit none
  private


  ! bin index map for mesh coordinates .........................................
  type, public :: bin_map
     ! number of bins in each macro-interval
     integer, allocatable :: nbins(:)

     ! index offset for each marco-interval
     integer, allocatable :: offset(:)

     contains
     procedure :: broadcast => broadcast_bin_map, collect, backward_transform
     procedure :: free => free_bin_map
  end type bin_map



  ! refined tpzmesh3d with index map for bins ..................................
  type, extends(tpzmesh3d), public :: refined_tpzmesh3d
     type(bin_map) :: umap(2)
     real(real64) :: dphi, ds
     integer :: checksum   ! checksum for original geometry

     contains
     procedure :: broadcast => broadcast_refined_tpzmesh3d
     procedure :: free => free_refined_tpzmesh3d
     procedure :: find_cell_index
     generic :: cell_index => find_cell_index
     procedure :: writenc
  end type refined_tpzmesh3d


  interface refined_tpzmesh3d
     procedure :: axisurf_tpzmesh3d, torosurf_tpzmesh3d
  end interface
  ! refined_tpzmesh3d ..........................................................



  ! fine mesh for census arrays on hypersurf3d .................................
  type, extends(r3grid), public :: hypermesh3d
     ! the block-structured domain
     type(block_structured), pointer :: block_structured

     ! the array of tpzmesh3d blocks
     type(refined_tpzmesh3d), pointer :: refined_tpzmesh3d(:)

     contains
     procedure :: broadcast => broadcast_hypermesh3d
     procedure :: free => free_hypermesh3d
     procedure :: io_metadata
     procedure :: ncells, cell_offset, cell_index, hypersurf3d_transform, area
  end type hypermesh3d


  interface hypermesh3d
     procedure :: init_rank1_resolution, init_scalar_resolution
  end interface
  ! hypermesh3d ................................................................



  public :: &
     tpzmesh3d_cell_area

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function refined_bins(nbins) result(this)
  !
  ! nbins: number of bins for each macro-interval
  !
  integer, intent(in) :: nbins(:)
  type(bin_map)       :: this

  integer :: i


  allocate (this%nbins(0:size(nbins)-1), source = nbins)
  allocate (this%offset(0:size(nbins)), source = 0)
  do i=1,size(nbins)
     this%offset(i) = this%offset(i-1) + nbins(i)
  enddo

  end function refined_bins
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function init_rank1_resolution(H, dphi, ds) result(this)
  !
  ! create refined mesh from hypersurf
  !
  use moose_units, only: METER, DEGREE
  use moose_r3grid, only: cylindrical_r3grid
  use moose_hypersurface, only: hypersurf3d
  use moose_axisurf
  use moose_torosurf
  class(hypersurf3d), intent(in) :: H
  real(real64),       intent(in) :: dphi(H%nsurfaces), ds(H%nsurfaces)
  type(hypermesh3d)              :: this

  integer :: i


  allocate (this%block_structured, source = block_structured(H%nsurfaces, 2, 3))
  allocate (this%refined_tpzmesh3d(H%nsurfaces))
  do i=1,H%nsurfaces
     select type(G => H%surfaces(i)%geometry)
     type is (axisurf)
        this%refined_tpzmesh3d(i) = axisurf_tpzmesh3d(G, dphi(i), ds(i))
     type is (torosurf)
        this%refined_tpzmesh3d(i) = torosurf_tpzmesh3d(G, dphi(i), ds(i))
     end select
     this%block_structured%blocks(i)%grid => this%refined_tpzmesh3d(i)
     this%block_structured%blocks(i)%key = H%surfaces(i)%key
  enddo
  call this%block_structured%init_index_offsets()
  this%r3grid = cylindrical_r3grid(this%block_structured, METER, DEGREE, 1, 2, 3, .false.)

  end function init_rank1_resolution
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function init_scalar_resolution(H, dphi, ds) result(this)
  !
  ! the same dphi and ds is applied to all surface patches
  !
  use moose_hypersurface, only: hypersurf3d
  class(hypersurf3d), intent(in) :: H
  real(real64),       intent(in) :: dphi, ds
  type(hypermesh3d)              :: this

  real(real64) :: dphi_n(H%nsurfaces), ds_n(H%nsurfaces)


  dphi_n = dphi
  ds_n = ds
  this = hypermesh3d(H, dphi_n, ds_n)

  end function init_scalar_resolution
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function axisurf_tpzmesh3d(A, dphi, ds) result(this)
  use moose_math, only: diff, linspace, pi, zero_cumsum
  use moose_axisurf, only: axisurf
  class(axisurf), intent(in) :: A
  real(real64),   intent(in) :: dphi, ds
  type(refined_tpzmesh3d)    :: this

  real(real64), allocatable :: v(:), d(:)
  real(real64) :: p1(2), p2(2)
  integer :: i, k, nphi, nv_sub


  ! toroidal resolution
  nphi = nint(A%dphi / dphi * 180 / pi)
  this%umap(1) = refined_bins([nphi])
  this%dphi = dphi
  this%ds = ds
  this%checksum = A%checksum()


  ! resolution for v-direction
  if (ds > 0.d0) then
     allocate (d, source = norm2(diff(A%P%nodes(), 1), 2))
     this%umap(2) = refined_bins(max(nint(d / ds), 1))
  else
     this%umap(2) = refined_bins([(1, i=1,A%nv)])
  endif
  if (A%vdef) then
     allocate (v(1:A%nv+1), source = A%v)
  else
     allocate (v, source = A%vfallback())
  endif


  ! construct refined mesh
  this%tpzmesh3d = tpzmesh3d(nphi+1, sum(this%umap(2)%nbins)+1, &
     "Toroidal angle [deg]", A%vlabel(), "r ["//A%units()//"]", "z ["//A%units()//"]", A%description())
  this%domain%u = linspace(A%phi0, A%phi0 + A%dphi, nphi+1) * 180 / pi
  i = 0
  p1 = A%P%node(0)
  do k=1,size(this%umap(2)%nbins)
     nv_sub = this%umap(2)%nbins(k-1)
     p2 = A%P%node(k)
     this%domain%v(:,i:i+nv_sub) = transpose(spread(linspace(v(k), v(k+1), nv_sub+1), 2, nphi+1))
     this%x(1,:,i:i+nv_sub) = transpose(spread(linspace(p1(1), p2(1), nv_sub+1), 2, nphi+1))
     this%x(2,:,i:i+nv_sub) = transpose(spread(linspace(p1(2), p2(2), nv_sub+1), 2, nphi+1))

     i = i + nv_sub
     p1 = p2
  enddo

  end function axisurf_tpzmesh3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function torosurf_tpzmesh3d(T, dphi, ds) result(this)
  use moose_math, only: diff, linspace, bilinspace, pi
  use moose_torosurf, only: torosurf
  class(torosurf), intent(in) :: T
  real(real64),    intent(in) :: dphi, ds
  type(refined_tpzmesh3d)     :: this

  real(real64), allocatable :: d(:,:)
  integer :: i, i1, i2, j, j1, j2, nu, nu_sub, nv, nv_sub


  this%checksum = T%checksum()
  nu = T%nu
  nv = T%nv
  associate (rz => T%rz, v => T%v)
  allocate (d(nv, nu))
  d = norm2(rz(:,1:nv,0:nu-1) + rz(:,1:nv,1:nu) - rz(:,0:nv-1,0:nu-1) - rz(:,0:nv-1,1:nu), 1) / 2


  ! toroidal resolution
  if (dphi > 0.d0) then
     this%umap(1) = refined_bins(max(nint(diff(T%phi) / dphi / pi * 180.d0), 1))
  else
     this%umap(1) = refined_bins([(1, i=1,nu)])
  endif
  this%dphi = dphi


  ! poloidal resolution
  if (ds > 0.d0) then
     this%umap(2) = refined_bins(max(nint(sum(d, 2) / nu / ds), 1))
  else
     this%umap(2) = refined_bins([(1, j=1,nv)])
  endif
  this%ds = ds


  ! construct refined mesh
  this%tpzmesh3d = tpzmesh3d(sum(this%umap(1)%nbins)+1, sum(this%umap(2)%nbins)+1, &
     "Toroidal angle [deg]", T%vlabel(), "r ["//T%units()//"]", "z ["//T%units()//"]", T%description())
  i1 = 0
  do i=0,size(this%umap(1)%nbins)-1
     nu_sub = this%umap(1)%nbins(i)
     i2 = i1 + nu_sub
     this%domain%u(i1:i2) = linspace(T%phi(i), T%phi(i+1), nu_sub+1)

     j1 = 0
     do j=0,size(this%umap(2)%nbins)-1
        nv_sub = this%umap(2)%nbins(j)
        j2 = j1 + nv_sub
        this%domain%v(i1:i2,j1:j2) = transpose(bilinspace(v(j:j+1,i:i+1), nv_sub+1, nu_sub+1))
        this%x(1,i1:i2,j1:j2) = transpose(bilinspace(rz(1,j:j+1,i:i+1), nv_sub+1, nu_sub+1))
        this%x(2,i1:i2,j1:j2) = transpose(bilinspace(rz(2,j:j+1,i:i+1), nv_sub+1, nu_sub+1))
        j1 = j2
     enddo
     i1 = i2
  enddo
  end associate
  this%domain%u  = this%domain%u / pi * 180.d0

  end function torosurf_tpzmesh3d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  pure function collect(this, u) result(bin_index)
  !
  ! find bin index for data point *u*
  !
  class(bin_map), intent(in) :: this
  real(real64),   intent(in) :: u
  integer                    :: bin_index

  real(real64) :: r
  integer :: i


  bin_index = -1
  if (u < 0.d0  .or.  u > size(this%nbins)) return


  i = int(u)
  if (i == size(this%nbins)) then
     bin_index = this%offset(size(this%nbins)) - 1
  else
     r = u - i
     bin_index = this%offset(i) + int(r * this%nbins(i))
  endif

  end function collect
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine backward_transform(this, uindex, u)
  !
  ! transform (uindex, u) from refined bins to macro-interval
  !
  use moose_algorithms, only: binary_search_R
  class(bin_map), intent(in   ) :: this
  integer,        intent(inout) :: uindex
  real(real64),   intent(inout) :: u

  integer :: k


  k = binary_search_R(this%offset, uindex) - 1
  u = (uindex - this%offset(k) + u) / this%nbins(k)
  uindex = k

  end subroutine backward_transform
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine broadcast_bin_map(this)
  use moose_mpi
  class(bin_map), intent(inout) :: this


  call proc(0)%broadcast_allocatable(this%nbins)
  call proc(0)%broadcast_allocatable(this%offset)

  end subroutine broadcast_bin_map
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_bin_map(this)
  class(bin_map), intent(inout) :: this


  deallocate (this%nbins, this%offset)

  end subroutine free_bin_map
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function find_cell_index(this, u) result(cell_index)
  class(refined_tpzmesh3d), intent(in) :: this
  real(real64),             intent(in) :: u(2)
  integer                              :: cell_index

  integer :: k(2)


  k(1) = this%umap(1)%collect(u(1))
  k(2) = this%umap(2)%collect(u(2))
  if (k(1) == -1  .or.  k(2) == -1) then
     cell_index = -1
     return
  endif
  cell_index = this%cell_index(k)

  end function find_cell_index
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine broadcast_refined_tpzmesh3d(this)
  use moose_mpi
  class(refined_tpzmesh3d), intent(inout) :: this


  call this%tpzmesh3d%broadcast()
  call this%umap(1)%broadcast()
  call this%umap(2)%broadcast()

  end subroutine broadcast_refined_tpzmesh3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_refined_tpzmesh3d(this)
  class(refined_tpzmesh3d), intent(inout) :: this


  call this%tpzmesh3d%free()
  call this%umap(1)%free()
  call this%umap(2)%free()

  end subroutine free_refined_tpzmesh3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  use moose_utils
  class(refined_tpzmesh3d), intent(in) :: this
  type(netcdf_dataset),     intent(in) :: N


  call N%put_att("dphi", this%dphi)
  call N%put_att("ds", this%ds)
  call N%put_att("checksum", this%checksum)
  call this%tpzmesh3d%writenc(N)

  end subroutine writenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine broadcast_hypermesh3d(this)
  use moose_mpi
  class(hypermesh3d), intent(inout) :: this

  integer :: i, n


  if (rank > 0) allocate (this%block_structured)
  this%domain => this%block_structured
  call this%r3grid%broadcast()
  n = size(this%block_structured%blocks)

  if (rank > 0) allocate (this%refined_tpzmesh3d(n))
  do i=1,n
     call this%refined_tpzmesh3d(i)%broadcast()
     this%block_structured%blocks(i)%grid => this%refined_tpzmesh3d(i)
  enddo

  end subroutine broadcast_hypermesh3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_hypermesh3d(this)
  class(hypermesh3d), intent(inout) :: this

  integer :: i


  do i=1,size(this%refined_tpzmesh3d)
     call this%refined_tpzmesh3d(i)%free()
  enddo
  deallocate (this%block_structured, this%refined_tpzmesh3d)
  call this%r3grid%free()

  end subroutine free_hypermesh3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function io_metadata(this)
  use moose_dict
  class(hypermesh3d), intent(in) :: this
  type(dict)                     :: io_metadata


  io_metadata = this%block_structured%io_metadata()

  end function io_metadata
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function ncells(this)
  class(hypermesh3d), intent(in) :: this
  integer                        :: ncells


  ncells = this%block_structured%ncells()

  end function ncells
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function cell_offset(this, n)
  class(hypermesh3d), intent(in) :: this
  integer,            intent(in) :: n
  integer                        :: cell_offset


  cell_offset = this%block_structured%cell_offset(n)

  end function cell_offset
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function cell_index(this, n, u)
  !
  ! index for cell that corresponds to surface coordinates (n, u)
  !
  class(hypermesh3d), intent(in) :: this
  integer,            intent(in) :: n
  real(real64),       intent(in) :: u(2)
  integer                        :: cell_index


  cell_index = -1
  if (n < lbound(this%refined_tpzmesh3d,1)  .or.  n > ubound(this%refined_tpzmesh3d,1)) return

  cell_index = this%refined_tpzmesh3d(n)%cell_index(u)
  if (cell_index == -1) return
  cell_index = cell_index + this%cell_offset(n)

  end function cell_index
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function hypersurf3d_transform(this, linear_index, u) result(x)
  !
  ! transform coordinates from refined mesh to hypersurf3d
  !
  use moose_hypersurface, only: hypersurf3d_coords
  class(hypermesh3d), intent(in) :: this
  integer,            intent(in) :: linear_index
  real(real64),       intent(in) :: u(2)
  type(hypersurf3d_coords)       :: x

  integer :: k(0:2)


  k = this%block_structured%cell_index(linear_index)
  x%surface_index = k(0)
  x%uindex = k(1)
  x%vindex = k(2)
  x%u = u(1)
  x%v = u(2)
  call this%refined_tpzmesh3d(x%surface_index)%umap(1)%backward_transform(x%uindex, x%u)
  call this%refined_tpzmesh3d(x%surface_index)%umap(2)%backward_transform(x%vindex, x%v)

  end function hypersurf3d_transform
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function area(this)
  use moose_math, only: pi
  class(hypermesh3d), intent(in) :: this
  real(real64)                   :: area(0:this%ncells()-1)

  integer :: i1, i2, k


  do k=1,size(this%refined_tpzmesh3d)
     i1 = this%cell_offset(k)
     i2 = this%cell_offset(k+1) - 1
     area(i1:i2) = pack(tpzmesh3d_cell_area(this%refined_tpzmesh3d(k)), .true.)
  enddo

  end function area
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  pure function tpzmesh3d_cell_area(this) result(a)
  !
  ! this procedure computes the (approximate!) area of the tpzmesh3d cells for
  ! the special case x1 = r, x2 = z, u = phi
  !
  use moose_math, only: pi, cross_product
  class(tpzmesh3d), intent(in   ) :: this
  real(real64)                    :: a(size(this%x,2)-1, size(this%x,3)-1)

  real(real64), allocatable :: x(:,:,:)
  integer :: i, j, nu, nv


  nu = size(this%x,2) - 1
  nv = size(this%x,3) - 1
  allocate (x(3, 0:nu, 0:nv))
  do i=0,nu
     x(1,i,:) = cos(this%domain%u(i) / 180 * pi) * this%x(1,i,:)
     x(2,i,:) = sin(this%domain%u(i) / 180 * pi) * this%x(1,i,:)
     x(3,i,:) = this%x(2,i,:)
  enddo


  do i=1,nu
  do j=1,nv
     a(i,j) = (norm2(cross_product(x(:,i,j-1)-x(:,i-1,j-1), x(:,i-1,j)-x(:,i-1,j-1))) &
             + norm2(cross_product(x(:,i,j-1)-x(:,i  ,j  ), x(:,i-1,j)-x(:,i  ,j  )))) / 2

  enddo
  enddo

  end function tpzmesh3d_cell_area
  !-----------------------------------------------------------------------------

end module moose_hypermesh3d
