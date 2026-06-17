! The geometry of the problem is specified here. This includes:
!
!    - the toroidal symmetry
!
!    - the edge-core interface
!
!    - the first wall, divertor targets and pumping surfaces
!
module firefly_geometry
  use iso_fortran_env
  use moose_geometry
  use flare_mmesh_unstructured_mmesh
  use flare_mmesh_unstructured_tracing
  implicit none


  ! toroidal resolution [deg] to use for triangulation of axisurf
  real(real64) :: axisurf_resolution = 1.d0

  ! poloidal increment for triangular approximation of core boundary
  integer :: core_boundary_incr = 4


  ! magnetic mesh
  type(mmesh), target :: mesh
  type(tracing_workspace) :: tracing
  logical :: mmesh_initialized = .false.

  ! toroidal symmetry of domain
  integer :: symmetry
  logical :: half_period
  real(real64) :: phi_bounds(2)   ! lower and upper toroidal bound [rad]
  real(real64) :: nvec_bounds(2,2)   ! R,Z components of normal vectors at toroidal bounds


  ! plasma facing components (incl. pumping surfaces)
  type(hypersurf3d), target :: pfc
  type(polygon2d) :: rzslice_casing(2)   ! R-Z slice of casing at lower and upper toroidal boundary
  integer, allocatable :: material_index(:), plasma_side(:)
  integer :: icasing, nsurfaces
  logical :: pfc_initialized = .false.

  ! triangulation of axisurf/torosurf geometry for tracing of neutral particles (intersection checks)
  type(trisurf), allocatable :: triangulated_surfaces(:)

  ! volume in core and edge/SOL reservoirs
  real(real64) :: volume(2)


  ! approximation of core boundary for tracing of neutral particles and divertor geometry validation
  type(trisurf) :: core_boundary
  type(torosurf) :: divertor_exclusion_boundary

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init_workspace(filename, lcfs, offset, seed)
  !
  ! Load magnetic mesh for fast reconstruction of field lines and initialize domain.
  !
  ! Optional:
  !    - user defined LCFS (default: inner mesh boundary)
  !    - offset from LCFS for divertor geometry validation [m]
  !    - seed for random number generator (used for field line diffusion)
  !
  use moose_error
  use moose_mpi
  use moose_math, only: pi2
  use moose_utils, only: user_option, random_seed1
  use moose_geometry, only: loadnc_torosurf
  character(len=*), intent(in) :: filename
  character(len=*), intent(in), optional :: lcfs
  real(real64),     intent(in), optional :: offset
  integer,          intent(in), optional :: seed

  real(real64) :: phiU
  integer :: ierr


  ! clean up first, if necessary
  if (pfc_initialized) then
     call tracing%free()
     call pfc%free()
     deallocate (triangulated_surfaces, material_index, plasma_side)
  endif
  if (mmesh_initialized) then
     call mesh%free()
  endif


  ! load magnetic mesh
  if (rank == 0) then
     print *, "Loading mmesh ..."
     mesh = loadnc_mmesh(filename)
     mesh%phi = mesh%phi / 360.d0 * pi2
     print *, "... done"
     print *
  endif
  if (nproc > 1) call mesh%broadcast()
  mmesh_initialized = .true.


  ! set symmetry of toroidal domain and related parameters
  symmetry = abs(mesh%symmetry)
  half_period = mesh%symmetry < 0

  phiU = pi2 / symmetry
  if (half_period) phiU = phiU / 2

  phi_bounds(1) = 0.d0
  phi_bounds(2) = phiU
  nvec_bounds(:,1) = [0.d0, 1.d0]
  nvec_bounds(:,2) = [-sin(phiU), cos(phiU)]


  ! set core boundary for tracing of neutral particles and divertor geometry validation
  if (user_option("", lcfs) /= "") then
     divertor_exclusion_boundary = loadnc_torosurf(lcfs, convert_units="m")
  else
     divertor_exclusion_boundary = mesh%core_boundary(core_boundary_incr)
  endif
  core_boundary = trisurf(divertor_exclusion_boundary)
  volume(1) = trisurf_volume(core_boundary)
  if (user_option(0.d0, offset) /= 0.d0) then
     call divertor_exclusion_boundary%lhshift(-offset, ierr)
     if (ierr /= 0) call ERROR("applying offset to LCFS failed", "init_workspace", ierr)
  endif


  ! initialize random seed
  if (user_option(0, seed) /= 0) call random_seed1(seed * (rank + 1))

  end subroutine init_workspace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_pfc(filename)
  !
  ! Load boundary geometry (first wall, divertor targets and baffles) for
  ! interception of magnetic field lines. This should be called after the
  ! workspace (magnetic mesh) has been initialized, and can be called again
  ! later to update the geometry.
  !
  use moose_error
  use moose_mpi
  use moose_geometry, only: axisurf, torosurf
  use flare_control, only: report
  character(len=*), intent(in) :: filename

  type(torosurf) :: T
  integer :: i, nsym, nphi


  ! mmesh must be initialized first
  if (.not.mmesh_initialized) call ERROR("magnetic mesh is not initialzied")


  ! clean up, if necessary
  if (pfc_initialized) then
     call tracing%free()
     call pfc%free()
     deallocate (triangulated_surfaces, material_index, plasma_side)
  endif


  ! load boundary geometry and initialize tracing workspace
  icasing = 1   ! TODO: verify that this is the correct index
  if (rank == 0) then
     if (report) print *, "Loading plasma facing components ..."
     pfc = loadnc_hypersurf3d(filename, convert_units="m")

     ! set phi0, dphi for axisurf from mesh
     do i=1,pfc%nsurfaces
        select type (S => pfc%surfaces(i)%geometry)
        type is (axisurf)
           S%phi0 = mesh%phi(0)
           S%dphi = mesh%phi(mesh%nphi-1) -  mesh%phi(0)
           S%symmetry = mesh%symmetry
        end select
     enddo
     if (report) then
        print *, "... done"
        print *
     endif
  endif
  call pfc%broadcast()
  tracing = tracing_workspace(mesh, pfc, report)
  nsurfaces = pfc%nsurfaces
  pfc_initialized = .true.


  ! set up triangulation of boundary for tracing of neutral particls
  allocate (triangulated_surfaces(nsurfaces), material_index(nsurfaces), plasma_side(nsurfaces))
  do i=1,nsurfaces
     select type(S => pfc%surfaces(i)%geometry)
     class is (torosurf)
        triangulated_surfaces(i) = trisurf(S)
        if (i == icasing) call set_rzslice_casing(S)

     class is (axisurf)
        nsym = symmetry;   if (half_period) nsym = 2 * nsym
        nphi = int(360.d0 / axisurf_resolution / nsym)
        T = torosurf(S, nsym, nphi)
        triangulated_surfaces(i) = trisurf(T)
        if (i == icasing) call set_rzslice_casing(T)
        call T%free()

     end select

     material_index(i) = pfc%surfaces(i)%geometry%metadata%getint("material_index", 2)
     plasma_side(i) = pfc%surfaces(i)%geometry%metadata%getint("plasma_side", 1)
     if (material_index(i) <= 0) call ERROR("material_index > 0 required for surface " // pfc%surfaces(i)%key)
     if (abs(plasma_side(i)) /= 1) call ERROR("abs(plasma_side) = 1 required for surface " // pfc%surfaces(i)%key)
  enddo


  ! compute volume in edge/SOL reservoir
  volume(2) = trisurf_volume(triangulated_surfaces(icasing)) - volume(1)

  contains
  !.............................................................................
  subroutine set_rzslice_casing(C)
  type(torosurf), intent(in) :: C


  rzslice_casing(1) = polygon2d(C%rz(:,:,0))
  rzslice_casing(2) = polygon2d(C%rz(:,:,C%nu))

  end subroutine set_rzslice_casing
  !.............................................................................
  end subroutine set_pfc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine assert_mmesh_and_pfc(procedure_name)
  use moose_error
  character(len=*), intent(in) :: procedure_name


  if (.not.mmesh_initialized) call ERROR("magnetic mesh is not initialized", procedure_name)
  if (.not.pfc_initialized) call ERROR("plasma facing components not initialized", procedure_name)

  end subroutine assert_mmesh_and_pfc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function validate_torosurf(T)
  !
  ! Confirm that *T* does not intersect core boundary.
  !
  use moose_error
  use moose_mpi
  use moose_algorithms, only: xsegments
  use moose_math,       only: pi
  use moose_polygon2d,  only: winding_number
  type(torosurf), intent(in) :: T
  logical                    :: validate_torosurf

  real(real64) :: p1(3), p2(3), x(3), s, u(2)
  integer :: i1, i2, j1, j2, stride, i2range(0:1), parallel_scan(0:nproc-1)


  ! 0. check input
  ! mmesh & lcfs must be initialized first
  if (.not.mmesh_initialized) call ERROR("magnetic mesh is not initialzied", "validate_torosurf")
  associate (lcfs => divertor_exclusion_boundary)

  ! check toroidal direction of torosurf
  stride = 1;   if (T%phi(T%nu) < T%phi(0)) stride = -1
  ! set i2-range
  i2range = [0, T%nu];   if (stride == -1) i2range = [T%nu, 0]
  ! verify toroidal positions
  if (T%nu > lcfs%nu) call ERROR("invalid toroidal resolution")
  if (any(abs(T%phi - lcfs%phi(i2range(0):i2range(1):stride)) > 1.d-12)) then
     do i1=0,T%nu
        print *, i1, [T%phi(i1), lcfs%phi(i2range(0)+stride*i1)] * 180 / pi
     enddo
     call ERROR("mismatching toroidal positions")
  endif


  ! 1. scan for intersection between poloidal segments
  parallel_scan = 0
  scan1: do i1=0,T%nu
     i2 = i2range(0) + stride * i1
     do j1=0,T%nv-1
     do j2=rank,lcfs%nv-1,nproc
        if (xsegments(T%rz(:,j1,i1), T%rz(:,j1+1,i1), lcfs%rz(:,j2,i2), lcfs%rz(:,j2+1,i2))) then
           parallel_scan(rank) = 1
           exit scan1
        endif
     enddo
     enddo
  enddo scan1
  call moose_mpi_sum(parallel_scan)
  validate_torosurf = all(parallel_scan == 0)
  if (.not.validate_torosurf) return


  ! 2. scan for intersection along toroidal direction
  parallel_scan = 0
  scan2: do i1=rank,T%nu-1,nproc
     do j1=0,T%nv
        p1 = [T%rz(:,j1,i1), T%phi(i1)]
        p2 = [T%rz(:,j1,i1+1), T%phi(i1+1)]
        if (lcfs%intersect(p1, p2, x, s, u)) then
           parallel_scan(rank) = 1
           exit scan2
        endif
     enddo
  enddo scan2
  call moose_mpi_sum(parallel_scan)
  validate_torosurf = all(parallel_scan == 0)
  if (.not.validate_torosurf) return


  ! 3. T must not be completely inside LCFS
  validate_torosurf = lcfs%winding_number(T%rz(:,0,0), i2range(0)) == 0
  end associate

  end function validate_torosurf
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function get_cell_number(x) result(icell)
  !
  ! identify cell index for a point given by its Cartesian coordinates [m].
  !
  use moose_math, only: cart_to_cyl
  real(real64), intent(in) :: x(3)
  integer                  :: icell

  integer :: wn


  icell = 1
  wn = core_boundary%winding_number(cart_to_cyl(x))
  if (wn == 0) icell = 2

  end function get_cell_number
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function trisurf_point(S) result(x)
  !
  ! evaluate Cartesian coordinates for point *S* on triangulated surface
  !
  use moose_hypersurface, only: hypersurf3d_coords
  type(hypersurf3d_coords), intent(in) :: S
  real(real64)                         :: x(3)


  x = triangulated_surfaces(S%surface_index)%interp(S)

  end function trisurf_point
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function trisurf_volume(this) result(V)
  class(trisurf), intent(in) :: this
  real(real64)               :: V

  real(real64) :: area, x0(3)
  integer :: i, j


  V = 0
  do i=0,this%nu-1
  do j=0,this%nv-1
     ! lower triangle
     x0 = this%x(:, j, i) + 0.5d0 * (this%tedge(:, 0, j, i) + this%pedge(:, 0, j, i))
     V = V + sum(x0 * this%nvect(:, 0, j, i)) * this%area(0, j, i)

     ! upper triangle
     x0 = this%x(:, j+1, i+1) + 0.5d0 * (this%tedge(:, 1, j, i) + this%pedge(:, 1, j, i))
     V = V + sum(x0 * this%nvect(:, 1, j, i)) * this%area(1, j, i)
  enddo
  enddo
  V = V / 3

  end function trisurf_volume
  !-----------------------------------------------------------------------------

end module firefly_geometry
