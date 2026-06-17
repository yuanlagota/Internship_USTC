module tools
  use kinds
  implicit none


  ! workspace for qmesh generators
  real(real64), allocatable :: qmesh_x(:,:,:)

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine bspline3d_slice(filename, idim1, x1, n2, n3, x2, x3, values)
  use moose_analysis, only: bspline3d, loadnc_bspline3d
  character(len=*), intent(in   ) :: filename
  integer,          intent(in   ) :: idim1, n2, n3
  real(real64),     intent(in   ) :: x1
  real(real64),     intent(  out) :: x2(n2), x3(n3), values(n2, n3)

  type(bspline3d) :: bspline
  real(real64) :: x(3), d2, d3
  integer :: i, idim2, idim3, j


  bspline = loadnc_bspline3d(filename)
  idim2 = mod(idim1, 3) + 1
  idim3 = mod(idim2, 3) + 1
  d2 = bspline%ub(idim2) - bspline%lb(idim2)
  d3 = bspline%ub(idim3) - bspline%lb(idim3)

  x(idim1) = x1
  do i=1,n2
     x(idim2) = bspline%lb(idim2) + d2 * (i - 1) / (n2 - 1)
     x2(i) = x(idim2)
     do j=1,n3
        x(idim3) = bspline%lb(idim3) + d3 * (j - 1) / (n3 - 1)
        x3(j) = x(idim3)
        values(i,j) = bspline%eval(x)
     enddo
  enddo

  end subroutine bspline3d_slice
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine quasi_orthogonal_qmesh(filename1, filename2, nu, nv, idir)
  use moose_grids,           only: qmesh
  use moose_geometry,        only: curve, loadtxt_curve
  use moose_qmesh_generator, only: make_qmesh => quasi_orthogonal_qmesh
  character(len=*), intent(in) :: filename1, filename2
  integer,          intent(in) :: nu, nv, idir

  class(curve), allocatable :: C1, C2
  type(qmesh) :: this


  if (allocated(qmesh_x)) deallocate(qmesh_x)
  allocate (C1, source=loadtxt_curve(filename1))
  allocate (C2, source=loadtxt_curve(filename2))
  this = make_qmesh(C1, C2, nu, nv, idir)
  allocate (qmesh_x, source=this%x)

  end subroutine quasi_orthogonal_qmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine nray_blocks_qmesh(filename1, filename2, nu, nv, nblocks)
  use moose_grids,           only: qmesh
  use moose_geometry,        only: curve, loadtxt_curve
  use moose_qmesh_generator, only: make_qmesh => nray_blocks_qmesh
  character(len=*), intent(in) :: filename1, filename2
  integer,          intent(in) :: nu, nv, nblocks

  class(curve), allocatable :: C1, C2
  type(qmesh) :: this


  if (allocated(qmesh_x)) deallocate(qmesh_x)
  allocate (C1, source=loadtxt_curve(filename1))
  allocate (C2, source=loadtxt_curve(filename2))
  this = make_qmesh(C1, C2, nu, nv, nblocks)
  allocate (qmesh_x, source=this%x)

  end subroutine nray_blocks_qmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine footpoints_qmesh(filename1, filename2, nu, nv, nblocks)
  use moose_grids,           only: qmesh
  use moose_geometry,        only: curve, loadtxt_curve
  use moose_qmesh_generator, only: make_qmesh => footpoints_qmesh
  character(len=*), intent(in) :: filename1, filename2
  integer,          intent(in) :: nu, nv, nblocks

  class(curve), allocatable :: C1, C2
  type(qmesh) :: this


  if (allocated(qmesh_x)) deallocate(qmesh_x)
  allocate (C1, source=loadtxt_curve(filename1))
  allocate (C2, source=loadtxt_curve(filename2))
  this = make_qmesh(C1, C2, nu, nv, nblocks)
  allocate (qmesh_x, source=this%x)

  end subroutine footpoints_qmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_qmesh_distance_contour_generator(u, v, i0, i1, P0, P1)
  use moose_grids, only: qmesh
  use moose_geometry, only: polygon2d
  use moose_qmesh_generator, only: aux_qmesh_distance_contour_generator_poly2d
  real(real64), intent(inout) :: u(:,:), v(:,:)
  integer,      intent(in   ) :: i0, i1
  real(real64), intent(in   ) :: P0(:,:), P1(:,:)

  type(qmesh) :: mesh
  type(polygon2d) :: P(0:1)


  mesh = qmesh(u, v)
  P(0) = polygon2d(P0)
  P(1) = polygon2d(P1)
  call aux_qmesh_distance_contour_generator_poly2d(mesh, i0, i1, P)
  u = mesh%u
  v = mesh%v

  end subroutine aux_qmesh_distance_contour_generator
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine trisurf_rzslice(torosurf_phi, torosurf_v, torosurf_rz, symmetry, phi, rzslice)
  use moose_math, only: pi
  use moose_geometry, only: torosurf, trisurf
  real(real64), intent(in   ) :: torosurf_phi(:), torosurf_v(:,:), torosurf_rz(:,:,:), phi
  integer,      intent(in   ) :: symmetry
  real(real64), intent(  out) :: rzslice(2, 0:2*(size(torosurf_rz, 2)-1))

  type(torosurf) :: this
  type(trisurf) :: T
  integer :: nu, nv, phi_order


  nu = size(torosurf_phi) - 1
  nv = size(torosurf_rz, 2) - 1
  phi_order = 1;   if (torosurf_phi(nu+1) < torosurf_phi(1)) phi_order = -1
  this = torosurf(nu, nv, symmetry, phi_order)
  this%phi = torosurf_phi / 180 * pi
  this%v = torosurf_v
  this%rz = torosurf_rz
  if (phi_order == -1) call this%reverse_phi_order()


  T = trisurf(this)
  rzslice = T%rzslice(phi / 180 * pi)

  end subroutine trisurf_rzslice
  !-----------------------------------------------------------------------------

end module tools
