module mmesh
  use kinds
  use flare_mmesh_unstructured_mmesh, unstructured_mmesh => mmesh
  implicit none


  ! workspace for unstructured mmesh
  type(unstructured_mmesh), private :: unstructured_mmesh_workspace
  logical :: unstructured_mmesh_workspace_allocated = .false.

  ! workspace rzmesh
  type(rzmesh), private :: rzmesh_workspace
  logical :: rzmesh_workspace_allocated = .false.

  ! workspace for rzbuffer
  real(real64), allocatable :: rzbuffer_x(:,:)

  contains
  !-----------------------------------------------------------------------------

include "_mmesh.f90"


! unstructured mmesh:
  !-----------------------------------------------------------------------------
  subroutine import_unstructured_mmesh(symmetry, nphi, phi, nzones, iphi_zone, nnodes, x, g, b, &
     nlines, izone_line, nbsect, inode_offset_nbsect, inode_offset_nlines, bsect, &
     ntubes, corner, next_tube, izone_tube, &
     rparam_tmap, iparam_tmap, nxmaps, iparam_xmap)
  integer,      intent(in) :: symmetry, nphi, nzones, nnodes, nlines, ntubes, nbsect, nxmaps
  real(real64), intent(in) :: phi(nphi)
  integer,      intent(in) :: iphi_zone(2, nzones)
  real(real64), intent(in) :: x(2, nnodes), g(2, nnodes), b(nnodes)
  integer,      intent(in) :: izone_line(nlines)
  integer,      intent(in) :: inode_offset_nbsect(nbsect)
  integer,      intent(in) :: inode_offset_nlines(nlines)
  integer,      intent(in) :: bsect(2, nbsect)
  integer,      intent(in) :: corner(4, ntubes)
  integer,      intent(in) :: next_tube(2, 4, ntubes)
  integer,      intent(in) :: izone_tube(ntubes)
  real(real64), intent(in) :: rparam_tmap(16, 2, ntubes)
  integer,      intent(in) :: iparam_tmap(2, 2, ntubes)
  integer,      intent(in) :: iparam_xmap(2, nxmaps)


  if (unstructured_mmesh_workspace_allocated) call free_unstructured_mmesh()
  associate (this => unstructured_mmesh_workspace)

  ! allocate arrays
  this = new_mmesh(symmetry, nphi, nzones, nnodes, nlines, ntubes, nbsect, nxmaps)
  unstructured_mmesh_workspace_allocated = .true.

  ! toroidal domain
  this%phi = phi
  this%iphi_zone = iphi_zone

  ! field lines
  this%x = x
  this%g = g
  this%b = b
  this%izone_line = izone_line
  this%inode_offset(-nbsect:-1) = inode_offset_nbsect
  this%inode_offset(0:nlines) = inode_offset_nlines
  if (nbsect > 0) this%bsect = bsect

  ! flux tubes
  this%corner = corner
  this%next_tube = next_tube
  this%izone_tube = izone_tube
  this%rparam_tmap = rparam_tmap
  this%iparam_tmap = iparam_tmap
  if (nxmaps > 0) this%iparam_xmap = iparam_xmap
  end associate

  end subroutine import_unstructured_mmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_unstructured_mmesh()


  call unstructured_mmesh_workspace%free()
  unstructured_mmesh_workspace_allocated = .false.

  end subroutine free_unstructured_mmesh
  !-----------------------------------------------------------------------------


! rzmesh:
  !-----------------------------------------------------------------------------
  subroutine mmesh_rzmesh(iphi, iside, nnodes, ncells, nbsect)
  integer, intent(in   ) :: iphi, iside
  integer, intent(  out) :: nnodes, ncells, nbsect


  if (rzmesh_workspace_allocated) call free_rzmesh()
  rzmesh_workspace = unstructured_mmesh_workspace%rzmesh(iphi, iside)
  nnodes = rzmesh_workspace%nnodes()
  ncells = rzmesh_workspace%ncells()
  nbsect = rzmesh_workspace%naux()
  rzmesh_workspace_allocated = .true.

  end subroutine mmesh_rzmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine get_rzmesh(nnodes, ncells, naux, x, quads, next_cell, aux_nodes, iline, itube)
  integer,      intent(in   ) :: nnodes, ncells, naux
  real(real64), intent(  out) :: x(2, nnodes)
  integer,      intent(  out) :: quads(4, ncells), next_cell(2, 4, ncells), aux_nodes(2, naux), iline(naux + nnodes), itube(ncells)


  x = rzmesh_workspace%x
  quads = rzmesh_workspace%quads
  next_cell = rzmesh_workspace%next_cell
  if (naux > 0) aux_nodes = rzmesh_workspace%aux_nodes
  iline = rzmesh_workspace%iline
  itube = rzmesh_workspace%itube

  end subroutine get_rzmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_rzmesh()


  call rzmesh_workspace%free()
  rzmesh_workspace_allocated = .false.

  end subroutine free_rzmesh
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  subroutine construct_flux_tubes(base_mesh, it_base, phi, filename)
  use moose_grids, only: qmesh
  use flare_mmesh, only: type_mmesh => mmesh, backend => construct_flux_tubes
  real(real64),     intent(in) :: base_mesh(:,:,:), phi(:)
  integer,          intent(in) :: it_base
  character(len=*), intent(in) :: filename

  type(qmesh) :: B
  type(type_mmesh) :: M


  B = qmesh(transpose(base_mesh(:,:,1)), transpose(base_mesh(:,:,2)))
  M = backend(B, it_base, phi)
  call M%savetxt(filename)

  end subroutine construct_flux_tubes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine trace_vertices(n, x0, phi, it0, report, x, iend)
  use flare_fieldline, only: fdriver
  use flare_mmesh_utils
  integer,      intent(in   ) :: n, it0
  real(real64), intent(in   ) :: x0(2, *), phi(:)
  logical,      intent(in   ) :: report
  real(real64), intent(  out) :: x(5, size(phi), n)
  integer,      intent(  out) :: iend(2, n)

  type(mdriver) :: M


  M%fdriver = fdriver(stop_at_boundary=.false.)
  M%report = report
  call M%trace_vertices(n, x0, phi, it0, x, iend)
  call M%free()

  end subroutine trace_vertices
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine rzbuffer(vacuum_boundary, plasma_boundary, width)
  use moose_geometry,     only: polygon2d
  use flare_mmesh_vacuum, only: make_buffer_zone
  real(real64), intent(in) :: vacuum_boundary(:,:), plasma_boundary(:,:), width

  type(polygon2d) :: V, P


  if (allocated(rzbuffer_x)) deallocate(rzbuffer_x)
  V = polygon2d(vacuum_boundary)
  P = polygon2d(plasma_boundary)


  call make_buffer_zone(V, P, width)
  allocate (rzbuffer_x, source=V%nodes())

  end subroutine rzbuffer
  !-----------------------------------------------------------------------------

end module mmesh
