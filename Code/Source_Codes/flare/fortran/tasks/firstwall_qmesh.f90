subroutine flare_task_firstwall_qmesh(filename, phi, nrad, npol, eps, output)
  !
  ! Construct mesh in R-Z plane between poloidally closed first-wall contour and user-defined inner boundary contour.
  !
  ! **Parameters:**
  !
  ! :filename:    Geometry file for inner mesh boundary.
  !
  ! :phi:         Toroidal angle [deg].
  !
  ! :nrad:        Number of radial points in mesh.
  !
  ! :npol:        Number of poloidal points in mesh.
  !
  ! :eps:         Offset [mm] from first wall for final mesh contour.
  !
  ! :output:      Name for output of mesh.
  !
  use iso_fortran_env
  use moose_error
  use moose_mpi
  use moose_utils, only: ordinal
  use moose_math,  only: pi, CYLINDRICAL_COORDINATES
  use moose_r3grid
  use moose_units
  use moose_geometry
  use moose_qmesh_generator
  use flare_model
  use flare_control
  use flare_boundary
  use flare_tasks
  implicit none
  character(len=*), intent(in) :: filename, output
  real(real64),     intent(in) :: phi, eps
  integer,          intent(in) :: nrad, npol

  type(polygon2d) :: P
  type(qmesh)     :: mesh
  type(r3grid)    :: G
  class(curve), allocatable :: C
  integer :: ierr


  call begin_task()
  if (rank == 0) then
     print 1000
     print *

     ! load user defined contour
     print 2000, trim(filename)
     print *
     allocate (C, source=loadtxt_curve(filename))

     ! slicing first wall at phi
     print 3000, phi
     print *
     P = firstwall_rzslice(phi / 180.d0 * pi)
     if (very_verbose) call P%savetxt("FIRSTWALL")

     ! apply offset
     print 4000, eps
     print *
     call P%shift(abs(eps*1.d-3) * P%orientation(), ierr)
     if (ierr /= 0) call ERROR("construction of offset surface from first wall failed")
     if (very_verbose) call P%savetxt("FIRSTWALL_OFFSET")

     ! begin mesh construction
     print 5000
     mesh = quasi_orthogonal_qmesh(C, bspline_polygon(P), nrad, npol, 1)

     ! r3grid
     G = cylindrical_r3grid(mesh, METER, DEGREE, 1, 2, phi, .true.)
     call G%savetxt(output)
  endif
 1000 format(1x,"Constructing mesh in annular domain between user defined contour and first wall ...")
 2000 format(3x,"- Loading user defined contour from ",a)
 3000 format(3x,"- Slicing first wall at ",f0.3," deg")
 4000 format(3x,"- Applying offset to first wall: ",f0.3," mm")
 5000 format(3x,"- Beginning mesh construction ...")


  call finalize_task()

end subroutine flare_task_firstwall_qmesh
