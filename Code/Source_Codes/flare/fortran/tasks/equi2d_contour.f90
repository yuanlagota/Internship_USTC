subroutine flare_task_equi2d_contour(r0, z0, direction, param, output)
  !
  ! Generate contour of toroidally symmetric equilibrium flux surface.
  !
  ! **Parameters:**
  !
  ! :r0:         R-coordinate [m] for reference point on flux surface.
  !
  ! :z0:         Z-coordinate [m] for reference point on flux surface.
  !
  ! :direction:  Trace direction:
  !
  !              :fwd, forward:            Trace in direction of poloidal field.
  !
  !              :bwd, backward:           Trace against direction of poloidal field.
  !
  !              :cw, clockwise:           Trace in clockwise direction.
  !
  !              :ccw, counter-clockwise:  Trace in counter-clockwise direction.
  !
  ! :param:      Parametrization of contour:
  !
  !              :arclength:               by arc length [m]
  !
  !              :magnetic_angle:          by (straight field line) poloidal angle [rad]
  !
  !              :geometric_angle:         by geometric poloidal angle [rad]
  !
  ! :output:     Filename for output of data set.
  !
  use iso_fortran_env
  use moose_mpi
  use flare_control
  use flare_equi2d
  use flare_model,    only: bfield, boundary2d
  use flare_fluxsurf2d
  use flare_tasks
  implicit none
  real(real64),     intent(in) :: r0, z0
  character(len=*), intent(in) :: direction, param, output

  type(fluxsurf2d) :: F
  integer :: idir


  call begin_task()
  if (rank /= 0) return
  print 1000
  print 1001, r0, z0
  idir = make_idir(direction, -bfield%equi%Bp_sign)
  print *
 1000 format("Tracing 2-D contour of (toroidally symmetric) equilibrium flux surface")
 1001 format(3x,"- Initial location: (",f0.3,", ",f0.3,") m"/)


  ! construct flux surface contour
  F = fluxsurf2d([r0, z0], idir, param, boundary2d)
  call F%savetxt(output)
  call finalize_task()

end subroutine flare_task_equi2d_contour
