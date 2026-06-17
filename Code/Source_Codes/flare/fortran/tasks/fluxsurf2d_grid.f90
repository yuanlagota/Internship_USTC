subroutine flare_task_fluxsurf2d_grid(psiN, ntheta, endpoint, orientation, param, output)
  !
  ! Construct grid on contour on flux surface of toroidally symmetric equilibrium.
  !
  ! **Parameters:**
  !
  ! :psiN:         Radial position of flux surface [normalized poloidal flux].
  !
  ! :ntheta:       Number of equidistant steps in poloidal direction.
  !
  ! :endpoint:     Include/exclude periodic point in grid.
  !
  ! :orientation:  Orientation of flux surface contour:
  !
  !                :fwd, forward:            in direction of poloidal field.
  !
  !                :bwd, backward:           against direction of poloidal field.
  !
  !                :cw, clockwise:           in clockwise direction.
  !
  !                :ccw, counter-clockwise:  in counter-clockwise direction.
  !
  ! :param:        Parametrization of flux surface contour, one of `arclength` or `magnetic_angle`.
  !
  ! :output:       Filename for output of grid.
  !
  use iso_fortran_env
  use moose_mpi
  use moose_math
  use moose_grids, only: r3grid, cylindrical_r3grid, cgrid
  use moose_units, only: METER, DEGREE
  use flare_control
  use flare_model
  use flare_fluxsurf2d
  use flare_tasks
  implicit none
  real(real64),     intent(in) :: psiN
  integer,          intent(in) :: ntheta
  logical,          intent(in) :: endpoint
  character(len=*), intent(in) :: orientation, param, output

  type(cgrid), target :: F
  type(r3grid) :: G
  integer      :: idir


  call begin_task()
  if (rank == 0) then
     print *, "Constructing grid on contour on flux surface of toroidally symmetric equilibrium ..."
     print *
     idir = make_idir(orientation, -bfield%equi%Bp_sign)
  endif
  if (rank /= 0) return


  call assert_equi2d("fluxsurf2d_grid")
  F = fluxsurf2d_grid(equi2d%rzcoords(psiN, 0.d0), ntheta, endpoint, idir, param, boundary2d)
  G = cylindrical_r3grid(F, METER, DEGREE, 1, 2, 0.d0, .false.)


  call G%savetxt(output)
  call finalize_task()

end subroutine flare_task_fluxsurf2d_grid
