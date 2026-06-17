subroutine flare_task_poincare_map_Rlinspace(R_start, R_end, nR, Z0, phi0, direction, phi_section, &
  nsym, npoints, nsections, bounded, output)
  !
  ! Create Poincare maps for equidistant initial points along radial direction.
  !
  ! **Parameters:**
  !
  ! :R_start:  Lower bound for major radius [m].
  !
  ! :R_end:    Upper bound for major radius [m].
  !
  ! :nR:       Number of equidistant points between ``R_start`` and ``R_end``.
  !
  ! :Z0:       Vertical position [m] for initial points.
  !
  ! :phi0:     Toroidal position [deg] for initial points.
  !
  ! See :ref:`poincare_map_grid` for other parameters.
  !
  use iso_fortran_env
  use moose_mpi
  use moose_math,   only: linspace, CYLINDRICAL_COORDINATES
  use moose_grids,  only: mesh1d, r3grid, cylindrical_r3grid
  use moose_units,  only: METER, DEGREE
  use flare_model,  only: equi2d, assert_equi2d
  implicit none
  real(real64),     intent(in) :: R_start, R_end, Z0, phi0, phi_section
  integer,          intent(in) :: nR, nsym, npoints, nsections
  character(len=*), intent(in) :: direction, output
  logical,          intent(in) :: bounded

  character(len=256) :: label
  type(mesh1d) :: G
  type(r3grid) :: R3


  G  = mesh1d(linspace(R_start, R_end, nR))
  R3 = cylindrical_r3grid(G, METER, DEGREE, 1, Z0, phi0, .true.)
  write (label, 1003) R_start, R_end, nR, Z0, phi0
 1003 format("Radial domain: R = ",f0.4," -> ",f0.4," m with ",i0," points at Z0 = ",f0.3," m and phi0 = ",f0.2," deg")


  call flare_task_poincare_map_grid_implementation(R3, label, direction, phi_section, nsym, npoints, nsections, bounded, output)

end subroutine flare_task_poincare_map_Rlinspace
