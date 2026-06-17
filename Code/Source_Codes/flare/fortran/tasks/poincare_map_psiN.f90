subroutine flare_task_poincare_map_psiN(psiN_start, psiN_end, npsiN, theta0, phi0, &
   direction, phi_section, nsym, npoints, nsections, bounded, output)
  !
  ! Create Poincare maps for equidistant initial points along radial direction.
  !
  ! **Parameters:**
  !
  ! :psiN_start:  Lower bound [normalized poloidal flux] for radial direction.
  !
  ! :psiN_end:    Upper bound [normalized poloidal flux] for radial direction.
  !
  ! :npsiN:       Number of equidistant points between ``psiN_start`` and ``psiN_end``.
  !
  ! :theta0:      Geometric poloidal angle [deg] for initial points.
  !
  ! :phi0:        Toroidal angle [deg] for initial points.
  !
  ! See :ref:`poincare_map_grid` for other parameters.
  !
  use iso_fortran_env
  use moose_mpi
  use moose_math,   only: pi, linspace, CYLINDRICAL_COORDINATES
  use moose_grids,  only: ugrid, r3grid, cylindrical_r3grid
  use moose_units,  only: METER, DEGREE
  use flare_model,  only: equi2d, assert_equi2d
  implicit none
  real(real64),     intent(in) :: psiN_start, psiN_end, theta0, phi0, phi_section
  integer,          intent(in) :: npsiN, nsym, npoints, nsections
  character(len=*), intent(in) :: direction, output
  logical,          intent(in) :: bounded

  character(len=256) :: label
  type(ugrid)  :: G
  type(r3grid) :: R3
  real(real64) :: psiN(0:npsiN-1), x(2), trad
  integer      :: i


  call assert_equi2d("poincare_map_psiN")
  trad = theta0 / 180.d0 * pi
  psiN = linspace(psiN_start, psiN_end, npsiN)
  G    = ugrid(npsiN, 2)
  do i=0,npsiN-1
     x = equi2d%rzcoords(psiN(i), trad)
     call G%set_node(i, x)
  enddo

  R3 = cylindrical_r3grid(G, METER, DEGREE, 1, 2, phi0, .true.)
  write (label, 1003) psiN_start, psiN_end, npsiN, theta0, phi0
 1003 format("Radial domain: psiN = ",f0.3," -> ",f0.3," with ",i0," points at theta = ",f0.3," deg and phi = ",f0.3," deg")


  call flare_task_poincare_map_grid_implementation(R3, label, direction, phi_section, nsym, npoints, nsections, bounded, output)

end subroutine flare_task_poincare_map_psiN
