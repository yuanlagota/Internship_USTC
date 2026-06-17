subroutine flare_task_equi2d_poloidal_angle(psiN, ntheta, output)
  !
  ! Generate mapping between arclength and straight field line poloidal angle on toroidally symmetric equilibrium flux surface.
  !
  ! **Parameters:**
  !
  ! :psiN:    Radial location of flux surface [normalized poloidal flux].
  !
  ! :ntheta:  Number of steps in poloidal direction.
  !
  ! :output:  Filename for output of data set.
  !
  use iso_fortran_env
  use moose_mpi
  use moose_analysis, only: interp, pchip
  use flare_control
  use flare_model
  use flare_fluxsurf2d
  use flare_tasks
  implicit none
  real(real64),     intent(in) :: psiN
  integer,          intent(in) :: ntheta
  character(len=*), intent(in) :: output

  type(fluxsurf2d) :: F
  type(interp)     :: map


  call begin_task()
  if (rank == 0) then
     print *, "Constructing mapping for straight field line angular coordinate"
     print *
     print 1000, psiN, ntheta
     print 1001, trim(output)
     print *
  endif
 1000 format(3x,"- Flux surface psiN = ",f0.3," with ",i0," support points")
 1001 format(8x,"saving output to ",a)


  F = fluxsurf2d(equi2d%rzcoords(psiN, 0.d0))

  map = pchip(F%arcl, F%theta)
  call map%plot(output, nsample=ntheta)

  call finalize_task()

end subroutine flare_task_equi2d_poloidal_angle
