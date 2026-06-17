subroutine flare_task_connection_histogram(psiN_start, psiN_end, nsteps, nturns, nsym, nphi, ntheta, param, output)
  !
  ! Compute connection length histogram for field lines starting on equilibrium flux surfaces.
  !
  ! **Parameters:**
  !
  ! :psiN_start:  Lower boundary [normalized poloidal flux] of radial domain.
  !
  ! :psiN_end:    Upper boundary [normalized poloidal flux] of radial domain.
  !
  ! :nsteps:      Number of steps in radial direction.
  !
  ! :nturns:      Max. number of toroidal turns for field line tracing.
  !
  ! :nsym:        Toroidal symmetry (360 deg / nsym).
  !
  ! :nphi:        Number of slices (equidistant steps in toroidal direction).
  !
  ! :ntheta:      Number of equidistant steps in poloidal direction.
  !
  ! :param:       Parametrization of poloidal direction:
  !
  !               :arclength:       Arc length along flux surface contour.
  !
  !               :magnetic_angle:  Straight field line coordinate.
  !
  ! :output:      Filename for output of field line statistics.
  !
  use iso_fortran_env
  use moose_math
  use flare_control
  use flare_model, only: bfield
  use flare_equi2d
  use flare_equi3d
  use flare_connection_length
  use flare_tasks
  implicit none
  real(real64),     intent(in) :: psiN_start, psiN_end
  integer,          intent(in) :: nsteps, nturns, nsym, nphi, ntheta
  character(len=*), intent(in) :: param, output

  type(connection_histogram) :: H
  real(real64) :: psiN(nsteps+1), turns(nturns)


  ! greeting
  call begin_task()
  if (rank == 0) then
     print *, "Computing connection length histogram ..."
     print *

     ! radial domain
     if (nsteps == 0) then
        print 1001, psiN_start
     elseif (nsteps > 0) then
        print 1002, psiN_start, psiN_end, nsteps
     else
        print 9001
        stop
     endif
     ! max. number of turns
     print 1003, nturns

     ! sample points on flux surface
     print 1004, nsym, nphi, ntheta
  endif
 1001 format(3x,"- At radial position [psiN]: ",f0.3)
 1002 format(3x,"- In radial domain [psiN]: ",f0.3," -> ",f0.3," with ",i0," steps")
 1003 format(8x,"tracing field lines up to ",i0," toroidal turns"/)
 1004 format(3x,"- Symmetry: ",i0,", toroidal x poloidal resolution: ",i0," x ",i0/)
 9001 format("ERROR in fieldline_loss: nsteps >= 0 required")


  ! execute task
  psiN  = linspace(psiN_start, psiN_end, nsteps+1)
  turns = geomspace(1.d0, 1.d0*nturns, nturns)
  select type(equi => bfield%equi)
  class is(equi2d)
     H = compute_connection_histogram(psiN, turns, nsym, nphi, ntheta, param)

  class is(equi3d)
     print *, "ERROR: loss_histogram for 3D equilibrium not implemented yet!"
     stop
  end select


  ! write output
  if (rank == 0) call H%savetxt(output)
  call H%free()
  call finalize_task()

end subroutine flare_task_connection_histogram
