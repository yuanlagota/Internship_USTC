subroutine flare_task_safety_factor(psiN_start, psiN_end, nsteps, output)
  !
  ! Compute safety factor profile for toroidally symmetric equilibrium.
  !
  ! **Parameters:**
  !
  ! :psiN_start:  Lower boundary [normalized poloidal flux] for radial domain.
  !
  ! :psiN_end:    Upper boundary [normalized poloidal flux] for radial domain.
  !
  ! :nsteps:      Number of equidistant steps between ``psiN_start`` and ``psiN_end``.
  !
  ! :output:      Filename for output of safety factor profile.
  !
  use iso_fortran_env
  use moose_mpi
  use moose_math,  only: linspace
  use moose_dataset
  use flare_control
  use flare_model
  use flare_fluxsurf2d
  use flare_tasks
  implicit none
  real(real64),     intent(in) :: psiN_start, psiN_end
  integer,          intent(in) :: nsteps
  character(len=*), intent(in) :: output

  real(real64), pointer :: element(:)
  type(fluxsurf2d) :: F
  type(dataset)    :: D
  real(real64) :: psiN(0:nsteps), q, L, area, Vprime, current, x(2)
  integer :: i


  ! greeting & input verification
  call begin_task()
  if (rank == 0) then
     print *, "Computing safety factor for toroidally symmetric equilibrium."
     if (psiN_end < psiN_start) then
        print 9001;   stop
     endif
!     if (psiN_end >= 1.d0) then
!        print 9002;   stop
!     endif
     print 1000, psiN_start, psiN_end, nsteps
     print *
  endif
  call assert_equi2d("safety_factor")
 1000 format(3x,"- Between ",f0.3," and ",f0.3," with ",i0," steps")
 9001 format("ERROR: psiN_start <= psiN_end required")
 9002 format("ERROR: psiN_end < 1 required")


  ! construct flux surfaces and compute safety factor
  D = dataset(7, nsteps+1, geometry="mesh1d(psiN)")
  call D%set_metadata(1, "psiN",   "Normalized Poloidal Flux")
  call D%set_metadata(2, "q",      "Safety Factor")
  call D%set_metadata(3, "R",      "Major radius", "m")
  call D%set_metadata(4, "L",      "Circumference", "m")
  call D%set_metadata(5, "A",      "Surface area", "m**2")
  call D%set_metadata(6, "Vprime", "Derivative of volume with poloidal flux")
  call D%set_metadata(7, "I",      "Current", "MA")
  psiN = linspace(psiN_start, psiN_end, nsteps+1)
  call progress_bar(0, nsteps+1)
  do i=rank,nsteps,nproc
     if (verbose) print *, i, psiN(i)

     x = equi2d%rzcoords(psiN(i), 0.d0)
     call fluxsurf2d_parameters(x, q, L, area, Vprime, current)
     element => D%element(i);   element = [psiN(i), q, x(1), L, area, Vprime, current]

     call progress_bar(i+1, nsteps+1)
  enddo
  call finalize_progress_bar()
  call D%allreduce()


  ! save output
  if (rank == 0) then
     call D%savetxt(output)
  endif
  call D%free()
  call finalize_task()

end subroutine flare_task_safety_factor
