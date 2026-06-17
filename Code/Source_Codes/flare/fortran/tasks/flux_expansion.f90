subroutine flare_task_flux_expansion(side, dr, nr, output)
  !
  ! Evaluate integral and local flux expansion on divertor target.
  !
  ! **Parameters:**
  !
  ! :side:     Side of strike point (fwd, bwd, cw, ccw).
  !
  ! :dr:       Radial interval [mm].
  !
  ! :nr:       Number of steps.
  !
  ! :output:   Name of output file.
  !
  use iso_fortran_env
  use moose_math
  use moose_dataset
  use flare_model
  use flare_control
  use flare_tasks
  use flare_fluxsurf2d
  implicit none
  character(len=*), intent(in) :: side, output
  real(real64),     intent(in) :: dr
  integer,          intent(in) :: nr

  real(real64), pointer :: column(:)
  type(dataset) :: D
  real(real64)  :: fX(2,nr)
  integer :: iside


  call begin_task()
  if (rank == 0) then
     print *, "Evaluating flux expansion on divertor target"
     print *
     print 1001, dr, nr
     iside = make_idir(side, -bfield%equi%Bp_sign)
  endif
  call proc(0)%broadcast(iside)
 1001 format(3x,"- Scanning radial interval of ",f0.1," mm at outboard midplane with ",i0," steps",/)


  fX = flux_expansion(iside, dr, nr)
  D = dataset(3, nr, geometry="grid1d(r)")
  call D%set_metadata(1, "r",   "R - R_sep", "mm")
  call D%set_metadata(2, "fX",  "integral flux expansion")
  call D%set_metadata(3, "dfX", "local flux expansion")
  column => D%column(1);   column = linspace(1.d0, 1.d0*nr, nr) / nr * dr
  column => D%column(2);   column = fX(1,:)
  column => D%column(3);   column = fX(2,:)


  if (rank == 0) call D%savetxt(output)
  call finalize_task()

end subroutine flare_task_flux_expansion
