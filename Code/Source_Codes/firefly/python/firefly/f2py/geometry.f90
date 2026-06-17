module geometry
  use kinds
  implicit none


  contains
  !-----------------------------------------------------------------------------

include "_geometry.f90"


  !-----------------------------------------------------------------------------
  subroutine init()
  use moose_mpi
  use flare_control, only: flare_init
  use firefly_version


  call flare_init(standalone=.false., greeting=.true.)

  if (rank == 0) then
     print *
     print *, "FIREFLY (" // achar(27) // "[94m" // version // achar(27) // "[0m)"
     print 1002
     print *
  endif
 1002 format(80("="))

  end subroutine init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function validate_torosurf(phi, rz, nsym, units)
  !
  ! explicit interface for validate_torosurf (pass arrays for torosurf components)
  !
  use moose_units,      only: length_scale
  use moose_math,       only: pi
  use moose_geometry,   only: torosurf, setup_torosurf
  use firefly_geometry, only: backend => validate_torosurf
  real(real64),     intent(in) :: phi(:), rz(:,:,:)
  integer,          intent(in) :: nsym
  character(len=*), intent(in) :: units
  logical                      :: validate_torosurf

  type(torosurf) :: T


  T = torosurf(size(phi)-1, size(rz, 2)-1, nsym)
  T%phi = phi / 180 * pi
  T%rz  = rz * length_scale(units)
  call setup_torosurf(T)
  validate_torosurf = backend(T)

  end function validate_torosurf
  !-----------------------------------------------------------------------------

end module geometry
