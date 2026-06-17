module cli
  use kinds
  implicit none

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine equi2d_autoconf(filename, dtype, nsample, rrange, zrange)
  use moose_utils,    only: str
  use moose_geometry, only: hypersurf2d, polygon2d_hypersurf
  use flare_control,  only: report
  use flare_equi2d
  use flare_jorek,    only: load_jorek_equi2d
  use flare_m3dc1,    only: load_m3dc1_equi2d
  use flare_model,    only: boundary
  character(len=*), intent(in) :: filename, dtype
  integer,          intent(in) :: nsample
  real(dp),         intent(in) :: rrange(2), zrange(2)

  class(equi2d), allocatable :: E
  type(hypersurf2d) :: Bslice
  logical :: ex
  integer :: iu, n(2)


  print 1000
  print *, "Automatic configuration of toroidally symmetric equilibrium ..."
  print *
  report = .true.
 1000 format(80("="))


  ! load equilibrium
  select case(dtype)
  case("geqdsk")
     allocate (E, source=load_geqdsk_equi2d(filename))

  case("sonnet")
     allocate (E, source=load_sonnet_equi2d(filename))

  case("jorek")
     allocate (E, source=load_jorek_equi2d(filename))

  case("m3dc1")
     allocate (E, source=load_m3dc1_equi2d(filename))

  case default
  end select


  ! scan equilibrium domain for X-points
  n = nsample
  Bslice = boundary%rzslice(0.d0)
  call E%find_xpoints(Bslice, ".", n, rrange, zrange)


  ! save configuration
  inquire (file=".bfield", exist=ex)
  if (ex) then
     print 2000
     print 2001, dtype, filename
  else
     open  (newunit=iu, file=".bfield")
     write (iu, 2001) dtype, filename
     close (iu)
  endif
 2000 format(1x,"You may need to update .bfield to include the following configuration:")
 2001 format("[equi2d_",a,"]",/,"filename:   ",a)

  end subroutine equi2d_autoconf
  !-----------------------------------------------------------------------------

end module cli
