module geqdsk
  use kinds
  implicit none


  integer      :: nr, nz, nbbbs, limitr
  real(real64) :: Rdim, Zdim, Rcentr, Rleft, Zmid, Rmaxis, Zmaxis, Simag, &
              Sibry, Bcentr, Current
  real(real64), allocatable :: fpol(:), pres(:), ffprim(:), pprime(:), psirz(:,:), qpsi(:)
  real(real64), allocatable :: rbbbs(:), zbbbs(:), rlim(:), zlim(:)


  logical, private :: first_run = .true.

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine load(filename)
  use flare_equi2d, only: geqdsk_dataset => geqdsk, load_geqdsk
  character(len=*), intent(in) :: filename

  type(geqdsk_dataset) :: g


  if (first_run) then
     first_run = .false.
  else
     deallocate (fpol, pres, ffprim, pprime, psirz, qpsi, rbbbs, zbbbs, rlim, zlim)
  endif

  g = load_geqdsk(filename)
  ! metadata
  nr      = g%nr
  nz      = g%nz
  Rdim    = g%Rdim
  Zdim    = g%Zdim
  Rcentr  = g%Rcentr
  Rleft   = g%Rleft
  Zmid    = g%Zmid
  Rmaxis  = g%Rmaxis
  Zmaxis  = g%Zmaxis
  Simag   = g%Simag
  Sibry   = g%Sibry
  Bcentr  = g%Bcentr
  Current = g%Current

  ! equilibrium data
  allocate (fpol,   source=g%fpol)
  allocate (pres,   source=g%pres)
  allocate (ffprim, source=g%ffprim)
  allocate (pprime, source=g%pprime)
  allocate (psirz,  source=transpose(g%psirz))
  allocate (qpsi,   source=g%qpsi)

  ! device and plasma boundary
  nbbbs   = g%nbbbs
  limitr  = g%limitr
  allocate (rbbbs, source=g%rbbbs)
  allocate (zbbbs, source=g%zbbbs)
  allocate (rlim,  source=g%rlim)
  allocate (zlim,  source=g%zlim)

  end subroutine load
  !-----------------------------------------------------------------------------

end module geqdsk
