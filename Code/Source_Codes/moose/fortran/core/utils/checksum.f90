module moose_checksum
  use iso_fortran_env
  implicit none
  private


  ! straight forward implementation of Adler-32 checksum
  interface checksum
     procedure :: checksum_real64_dim1, checksum_real64_dim2, checksum_real64_dim3, &
        checksum_real64_dim4, checksum_int32_dim1, checksum_int32_dim2, &
        checksum_int32_dim3, checksum_int32_dim4
  end interface checksum

  public :: checksum

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function checksum_real64_dim1(array) result(checksum)
  real(real64), intent(in) :: array(:)
  integer                  :: checksum


  checksum = checksum_real64(size(array), array)

  end function checksum_real64_dim1
  !-----------------------------------------------------------------------------
  function checksum_real64_dim2(array) result(checksum)
  real(real64), intent(in) :: array(:,:)
  integer                  :: checksum


  checksum = checksum_real64(size(array), array)

  end function checksum_real64_dim2
  !-----------------------------------------------------------------------------
  function checksum_real64_dim3(array) result(checksum)
  real(real64), intent(in) :: array(:,:,:)
  integer                  :: checksum


  checksum = checksum_real64(size(array), array)

  end function checksum_real64_dim3
  !-----------------------------------------------------------------------------
  function checksum_real64_dim4(array) result(checksum)
  real(real64), intent(in) :: array(:,:,:,:)
  integer                  :: checksum


  checksum = checksum_real64(size(array), array)

  end function checksum_real64_dim4
  !-----------------------------------------------------------------------------
  function checksum_real64(n, array) result(checksum)
  integer,      intent(in) :: n
  real(real64), intent(in) :: array(*)
  integer                  :: checksum


  checksum = checksum_int32(n, transfer(array(1:n), checksum, 2*n))

  end function checksum_real64
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function checksum_int32_dim1(array) result(checksum)
  integer, intent(in) :: array(:)
  integer             :: checksum


  checksum = checksum_int32(size(array), array)

  end function checksum_int32_dim1
  !-----------------------------------------------------------------------------
  function checksum_int32_dim2(array) result(checksum)
  integer, intent(in) :: array(:,:)
  integer             :: checksum


  checksum = checksum_int32(size(array), array)

  end function checksum_int32_dim2
  !-----------------------------------------------------------------------------
  function checksum_int32_dim3(array) result(checksum)
  integer, intent(in) :: array(:,:,:)
  integer             :: checksum


  checksum = checksum_int32(size(array), array)

  end function checksum_int32_dim3
  !-----------------------------------------------------------------------------
  function checksum_int32_dim4(array) result(checksum)
  integer, intent(in) :: array(:,:,:,:)
  integer             :: checksum


  checksum = checksum_int32(size(array), array)

  end function checksum_int32_dim4
  !-----------------------------------------------------------------------------
  function checksum_int32(n, array) result(checksum)
  integer, intent(in) :: n
  integer, intent(in) :: array(*)
  integer             :: checksum

  integer, parameter :: mod_adler32 = 65521

  integer :: a, b, i


  a = 1
  b = 0
  do i=1,n
     a = mod(a + array(i), mod_adler32)
     b = mod(b + a, mod_adler32)
  enddo
  checksum = ior(b * 65536, a)

  end function checksum_int32
  !-----------------------------------------------------------------------------

end module moose_checksum
