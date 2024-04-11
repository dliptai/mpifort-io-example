module io_read
  use mpi
  implicit none

  integer(kind=MPI_OFFSET_KIND) :: offset = 0

contains

  subroutine update_offset(fh, dtype)
    integer, intent(in) :: fh, dtype
    integer :: ierr
    integer(kind=MPI_OFFSET_KIND) :: bytes

    ! get the number of bytes read
    call mpi_file_get_type_extent(fh, dtype, bytes, ierr)
    offset = offset + bytes
    call mpi_barrier(MPI_COMM_WORLD, ierr)

  end subroutine update_offset

  subroutine read_fake_recordmarker(fh)
    integer, intent(in) :: fh
    integer :: rm
    integer :: ierr

    call mpi_file_set_view(fh, offset, MPI_INTEGER4, MPI_INTEGER4, 'native', MPI_INFO_NULL, ierr)
    call mpi_file_read_all(fh, rm, 1, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
    call update_offset(fh, MPI_INTEGER4)

  end subroutine read_fake_recordmarker

  subroutine read_header(fh, header1, header2, legacy)
    integer, intent(in) :: fh
    integer, intent(out) :: header1
    real(8), intent(out) :: header2
    logical, intent(in), optional :: legacy
    logical :: islegacy = .false.
    integer :: ierr
    if (present(legacy)) islegacy = legacy

    if (islegacy) call read_fake_recordmarker(fh)

    call mpi_file_set_view(fh, offset, MPI_INTEGER4, MPI_INTEGER4, 'native', MPI_INFO_NULL, ierr)
    call mpi_file_read_all(fh, header1, 1, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
    call update_offset(fh, MPI_INTEGER4)

    call mpi_file_set_view(fh, offset, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ierr)
    call mpi_file_read_all(fh, header2, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
    call update_offset(fh, MPI_REAL8)

    if (islegacy) call read_fake_recordmarker(fh)

    call mpi_barrier(MPI_COMM_WORLD, ierr)

  end subroutine read_header

  subroutine read_array_mpi(fh, buffer, type_mpi_array, legacy)
    use utils, only: total_elements
    integer, intent(in) :: fh
    real(8), intent(out) :: buffer(:,:,:)
    integer, intent(in) :: type_mpi_array
    logical, intent(in), optional :: legacy
    logical :: islegacy = .false.
    integer :: status(MPI_STATUS_SIZE), ierr
    if (present(legacy)) islegacy = legacy

    ! "read" fake record marker
    if (islegacy) call read_fake_recordmarker(fh)

    call mpi_file_set_view(fh, offset, MPI_DOUBLE_PRECISION, type_mpi_array, 'native', MPI_INFO_NULL, ierr)
    call mpi_file_read_all(fh, buffer, total_elements(buffer), MPI_DOUBLE_PRECISION, status, ierr)
    call update_offset(fh, type_mpi_array)

    ! "read" fake record marker
    if (islegacy) call read_fake_recordmarker(fh)

  end subroutine read_array_mpi

end module io_read
