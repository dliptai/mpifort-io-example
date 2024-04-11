module io_write
  use mpi
  implicit none

contains

  subroutine write_fake_recordmarker(fh)
    use utils, only: get_file_end
    integer, intent(in) :: fh
    integer :: rm = -1 ! Assmuing record marker is the size of an integer
    integer :: ierr
    integer(kind=MPI_OFFSET_KIND) :: end_bytes

    call get_file_end(fh, end_bytes)
    call mpi_file_set_view(fh, end_bytes, MPI_INTEGER4, MPI_INTEGER4, 'native', MPI_INFO_NULL, ierr)
    call mpi_file_write_all(fh, rm, 1, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)

  end subroutine write_fake_recordmarker

  subroutine write_header(fh, header1, header2, legacy)
    use utils, only: total_elements, get_file_end
    integer, intent(in) :: fh
    integer, intent(in) :: header1
    real(8), intent(in) :: header2
    logical, intent(in), optional :: legacy
    logical :: islegacy = .false.
    integer :: ierr
    integer(kind=MPI_OFFSET_KIND) :: end_bytes
    if (present(legacy)) islegacy = legacy

    if (islegacy) call write_fake_recordmarker(fh)

    call get_file_end(fh, end_bytes)
    call mpi_file_set_view(fh, end_bytes, MPI_INTEGER4, MPI_INTEGER4, 'native', MPI_INFO_NULL, ierr)
    call mpi_file_write_all(fh, header1, 1, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)

    call get_file_end(fh, end_bytes)
    call mpi_file_set_view(fh, end_bytes, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ierr)
    call mpi_file_write_all(fh, header2, 1, MPI_REAL8, MPI_STATUS_IGNORE, ierr)

    if (islegacy) call write_fake_recordmarker(fh)

    call mpi_barrier(MPI_COMM_WORLD, ierr)

  end subroutine write_header

  subroutine write_array_mpi(fh, buffer, type_mpi_array, legacy)
    use utils, only: total_elements, get_file_end
    integer, intent(in) :: fh
    real(8), intent(in), allocatable :: buffer(:,:,:)
    integer, intent(in) :: type_mpi_array
    logical, intent(in), optional :: legacy
    logical :: islegacy = .false.
    integer :: status(MPI_STATUS_SIZE), ierr
    integer(kind=MPI_OFFSET_KIND) :: end_bytes
    if (present(legacy)) islegacy = legacy

    if (islegacy) call write_fake_recordmarker(fh)

    call get_file_end(fh, end_bytes)
    call mpi_file_set_view(fh, end_bytes, MPI_REAL8, type_mpi_array, 'native', MPI_INFO_NULL, ierr)
    call mpi_file_write_all(fh, buffer, total_elements(buffer), MPI_REAL8, status, ierr)

    if (islegacy) call write_fake_recordmarker(fh)

  end subroutine write_array_mpi

end module io_write
