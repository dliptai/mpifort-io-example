program Main
  use mpi
  use utils
  use io_write
  use io_read

  implicit none

  integer :: ierr, type_mpi_array
  integer :: mu,i,j,k, ranki, diff
  character(len=*), parameter :: filename = 'output.dat'
  integer, allocatable, dimension(:,:,:,:) :: arr
  integer, external :: xfunc
  logical :: legacy = .false.
  integer :: fh, error

  call init_mpi

  is_global = 1 ; ie_global = 5
  js_global = 1 ; je_global = 5
  ks_global = 1 ; ke_global = 5
  mus       = 1 ; mue       = 5
  nmu = mue - mus + 1

  call decompose_domain
  call create_type_mpi_array(type_mpi_array)

  ! allocate and initialise array
  allocate(arr(mus:mue,is:ie,js:je,ks:ke))

  do mu = mus,mue
    do i = is,ie
      do j = js,je
        do k = ks,ke
          arr(mu,i,j,k) = xfunc(mu,i,j,k)
        end do
      end do
    end do
  end do

  ! create empty file on master task
  if (myrank == 0) then
    open(newunit=fh, file=filename, status='replace', form='unformatted', action='write', access='stream')
    close(unit=fh)
  end if
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  call open_file_write(filename, fh)
  call write_int4_array_4d(fh, arr, type_mpi_array)
  call close_file(fh)
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! put junk in array
  arr = -1.0

  ! -------------- Now read back the same file -----------------------------------------

  ! call is_file_legacy(filename, legacy) ! this doesn't work now, because the first thing in the header is an integer...

  call open_file_read(filename, fh)
  call read_int4_array_4d(fh, arr, type_mpi_array)
  call close_file(fh)

  do ranki=0,nprocs-1
    if (myrank == ranki) then
      do mu = mus,mue
        do i = is,ie
          do j = js,je
            do k = ks,ke
              print*, 'Rank ', myrank, mu, i, j, k, arr(mu,i,j,k)
            end do
          end do
        end do
      end do
    endif
    call mpi_barrier(MPI_COMM_WORLD, ierr)
  end do

  error = 0

  if (myrank ==0 ) write(*,'(A)') '>>> Checking arrays...'
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  do ranki = 0, nprocs-1
    if (myrank == ranki) then
      do mu = mus,mue
        do i = is,ie
          do j = js,je
            do k = ks,ke
              diff = abs(arr(mu,i,j,k) - xfunc(mu,i,j,k))
              if (diff > 0) then
                print*, 'Rank ', myrank, 'has a problem at ', mu, i, j, k, arr(mu,i,j,k), diff
                error = error + 1
              end if
            end do
          end do
        end do
      end do
    endif
    call mpi_barrier(MPI_COMM_WORLD, ierr)
  end do

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  if (error == 0) then
    print*, 'Rank ', myrank, ' All tests passed'
  else
    print*, 'Rank ', myrank, ' N tests failed:', error
  end if

  call MPI_FINALIZE(ierr)

end program Main

integer function xfunc(i, j, k, l)
  implicit none
  integer, intent(in) :: i, j, k, l

  xfunc = 1000*i + 100*j + 10*k + 1*l

end function xfunc
