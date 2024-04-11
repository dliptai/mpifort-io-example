program Main
  use mpi
  use vars
  use utils
  use io

  implicit none

  integer :: ierr
  integer :: i,j,k
  character(len=100) :: filename
  real(8), external :: xfunc
  integer :: fh

  call init_mpi
  ! Get filename from command-line argument on master rank
  if (myrank == 0) then
    call getarg(1, filename)
    if (trim(filename) == "") then
      print*, "No filename provided. Program will terminate."
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
  end if

  ! Broadcast filename to all ranks
  call MPI_BCAST(filename, len(filename), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  ! Set bounds of full domain
  is_global=1; ie_global=500
  js_global=1; je_global=1
  ks_global=1; ke_global=1

  call decompose_domain
  call create_type_mpi_array

  allocate(d(is:ie,js:je,ks:ke))
  allocate(v1,v2,v3,e,mold=d)

  call read_hormone_dump(trim(filename))

  do i = is, ie
    write(myrank+100,*) i, d(i,:,:)
  enddo

  call MPI_FINALIZE(ierr)

end program Main

real(8) function xfunc(i, j, k)
  implicit none
  integer, intent(in) :: i, j, k

  xfunc = real(100*i + 10*j + k, kind=8)

end function xfunc
