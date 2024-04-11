program Main
  use mpi
  use vars
  use utils
  use io

  implicit none

  integer :: ierr
  integer :: i,j,k
  character(len=*), parameter :: filename = 'sodshock_x.bin'
  real(8), external :: xfunc
  integer :: fh

  call init_mpi

  ! Set bounds of full domain
  is_global=1; ie_global=500
  js_global=1; je_global=1
  ks_global=1; ke_global=1

  call decompose_domain
  call create_type_mpi_array

  allocate(d(is:ie,js:je,ks:ke))
  allocate(v1,v2,v3,e,mold=d)

  call read_hormone_dump(filename)

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
