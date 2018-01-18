program image_recon
use mpi 

implicit none
integer m,n,k,numit,nn
parameter (m = 239)
parameter (n = 432)
! change numit to 0 to see the original picture
parameter (numit = 10000)
character (len=32) :: filename

integer MPI_COMM_CART
integer, dimension(1) :: dims
logical, dimension(1) :: periods
integer :: ierror, up, down, my_rank, comm_size, tag
double precision :: time
real, allocatable, dimension(:,:) :: im, old, new, buf

call MPI_Init(ierror)
call MPI_Comm_size(MPI_COMM_WORLD, comm_size, ierror)
dims(1) = comm_size
periods(1) = .FALSE.

call MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, .TRUE., MPI_COMM_CART, ierror)
call MPI_Comm_rank(MPI_COMM_CART, my_rank, ierror)
call MPI_Cart_shift(MPI_COMM_CART, 0, 1, up, down, ierror)
tag = 1

!Allocate working arrays
allocate( im(0:m+1,0:n+1))
allocate(old(0:m+1,0:n+1))
allocate(new(0:m+1,0:n+1))
allocate(buf(m,n))

im(:,  :) = 0.d0
old(:, :) = 0.d0

if (my_rank.eq.0) then
  write(filename,'(a,I3,a,I3,a)'),"edge",m,"x",n,".dat"
  call datread(filename, buf, m, n)
end if

! Number of processes should divide n
nn = n / comm_size

! Timing
time = MPI_Wtime()

! Dispatch the data across processes
write(0, *) 'Scatter begin'
if (my_rank.eq.0) then
  call MPI_Scatter(buf, nn * m, MPI_REAL, MPI_IN_PLACE, nn * m, MPI_REAL, 0, MPI_COMM_CART, ierror)
else
  call MPI_Scatter(buf, nn * m, MPI_REAL,          buf, nn * m, MPI_REAL, 0, MPI_COMM_CART, ierror)
end if
write(0, *) 'Scatter completed'
im(1:m,1:nn) = buf(1:m,1:nn)
old(1:m,1:nn) = buf(1:m,1:nn)

! Main iterative loop
do k = 1, numit
! Exchange the hallo's
! send top hallo and receive from bottom 
!  call MPI_Sendrecv(old(1:m,1 ), m, MPI_REAL, up  , tag, old(1:m,nn+1), m, &
!                    MPI_REAL, down, tag, MPI_COMM_CART, MPI_STATUS_IGNORE, ierror)
! right hallo
!  call MPI_Sendrecv(old(1:m,nn), m, MPI_REAL, down, tag, old(1:m,0   ), m, &
!                    MPI_REAL, up  , tag, MPI_COMM_CART, MPI_STATUS_IGNORE, ierror)

  new(1:m,1:nn) = 0.25e0*(old(0:m-1,1:nn) + old(2:m+1,1:nn) + old(1:m,0:nn-1) + old(1:m,2:nn+1) - im(1:m,1:nn))
  old(1:m,1:nn)=new(1:m,1:nn)

  if ( mod(k,1000) == 0 .AND. my_rank == 0 ) then
    write(*,*) 'Iteration # ',k
  end if
end do

! Copy result into communication buffer
buf(1:m,1:nn)=old(1:m,1:nn)

! Collect the data and write into file
if (my_rank.eq.0) then
  call MPI_Gather(MPI_IN_PLACE, nn * m, MPI_REAL, buf, nn * m, MPI_REAL, 0, MPI_COMM_CART, ierror)
  time = MPI_Wtime() - time
  write(*,*) 'Execution time = ', time, ' s' 
  call pgmwrite('image.pgm', buf , m, n)
else
  call MPI_Gather(         buf, nn * m, MPI_REAL, buf, nn * m, MPI_REAL, 0, MPI_COMM_CART, ierror)
end if

!DeAllocate working arrays
deallocate( im)
deallocate(old)
deallocate(new)
deallocate(buf)

call MPI_Finalize(ierror)

end

