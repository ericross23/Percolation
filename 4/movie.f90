program eden
  implicit none
  integer,allocatable,dimension(:,:) :: site
  integer :: n_columns, n_rows, cluster_size, i, j, t_col, t_row, loops
  real :: randnum
  character(len=20) :: filename
  character(len=5) :: sizechar
  character(len=4) :: csv
  call init_random_seed()

  csv = '.csv'
  n_columns = 200
  n_rows = 200
  cluster_size = 10000
  allocate(site(n_columns+1,n_rows+1))

  site = 0

  ! Place the seed in the middle of the lattice.
  site(n_columns/2,n_rows/2) = 1
  do loops = 1, cluster_size-1

    ! Randomly select an occupied site
    100 call random_number(randnum)
    t_col = ceiling(randnum*n_columns) + 1
    call random_number(randnum)
    t_row = ceiling(randnum*n_rows) + 1
    if(site(t_col,t_row) == 0) go to 100
    if(site(t_col+1,t_row) == 1 .and. site(t_col-1,t_row) == 1 .and. &
       site(t_col,t_row+1) == 1 .and. site(t_col,t_row-1) == 1) go to 100

    ! Randomly populate a nearest neighbor of the randomly selected site
    110 call random_number(randnum)
    if (randnum >= 0 .and. randnum < 0.25) then
      if (site(t_col + 1, t_row) == 0 .and. t_col <= n_columns) then
        site(t_col + 1, t_row) = 1
      else
        go to 110
      end if
    else if (randnum >= 0.25 .and. randnum < 0.50) then
      if (site(t_col - 1, t_row) == 0 .and. t_col > 0) then
        site(t_col - 1, t_row) = 1
      else
        go to 110
      end if
    else if (randnum >= 0.50 .and. randnum < 0.75) then
      if (site(t_col, t_row + 1) == 0 .and. t_row <= n_rows) then
        site(t_col, t_row + 1) = 1
      else
        go to 110
      end if
    else
      if (site(t_col, t_row - 1) == 0 .and. t_row > 0) then
        site(t_col, t_row - 1) = 1
      else
        go to 110
      end if
    end if
    print*, loops

    write(sizechar,'(I5)') loops
    filename = trim(adjustl(sizechar)) // trim(adjustl(csv))
    open(10,file=filename)
    do i = 1, n_columns
      do j = 1, n_rows
        write(10,'(I6)',advance='NO') site(i,j)
      end do
      write(10,*)
    end do
    close(10)
  end do
end program eden

subroutine init_random_seed()
    implicit none
    integer                          :: i, n, clock
    integer,allocatable,dimension(:) :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i-1, i = 1, n) /)
    call random_seed(put = seed)

end subroutine
subroutine generate_seed(seed)
    implicit none
    integer :: seed, clock

    call system_clock(count=clock)
    seed = clock + 37

end subroutine
