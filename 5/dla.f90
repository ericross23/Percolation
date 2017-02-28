program dla

  implicit none
  integer,allocatable,dimension(:,:) :: site
  integer :: n_columns, n_rows, cluster_size, i, j, t_col, t_row
  integer :: n_occ
  real :: randnum, radius, dist
  real, parameter :: pi = 4.0 * atan(1.0)
  character(len=20) :: filename
  character(len=5) :: sizechar
  character(len=4) :: csv

  csv = '.csv'
  call init_random_seed()
  open(10,file='radius.csv')

  n_columns = 200
  n_rows = 200
  cluster_size = 10000
  allocate(site(n_columns+1,n_rows+1))

  radius = 10
  site = 0
  n_occ = 0
  ! Place the seed in the middle of the lattice.
  site(n_columns/2,n_rows/2) = 1
  do
    ! Place particle randomly on a circle of radius
    100 call random_number(randnum)
    randnum = 2 * pi * randnum
    t_col = ceiling(radius * cos(randnum) + n_columns / 2)
    call random_number(randnum)
    randnum = 2 * pi * randnum
    t_row = ceiling(radius * sin(randnum) + n_rows / 2)

    if (site(t_col,t_row) == 1) go to 100

    do
      call random_number(randnum)
      if(randnum < 0.25) then
        t_col = t_col + 1
      else if(randnum < 0.50) then
        t_col = t_col - 1
      else if(randnum < 0.75) then
        t_row = t_row + 1
      else
        t_row = t_row - 1
      end if
      dist = sqrt(real(n_columns/2 - t_col) ** 2 + real(n_rows/2 - t_row) ** 2)
      if(dist >= min(real((n_columns-2)/2.0), radius + 25.0)) go to 100

      if(site(t_col + 1, t_row) == 1 .or. site(t_col - 1, t_row) == 1 .or. &
      site(t_col, t_row + 1) == 1 .or. site(t_col, t_row - 1) == 1) then
        call random_number(randnum)
        if (randnum < 0.1) then
          site(t_col,t_row) = 1
          n_occ = n_occ + 1
          if(dist > radius) radius = dist
          write(10,*), n_occ, radius
          exit
        end if
      end if
    end do

    if(n_occ == 10000) then
    write(sizechar,'(I5)') n_occ
    filename = trim(adjustl(sizechar)) // trim(adjustl(csv))
    open(10,file=filename)
    do i = 1, n_columns
      do j = 1, n_rows
        write(10,'(I6)',advance='NO') site(i,j)
      end do
      write(10,*)
    end do
    close(10)
    end if
    print*, n_occ
    if(n_occ > cluster_size) exit
  end do

end program


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
