program label_main
    integer,allocatable,dimension(:,:) :: site
    integer,allocatable,dimension(:)   :: cluster_size
    integer :: i, j, label,iprint
    integer :: iter, maxiter, occ_sites, span_label
    real    :: randnum, p_c, p_inf, p_avg
    logical :: spanning


    call init_random_seed()


    n_columns = 10
    n_rows = 10
    iprint = 1
    allocate(site(n_columns,n_rows))
    allocate(cluster_size(n_columns*100))

    open(30,file='p_inf.csv')
    maxiter = 1

    p_c = 0.5

      p_avg = 0
      do iter = 1, maxiter
        label = 0
        cluster_size = 0
        occ_sites = 0
        do i = 1, n_columns
            do j = 1, n_rows
                call random_number(randnum)
                if (randnum < p_c) then
                    site(i,j) = -1
                    occ_sites = occ_sites + 1
                else
                    site(i,j) = 0
                end if
            end do
        end do

        do i = 1, n_columns
            do j = 1, n_rows
                if(site(i,j) == -1) then
                    label = label + 1
                    call label_site(i,j,label,site,n_columns,n_rows)
                end if
            end do
        end do

        spanning = .FALSE.
        do i = 1, n_columns
          do j = 1, n_columns
            if(site(i,1) == site(j,n_rows) .and. site(i,1) /= 0) then
              span_label = site(i,1)
              spanning = .TRUE.
            end if
            if(spanning) exit
          end do
          if(spanning) exit
        end do

        if(.NOT. spanning) then
          do i = 1, n_rows
            do j = 1, n_rows
              if(site(1,i) == site(n_columns,j) .and. site(1,i) /=0) then
                span_label = site(1,i)
                spanning = .TRUE.
              end if
              if(spanning) exit
            end do
            if(spanning) exit
          end do
        end if

        if(.NOT. spanning) then
          span_label = 0
          p_inf = 0.d0
        else
          p_inf = real(cluster_size(span_label)) / real(occ_sites)
        end if

        p_avg = p_avg + p_inf
      end do
      p_avg = p_avg / real(maxiter)
      print*, p_c, p_avg
      write(30,*) p_c, p_avg



    contains
        recursive subroutine label_site(i,j,label,site,n_columns,n_rows)
            implicit none
            integer,dimension(300,300) :: site
            integer :: i, j, m, n, label, n_columns, n_rows
            site(i,j) = label
            cluster_size(label) = cluster_size(label) + 1

            m = i - 1
            n = j
            if(m /= 0 .and. m <= n_columns) then
              if (site(m,n) == -1) call label_site(m,n,label,site,n_columns,n_rows)
            end if

            m = i + 1
            n = j
            if(m /= 0 .and. m <= n_columns) then
              if (site(m,n) == -1) call label_site(m,n,label,site,n_columns,n_rows)
            end if

            m = i
            n = j - 1
            if(n /= 0 .and. n <= n_rows) then
              if (site(m,n) == -1) call label_site(m,n,label,site,n_columns,n_rows)
            end if

            m = i
            n = j + 1
            if(n /= 0 .and. n <= n_rows) then
              if (site(m,n) == -1) call label_site(m,n,label,site,n_columns,n_rows)
            end if

        end subroutine
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
