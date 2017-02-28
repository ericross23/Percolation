program label_main
    integer,allocatable,dimension(:,:) :: site, spanning_cluster
    integer,allocatable,dimension(:)   :: cluster_size
    integer :: i, j, label, transv_ave, lon_ave, iprint, k, l
    integer ::lon_max, transv_max, lon, transv, i_min, i_max, j_min, j_max
    integer :: iter, maxiter, occ_sites, span_label, box_size, box_count, n_boxes
    real    :: avg_boxes
    real    :: randnum, p_c, start, finish, p_inf, p_avg
    logical :: spanning, occ_box

    call cpu_time(start)
    call init_random_seed()


    n_columns = 75
    n_rows = 75
    allocate(site(n_columns,n_rows))
    allocate(spanning_cluster(n_columns,n_rows))
    allocate(cluster_size(n_columns*n_rows))
    iprint = 0
    if(iprint == 2) open(10,file='cluster.csv')
    if(iprint == 2) open(20,file='clusterlabel.csv')
    open(30,file='p_inf.csv')
    maxiter = 100
    box_size = 3



    p_c = 0.592746

      p_avg = 0
      do iter = 1, maxiter
        75 lon_max = 0
        transv_max = 0
        lon_ave = 0
        transv_ave = 0
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
                if(iprint == 1) write(*,'(I6)',advance='NO') site(i,j)
                if(iprint == 2) write(10,'(I6)',advance='NO') site(i,j)
            end do
            if(iprint == 1) print*,
            if(iprint == 2) write(10,*)
        end do

        do i = 1, n_columns
            do j = 1, n_rows
                if(site(i,j) == -1) then
                    label = label + 1
                    lon = 0
                    transv = 0
                    i_min = i
                    i_max = i
                    j_min = j
                    j_max = j
                    call label_site(i,j,label,site,n_columns,n_rows,i_min,i_max,j_min,j_max)
                    transv = i_max - i_min + 1
                    if(transv_max < transv) transv_max = transv
                    transv_ave = transv_ave + transv
                    lon = j_max - j_min + 1
                    if(lon_max < lon) lon_max = lon
                    lon_ave = lon_ave + lon
                end if
            end do
        end do

        if(iprint==1 .or. iprint == 2) then
        print*, '-----------------------------------------------------------------'
        print*, 'Result: '
            do i = 1, n_columns
                do j = 1, n_rows
                    if(iprint == 1) write(*,'(I6)',advance='NO') site(i,j)
                    if(iprint == 2) write(20,'(I6)',advance='NO') site(i,j)
                end do
                if(iprint == 1) print*,
                if(iprint == 2) write(20,*)
            end do
        print*, '----------------------------------------------------------------'
        end if
        if(label /= 0) lon_ave = lon_ave / label
        if(label /= 0) transv_ave = transv_ave / label

        if(iprint == 1) then
          print*, 'Maximum cluster width: ', lon_max, 'Average cluster width: ', &
                  lon_ave, 'Lattice width: ', n_rows
          print*, 'Maximum cluster length: ', transv_max, 'Average cluster length: ',&
                  transv_ave, 'Lattice length: ', n_columns
        end if

        if(iprint == 1) then
          if(lon_max == n_rows .or. transv_max == n_columns) then
              print*, 'Cluster is spanning.'
          else
              print*, 'Cluster is not spanning.'
          end if
        end if
        call cpu_time(finish)

        if(iprint == 1) print*,  'Time = ', finish-start, 'seconds.'

        if(iprint == 1) then
          do i = 1, n_columns
            print*, i, cluster_size(i)
          end do
          print*, 'Occupied sites: ', occ_sites, sum(cluster_size)
        end if

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
        if(iprint==1) print*, "Spanning cluster: ", span_label, p_inf
        p_avg = p_avg + p_inf
        if(iprint == 2) close(10)
        if(iprint == 2) close(20)

        if(.NOT. spanning) go to 75
        ! Place spanning cluster into new object
        do i = 1, n_columns
          do j = 1, n_rows
            if(site(i,j) == span_label) then
              spanning_cluster(i,j) = 1
            else
              spanning_cluster(i,j) = 0
            end if
          end do
        end do

        if(iprint == 1) print*, 'Spanning cluster only: '
        do i = 1, n_columns
            do j = 1, n_rows
                if(iprint == 1) write(*,'(I6)',advance='NO') spanning_cluster(i,j)
            end do
            if(iprint == 1) print*,
        end do

        n_boxes = n_columns / box_size
        box_count = 0
        do i = 1, n_boxes
          do j = 1, n_boxes
            occ_box = .FALSE.
            do k = (i-1)*box_size + 1, (i-1)*box_size + box_size
              do l = (j-1)*box_size + 1, (j-1) * box_size + box_size
                if(spanning_cluster(k,l) == 1) then
                  occ_box = .TRUE.
                end if
              end do
            end do
            if(occ_box) box_count = box_count + 1
          end do
        end do

        print*, 'Individual box count: ', box_count
        avg_boxes = avg_boxes + real(box_count)
      end do
      p_avg = p_avg / maxiter
      avg_boxes = avg_boxes / maxiter
      print*, 'p_c ',p_c, 'p_inf ',p_avg
      write(30,*) p_c, p_avg
      print*, 'Average occupied boxes: ', avg_boxes


    contains
        recursive subroutine label_site(i,j,label,site,n_columns,n_rows,i_min,i_max,j_min,j_max)
            implicit none
            integer,allocatable,dimension(:,:) :: site
            integer :: i, j, m, n, label, n_columns, n_rows, i_min, i_max, j_min, j_max
            site(i,j) = label
            cluster_size(label) = cluster_size(label) + 1
            if (i_min > i) i_min = i
            if (i_max < i) i_max = i
            if (j_min > j) j_min = j
            if (j_max < j) j_max = j

            m = i - 1
            n = j
            if(m /= 0 .and. m <= n_columns) then
              if (site(m,n) == -1) call label_site(m,n,label,site,n_columns,n_rows,i_min,i_max,j_min,j_max)
            end if

            m = i + 1
            n = j
            if(m /= 0 .and. m <= n_columns) then
              if (site(m,n) == -1) call label_site(m,n,label,site,n_columns,n_rows,i_min,i_max,j_min,j_max)
            end if

            m = i
            n = j - 1
            if(n /= 0 .and. n <= n_rows) then
              if (site(m,n) == -1) call label_site(m,n,label,site,n_columns,n_rows,i_min,i_max,j_min,j_max)
            end if

            m = i
            n = j + 1
            if(n /= 0 .and. n <= n_rows) then
              if (site(m,n) == -1) call label_site(m,n,label,site,n_columns,n_rows,i_min,i_max,j_min,j_max)
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
