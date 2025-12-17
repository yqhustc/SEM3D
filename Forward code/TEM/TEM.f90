
    program TEM
    use fwd_modeling
    use file_option
    use mkl_service
    use omp_lib
    !use system_lib
    USE IFPORT
    implicit none
    type(model) :: fwd_mod

    real(8) :: t1,t2,t3,t4,t5,t6
    t1=timef()
    call read_mesh(fwd_mod)
    t2=timef()
    print*, "Time passes:",(t2-t1),"s"
    call read_ctl(fwd_mod)
    t3=timef()
    print*, "Time passes:",(t3-t1),"s"
    call compute_boundary(fwd_mod)
    t4=timef()
    print*, "Time passes:",(t4-t1),"s"
    call compute_matrix(fwd_mod)
    t5=timef()
    print*, "Time passes:",(t5-t1),"s"
    call time_loop_1d(fwd_mod)
    t6=timef()
    print*, "Time passes:",(t6-t1),"s"
    call clear_model(fwd_mod)

    end program TEM

