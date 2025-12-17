    module fwd_modeling
    USE IFPORT
    implicit none
    private
    integer :: ctl_unit = 101, fed_unit = 102, opt_unit = 11
    real(8),parameter :: pi  = 3.141592653589793238462643383279502884197_8
    real(8),parameter :: mu = pi * 4.d-7
    type :: integer_1d_ptr
        integer, pointer :: val(:)
    end type
    type :: integer_2d_ptr
        integer, pointer :: val(:,:)
    end type
    type :: real1d_ptr
        real(8), pointer :: val(:)
    end type
    type :: real2d_ptr
        real(8), pointer :: val(:,:)
    end type
    type :: real4d_ptr
        real(8), pointer :: val(:,:,:,:)
    end type
    type model
        !model information
        integer :: rec_num,t_inp_len,t_res_len,fault_num
        integer,allocatable :: x_num(:),y_num(:),src_num(:)
        !time series
        real(8),allocatable :: t_inp(:),t_sgl(:)
        !source parameter
        type(real1d_ptr),allocatable :: pulse_current(:)
        integer,allocatable :: degree_val(:)
        integer :: degree_range
        type(integer_2d_ptr), allocatable :: src_centre(:)
        !temperatrue
        real(8) :: T_int,ana_cof,T_end
        integer :: move
        !field result
        type(real4d_ptr), allocatable :: fed_all(:)
        type(real1d_ptr), allocatable :: Normal_Factor(:)
        !objective function weight
        real(8),allocatable :: weight(:)
        real(8),allocatable :: w1(:),w2(:),w3(:),w4(:),w5(:),w6(:)
        !file
        character(80) :: file_root,inv_data
        character(80),allocatable :: forward_name(:)
        real(8) :: intensity_factor,smooth_weight_a,smooth_weight_t,entropy_weight_a,entropy_weight_t,L1_weight,single_weight
    end type model
    public :: model,integer_1d_ptr,integer_2d_ptr,real1d_ptr,real4d_ptr,real2d_ptr
    public :: read_ctl,clear_model,fed_interp_parallel,norm,del2,loss_func,simulated_annealing_t,normal,rand_par
    contains

    subroutine read_ctl(cmod)
    implicit none
    type(model),intent(out) :: cmod
    !file
    character(80) :: ctl_file,fed_file,fed_num,deg_num
    character(180) :: sLine, sCode, sValue
    logical :: bComment
    real(8) :: time1,time2
    real(8) :: t1,dt,t2,t_res_min,t_res
    real(8) :: d_dg,dmin,dmax
    type(real2d_ptr),allocatable :: all_Normal_Factor(:)
    integer :: i,j,k,ierr,step,nt,nt1,fn,m
    cmod%rec_num=0
    cmod%T_int=0.d0
    cmod%ana_cof=0.d0
    cmod%T_end=0.d0
    t_res_min=0.d0
    call getarg(2,cmod%file_root)
    ctl_file = "inv_file/"//trim(cmod%file_root)//".inv"
    print*, "==================================building ctl file..."
    print*, "--- reading control file ..."
    open(ctl_unit,file=trim(ctl_file))
    do while (.true.)
        read( ctl_unit, '(A180)', iostat = ierr ) sLine
        if (ierr /= 0) exit
        call parseCode( len(sLine), sLine, sCode, sValue, bComment )
        if( bComment ) cycle
        select case (trim(sCode))
        case('fault_num')
            read(sValue,*) cmod%fault_num
        case('forward_data')
            allocate(cmod%forward_name(cmod%fault_num))
            do i=1,cmod%fault_num
                read(ctl_unit,*) cmod%forward_name(i)
            end do
        case('inv_data')
            read(ctl_unit,*) cmod%inv_data
        case('src_number')
            allocate(cmod%x_num(cmod%fault_num),cmod%y_num(cmod%fault_num),cmod%src_num(cmod%fault_num))
            allocate(cmod%src_centre(cmod%fault_num))
            do fn=1,cmod%fault_num
                read(ctl_unit,*) cmod%x_num(fn),cmod%y_num(fn)
                print*,cmod%x_num(fn),cmod%y_num(fn)
                cmod%src_num(fn)=cmod%x_num(fn)*cmod%y_num(fn)
                allocate(cmod%src_centre(fn)%val(cmod%src_num(fn),2))
                k=1
                do i=1,cmod%x_num(fn)
                    do j=1,cmod%y_num(fn)
                        cmod%src_centre(fn)%val(k,1) = j
                        cmod%src_centre(fn)%val(k,2) = i
                        k=k+1
                    end do
                end do
                if (cmod%src_num(fn)/=k-1) print*,"source number error!!"
            end do
        case('degree_range')
            read(sValue,*) cmod%degree_range
            allocate(cmod%degree_val(cmod%degree_range))
            read(ctl_unit,*) dmin,d_dg,dmax
            do i=1,cmod%degree_range
                cmod%degree_val(i)=dmin+(i-1)*d_dg
            end do
            if (cmod%degree_val(cmod%degree_range)/=dmax) then
                stop 'degree error!'
            end if
        case('rec_number')
            read(sValue,*) cmod%rec_num
        case('result_time_length')
            read(sValue,*) cmod%t_res_len
        case('interpolation_time')
            read(ctl_unit,*) t1,dt,t2
            cmod%t_inp_len=int((t2-t1)/dt+1)
            allocate(cmod%t_inp(cmod%t_inp_len))
            cmod%t_inp= [((i-1)*dt,i=1,cmod%t_inp_len)]
        case('initial_temperatrue')
            read(sValue,*) cmod%T_int
        case('annealing_coefficient')
            read(sValue,*) cmod%ana_cof
        case('termination_temperature')
            read(sValue,*) cmod%T_end
        case('move')
            read(sValue,*) cmod%move
        case('intensity_range_factor')
            read(sValue,*) cmod%intensity_factor
        case('smooth_constraint_weight_a')
            read(sValue,*) cmod%smooth_weight_a
        case('smooth_constraint_weight_t')
            read(sValue,*) cmod%smooth_weight_t
        case('entropy_constraint_weight_a')
            read(sValue,*) cmod%entropy_weight_a
        case('entropy_constraint_weight_t')
            read(sValue,*) cmod%entropy_weight_t
        case('single_constraint_weight')
            read(sValue,*) cmod%single_weight
        case('sparse_constraint_weight')
            read(sValue,*) cmod%L1_weight
        end select
    end do
    close(ctl_unit)
    print*, "--- calculation parameters ..."
    print*, "number of receivers      :", cmod%rec_num
    print*, "number of sources        :", cmod%src_num
    print*, "time frame               :", cmod%t_inp(1),cmod%t_inp(cmod%t_inp_len)
    print*, "temperatrue              :", cmod%T_int,cmod%T_end
    print*, "annealing coefficient    :", cmod%ana_cof
    print*, "intensity range factor    :", cmod%intensity_factor
    print*, "smooth_constraint_weight    :", cmod%smooth_weight_a,cmod%smooth_weight_t
    print*, "entropy_constraint_weight    :", cmod%entropy_weight_a,cmod%entropy_weight_t
    print*, "single_constraint_weight    :", cmod%single_weight
    print*, "L1_constraint_weight    :", cmod%L1_weight
    print*, "==================================building field file..."
    time1=timef()
    allocate(cmod%t_sgl(cmod%t_res_len),cmod%fed_all(cmod%fault_num))
    allocate(all_Normal_Factor(cmod%fault_num),cmod%Normal_Factor(cmod%fault_num))
    allocate(cmod%pulse_current(cmod%fault_num))
    do fn=1,cmod%fault_num
        allocate(cmod%fed_all(fn)%val(cmod%rec_num*3,cmod%t_res_len,cmod%src_num(fn),cmod%degree_range))
        do m=1,cmod%degree_range
            if (m<10) then
                 write(deg_num,'(I1)') m
            else if (m<100) then
                 write(deg_num,'(I2)') m
            end if
            do i=1,cmod%src_num(fn)
                if (i<10) then
                    write(fed_num,'(I1)') i
                else if (i<100) then
                    write(fed_num,'(I2)') i
                else
                    write(fed_num,'(I3)') i
                end if
                fed_file = "forward_data/"//trim(cmod%forward_name(fn))//"/sph"//trim(deg_num)//"/sph"//trim(deg_num)//"_"//trim(fed_num)//".field"
                open(fed_unit,file=trim(fed_file))
                read(fed_unit,*) step
                do j=1,cmod%t_res_len
                    read(fed_unit,*) step,cmod%t_sgl(j),(cmod%fed_all(fn)%val(k,j,i,m),k=1,cmod%rec_num*3)
                end do
                close(fed_unit)
                fed_unit=fed_unit+1
            end do
        end do
        print*,"---------fault:",fn
        do i=1,cmod%degree_range
            print*,"feild",i,"info:"
            write(*,"(2es15.6)") cmod%fed_all(fn)%val(1,1,1,i),cmod%fed_all(fn)%val(cmod%rec_num*3,cmod%t_res_len,cmod%src_num(fn),i)
        end do
        !cmod%fed_all
        allocate(all_Normal_Factor(fn)%val(cmod%degree_range,cmod%src_num(fn)),cmod%Normal_Factor(fn)%val(cmod%src_num(fn)))
        do m=1,cmod%degree_range
            do i=1,cmod%src_num(fn)
                all_Normal_Factor(fn)%val(m,i)=1.d-20/maxval(abs(cmod%fed_all(fn)%val(1,:,i,m)))
            end do
        end do
        do i=1,cmod%src_num(fn)
            cmod%Normal_Factor(fn)%val(i)=sum(all_Normal_Factor(fn)%val(:,i))/cmod%degree_range
        end do
        do i=1,cmod%x_num(fn)
            cmod%Normal_Factor(fn)%val(1+(i-1)*cmod%y_num(fn):i*cmod%y_num(fn))=sum(cmod%Normal_Factor(fn)%val(1+(i-1)*cmod%y_num(fn):i*cmod%y_num(fn)))/cmod%y_num(fn)
        end do
        do i=1,cmod%src_num(fn)
            cmod%fed_all(fn)%val(:,:,i,:)=cmod%fed_all(fn)%val(:,:,i,:)*cmod%Normal_Factor(fn)%val(i)
        end do
        allocate(cmod%pulse_current(fn)%val(cmod%rec_num*3*cmod%src_num(fn)*cmod%degree_range))
        k=1
        do m=1,cmod%degree_range
            do i=1,cmod%src_num(fn)
                do j=1,cmod%rec_num*3
                    cmod%pulse_current(fn)%val(k)=maxval(abs(cmod%fed_all(fn)%val(j,:,i,m)))
                    k=k+1
                end do
            end do
        end do
    end do
    time2=timef()
    print*, "time info",cmod%t_sgl(1),cmod%t_sgl(cmod%t_res_len)
    print*, "Time passes:",(time2-time1),"s"
    do i = 1, cmod%fault_num
      if (associated(all_Normal_Factor(i)%val)) then
        deallocate(all_Normal_Factor(i)%val)
      end if
    end do
    if (allocated(all_Normal_Factor)) deallocate(all_Normal_Factor)
    
    end subroutine read_ctl

    
    subroutine simulated_annealing_t(cmod)
    implicit none
    type(model),intent(inout) :: cmod
    real(8) :: fed_true(cmod%rec_num*3,cmod%t_inp_len),fed_ini(cmod%rec_num*3,cmod%t_inp_len),fed_new(cmod%rec_num*3,cmod%t_inp_len),&
        &fed_mid(cmod%rec_num*3,cmod%t_inp_len),res(cmod%rec_num*3,cmod%t_inp_len)
    real(8) :: fed_interp(cmod%rec_num*3,cmod%t_inp_len,cmod%fault_num),fed_interp_new(cmod%rec_num*3,cmod%t_inp_len,cmod%fault_num)
    type(real1d_ptr),allocatable :: A_ini(:),time_ini(:),A_new(:),time_new(:),t_best(:),A_best(:)
    type(integer_1d_ptr),allocatable :: degree_ini(:),degree_new(:),degree_best(:)
    real(8) :: E_best,E_cur,E_new,E_mid,temp,u,A_low,A_up(cmod%fault_num),t_low,t_up,e1,e2,error
    integer :: i,j,stop_num,d1,fn,k,sn
    real(8) :: time1,time2
    character(40) :: step_file,step_num
    logical(kind=4) :: ierr1,errnum
    integer(kind=4) :: istatus1
    character(len=256) :: cmd
    print*, "==================================creating a new output folder..."
    inquire(DIRECTORY="inv_result/"//trim(cmod%file_root)//"_result", EXIST=ierr1)
    if(ierr1) then
       print*,'The directory have existed and not been needed to create'
    else
       print*,'The directory not exist and creat it'
       cmd = "mkdir -p inv_result/" // trim(cmod%file_root) // "_result"
       print*,trim(cmd)
       istatus1=system(trim(cmd))
       if(istatus1==-1) then
           errnum=ierrno()
           print*,'Error=',errnum,' inquire the Intel Visual Fortran help document'
           stop ' Folder creating is fail'
       end if
    end if
    open(92,file="forward_data/"//trim(cmod%inv_data)//"/field_true.field")
    do i=1,cmod%t_inp_len
        read(92,*) (fed_true(j,i),j=1,cmod%rec_num*3)
    end do
    close(92)
    print*, "--- file information ..."
    print*, "Reverse event:", "forward_data/"//trim(cmod%inv_data)//"/field_true.field"
    !print*,fed_true(1:3,cmod%t_inp_len)
    !search range
    allocate(A_ini(cmod%fault_num),time_ini(cmod%fault_num),A_new(cmod%fault_num),time_new(cmod%fault_num),&
        &t_best(cmod%fault_num),A_best(cmod%fault_num),degree_ini(cmod%fault_num),degree_new(cmod%fault_num),degree_best(cmod%fault_num)) 
    A_low=0.d0
    do fn=1,cmod%fault_num
        A_up(fn)=maxval(abs(fed_true))/minval(cmod%pulse_current(fn)%val)/cmod%intensity_factor
    end do
    t_low=0.d0
    t_up=maxval(cmod%t_inp)
    call random_seed ()
    do fn=1,cmod%fault_num
        allocate(A_ini(fn)%val(cmod%src_num(fn)),time_ini(fn)%val(cmod%src_num(fn)))
        allocate(A_new(fn)%val(cmod%src_num(fn)),time_new(fn)%val(cmod%src_num(fn)))
        allocate(degree_ini(fn)%val(cmod%src_num(fn)),degree_new(fn)%val(cmod%src_num(fn)),degree_best(fn)%val(cmod%src_num(fn)))
        allocate(t_best(fn)%val(cmod%src_num(fn)),A_best(fn)%val(cmod%src_num(fn)))
        call random_number (A_ini(fn)%val)
        call random_number (time_ini(fn)%val)
        do i=1,cmod%src_num(fn)
            degree_ini(fn)%val(i)=rand_int(cmod%degree_range)
        end do
        A_ini(fn)%val=(A_up(fn)-A_low)*A_ini(fn)%val+A_low
        time_ini(fn)%val=(t_up-t_low)*time_ini(fn)%val+t_low
        !write(cmd,'(I1)') fn
        !open(87,file="inv_result/"//trim(cmod%file_root)//"_result/src_init"//trim(cmd)//".txt")
        !do i=1,cmod%src_num(fn)
        !    write(87,"(2es15.6,i6)") a_ini(fn)%val(i)*cmod%Normal_Factor(fn)%val(i),time_ini(fn)%val(i),degree_ini(fn)%val(i)
        !end do
        !close(87)
        !open(777,file="inv_result/"//trim(cmod%file_root)//"_result/factor"//trim(cmd)//".txt")
        !do i=1,cmod%src_num(fn)
        !    write(777,"(es15.6)") cmod%Normal_Factor(fn)%val(i)*A_up(fn)
        !end do
        !close(777)
    end do
    allocate(cmod%w1(cmod%fault_num),cmod%w2(cmod%fault_num),cmod%w3(cmod%fault_num),&
        &cmod%w4(cmod%fault_num),cmod%w5(cmod%fault_num),cmod%w6(cmod%fault_num))
    fed_ini=0.d0
    do fn=1,cmod%fault_num
        call fed_interp_parallel(cmod,A_ini(fn)%val,time_ini(fn)%val,degree_ini(fn)%val,fed_mid,fn)
        fed_ini=fed_ini+fed_mid
        fed_interp(:,:,fn)=fed_mid
    end do
    !open(87,file="inv_result/"//trim(cmod%file_root)//"_result/field_init.field")
    !do i=1,cmod%t_inp_len
    !    write(87,"("//itoa(cmod%rec_num*3)//"e15.5)") (fed_ini(j,i),j=1,cmod%rec_num*3)
    !end do
    !close(87)
    res=fed_ini-fed_true
    E_cur=0.d0
    do fn=1,cmod%fault_num
        call loss_func(cmod,res,A_ini(fn)%val,time_ini(fn)%val,degree_ini(fn)%val,E_mid,.true.,fn)
        E_cur=E_cur+e_mid
    end do
    E_cur=E_cur/cmod%fault_num
    print*, "==================================begining simulated annealing (time)..."
    print*, 'the searching range of A:',A_low,'-(',A_up,')'
    print*, 'the searching range of time:',t_low,'-',t_up
    print*, 'initial error:',E_cur
    temp=cmod%T_int
    e1=1.d0
    e2=1.d0
    error=0.d0
    stop_num=0
    k=1
    e_best=e_cur
    open(88,file="inv_result/"//trim(cmod%file_root)//"_result/error.txt")
    do while (temp>cmod%T_end)
        time1=timef()
        do i=1,sum(cmod%src_num)
            do j=1,cmod%move
                fed_interp_new=fed_interp
                if (i<=cmod%src_num(1)) then
                    do fn=1,cmod%fault_num
                        degree_new(fn)%val=degree_ini(fn)%val
                        a_new(fn)%val=a_ini(fn)%val
                        time_new(fn)%val=time_ini(fn)%val
                    end do
                    call rand_par(A_ini(1)%val,A_new(1)%val,i,A_low,A_up(1))
                    call rand_par(time_ini(1)%val,time_new(1)%val,i,t_low,t_up) 
                    d1=rand_int(cmod%degree_range)
                    do while (degree_new(1)%val(i)==d1)
                        d1=rand_int(cmod%degree_range)
                    end do
                    degree_new(1)%val(i)=d1
                    call fed_interp_parallel(cmod,A_new(1)%val,time_new(1)%val,degree_new(1)%val,fed_mid,1)
                    fed_interp_new(:,:,1)=fed_mid
                else
                    do fn=1,cmod%fault_num
                        degree_new(fn)%val=degree_ini(fn)%val
                        a_new(fn)%val=a_ini(fn)%val
                        time_new(fn)%val=time_ini(fn)%val
                    end do
                    call rand_par(A_ini(cmod%fault_num)%val,A_new(cmod%fault_num)%val,i-cmod%src_num(1),A_low,A_up(cmod%fault_num))
                    call rand_par(time_ini(cmod%fault_num)%val,time_new(cmod%fault_num)%val,i-cmod%src_num(1),t_low,t_up)
                    d1=rand_int(cmod%degree_range)
                    do while (degree_new(cmod%fault_num)%val(i-cmod%src_num(1))==d1)
                        d1=rand_int(cmod%degree_range)
                    end do
                    degree_new(cmod%fault_num)%val(i-cmod%src_num(1))=d1
                    call fed_interp_parallel(cmod,A_new(cmod%fault_num)%val,time_new(cmod%fault_num)%val,degree_new(cmod%fault_num)%val,fed_mid,cmod%fault_num)
                    fed_interp_new(:,:,cmod%fault_num)=fed_mid
                end if
                !caculation loss function
                fed_new=0.d0
                do fn=1,cmod%fault_num
                    fed_new=fed_new+fed_interp_new(:,:,fn)
                end do
                res=fed_new-fed_true
                E_new=0.d0
                do fn=1,cmod%fault_num
                    call loss_func(cmod,res,A_new(fn)%val,time_new(fn)%val,degree_new(fn)%val,E_mid,.false.,fn)
                    E_new=E_new+e_mid
                end do
                E_new=E_new/cmod%fault_num
                u=rand0()
                if (E_cur>E_new) then
                    do fn=1,cmod%fault_num
                        A_ini(fn)%val=A_new(fn)%val
                        time_ini(fn)%val=time_new(fn)%val
                        degree_ini(fn)%val=degree_new(fn)%val
                    end do
                    fed_interp=fed_interp_new
                    E_cur=E_new
                else if (exp(-(E_new-E_cur)/temp)>u) then
                    do fn=1,cmod%fault_num
                        A_ini(fn)%val=A_new(fn)%val
                        time_ini(fn)%val=time_new(fn)%val
                        degree_ini(fn)%val=degree_new(fn)%val
                    end do
                    fed_interp=fed_interp_new
                    E_cur=E_new
                end if
                if (e_best>e_cur) then
                    do fn=1,cmod%fault_num
                        a_best(fn)%val=A_ini(fn)%val
                        t_best(fn)%val=time_ini(fn)%val
                        degree_best(fn)%val=degree_ini(fn)%val
                    end do
                    e_best=e_cur
                end if
            end do
        end do
        time2=timef()
        temp=temp*cmod%ana_cof
        write(88,"(2es15.6)") temp,e_cur
        print*, 'step:',k,e_cur,e_best,'time:',time2-time1,'temp:',temp
        if (mod(k,50)==0) then
            if (k<100) then
                write(step_num,'(I2)') k
            else if (k<1000) then
                write(step_num,'(I3)') k
            else
                write(step_num,'(I4)') k
            end if
            fed_new=0.d0
            do fn=1,cmod%fault_num
                write(cmd,'(I1)') fn
                step_file="inv_result/"//trim(cmod%file_root)//"_result/step"//trim(step_num)//"_f"//trim(cmd)//".txt"
                open(90,file=trim(step_file))
                do i=1,cmod%src_num(fn)
                    write(90,"(2es15.6,i6)") a_best(fn)%val(i)*cmod%Normal_Factor(fn)%val(i),t_best(fn)%val(i),degree_best(fn)%val(i)
                end do
                close(90)
                call fed_interp_parallel(cmod,a_best(fn)%val,t_best(fn)%val,degree_best(fn)%val,fed_mid,fn)
                fed_new=fed_new+fed_mid
            end do
            !output field  
            step_file="inv_result/"//trim(cmod%file_root)//"_result/step"//trim(step_num)//".field"
            open(999,file=trim(step_file))
            do i=1,cmod%t_inp_len
                write(999,"("//itoa(cmod%rec_num*3)//"e15.5)") (fed_new(j,i),j=1,cmod%rec_num*3)
            end do
            close(999)
        end if
        k=k+1
        !Error judgment
        e2=e1
        e1=e_cur
        error=abs(e1-e2)/e2
        if (error<0.005) then
            stop_num=stop_num+1
        else
            stop_num=0
        end if
        if (stop_num>50) exit
    end do
    close(88)
    print*, "fin_error",e_best
    fed_new=0.d0
    do fn=1,cmod%fault_num
        write(cmd,'(I1)') fn
        open(89,file="inv_result/"//trim(cmod%file_root)//"_result/src_inv_f"//trim(cmd)//".txt")
        do i=1,cmod%src_num(fn)
            write(89,"(2es15.6,i6)") a_best(fn)%val(i)*cmod%Normal_Factor(fn)%val(i),t_best(fn)%val(i),degree_best(fn)%val(i)
        end do
        close(89)
        call fed_interp_parallel(cmod,a_best(fn)%val,t_best(fn)%val,degree_best(fn)%val,fed_mid,fn)
        fed_new=fed_new+fed_mid
    end do
    open(91,file="inv_result/"//trim(cmod%file_root)//"_result/field_inv.field")
    do i=1,cmod%t_inp_len
        write(91,"("//itoa(cmod%rec_num*3)//"e15.5)") (fed_new(j,i),j=1,cmod%rec_num*3)
    end do
    close(91)
    do i = 1, cmod%fault_num
      if (associated(A_ini(i)%val)) then
        deallocate(A_ini(i)%val)
      end if
      if (associated(A_new(i)%val)) then
        deallocate(A_new(i)%val)
      end if
      if (associated(A_best(i)%val)) then
        deallocate(A_best(i)%val)
      end if
      if (associated(time_ini(i)%val)) then
        deallocate(time_ini(i)%val)
      end if
      if (associated(time_new(i)%val)) then
        deallocate(time_new(i)%val)
      end if
      if (associated(t_best(i)%val)) then
        deallocate(t_best(i)%val)
      end if
      if (associated(degree_ini(i)%val)) then
        deallocate(degree_ini(i)%val)
      end if
      if (associated(degree_new(i)%val)) then
        deallocate(degree_new(i)%val)
      end if
      if (associated(degree_best(i)%val)) then
        deallocate(degree_best(i)%val)
      end if
    end do
    if (allocated(A_ini)) deallocate(A_ini)
    if (allocated(A_new)) deallocate(A_new)
    if (allocated(A_best)) deallocate(A_best)
    if (allocated(time_ini)) deallocate(time_ini)
    if (allocated(time_new)) deallocate(time_new)
    if (allocated(t_best)) deallocate(t_best)
    if (allocated(degree_ini)) deallocate(degree_ini)
    if (allocated(degree_new)) deallocate(degree_new)
    if (allocated(degree_best)) deallocate(degree_best)
    end subroutine simulated_annealing_t
    

    subroutine loss_func(cmod,res,A,time,degree,loss,cal_weight,fault_num)
    implicit none
    type(model),intent(inout) :: cmod
    real(8),intent(in) :: A(:),time(:),res(:,:)
    integer,intent(in) :: degree(:),fault_num
    logical,intent(in) :: cal_weight
    real(8),intent(out) :: loss
    real(8) :: A_2d(cmod%y_num(fault_num),cmod%x_num(fault_num)),time_2d(cmod%y_num(fault_num),cmod%x_num(fault_num)),&
        &Ax_del(cmod%y_num(fault_num)*cmod%x_num(fault_num)),Ay_del(cmod%y_num(fault_num)*cmod%x_num(fault_num)),&
        &time_del(cmod%y_num(fault_num)*cmod%x_num(fault_num)),n_res(cmod%rec_num*3)
    real(8) :: Ax(cmod%y_num(fault_num)*cmod%x_num(fault_num)),Ay(cmod%y_num(fault_num)*cmod%x_num(fault_num)),&
        &Ax_2d(cmod%y_num(fault_num),cmod%x_num(fault_num)),Ay_2d(cmod%y_num(fault_num),cmod%x_num(fault_num)),&
        &A_true(cmod%y_num(fault_num)*cmod%x_num(fault_num))
    real(8) :: n1,n2,n_a,n_time,t_len,t_lq,n_ax,n_ay,m_sum,m_single
    real(8) :: aa(cmod%y_num(fault_num)*cmod%x_num(fault_num),2),cluster_centre(2),nn,delta,Q_a,Q_t,S_me0_a,S_me0_t
    integer :: i
    loss=0.d0
    t_len=size(cmod%t_inp)
    t_lq=sqrt(t_len)
    if (.not. allocated(cmod%weight)) allocate(cmod%weight(cmod%rec_num*3))
    do i=1,cmod%y_num(fault_num)*cmod%x_num(fault_num)
        Ax(i)=A(i)*cos(cmod%degree_val(degree(i))*pi/180)*cmod%Normal_Factor(fault_num)%val(i)
        Ay(i)=A(i)*sin(cmod%degree_val(degree(i))*pi/180)*cmod%Normal_Factor(fault_num)%val(i)
    end do
    Ax_2d=reshape(Ax,[cmod%y_num(fault_num),cmod%x_num(fault_num)])
    Ay_2d=reshape(Ay,[cmod%y_num(fault_num),cmod%x_num(fault_num)])
    time_2d=reshape(time,[cmod%y_num(fault_num),cmod%x_num(fault_num)])
    call del2(Ax_2d,Ax_del,cmod%y_num(fault_num),cmod%x_num(fault_num))
    call del2(Ay_2d,Ay_del,cmod%y_num(fault_num),cmod%x_num(fault_num))
    call del2(time_2d,time_del,cmod%y_num(fault_num),cmod%x_num(fault_num))
    do i=1,cmod%rec_num*3
        call norm(res(i,:),n1,1)
        call norm(res(i,:),n2,2)
        n_res(i)=n1/cmod%t_inp_len+n2/t_lq
        if (cal_weight) then
            cmod%weight(i)=1.d0/n_res(i)
        end if
        loss=loss+cmod%weight(i)*n_res(i)
    end do
    loss=loss/(cmod%rec_num*3)
    call norm(Ax_del,n_Ax,4)
    call norm(Ay_del,n_Ay,4)
    n_a=sqrt((n_ax+n_ay)/(cmod%y_num(fault_num)*cmod%x_num(fault_num)))
    call norm(time_del,n_time,3)
    if (cal_weight) cmod%w1(fault_num)=loss*cmod%smooth_weight_a/(n_A)
    if (cal_weight) cmod%w2(fault_num)=loss*cmod%smooth_weight_t/(n_time)
    !Minimum entropy constraint
    A_true=A*cmod%Normal_Factor(fault_num)%val
    delta=1.d-7
    Q_a=sum(abs(A_true))+cmod%src_num(fault_num)*delta
    S_me0_a=0.d0
    do i=1,cmod%src_num(fault_num)
        S_me0_a=S_me0_a+(abs(A_true(i))+delta)/Q_a*log((abs(A_true(i))+delta)/Q_a)
    end do
    if (cal_weight) cmod%w3(fault_num)=loss*cmod%entropy_weight_a/(S_me0_a)
    Q_t=sum(time)+cmod%src_num(fault_num)*delta
    S_me0_t=0.d0
    do i=1,cmod%src_num(fault_num)
        S_me0_t=S_me0_t+(time(i)+delta)/Q_t*log((time(i)+delta)/Q_t)
    end do
    if (cal_weight) cmod%w6(fault_num)=loss*cmod%entropy_weight_t/(S_me0_t)
    !L1 constraint
    m_sum=sum(abs(A_true))
    if (cal_weight) cmod%w4(fault_num)=loss*cmod%L1_weight/(m_sum)
    !single-cluster constraint
    aa(:,1)=abs(A_true)
    aa(:,2)=abs(A_true)
    cluster_centre=SUM(cmod%src_centre(fault_num)%val*aa, dim=1)/sum(abs(A_true))
    m_single=0
    do i=1,cmod%src_num(fault_num)
        call norm(cmod%src_centre(fault_num)%val(i,:)-cluster_centre,nn,2)
        m_single=m_single+abs(A_true(i))*(nn**2)
    end do
    if (cal_weight) cmod%w5(fault_num)=loss*cmod%single_weight/(m_single)
    loss=loss+cmod%w1(fault_num)*n_A+cmod%w2(fault_num)*n_time+cmod%w3(fault_num)*S_me0_a+cmod%w6(fault_num)*S_me0_t+cmod%w4(fault_num)*m_sum+cmod%w5(fault_num)*m_single
    !print*,n_A,n_time,S_me0_a,S_me0_t,m_sum,m_single
    end subroutine loss_func

    !subroutine fed_interp(cmod,A,time,degree,fed_out,fault_num)
    !USE IFPORT
    !implicit none
    !type(model),intent(in) :: cmod
    !real(8),intent(in) :: A(:),time(:)
    !integer,intent(in) :: degree(:),fault_num
    !real(8),intent(out) :: fed_out(:,:)
    !real(8) :: y_out(3*cmod%rec_num,cmod%t_inp_len),t(size(cmod%t_sgl))
    !real(8) :: w
    !integer :: i,isrc,j,IND
    !real(8) :: x_log1(size(cmod%t_sgl)),x_log2(size(cmod%t_inp))
    !x_log2 = log10(cmod%t_inp)
    !fed_out = 0.d0
    !!!$OMP PARALLEL DO private(i,isrc,IND,W,t,x_log1)
    !do isrc=1,size(A)
    !    t = cmod%t_sgl + time(isrc)
    !    x_log1 = log10(t)
    !    DO I = 1,size(cmod%t_inp)
    !        CALL locate (x_log1,size(t),x_log2(I),IND)
    !        IF (IND==0 .or. ind==size(t)) THEN
    !            y_out(:,i)=0.d0
    !        ELSE
    !            w = (x_log2(i)-x_log1(IND)) / (x_log1(IND+1)-x_log1(IND))
    !            y_out(:,i) = ( (1.D0-w) * cmod%fed_all(fault_num)%val(:,IND,isrc,degree(isrc)) + w * cmod%fed_all(fault_num)%val(:,IND+1,isrc,degree(isrc)) ) * A(isrc)
    !        END IF
    !    END DO
    !    fed_out = fed_out + y_out
    !end do
    !!!$OMP END PARALLEL DO
    !end subroutine fed_interp
    
    subroutine fed_interp_parallel(cmod,A,time,degree,fed_out,fault_num)
    USE IFPORT
    implicit none
    type(model),intent(in) :: cmod
    real(8),intent(in) :: A(:),time(:)
    integer,intent(in) :: degree(:),fault_num
    real(8),intent(out) :: fed_out(:,:)
    real(8) :: y_out(3*cmod%rec_num,cmod%t_inp_len),t(size(cmod%t_sgl))
    real(8) :: w
    integer :: i,isrc,j,IND
    real(8) :: x_log1(size(cmod%t_sgl)),x_log2(size(cmod%t_inp))
    x_log2 = log10(cmod%t_inp)
    fed_out = 0.d0
    !!$OMP PARALLEL DO private(i,isrc,IND,W,t,x_log1)
    do isrc=1,size(A)
        t = cmod%t_sgl + time(isrc)
        x_log1 = log10(t)
        !$OMP PARALLEL DO private(i,IND,W)
        DO I = 1,size(cmod%t_inp)
            CALL locate (x_log1,size(t),x_log2(I),IND)
            IF (IND==0 .or. ind==size(t)) THEN
                y_out(:,i)=0.d0
            ELSE
                w = (x_log2(i)-x_log1(IND)) / (x_log1(IND+1)-x_log1(IND))
                y_out(:,i) = ( (1.D0-w) * cmod%fed_all(fault_num)%val(:,IND,isrc,degree(isrc)) + w * cmod%fed_all(fault_num)%val(:,IND+1,isrc,degree(isrc)) ) * A(isrc)
            END IF
            fed_out(:,i) = fed_out(:,i) + y_out(:,i)
        END DO
        !$OMP END PARALLEL DO
    end do
    !!$OMP END PARALLEL DO
    end subroutine fed_interp_parallel



    SUBROUTINE locate(xx,n,x,j)
    INTEGER :: j,n
    REAL(8) :: x,xx(n)
    !Given an array xx(1:n), and given a value x, returns a value j such that x is between
    !xx(j) and xx(j+1). xx(1:n) must be monotonic, either increasing or decreasing. j=0
    !or j=n is returned to indicate that x is out of range.
    INTEGER :: jl,jm,ju
    jl=0 !Initialize lower
    ju=n+1 !and upper limits.
10  if(ju-jl.gt.1)then !If we are not yet done,
        jm=(ju+jl)/2 !compute a midpoint,
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
            jl=jm !and replace either the lower limit
        else
            ju=jm !or the upper limit, as appropriate.
        endif
        goto 10 !Repeat until
    end if !the test condition 10 is satisfied.
    if(x.eq.xx(1))then !Then set the output
        j=1
    else if(x.eq.xx(n))then
        j=n-1
    else
        j=jl
    endif
    return
    END

    subroutine del2(x,y_out,row,column)
    implicit none
    real(8),intent(in) :: x(:,:)
    integer,intent(in) :: row,column
    real(8),intent(out) :: y_out(:)
    real(8),allocatable :: x_new(:,:),y(:,:)
    integer :: i,j
    y_out=0.d0
    allocate(x_new(row+2,column+2),y(row,column))
    x_new=0.d0
    y=0.d0
    do i=1,row
        do j=1,column
            x_new(i+1,j+1)=x(i,j)
        end do
    end do
    do i=2,row+1
        do j=2,column+1
            y(i-1,j-1)=x_new(i,j)-0.25*(x_new(i-1,j)+x_new(i+1,j)+x_new(i,j-1)+x_new(i,j+1))
        end do
    end do
    y_out=reshape(y,[row*column])
    end subroutine del2

    function rand0()   !returns random number between 0 - 1
    implicit none
    integer , save :: flag=0
    double precision :: rand0
    !print*,flag
    if(flag==0) then
        call random_seed()
        flag = 1
    endif
    call random_number(rand0)     ! built in fortran 90 random number function
    end function rand0

    subroutine normal(mean,sigma,randn)
    implicit none
    real(8),intent(in) :: mean,sigma
    real(8),intent(out) :: randn
    !integer ,save :: flag
    real(8),parameter :: pi  = 3.141592653589793238462643383279502884197_8
    real(8) :: x1,x2,y
    x1=rand0();x2=rand0();
    y=sqrt(-2.0d0*log(x1))*cos(2.0d0*pi*x2)
    randn=mean+sigma*y
    end subroutine normal

    function rand_int(i)!Generate random integers of 1-i
    implicit none
    integer :: rand_int,i
    real(8) :: r
    r=rand0()
    rand_int=1+int(r*i)
    end function rand_int

    subroutine norm(x,y,p)
    implicit none
    real(8),intent(in) :: x(:)
    integer,intent(in) :: p
    real(8),intent(out) :: y
    if (p==1) then
        y=sum(abs(x))
    else if (p==2) then
        y=sqrt(dot_product(x,x))
    else if (p==3) then
        y=sqrt(dot_product(x,x)/size(x))
    else if (p==4) then
        y=dot_product(x,x)
    end if
    end subroutine norm

    subroutine rand_par(par,par_new,num,low,up)
    implicit none
    real(8),intent(in) :: par(:),low,up
    integer,intent(in) :: num
    real(8),intent(out) :: par_new(:)
    real(8) :: r,u
    par_new=par
    call normal(par(num),(up-low)/2,par_new(num))
    if (par_new(num)>up) then
        u=rand0()
        r=u*(up-par(num))
        par_new(num)=par(num)+r
    elseif (par_new(num)<low) then
        u=rand0()
        r=u*(par(num)-low)
        par_new(num)=par(num)-r
    end if
    end subroutine rand_par

    subroutine parseCode( nLen, sLine, sCode, sValue, bComment )
    ! parse a line read from a file into a code & value. Force the code to be all lowercase with no ending colon. Terminate the line at a '%' or '!' sign (these allow for user comments!)
    implicit none
    ! Args
    integer, intent(in)   :: nLen
    character(nLen)       :: sLine
    character(nLen), intent(out) :: sCode, sValue
    logical, intent(out)    :: bComment
    ! Local vars
    integer :: iFrom, iTo
    ! Init returns
    bComment = .false.
    sCode = ' '
    sValue = ' '
    ! Convert all tab characters to spaces
    forall( iTo = 1:nLen, ichar(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '
    ! Skip any beginning blanks
    do iFrom = 1,nLen
        if (sLine(iFrom:iFrom) .ne. ' ') exit
    enddo
    ! If the first char is a comment char, then the whole line is a comment. Also, if the line is blank, consider it a comment.
    if (iFrom >= nLen) then !KWK may 2009 pulled this out in from since sometimes iFrom > nLen and this kills (iFrom:iFrom) below
        bComment = .true.
        return
    endif
    if(  sLine(iFrom:iFrom) == '%' .or. sLine(iFrom:iFrom) == '!' ) then
        bComment = .true.
        return
    endif
    ! Pull off the code value. Cvt to lowercase as we go.
    iTo = index(sLine,':') - 1
    if (iTo < iFrom) then
        !write(*,*) 'Parsing Error: missing colon in line:',sLine
        return
    endif
    sCode = sLine(iFrom:iTo)
    call Lower(sCode)
    ! Skip spaces after the colon
    do iFrom = iTo+2,nLen
        if (sLine(iFrom:iFrom) .ne. ' ') exit
    enddo
    ! Get the rest, up to any comment
    sValue = sLine(iFrom:)
    iTo = len_trim(sValue)
    iFrom = index(sValue,'%')
    if (iFrom > 0 .and. iFrom < iTo) then
        sValue(iFrom:iTo) = ' '
    endif
    iFrom = index(sValue,'!')
    if (iFrom > 0 .and. iFrom < iTo) then
        sValue(iFrom:iTo) = ' '
    endif
    !call Lower(sValue)   ! No: Some values are filenames which are case-sensitive on UNIX!
    end subroutine parseCode

    subroutine Lower(s)
    ! convert string to lower case
    character(*), intent(out)  :: s
    integer :: i
    do  i=1,len_trim(s)
        if  ( s(i:i) >= 'A' .and. s(i:i) <= 'Z' ) then
            s(i:i) = char(ichar(s(i:i)) + 32)
        endif
    enddo
    end subroutine Lower

    function itoa(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i0)') i
    res = trim(tmp)
    end function

    subroutine clear_model (cmod)
    implicit none
    type(model),intent(inout) :: cmod
    integer :: i
    if (allocated(cmod%degree_val)) deallocate(cmod%degree_val)
    if (allocated(cmod%forward_name)) deallocate(cmod%forward_name)
    if (allocated(cmod%t_inp)) deallocate(cmod%t_inp)
    if (allocated(cmod%t_sgl)) deallocate(cmod%t_sgl)
    if (allocated(cmod%weight)) deallocate(cmod%weight)
    if (allocated(cmod%x_num)) deallocate(cmod%x_num)
    if (allocated(cmod%y_num)) deallocate(cmod%y_num)
    if (allocated(cmod%src_num)) deallocate(cmod%src_num)
    if (allocated(cmod%w1)) deallocate(cmod%w1)
    if (allocated(cmod%w2)) deallocate(cmod%w2)
    if (allocated(cmod%w3)) deallocate(cmod%w3)
    if (allocated(cmod%w4)) deallocate(cmod%w4)
    if (allocated(cmod%w5)) deallocate(cmod%w5)
    if (allocated(cmod%w6)) deallocate(cmod%w6)
    
    do i = 1, cmod%fault_num
      if (associated(cmod%src_centre(i)%val)) then
        deallocate(cmod%src_centre(i)%val)
      end if
      if (associated(cmod%fed_all(i)%val)) then
        deallocate(cmod%fed_all(i)%val)
      end if
      if (associated(cmod%Normal_Factor(i)%val)) then
        deallocate(cmod%Normal_Factor(i)%val)
      end if
      if (associated(cmod%pulse_current(i)%val)) then
        deallocate(cmod%pulse_current(i)%val)
      end if
    end do
    if (allocated(cmod%src_centre)) deallocate(cmod%src_centre)
    if (allocated(cmod%fed_all)) deallocate(cmod%fed_all)
    if (allocated(cmod%Normal_Factor)) deallocate(cmod%Normal_Factor)
    if (allocated(cmod%pulse_current)) deallocate(cmod%pulse_current)
    end subroutine clear_model

    end module fwd_modeling

    program TEM_yu
    USE IFPORT
    use fwd_modeling
    implicit none
    type(model) :: fwd_mod
    
    call read_ctl(fwd_mod)
    call simulated_annealing_t(fwd_mod)
    call clear_model(fwd_mod)

    end program TEM_yu




