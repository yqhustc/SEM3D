    include "mkl_spblas.f90"
    include "mkl_service.f90"
    include 'mkl_pardiso.f90'
    include 'mkl_sparse_handle.f90'
    !include 'fwd_modeling.f90'
    module file_option
    use fwd_modeling
    USE IFCORE, only : commitqq
    USE IFPOSIX, only : PXFGETCWD
    use MKL_SPBLAS
    USE mkl_pardiso
    implicit none
    private
    real(8),parameter :: pi  = 3.141592653589793238462643383279502884197_8
    real(8),parameter :: mu = pi * 4.d-7
    real(8),parameter :: epsilon = 8.8541878176d-12
    integer,parameter :: ctl_unit = 101, mesh_unit = 102, opt_unit=103
    integer,parameter :: ind_tet(6,2)=reshape([1,1,1,2,4,3,2,3,4,3,2,4],[6,2])
    integer,parameter :: pla_tet(4,3)=reshape([1,1,2,4,2,3,3,5,4,5,6,6],[4,3])
    type(MATRIX_DESCR) :: DESCR
    INTEGER,ALLOCATABLE :: OPTS_UNITS(:)
    character(80) :: file_root
    ! pardiso variables
    TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt_par(:)
    INTEGER, ALLOCATABLE :: iparm_par(:)
    INTEGER :: idum_par(1)
    REAL(8) :: ddum_par(1)
    INTEGER :: maxfct_par, mnum_par, mtype_par, phase_par, n_par, nrhs_par, error_par, msglvl_par, nnz_par
    ! time variables
    real(8) :: time_frac,time_stepping
    real(8) :: time_3,time_4,time_1,time_2
    !每个点对应的坐标（num_nod,3) :: point
    !每条棱边两个端点对应的全局编号（num_edg,2) :: edge
    !每个单元对应四个点的全局编号（num_ele,4）:: tet
    !每个单元对应六条棱边的全局编号,不分正负（num_ele,6)  :: tet_edge
    !场源上的棱边 :: src_edge,src_edg
    !xxx_phy 元素的物理属性    
    type model
        character(20) :: ctl_file,mesh_file,opt_file
        !read_mesh info
        real(8),allocatable :: point(:,:) 
        integer,allocatable :: tet(:,:),edge(:,:),tet_edge(:,:)
        integer,allocatable :: src_phy(:),tet_phy(:),tet_drc(:,:)
        integer,allocatable :: src_edge(:,:),src_edg(:)
        integer :: nod_num,seg_num,tet_num,edg_num
        ! output variables: ex ey ez dx dy dz bx by bz
        logical :: o_f(9)
        integer,allocatable :: o_v(:)
        ! receiver
        integer :: rec_num
        real(8),allocatable :: rec_cor(:,:)
        integer,allocatable :: rec_tet(:)
        ! domain character
        integer :: dom_num
        integer,allocatable :: dom_ind(:)
        real(8),allocatable :: sig_3d(:),per_3d(:),mag_3d(:)
        ! time variables
        integer :: step_mthd,step_dvsi,step_dble,step_dblesize
        real(8) :: time_maxm
        ! source info from the control file
        integer :: src_num
        integer,allocatable :: src_pulsetype(:)
        real(8),allocatable :: src_pulsewidth(:),src_current(:),src_delay(:)
        integer,allocatable :: src_tp(:),src_count_edges(:)
        real(8),allocatable :: src_end_1(:,:),src_end_2(:,:)
        integer,allocatable :: src_drc(:)
        ! edge not at the boundary
        integer,allocatable :: nobond_edge(:)
        integer :: nbeg_num
        ! system matrix
        type(SPARSE_MATRIX_T) :: a_mkl,b_mkl,c_mkl
        !solution output
        real(8),allocatable :: last_b_f(:,:,:)
        !source return
        real(8) :: add
        logical :: back
    end type model
    public :: model
    public :: read_mesh,read_ctl,compute_boundary,compute_matrix,clear_model,time_loop_1d
    
contains
    
    subroutine set_file(cmod)
    implicit none
    type(model),intent(out) :: cmod
    call getarg(2,file_root)
    ! set file names
    cmod%ctl_file = trim(file_root)//".cntl"
    cmod%mesh_file = trim(file_root)//".mesh" 
    end subroutine set_file
    
    subroutine read_mesh(cmod)
    implicit none
    type(model),intent(inout) :: cmod
    character(len=200) :: buffer
    character(len=5) :: word
    integer :: ierr,i,j,k,fin
    integer,allocatable :: edge_1(:),edge_2(:),edg_phy(:,:),indices(:),value(:),&
        &edge_01(:),edge_02(:)
    print*, "==================================building mesh..."
    print*, "--- reading mesh file ..."
    call set_file(cmod)
    open(mesh_unit,file=trim(cmod%mesh_file))
    do
        read (mesh_unit,"(a)",iostat=ierr) buffer
        if (ierr /= 0) exit
        read (buffer,*) word
        if (word == "Verti") then
            read (mesh_unit,*) cmod%nod_num
            allocate(cmod%point(cmod%nod_num,3))
            do i=1,cmod%nod_num
                read (mesh_unit,*) cmod%point(i,1),cmod%point(i,2),cmod%point(i,3)
				cmod%point(i,1) = cmod%point(i,1)*1000
				cmod%point(i,2) = cmod%point(i,2)*1000
				cmod%point(i,3) = cmod%point(i,3)*1000
            end do
        else if (word == "Edges") then
            read (mesh_unit,*) cmod%seg_num
            allocate(cmod%src_edge(cmod%seg_num,2),cmod%src_phy(cmod%seg_num))
            do i=1,cmod%seg_num
                read (mesh_unit,*) cmod%src_edge(i,1),cmod%src_edge(i,2),cmod%src_phy(i)
            end do
        else if (word == "Tetra") then
            read (mesh_unit,*) cmod%tet_num
            allocate(cmod%tet(cmod%tet_num,4),cmod%tet_phy(cmod%tet_num))
            do i=1,cmod%tet_num
                read (mesh_unit,*) cmod%tet(i,1),cmod%tet(i,2),cmod%tet(i,3),cmod%tet(i,4),cmod%tet_phy(i)
            end do
        end if
    end do
    close(mesh_unit)
    print*, "--- rebuilding mesh edge ..."
    allocate(edge_1(cmod%tet_num*6),edge_2(cmod%tet_num*6),cmod%tet_drc(cmod%tet_num,6),&
        &edge_01(cmod%tet_num*6),edge_02(cmod%tet_num*6))
    k=0
    do i=1,cmod%tet_num
        do j=1,6
            k=k+1
            if (cmod%tet(i,ind_tet(j,1))>cmod%tet(i,ind_tet(j,2))) then
                cmod%tet_drc(i,j)=-1
                edge_1(k) = cmod%tet(i,ind_tet(j,2))
                edge_2(k) = cmod%tet(i,ind_tet(j,1))
            else
                cmod%tet_drc(i,j)=1
                edge_1(k) = cmod%tet(i,ind_tet(j,1))
                edge_2(k) = cmod%tet(i,ind_tet(j,2))
            end if  
        end do
    end do 
    edge_01=edge_1
    edge_02=edge_2
    call sort2d(edge_1,edge_2)
    call unique2d(edge_1,edge_2,size(edge_1))
    cmod%edg_num = size(edge_1)
    allocate(cmod%edge(cmod%edg_num,2))
    cmod%edge(:,1) = edge_1
    cmod%edge(:,2) = edge_2
    print*, "number of points:",cmod%nod_num
    print*, "number of edges :",cmod%edg_num
    print*, "number of tets  :",cmod%tet_num
    call unique1d(edge_1,indices,value)
    allocate(cmod%tet_edge(cmod%tet_num,6))
    !!$OMP PARALLEL DO private(i,j,k,fin)!变化的量
    do i=1,cmod%tet_num
        do j=1,6
            k=search(value,min(edge_01((i-1)*6+j),edge_02((i-1)*6+j)))
            if (k<size(value)) then
                fin=indices(k+1)-1
            else
                fin=size(edge_1)
            end if
           cmod%tet_edge(i,j)=indices(k)+search(edge_2(indices(k):fin),max(edge_01((i-1)*6+j),edge_02((i-1)*6+j)))-1
        end do
    end do
    !!$OMP end PARALLEL DO
    !search source edge
    allocate(cmod%src_edg(cmod%seg_num))
    do i=1,cmod%seg_num
        k=search(value,min(cmod%src_edge(i,1),cmod%src_edge(i,2)))
        if (k<size(value)) then
            fin=indices(k+1)-1
        else
            fin=size(edge_1)
        end if
        cmod%src_edg(i)=(indices(k)+search(edge_2(indices(k):fin),max(cmod%src_edge(i,1),cmod%src_edge(i,2)))-1)
    end do   
    !open(opt_unit,file=trim(cmod%opt_file))
    !do i=1,cmod%tet_num
    !    write(opt_unit,*) cmod%tet_edge(i,:)
    !end do
    !close(opt_unit)
    end subroutine read_mesh
    
    subroutine read_ctl(cmod)
    implicit none
    type(model),intent(inout) :: cmod
    real(8),allocatable :: sig_iso(:),per_iso(:),mag_iso(:)
    real(8),allocatable :: sig_ani(:,:),sig_ang(:,:),per_ani(:,:),per_ang(:,:),mag_ani(:,:),mag_ang(:,:)
    character(180) ::sLine, sCode, sValue
    logical ::bComment
    integer :: ierr,i,j,tmp(9)
    real(8) :: dx(3,3),dy(3,3),dz(3,3),direc(3)
    logical :: find_phy
    cmod%rec_num = 0
    cmod%dom_num = 0
    cmod%src_num = 0
    cmod%step_mthd=0
    cmod%step_dvsi = 0
    cmod%time_maxm =0.d0
    cmod%step_dble = 0
    cmod%o_f = .false.
    print*, "==================================building ctl file..."
    print*, "--- reading control file ..."
    open(ctl_unit,file=trim(cmod%ctl_file))
    do while (.true.)
        read( ctl_unit, '(A180)', iostat = ierr ) sLine
        if (ierr /= 0) exit
        call parseCode( len(sLine), sLine, sCode, sValue, bComment )
        if( bComment ) cycle
        select case (trim(sCode))
        case ('receivers')
            read(sValue,*) cmod%rec_num
            allocate(cmod%rec_cor(cmod%rec_num,3))
            do i=1,cmod%rec_num
                read (ctl_unit,*) cmod%rec_cor(i,1),cmod%rec_cor(i,2),cmod%rec_cor(i,3)
				cmod%rec_cor(i,1) = cmod%rec_cor(i,1)*1000
				cmod%rec_cor(i,2) = cmod%rec_cor(i,2)*1000
				cmod%rec_cor(i,3) = cmod%rec_cor(i,3)*1000
            end do
    
        case ('isotropic')
            read(sValue,*) cmod%dom_num
            allocate(cmod%dom_ind(cmod%dom_num))
            allocate(sig_iso(cmod%dom_num),per_iso(cmod%dom_num),mag_iso(cmod%dom_num))
            do i=1,cmod%dom_num
                read (ctl_unit,*) cmod%dom_ind(i),sig_iso(i),per_iso(i),mag_iso(i)
            end do
            allocate(cmod%sig_3d(cmod%tet_num),cmod%per_3d(cmod%tet_num),cmod%mag_3d(cmod%tet_num))
            do i=1,cmod%tet_num
                find_phy = .false.
                do j=1,cmod%dom_num
                    if (cmod%tet_phy(i)==cmod%dom_ind(j)) then
                        find_phy = .true.
                        cmod%sig_3d(i) = sig_iso(j)
                        cmod%per_3d(i) = per_iso(j) * epsilon
                        cmod%mag_3d(i) = mag_iso(j) * mu
                    end if
                end do
                if (.not.find_phy) stop("error - physical domain not found!")
            end do
        !case ('boundary condition')
        !    read(sValue,*) cmod%bound_cond
        case ('source')
            !source: ind source_type(1 electric, 2 magnetic) pulse_type (1 fenoglio 2 Gaussian pulse) t0 current
            !   1    1    1.0e-7    1.0
            !   51
            !   x_start y_start z_start
            !   x_end y_end z_end
            !   ... ...
            read(sValue,*) cmod%src_num

            allocate(cmod%src_tp(cmod%src_num),cmod%src_end_1(cmod%src_num,3),cmod%src_end_2(cmod%src_num,3),cmod%src_pulsewidth(cmod%src_num),cmod%src_current(cmod%src_num),cmod%src_delay(cmod%src_num),cmod%src_pulsetype(cmod%src_num))
            do i=1,cmod%src_num
                read (ctl_unit,*) cmod%src_tp(i),cmod%src_pulsetype(i),cmod%src_delay(i),cmod%src_pulsewidth(i),cmod%src_current(i),&
                    cmod%src_end_1(i,1),cmod%src_end_1(i,2),cmod%src_end_1(i,3),cmod%src_end_2(i,1),cmod%src_end_2(i,2),cmod%src_end_2(i,3)
            end do
            cmod%src_end_1=cmod%src_end_1*1000
            cmod%src_end_2=cmod%src_end_2*1000
            ! src_count_edges
            allocate(cmod%src_count_edges(cmod%src_num))
            do i=1,cmod%src_num
                cmod%src_count_edges(i) = count(cmod%src_phy==cmod%src_tp(i))
            end do
            ! src_drc
            print*, "--- rebuild source info ..."
            allocate(cmod%src_drc(cmod%seg_num))
            do i=1,cmod%seg_num
                direc(1) = cmod%point(cmod%src_edge(i,2),1)-cmod%point(cmod%src_edge(i,1),1)
                direc(2) = cmod%point(cmod%src_edge(i,2),2)-cmod%point(cmod%src_edge(i,1),2)
                direc(3) = cmod%point(cmod%src_edge(i,2),3)-cmod%point(cmod%src_edge(i,1),3)
                do j=1,cmod%src_num
                    if (cmod%src_phy(i)==cmod%src_tp(j)) then
                        if (dot_product((cmod%src_end_2(j,:)-cmod%src_end_1(j,:)),direc)>0) then
                            if (cmod%src_edge(i,1)<cmod%src_edge(i,2)) then
                                cmod%src_drc(i) = 1
                            else
                                cmod%src_drc(i) =-1
                            end if
                        else
                            if (cmod%src_edge(i,1)<cmod%src_edge(i,2)) then
                                cmod%src_drc(i) =-1
                            else
                                cmod%src_drc(i) = 1
                            end if
                        end if
                    end if
                end do
            end do
        case ('sp_mthd')
            read(sValue,*) cmod%step_mthd
        case ('sp_division')
            read(sValue,*) cmod%step_dvsi
        case ('time_maximum')
            read(sValue,*) cmod%time_maxm
        case ('sp_double')
            read(sValue,*) cmod%step_dble
        case ('sp_dblesize')
            read(sValue,*) cmod%step_dblesize
        case ('output')
            read (ctl_unit,*) tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6),tmp(7),tmp(8),tmp(9)
            do i=1,9
                if (tmp(i)==1) then
                    cmod%o_f(i)=.true.
                else
                    cmod%o_f(i)=.false.
                end if
            end do
            allocate(cmod%o_v(0))
            do i=1,9
                if(cmod%o_f(i)) cmod%o_v=[cmod%o_v,i]
            end do
            case default
        end select
    end do
    close(ctl_unit)
    print*, "--- rebuild receiver info ..."
    print*, "--- print out meshing info ..."
    write(*,"(a,f12.3,8x,a,f12.3)")," upper bound x:",maxval(cmod%point(:,1))," lower bound x:",minval(cmod%point(:,1))
    write(*,"(a,f12.3,8x,a,f12.3)")," upper bound y:",maxval(cmod%point(:,2))," lower bound y:",minval(cmod%point(:,2))
    write(*,"(a,f12.3,8x,a,f12.3)")," upper bound z:",maxval(cmod%point(:,3))," lower bound z:",minval(cmod%point(:,3))
    print*, "--- print out source info ..."
    print*,"number of source edges   :", cmod%seg_num
    do i=1,cmod%src_num
        write(*,"(a,i3,a,i4,a,i4)"), " Source ",i,": Phy num",cmod%src_tp(i),", Edg Num:",cmod%src_count_edges(i)
        write(*,"(a,f12.3,f12.3,f12.3)"), " coor1:",cmod%src_end_1(i,1),cmod%src_end_1(i,2),cmod%src_end_1(i,3)
        write(*,"(a,f12.3,f12.3,f12.3)"), " coor2:",cmod%src_end_2(i,1),cmod%src_end_2(i,2),cmod%src_end_2(i,3)
    end do
    allocate(cmod%rec_tet(cmod%rec_num))
    cmod%rec_tet = 0
    do i = 1,cmod%rec_num
        call p_in_tet (cmod,[cmod%rec_cor(i,1),cmod%rec_cor(i,2),cmod%rec_cor(i,3)],cmod%rec_tet(i))
    end do
    print*, "--- print out receiver info ..."
    print*, "number of receivers      :", cmod%rec_num
    do i=1,cmod%rec_num
         write(*,"(a,3f16.2,a,i10)") " Coor:",cmod%rec_cor(i,1),cmod%rec_cor(i,2),cmod%rec_cor(i,3), " in tet num:",cmod%rec_tet(i)
    end do
    write(*,*) "Following field components will be output:",cmod%o_f
    print*, "--- print out domain info ..."
    do i = 1,cmod%dom_num
         write(*,"(a,i4,a,i8,a)"), " Isotropic domain",cmod%dom_ind(i)," is discretized into", count(cmod%tet_phy==cmod%dom_ind(i)), " tets."
         do j=1,cmod%tet_num
             if (cmod%dom_ind(i)==cmod%tet_phy(j)) then
                 write(*,"(a,e10.3,3x,a,e10.3,3x,a,e10.3)"), " Conductivity:",cmod%sig_3d(j),"Permittivity:",cmod%per_3d(j), "Magnetic permeability:",cmod%mag_3d(j)
                 exit
             end if
         end do
    end do
    end subroutine read_ctl
    
    subroutine compute_boundary(cmod)
    implicit none
    type(model),intent(inout) :: cmod
    logical :: mask(cmod%edg_num)
    integer,allocatable :: indices(:),bond_edge3(:),bond_edge(:),all_edge(:)
    integer,allocatable :: tet_plane1(:),tet_plane2(:),tet_plane3(:) !four planes of all tets
    integer :: plane(3) !one plane
    integer :: k,i,j
    print*, "==================================boundary processing..."
    allocate(tet_plane1(4*cmod%tet_num),tet_plane2(4*cmod%tet_num),tet_plane3(4*cmod%tet_num))
    k=0
    do i=1,cmod%tet_num
        do j=1,4
           k=k+1
           plane=[cmod%tet_edge(i,pla_tet(j,1)),cmod%tet_edge(i,pla_tet(j,2)),cmod%tet_edge(i,pla_tet(j,3))]
           call sort1d(plane)
           tet_plane1(k)=plane(1)
           tet_plane2(k)=plane(2)
           tet_plane3(k)=plane(3)
        end do
    end do
    call sort3d(tet_plane1,tet_plane2,tet_plane3)
    call bound(tet_plane1,tet_plane2,tet_plane3,size(tet_plane1))
    allocate(bond_edge3(3*size(tet_plane1)))
    bond_edge3=[tet_plane1,tet_plane2,tet_plane3]
    call QsortC(bond_edge3)
    call unique1d(bond_edge3,indices,bond_edge)
    print*, "number of boundary edge:",size(bond_edge)
    allocate(all_edge(cmod%edg_num))
    all_edge=[1:cmod%edg_num]
    do i=1,cmod%edg_num
        if (search(bond_edge,all_edge(i))==0) then
            mask(i)=.true.
        else
            mask(i)=.false.
        end if
    end do
    allocate(cmod%nobond_edge,source=pack(all_edge,mask))
    cmod%nbeg_num=size(cmod%nobond_edge)
    print*, "number of no boundary edge:",cmod%nbeg_num!内部的边
    !open(opt_unit,file=trim(cmod%opt_file))
    !do i=1,size(cmod%nobond_edge)
    !    write(opt_unit,*) cmod%nobond_edge(i)
    !end do
    !close(opt_unit)
    
    end subroutine compute_boundary
    
    
    
    subroutine compute_matrix(cmod)
    USE IFPORT
    implicit none
    type(model),intent(inout) :: cmod
    type(SPARSE_MATRIX_T) :: at_mkl,bt_mkl,ct_mkl
    integer :: i,m,n,k,info,l1,l2
    real(8) :: E(6,6),F(6,6),a(6,6),b(6,6),c(6,6)
    integer,allocatable :: XY1(:),XY2(:)
    real(8),allocatable :: a1(:),b1(:),c1(:)
    print*, "==================================mesh processing..."
    allocate(XY1(6*6*cmod%tet_num),XY2(6*6*cmod%tet_num))
    allocate(a1(6*6*cmod%tet_num),b1(6*6*cmod%tet_num),c1(6*6*cmod%tet_num))
    k=0
    do i=1,cmod%tet_num
        call intval(E,F,cmod%point(cmod%tet(i,1),:),cmod%point(cmod%tet(i,2),:),&
            &cmod%point(cmod%tet(i,3),:),cmod%point(cmod%tet(i,4),:))
        a=cmod%per_3d(i)*F
        b=cmod%sig_3d(i)*F
        c=1.d0/cmod%mag_3d(i)*E
        do m=1,6
            do n=1,6
                if (cmod%tet_edge(i,m)<=cmod%tet_edge(i,n)) then
                    l1=search(cmod%nobond_edge,cmod%tet_edge(i,m))
                    l2=search(cmod%nobond_edge,cmod%tet_edge(i,n))   
                    if (l1/=0.and.l2/=0) then 
                        k=k+1
                        XY1(k)=l1
                        XY2(k)=l2
                        a1(k)=cmod%tet_drc(i,m)*cmod%tet_drc(i,n)*a(m,n)
                        b1(k)=cmod%tet_drc(i,m)*cmod%tet_drc(i,n)*b(m,n)
                        c1(k)=cmod%tet_drc(i,m)*cmod%tet_drc(i,n)*c(m,n)
                    end if
                end if
            end do
        end do           
    end do
    print*, "--- handle global matrix A ..."
    info = mkl_sparse_d_create_coo (at_mkl, SPARSE_INDEX_BASE_one, cmod%nbeg_num, cmod%nbeg_num, k, XY1(1:k), XY2(1:k), a1(1:k))
    info = mkl_sparse_convert_csr (at_mkl, SPARSE_OPERATION_NON_TRANSPOSE, cmod%a_mkl)
    info = MKL_SPARSE_OPTIMIZE(cmod%a_mkl)
    info = mkl_sparse_destroy(at_mkl)
    print*, "--- handle global matrix B ..."
    info = mkl_sparse_d_create_coo (bt_mkl, SPARSE_INDEX_BASE_one, cmod%nbeg_num, cmod%nbeg_num, k, XY1(1:k), XY2(1:k), b1(1:k))
    info = mkl_sparse_convert_csr (bt_mkl, SPARSE_OPERATION_NON_TRANSPOSE, cmod%b_mkl)
    info = MKL_SPARSE_OPTIMIZE(cmod%b_mkl)
    info = mkl_sparse_destroy(bt_mkl)
    print*, "--- handle global matrix C ..."
    info = mkl_sparse_d_create_coo (ct_mkl, SPARSE_INDEX_BASE_one, cmod%nbeg_num, cmod%nbeg_num, k, XY1(1:k), XY2(1:k), c1(1:k))
    info = mkl_sparse_convert_csr (ct_mkl, SPARSE_OPERATION_NON_TRANSPOSE, cmod%c_mkl)
    info = MKL_SPARSE_OPTIMIZE(cmod%c_mkl)
    info = mkl_sparse_destroy(ct_mkl)
    end subroutine compute_matrix
    
    subroutine time_loop_1d(cmod)
    USE IFPORT
    !use pardiso_multi_real
    use source_mod
    use loop_solution
    implicit none
    type(model),intent(inout) :: cmod
    integer :: nsys,nnzsys
    integer,allocatable :: isys(:),jsys(:)
    real(8),allocatable :: vsys(:)
    ! mkl matrix
    type(SPARSE_MATRIX_T) :: p_mkl,tmp_p_mkl
    ! loop variables
    integer :: info
    ! e temp vectors
    real(8),allocatable :: e1(:,:),e2(:,:),e3(:,:),e4(:,:),r_h(:,:),r_h_1(:,:),rh(:),e_tmp_1(:)
    ! time variables
    real(8) :: time_cur,step_cur,error,time_limit
    integer :: step,not_changed,i,isrc
    ! folder
    logical(kind=4) :: ierr1,errnum
    integer(kind=4) :: istatus1
    character(len=256) :: cmd
    print*, "==================================creating a new output folder..."
    inquire(DIRECTORY=trim(file_root), EXIST=ierr1)
    if(ierr1) then
       print*,'The directory have existed and not been needed to create'
    else
       print*,'The directory not exist and creat it'
       cmd = "mkdir -p " // trim(file_root)
       print*,trim(cmd)
       istatus1=system(trim(cmd))
       if(istatus1==-1) then
           errnum=ierrno()
           print*,'Error=',errnum,' inquire the Intel Visual Fortran help document'
           stop ' Folder creating is fail'
       end if
    end if
    
    print*, "==================================time loop..."
    print*, "--- init time variables ..."
    time_frac=0.d0
    step = 0
    time_cur = 0.d0
    step_cur = minval(cmod%src_pulsewidth) / cmod%step_dvsi
    ! how many time step should be performed
    DESCR%TYPE = SPARSE_MATRIX_TYPE_SYMMETRIC
    DESCR%mode = SPARSE_FILL_MODE_UPPER
    DESCR%diag = SPARSE_DIAG_NON_UNIT
    ALLOCATE(OPTS_UNITS(cmod%src_num))
    do isrc = 1,cmod%src_num
        OPTS_UNITS(ISRC) = 1000+ISRC
        open (OPTS_UNITS(ISRC),name=trim(file_root)//"/"//trim(file_root)//"_"//itoa(isrc)//".field")
    END DO
    print*, "--- start time loop ..."
    time_limit = cmod%time_maxm
    time_1 = timef()
    not_changed=0
    cmod%add=0.d0
    cmod%back=.false.
    do while (time_cur<time_limit)
        step = step + 1
        time_cur = time_cur + step_cur
        if (step==1) then
            ! SYSTEM MATRIX
            call mkl_add_3 (cmod%nbeg_num, 0.d0,cmod%a_mkl, 3.d0, cmod%b_mkl, 2.d0*step_cur, cmod%c_mkl, p_mkl)
            call mkl_to_csr (p_mkl,nsys,nnzsys,vsys,isys,jsys)
            ! e and r_h variables
            allocate(e1(nsys,cmod%src_NUM),e2(nsys,cmod%src_NUM),e3(nsys,cmod%src_NUM),e4(nsys,cmod%src_NUM),r_h(nsys,cmod%src_NUM),r_h_1(nsys,cmod%src_NUM),rh(nsys*cmod%src_NUM),e_tmp_1(nsys*cmod%src_NUM))
            e1=0.d0; e2=0.d0; e3=0.d0; e4=0.d0; r_h=0.d0; r_h_1=0.d0;
            call solution_output (cmod,0,time_cur-step_cur,step_cur,e2)
            ! PARDISO AND system matrix
            !call pardiso_alloc(2)
            time_3=timef()
            call pardiso_direct_1 (.true.,nsys,nnzsys,isys,jsys,vsys)
            time_4=timef()
            time_frac = time_frac+time_4-time_3
            deallocate(isys,jsys,vsys)
            write(*,"(a,i6,a,d10.3,a,d10.3)") "--- fractoratization finished at step: ",step,", current time: ",time_cur,", current step: ",step_cur
        end if

        if  (mod(step,cmod%step_dble)==0 .and. step>cmod%step_dvsi) then
            step_cur = step_cur * cmod%step_dblesize
            time_cur = time_cur + step_cur - step_cur/cmod%step_dblesize
            call mkl_add_3 (cmod%nbeg_num, 0.d0,cmod%a_mkl, 3.d0, cmod%b_mkl, 2.d0*step_cur, cmod%c_mkl, tmp_p_mkl)
            call mkl_to_csr (tmp_p_mkl,nsys,nnzsys,vsys,isys,jsys)
            call pardiso_direct_3
            time_3=timef()
            call pardiso_direct_1 (.true.,nsys,nnzsys,isys,jsys,vsys)
            time_4=timef()
            time_frac = time_frac+time_4-time_3
            call sys_right (step,time_cur,cmod,r_h)
            do i=1,cmod%src_NUM
            e_tmp_1(1+nsys*(i-1):nsys*i) = 4.d0*e2(:,i)-e4(:,i)
            end do
            info = mkl_sparse_d_mm (SPARSE_OPERATION_NON_TRANSPOSE, 1.d0, cmod%b_mkl, DESCR, SPARSE_LAYOUT_COLUMN_MAJOR, e_tmp_1, cmod%src_NUM, nsys, 0.d0, rh, nsys)
            r_h_1=reshape(rh,[nsys,cmod%src_NUM])
            r_H = r_h_1 - 2.d0 * step_cur * r_h
            call pardiso_direct_2 (r_h,e1)
            write(*,"(a,i6,a,d10.3,a,d10.3)") "--- fractoratization finished at step: ",step,", current time: ",time_cur,", current step: ",step_cur
            e4 = e4
            e3 = e2
            e2 = e1
            call solution_output (cmod,step,time_cur,step_cur,e2)
            cycle
        end if
        ! sys_right and solve
        call sys_right (step,time_cur,cmod,r_h)
        do i=1,cmod%src_NUM
            e_tmp_1(1+nsys*(i-1):nsys*i) = 4.d0*e2(:,i)-e3(:,i)
        end do
        info = mkl_sparse_d_mm (SPARSE_OPERATION_NON_TRANSPOSE, 1.d0, cmod%b_mkl, DESCR, SPARSE_LAYOUT_COLUMN_MAJOR, e_tmp_1, cmod%src_NUM, nsys, 0.d0, rh, nsys)
        r_h_1=reshape(rh,[nsys,cmod%src_NUM])
        r_H = r_h_1 - 2.d0 * step_cur * r_h
        call pardiso_direct_2 (r_h,e1)
        ! output
        call solution_output (cmod,step,time_cur,step_cur,e1)
        e4 = e3
        e3 = e2
        e2 = e1
    end do
    do isrc = 1,cmod%src_num
        CLOSE (OPTS_UNITS(ISRC))
    END DO
    call pardiso_direct_3
    !call pardiso_dealloc
    time_2 = timef()
    time_stepping = time_2-time_1-time_frac
    write(*,"(a,f12.2,a)")," time of fractorization: ",time_frac," seconds."
    write(*,"(a,f12.2,a)")," time of stepping: ",time_stepping," seconds."
    end subroutine time_loop_1d
    
    
    subroutine clear_model (cmod)
    implicit none
    type(model),intent(inout) :: cmod
    integer :: info
    ! node
    cmod%nod_num=0
    if (allocated(cmod%point)) deallocate(cmod%point)
    ! source
    cmod%seg_num=0
    if (allocated(cmod%src_edge)) deallocate(cmod%src_edge)
    if (allocated(cmod%src_phy)) deallocate(cmod%src_phy)
    if (allocated(cmod%src_edg)) deallocate(cmod%src_edg)
    if (allocated(cmod%src_drc)) deallocate(cmod%src_drc)
    ! source info from the control file
    cmod%src_num=0
    cmod%src_pulsetype=0
    cmod%src_pulsewidth=0.d0
    cmod%src_current=0.d0
    if (allocated(cmod%src_tp)) deallocate(cmod%src_tp)
    if (allocated(cmod%src_count_edges)) deallocate(cmod%src_count_edges)
    if (allocated(cmod%src_end_1)) deallocate(cmod%src_end_1)
    if (allocated(cmod%src_end_2)) deallocate(cmod%src_end_2)
    ! tet
    cmod%tet_num=0
    if (allocated(cmod%tet_edge)) deallocate(cmod%tet_edge)
    if (allocated(cmod%tet_phy)) deallocate(cmod%tet_phy)
    if (allocated(cmod%tet_drc)) deallocate(cmod%tet_drc)
    ! sur_edg
    !if (allocated(cmod%sur_edg_list)) deallocate(cmod%sur_edg_list)
    !if (allocated(cmod%sur_edg_index)) deallocate(cmod%sur_edg_index)
    ! edge
    cmod%edg_num=0
    if (allocated(cmod%edge)) deallocate(cmod%edge)
    ! receiver
    cmod%rec_num=0
    if (allocated(cmod%rec_cor)) deallocate(cmod%rec_cor)
    if (allocated(cmod%rec_tet)) deallocate(cmod%rec_tet)
    !if (allocated(cmod%n_arr)) deallocate(cmod%n_arr)
    !if (allocated(cmod%d_arr)) deallocate(cmod%d_arr)
    ! domain character
    cmod%dom_num=0
    if (allocated(cmod%dom_ind)) deallocate(cmod%dom_ind)
    if (allocated(cmod%sig_3d)) deallocate(cmod%sig_3d)
    if (allocated(cmod%per_3d)) deallocate(cmod%per_3d)
    if (allocated(cmod%mag_3d)) deallocate(cmod%mag_3d)
    ! system matrix
    info = mkl_sparse_destroy (cmod%a_mkl)
    info = mkl_sparse_destroy (cmod%b_mkl)
    info = mkl_sparse_destroy (cmod%c_mkl)
    ! surf_kind: 0 - PEC boundary ; 1 - 1st order ABC
    !cmod%bound_cond=0
    ! time variables
    cmod%step_mthd=0
    cmod%step_dvsi=0
    cmod%step_dble=0
    cmod%step_dblesize=0
    cmod%time_maxm=0.d0
    ! field variables
    if (allocated(cmod%last_b_f)) deallocate(cmod%last_b_f)
    !cmod%ind_f=0
    !cmod%init_size_F=0
    ! output variables: ex ey ez dx dy dz bx by bz
    !cmod%ind_f=0
    !cmod%init_size_F=0
    cmod%o_f=.false.
    if (allocated(cmod%o_v)) deallocate(cmod%o_v)
    ! vtk output variables
    !cmod%vtk_num = 0
    if (allocated(cmod%src_pulsetype))  deallocate(cmod%src_pulsetype)
    if (allocated(cmod%src_pulsewidth))  deallocate(cmod%src_pulsewidth)
    if (allocated(cmod%src_current))  deallocate(cmod%src_current)
    if (allocated(cmod%src_delay))  deallocate(cmod%src_delay)
    end subroutine clear_model
    
    subroutine p_in_tet (cmod,p,in_which)
    implicit none
    type(model),intent(in) :: cmod
    real(8),intent(in) :: p(:)
    integer,intent(out) :: in_which
    real(8) :: x(4),y(4),z(4),d(5)
    integer :: i,j
    ! process
    in_which = 0
    !!$OMP PARALLEL DO private(i,j,x,y,z,d)
    !in_num = 0
    do i = 1,cmod%tet_num
        do j=1,4
            x(j) = cmod%point(cmod%tet(i,j),1)
            y(j) = cmod%point(cmod%tet(i,j),2)
            z(j) = cmod%point(cmod%tet(i,j),3)
        end do
        d(1) = det4(reshape([x(1),x(2),x(3),x(4),y(1),y(2),y(3),y(4),z(1),z(2),z(3),z(4),1.d0,1.d0,1.d0,1.d0],[4,4]))
        d(2) = det4(reshape([p(1),x(2),x(3),x(4),p(2),y(2),y(3),y(4),p(3),z(2),z(3),z(4),1.d0,1.d0,1.d0,1.d0],[4,4]))
        d(3) = det4(reshape([x(1),p(1),x(3),x(4),y(1),p(2),y(3),y(4),z(1),p(3),z(3),z(4),1.d0,1.d0,1.d0,1.d0],[4,4]))
        d(4) = det4(reshape([x(1),x(2),p(1),x(4),y(1),y(2),p(2),y(4),z(1),z(2),p(3),z(4),1.d0,1.d0,1.d0,1.d0],[4,4]))
        d(5) = det4(reshape([x(1),x(2),x(3),p(1),y(1),y(2),y(3),p(2),z(1),z(2),z(3),p(3),1.d0,1.d0,1.d0,1.d0],[4,4]))
        if ( (d(1)>=0.d0.and.d(2)>=0.d0.and.d(3)>=0.d0.and.d(4)>=0.d0.and.d(5)>=0.d0).or.(d(1)<=0.d0.and.d(2)<=0.d0.and.d(3)<=0.d0.and.d(4)<=0.d0.and.d(5)<=0.d0) ) then
            in_which = i
            exit
        end if
    end do
    if (in_which==0) stop 'no rec tet'
    end subroutine p_in_tet
    
    subroutine solution_output (cmod,step,t,step_cur,res_e)
    use loop_solution
    implicit none
    type(model),intent(inout) :: cmod
    integer,intent(in) :: step
    real(8),intent(in) :: t,step_cur
    real(8),allocatable,intent(in) :: res_e(:,:)
    integer :: i,j,isrc
    real(8) :: a_f(cmod%src_num,cmod%rec_num,9),all_e(cmod%edg_num,cmod%src_num)
    real(8) :: N(6,3),curlN(6,3)
    a_f=0.d0
    all_e=0.d0
    do i=1,cmod%nbeg_num
        all_e(cmod%nobond_edge(i),:)=res_e(i,:)
    end do
    do isrc=1,cmod%src_num
        do i=1,cmod%rec_num
            call Ne(N,curlN,cmod%point(cmod%tet(cmod%rec_tet(i),1),:),cmod%point(cmod%tet(cmod%rec_tet(i),2),:),&
            &cmod%point(cmod%tet(cmod%rec_tet(i),3),:),cmod%point(cmod%tet(cmod%rec_tet(i),4),:),&
            &cmod%rec_cor(i,:))
            do j=1,6
               !e field
               a_f(isrc,i,1)=a_f(isrc,i,1)+N(j,1)*cmod%tet_drc(cmod%rec_tet(i),j)*all_e(cmod%tet_edge(cmod%rec_tet(i),j),isrc)
               a_f(isrc,i,2)=a_f(isrc,i,2)+N(j,2)*cmod%tet_drc(cmod%rec_tet(i),j)*all_e(cmod%tet_edge(cmod%rec_tet(i),j),isrc)
               a_f(isrc,i,3)=a_f(isrc,i,3)+N(j,3)*cmod%tet_drc(cmod%rec_tet(i),j)*all_e(cmod%tet_edge(cmod%rec_tet(i),j),isrc)
               !dB field
               a_f(isrc,i,4)=a_f(isrc,i,4)-curlN(j,1)*cmod%tet_drc(cmod%rec_tet(i),j)*all_e(cmod%tet_edge(cmod%rec_tet(i),j),isrc)
               a_f(isrc,i,5)=a_f(isrc,i,5)-curlN(j,2)*cmod%tet_drc(cmod%rec_tet(i),j)*all_e(cmod%tet_edge(cmod%rec_tet(i),j),isrc)
               a_f(isrc,i,6)=a_f(isrc,i,6)-curlN(j,3)*cmod%tet_drc(cmod%rec_tet(i),j)*all_e(cmod%tet_edge(cmod%rec_tet(i),j),isrc)
            end do
        end do
    end do
    if (.not.allocated(cmod%last_b_f)) allocate(cmod%last_b_f(cmod%src_num,cmod%rec_num,3))
    if (step==0) then
        cmod%last_b_f=0.d0
    else
        cmod%last_b_f=cmod%last_b_f+a_f(:,:,4:6)*step_cur
    end if
    !B field
    a_f(:,:,7:9)=cmod%last_b_f
    do isrc = 1,cmod%src_num
    write(optS_unitS(isrc),"(i6,e15.5,2x,"//itoa(cmod%rec_num)//"("//itoa(count(cmod%o_f))//"e15.5,2x))") step,t,(a_f(isrc,i,cmod%o_v),i=1,cmod%rec_num)
    end do
    end subroutine solution_output
    
    subroutine sys_right (step,t,cmod,r_h)
    use source_mod
    implicit none
    ! in/out variables
    integer,intent(in) :: step
    real(8),intent(in) :: t
    type(model),intent(inout) :: cmod
    real(8),intent(out),allocatable :: r_h(:,:)
    ! local variables
    real(8) :: didt
    integer :: i,j,k
    logical :: find_phy
    ! process
    allocate(r_h(cmod%nbeg_num,cmod%src_num))
    r_h = 0.d0
    if (cmod%src_pulsetype(1)==2) then
    if ( step >= cmod%step_dvsi ) return
    end if
    !if (cmod%src_pulsetype==4) then
    !if ( t > cmod%src_pulsewidth ) return
    !end if
   
    if (cmod%back) return

    do i=1,cmod%seg_num
        find_phy = .false.
        do j=1,cmod%src_num
            if (cmod%src_phy(i)==cmod%src_tp(j)) then
                find_phy = .true.
                call time_const_p (cmod%src_pulsetype(j),t,cmod%src_pulsewidth(j),cmod%src_current(j),didt)  
                didt = didt / cmod%src_count_edges(j)
                if ((cmod%add+didt)<=0) then
                    didt=-cmod%add
                    cmod%back=.true.
                end if 
                k=search(cmod%nobond_edge,cmod%src_edg(i))
                r_h(k,j) = didt * cmod%src_drc(i) !* cmod%edg_len(cmod%src_edg(i))
            end if
        end do
        if (.not.find_phy) stop("error - physical domain not found!") 
    end do
    cmod%add=cmod%add+didt
    r_h = r_h / cmod%src_pulsewidth(1) * 3.d0
    end subroutine sys_right
    
    subroutine pardiso_direct_1 (symm,n_in,nnz_in,ia_in,ja_in,a_in_double)
    IMPLICIT NONE
    logical,intent(in) :: symm
    integer,intent(in) :: n_in,nnz_in
    INTEGER,ALLOCATABLE,intent(in) :: ia_in(:),ja_in(:)
    REAL(8),ALLOCATABLE,intent(in) :: a_in_double(:)
    if (symm) then
        mtype_par = 2
    else
        mtype_par = 11
    end if
    n_par = n_in
    nnz_par = nnz_in
    ALLOCATE (pt_par(64),iparm_par(64))
    call pardisoinit (pt_par, mtype_par, iparm_par)
    iparm_par(2) = 3
    iparm_par(24) = 1
    nrhs_par=1; maxfct_par=1; mnum_par=1; error_par=0; msglvl_par=0
    phase_par = 12
    CALL pardiso (pt_par,maxfct_par,mnum_par,mtype_par,phase_par,n_par,a_in_double,ia_in,ja_in,idum_par,nrhs_par,iparm_par,msglvl_par,ddum_par,ddum_par,error_par)
    IF (error_par /= 0) THEN
        if (symm) then
            mtype_par = -2
            print*,"trying indefinite frac..."
            CALL pardiso (pt_par,maxfct_par,mnum_par,mtype_par,phase_par,n_par,a_in_double,ia_in,ja_in,idum_par,nrhs_par,iparm_par,msglvl_par,ddum_par,ddum_par,error_par)
            IF (error_par /= 0) then
                print*, 'ERROR - The following error was detected during analysis step of pardiso solver: ', error_par
                stop
            end if
        ELSE
            print*, 'ERROR - The following error was detected during analysis step of pardiso solver: ', error_par
            stop
        END IF
    END IF
    END subroutine pardiso_direct_1

    subroutine pardiso_direct_2 (b_in_double,x_out_double)
    implicit none
    real(8),allocatable :: b_in_double(:,:),x_out_double(:,:)
    phase_par = 33
    CALL pardiso (pt_par,maxfct_par,mnum_par,mtype_par,phase_par,n_par,ddum_par,idum_par,idum_par,idum_par,size(b_in_double,2),iparm_par,msglvl_par,b_in_double,x_out_double,error_par)
    IF (error_par/= 0) THEN
        print*, 'ERROR - The following error was detected during solve step of pardiso solver: ', error_par
        stop
    END IF
    END subroutine pardiso_direct_2

    subroutine pardiso_direct_3
    implicit none
    phase_par = -1
    CALL pardiso (pt_par,maxfct_par,mnum_par,mtype_par,phase_par,n_par,ddum_par,idum_par,idum_par,idum_par,nrhs_par,iparm_par,msglvl_par,ddum_par,ddum_par,error_par)
    IF(ALLOCATED(pt_par)) DEALLOCATE(pt_par)
    IF(ALLOCATED(iparm_par)) DEALLOCATE(iparm_par)
    IF (error_par /= 0) THEN
        print*, 'ERR - The following error was detected during terminal step of pardiso solver: ',error_par
        stop
    END IF
    END subroutine pardiso_direct_3
    
end module file_option