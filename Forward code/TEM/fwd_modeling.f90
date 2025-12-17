module fwd_modeling
    implicit none
    contains
    
    subroutine swap(a,b)
    implicit none
    integer :: a,b,k
    if (a>b) then
        k=a
        a=b
        b=k
    end if
    end subroutine swap
    
    subroutine sort1d(src)
    implicit none
    integer length
    integer min,i,j,tmp
    integer,INTENT(INOUT) :: src(:)
    length=size(src)
    do i=1,length-1
      min=i
      do j=i+1,length
          if (src(j) < src(min)) then
              min=j
          end if
      end do
      tmp = src(min)
      src(min) = src(i)
      src(i) = tmp 
    end do
    end subroutine sort1d
    
    subroutine sort2d(edge1,edge2)
    implicit none
    integer :: i,k,length
    integer, INTENT(INOUT) :: edge1(:),edge2(:)
    integer,allocatable :: indices(:)
    call QsortC(edge1,edge2)
    call unique1d(edge1,indices)
    do i=1,size(indices)
        if (i<size(indices)) then
            k=indices(i+1)-1
        else
            k=size(edge1)
        end if
        call sort1d(edge2(indices(i):k))
    end do
    deallocate(indices)
    end subroutine sort2d
    
    subroutine sort3d(edge1,edge2,edge3)
    implicit none
    integer :: i,k,length
    integer, INTENT(INOUT) :: edge1(:),edge2(:),edge3(:)
    integer,allocatable :: indices(:)
    call QsortC(edge1,edge2,edge3)
    call unique1d(edge1,indices)
    do i=1,size(indices)
        if (i<size(indices)) then
            k=indices(i+1)-1
        else
            k=size(edge1)
        end if
        call sort2d(edge2(indices(i):k),edge3(indices(i):k))
    end do
    deallocate(indices)
    end subroutine sort3d
    
    
    
    subroutine unique1d(edge1,indices,value)
    implicit none
    integer :: kx
    integer,intent(in) :: edge1(:)!,put(:)
    integer,allocatable,intent(out) :: indices(:)
    integer,allocatable,optional :: value(:)
    logical :: mask(size(edge1))
    mask = .true.
    do kx=size(edge1),2,-1
      mask(kx)= .not.(edge1(kx-1)==edge1(kx))
    end do
    allocate(indices,source=pack([(kx,kx=1,size(edge1))],mask))
    if (present(value)) allocate(value,source=pack(edge1,mask))
    !allocate(put(size(indices)))
    !do kx=1,size(indices)
    !    put(kx)=src(indices(kx))
    !end do
    !deallocate(src)
    !allocate(src(size(indices)))
    !src=put
    end subroutine unique1d
    
    subroutine unique2d(edge1,edge2,length)
    implicit none
    integer :: kx,length 
    integer,allocatable :: edge1(:),edge2(:)
    logical :: mask(length)
    integer :: put1(size(edge1)),put2(size(edge2))
    put1=edge1
    put2=edge2
    mask(1)=.true.
    do kx=length,2,-1
      mask(kx)= .not.((put1(kx-1)==put1(kx)).and.(put2(kx-1)==put2(kx)))
    end do
    deallocate(edge1,edge2)
    allocate(edge1,source=pack(put1,mask))
    allocate(edge2,source=pack(put2,mask))
    end subroutine unique2d
    
    subroutine repeat3d(edge1,edge2,edge3,length)
    implicit none
    integer :: kx,length 
    integer,allocatable :: edge1(:),edge2(:),edge3(:)
    logical :: mask(length)
    integer :: put1(size(edge1)),put2(size(edge2)),put3(size(edge3))
    put1=edge1
    put2=edge2
    put3=edge3
    if ((put1(1)==put1(2)).and.(put2(1)==put2(2))) then
        mask(1)=.true.
    else
        mask(1)=.false.
    end if
    do kx=2,length-1
      mask(kx)= ((put1(kx+1)==put1(kx)).and.(put2(kx+1)==put2(kx)))
    end do
    mask(length)=.false.
    deallocate(edge1,edge2,edge3)
    allocate(edge1,source=pack(put1,mask))
    allocate(edge2,source=pack(put2,mask))
    allocate(edge3,source=pack(put3,mask))
    end subroutine repeat3d
    
    subroutine bound(edge1,edge2,edge3,length)
    !occur one time surface
    implicit none
    integer :: kx,length 
    integer,allocatable :: edge1(:),edge2(:),edge3(:)
    logical :: mask(length)
    integer :: put1(size(edge1)),put2(size(edge2)),put3(size(edge3))
    put1=edge1
    put2=edge2
    put3=edge3
    if ((put1(1)==put1(2)).and.(put2(1)==put2(2))) then
        mask(1)=.false.
    else
        mask(1)=.true.
    end if
    if ((put1(length-1)==put1(length)).and.(put2(length-1)==put2(length))) then
        mask(length)=.false.
    else
        mask(length)=.true.
    end if
    do kx=2,length-1
      mask(kx)=.not.(((put1(kx+1)==put1(kx)).and.(put2(kx+1)==put2(kx))).or.((put1(kx-1)==put1(kx)).and.(put2(kx-1)==put2(kx))))
    end do
    
    deallocate(edge1,edge2,edge3)
    allocate(edge1,source=pack(put1,mask))
    allocate(edge2,source=pack(put2,mask))
    allocate(edge3,source=pack(put3,mask))
    end subroutine bound

  recursive subroutine QsortC(A,B,C)
   implicit none
    integer, intent(in out), dimension(:) :: A
    integer, intent(in out), dimension(:), optional :: C,B
    integer :: iq
    if(size(A) > 1) then
      if (present(C)) then
        call Partition(A,iq,B,C)
        call QsortC(A(:iq-1),B(:iq-1),C(:iq-1))
        call QsortC(A(iq:),B(iq:),C(iq:))
      elseif (present(B)) then
        call Partition(A,iq,B)
        call QsortC(A(:iq-1),B(:iq-1))
        call QsortC(A(iq:),B(iq:))
      else
        call Partition(A,iq)
        call QsortC(A(:iq-1))
        call QsortC(A(iq:)) 
      end if
    end if
  end subroutine QsortC

  subroutine Partition(A,marker,B,C)
  implicit none
    integer, intent(in out), dimension(:) :: A
    integer, intent(in out), dimension(:), optional :: C,B
    integer, intent(out) :: marker
    integer :: i, j
    real :: temp
    real :: x      
    x = A(1)
    i= 0
    j= size(A) + 1
    do
      j = j-1
      do
        if (A(j) <= x) exit
        j = j-1
      end do
      i = i+1
      do
        if (A(i) >= x) exit
        i = i+1
      end do
      if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
        if (present(B)) then
        temp = B(i)
        B(i) = B(j)
        B(j) = temp
        end if
        if (present(C)) then
            temp = C(i)
            C(i) = C(j)
            C(j) = temp
        end if
        
      elseif (i == j) then
        marker = i+1
        return
      else
        marker = i
        return
      endif
    end do
  end subroutine Partition

  
  integer function search(list,x) 
    implicit none
    integer,intent(in) :: x
    integer,INTENT(IN),target :: list(:)
    integer, pointer :: p(:)
    integer :: offset,mid
    p => list
    search = 0
    offset = 0
    do while (size(p) > 0)
        mid = size(p)/2 + 1
        if (p(mid) > x) then
            p => p(:mid-1)
        else if (p(mid) < x) then
            offset = offset + mid
            p => p(mid+1:)
        else
            search = offset + mid
            exit
        end if
    end do
  end function search
  
  
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
  
  subroutine Lower( s )
    ! convert string to lower case
  character(*), intent(out)  :: s
  integer :: i
  do  i=1,len_trim(s)
      if  ( s(i:i) >= 'A' .and. s(i:i) <= 'Z' ) then
          s(i:i) = char(ichar(s(i:i)) + 32)
      endif
  enddo
  end subroutine Lower
   
  real(8) function det4(A)
  IMPLICIT NONE
  real(8), INTENT(IN)  :: A(:,:)
  det4 =  A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+ &
      A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+ &
      A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
  RETURN
  END FUNCTION det4
  
  real(8) FUNCTION det3 (A)
  implicit none
  real(8), INTENT(IN)  :: A(:,:)
  det3 = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3)  + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)
  RETURN
  end function det3
  
  subroutine intval(E,F,nod1,nod2,nod3,nod4)
  implicit none
  real(8),intent(in) :: nod1(3),nod2(3),nod3(3),nod4(3)
  real(8),intent(inout) :: E(6,6),F(6,6)
  real(8) :: b(4),c(4),d(4),l(6),ff(4,4)
  real(8) :: V
  integer :: cor(6,2)=reshape([1,1,1,2,4,3,2,3,4,3,2,4],[6,2])
  integer :: i,j,m,n
  b(1)=-det3(reshape([1.d0,1.d0,1.d0,nod2(2),nod3(2),nod4(2),nod2(3),nod3(3),nod4(3)],[3,3]));
  b(2)=det3(reshape([1.d0,1.d0,1.d0,nod1(2),nod3(2),nod4(2),nod1(3),nod3(3),nod4(3)],[3,3]));
  b(3)=-det3(reshape([1.d0,1.d0,1.d0,nod1(2),nod2(2),nod4(2),nod1(3),nod2(3),nod4(3)],[3,3]));
  b(4)=det3(reshape([1.d0,1.d0,1.d0,nod1(2),nod2(2),nod3(2),nod1(3),nod2(3),nod3(3)],[3,3]));
  c(1)=det3(reshape([1.d0,1.d0,1.d0,nod2(1),nod3(1),nod4(1),nod2(3),nod3(3),nod4(3)],[3,3]));
  c(2)=-det3(reshape([1.d0,1.d0,1.d0,nod1(1),nod3(1),nod4(1),nod1(3),nod3(3),nod4(3)],[3,3]));
  c(3)=det3(reshape([1.d0,1.d0,1.d0,nod1(1),nod2(1),nod4(1),nod1(3),nod2(3),nod4(3)],[3,3]));
  c(4)=-det3(reshape([1.d0,1.d0,1.d0,nod1(1),nod2(1),nod3(1),nod1(3),nod2(3),nod3(3)],[3,3]));
  d(1)=-det3(reshape([1.d0,1.d0,1.d0,nod2(1),nod3(1),nod4(1),nod2(2),nod3(2),nod4(2)],[3,3]));
  d(2)=det3(reshape([1.d0,1.d0,1.d0,nod1(1),nod3(1),nod4(1),nod1(2),nod3(2),nod4(2)],[3,3]));
  d(3)=-det3(reshape([1.d0,1.d0,1.d0,nod1(1),nod2(1),nod4(1),nod1(2),nod2(2),nod4(2)],[3,3]));
  d(4)=det3(reshape([1.d0,1.d0,1.d0,nod1(1),nod2(1),nod3(1),nod1(2),nod2(2),nod3(2)],[3,3]));
  l(1)=norm(nod2-nod1)
  l(2)=norm(nod3-nod1)
  l(3)=norm(nod4-nod1)
  l(4)=norm(nod3-nod2)
  l(5)=norm(nod4-nod2)
  l(6)=norm(nod4-nod3)
  V=1.d0/6.d0*det4(reshape([1.d0,1.d0,1.d0,1.d0,nod1(1),nod2(1),nod3(1),nod4(1),&
      &nod1(2),nod2(2),nod3(2),nod4(2),nod1(3),nod2(3),nod3(3),nod4(3)],[4,4]))
  do i=1,4
     do j=1,4
         ff(i,j)=b(i)*b(j)+c(i)*c(j)+d(i)*d(j)
     end do 
  end do
  do i=1,6
      do j=1,6
         E(i,j)=(4.d0*l(i)*l(j)*V/(6.d0*V)**4)*((c(cor(i,1))*d(cor(i,2))-c(cor(i,2))*d(cor(i,1)))*(c(cor(j,1))*d(cor(j,2))-c(cor(j,2))*d(cor(j,1)))+&
             &(b(cor(i,2))*d(cor(i,1))-b(cor(i,1))*d(cor(i,2)))*(b(cor(j,2))*d(cor(j,1))-b(cor(j,1))*d(cor(j,2)))+&
             &(b(cor(i,1))*c(cor(i,2))-b(cor(i,2))*c(cor(i,1)))*(b(cor(j,1))*c(cor(j,2))-b(cor(j,2))*c(cor(j,1))))
         F(i,j)=(l(i)*l(j)/(720.d0*V))*(ff(cor(i,2),cor(j,2))*(1.d0+delt(cor(i,1),cor(j,1)))-ff(cor(i,1),cor(j,2))*(1.d0+delt(cor(i,2),cor(j,1)))&
             &-ff(cor(i,2),cor(j,1))*(1.d0+delt(cor(i,1),cor(j,2)))+ff(cor(i,1),cor(j,1))*(1.d0+delt(cor(i,2),cor(j,2)))) 
      end do
  end do
  end subroutine intval
  
  subroutine Ne(N,curlN,nod1,nod2,nod3,nod4,p)
  implicit none
  real(8),intent(in) :: nod1(3),nod2(3),nod3(3),nod4(3),p(3)
  real(8),intent(inout) :: N(6,3),curlN(6,3)
  real(8) :: a(4),b(4),c(4),d(4),l(6)
  real(8) :: V
  integer :: cor(6,2)=reshape([1,1,1,2,4,3,2,3,4,3,2,4],[6,2])
  integer :: i
  a(1)=det3(reshape([nod2(1),nod3(1),nod4(1),nod2(2),nod3(2),nod4(2),nod2(3),nod3(3),nod4(3)],[3,3]))
  a(2)=-det3(reshape([nod1(1),nod3(1),nod4(1),nod1(2),nod3(2),nod4(2),nod1(3),nod3(3),nod4(3)],[3,3]))
  a(3)=det3(reshape([nod1(1),nod2(1),nod4(1),nod1(2),nod2(2),nod4(2),nod1(3),nod2(3),nod4(3)],[3,3]))
  a(4)=-det3(reshape([nod1(1),nod2(1),nod3(1),nod1(2),nod2(2),nod3(2),nod1(3),nod2(3),nod3(3)],[3,3]))
  b(1)=-det3(reshape([1.d0,1.d0,1.d0,nod2(2),nod3(2),nod4(2),nod2(3),nod3(3),nod4(3)],[3,3]));
  b(2)=det3(reshape([1.d0,1.d0,1.d0,nod1(2),nod3(2),nod4(2),nod1(3),nod3(3),nod4(3)],[3,3]));
  b(3)=-det3(reshape([1.d0,1.d0,1.d0,nod1(2),nod2(2),nod4(2),nod1(3),nod2(3),nod4(3)],[3,3]));
  b(4)=det3(reshape([1.d0,1.d0,1.d0,nod1(2),nod2(2),nod3(2),nod1(3),nod2(3),nod3(3)],[3,3]));
  c(1)=det3(reshape([1.d0,1.d0,1.d0,nod2(1),nod3(1),nod4(1),nod2(3),nod3(3),nod4(3)],[3,3]));
  c(2)=-det3(reshape([1.d0,1.d0,1.d0,nod1(1),nod3(1),nod4(1),nod1(3),nod3(3),nod4(3)],[3,3]));
  c(3)=det3(reshape([1.d0,1.d0,1.d0,nod1(1),nod2(1),nod4(1),nod1(3),nod2(3),nod4(3)],[3,3]));
  c(4)=-det3(reshape([1.d0,1.d0,1.d0,nod1(1),nod2(1),nod3(1),nod1(3),nod2(3),nod3(3)],[3,3]));
  d(1)=-det3(reshape([1.d0,1.d0,1.d0,nod2(1),nod3(1),nod4(1),nod2(2),nod3(2),nod4(2)],[3,3]));
  d(2)=det3(reshape([1.d0,1.d0,1.d0,nod1(1),nod3(1),nod4(1),nod1(2),nod3(2),nod4(2)],[3,3]));
  d(3)=-det3(reshape([1.d0,1.d0,1.d0,nod1(1),nod2(1),nod4(1),nod1(2),nod2(2),nod4(2)],[3,3]));
  d(4)=det3(reshape([1.d0,1.d0,1.d0,nod1(1),nod2(1),nod3(1),nod1(2),nod2(2),nod3(2)],[3,3]));
  l(1)=norm(nod2-nod1)
  l(2)=norm(nod3-nod1)
  l(3)=norm(nod4-nod1)
  l(4)=norm(nod3-nod2)
  l(5)=norm(nod4-nod2)
  l(6)=norm(nod4-nod3)
  V=1.d0/6.d0*det4(reshape([1.d0,1.d0,1.d0,1.d0,nod1(1),nod2(1),nod3(1),nod4(1),&
      &nod1(2),nod2(2),nod3(2),nod4(2),nod1(3),nod2(3),nod3(3),nod4(3)],[4,4]))
  do i=1,6
      N(i,1)=(1.d0/(6.d0*V)*(a(cor(i,1))+b(cor(i,1))*p(1)+c(cor(i,1))*p(2)+d(cor(i,1))*p(3))*1.d0/(6.d0*V)*b(cor(i,2))-&
          &1.d0/(6.d0*V)*(a(cor(i,2))+b(cor(i,2))*p(1)+c(cor(i,2))*p(2)+d(cor(i,2))*p(3))*1.d0/(6.d0*V)*b(cor(i,1)))*l(i)
      N(i,2)=(1.d0/(6.d0*V)*(a(cor(i,1))+b(cor(i,1))*p(1)+c(cor(i,1))*p(2)+d(cor(i,1))*p(3))*1.d0/(6.d0*V)*c(cor(i,2))-&
          &1.d0/(6.d0*V)*(a(cor(i,2))+b(cor(i,2))*p(1)+c(cor(i,2))*p(2)+d(cor(i,2))*p(3))*1.d0/(6.d0*V)*c(cor(i,1)))*l(i)
      N(i,3)=(1.d0/(6.d0*V)*(a(cor(i,1))+b(cor(i,1))*p(1)+c(cor(i,1))*p(2)+d(cor(i,1))*p(3))*1.d0/(6.d0*V)*d(cor(i,2))-&
          &1.d0/(6.d0*V)*(a(cor(i,2))+b(cor(i,2))*p(1)+c(cor(i,2))*p(2)+d(cor(i,2))*p(3))*1.d0/(6.d0*V)*d(cor(i,1)))*l(i)
      curlN(i,1)=2.d0*l(i)/(6.d0*V)**2*(c(cor(i,1))*d(cor(i,2))-c(cor(i,2))*d(cor(i,1)))
      curlN(i,2)=2.d0*l(i)/(6.d0*V)**2*(b(cor(i,2))*d(cor(i,1))-b(cor(i,1))*d(cor(i,2)))
      curlN(i,3)=2.d0*l(i)/(6.d0*V)**2*(b(cor(i,1))*c(cor(i,2))-b(cor(i,2))*c(cor(i,1)))
  end do
  end subroutine Ne
  
  real(8) function norm(x)
  implicit none
  real(8) :: x(3)
  norm=sqrt(x(1)**2+x(2)**2+x(3)**2)
  end function norm
  
  real(8) function delt(i,j)
  integer :: i,j
  if (i==j) then
      delt=1.d0
  else
      delt=0.d0
  end if 
  end function delt
  
end module fwd_modeling
    
! module pardiso_multi_real
!    ! a module calling MKL Pardiso to solve symmetric/nonsymmetric real linear equations
!    ! ----------------------------------------------------------
!    ! public function protos for calling
!    !pardiso_alloc (max_matrix)
!    !pardiso_direct_1 (num_m,upper,n_in,nnz_in,ia,ja,va)
!    !pardiso_direct_2 (num_m,b_in,x_out,ia,ja,va)
!    !pardiso_direct_3 (num_m)
!    !pardiso_dealloc
!    ! ----------------------------------------------------------
!    USE mkl_pardiso
!    implicit none
!    private
!    ! Internal solver memory pointer pt and iparm vector
!    type(MKL_PARDISO_HANDLE),allocatable :: pt(:,:)
!    INTEGER, ALLOCATABLE :: iparm(:,:)
!    ! temporal variables
!    INTEGER :: idum(1)
!    REAL(8) :: ddum(1)
!    ! other variables
!    INTEGER :: maxfct=1, mnum=1, nrhs=1, msglvl=0, error, phase
!    integer,allocatable :: mtype(:), n(:), nnz(:)
!    ! public subroutines
!    public :: pardiso_alloc,pardiso_direct_1,pardiso_direct_2,pardiso_direct_3,pardiso_dealloc
!    ! module subroutine
!    contains
!    ! pardiso_alloc
!    subroutine pardiso_alloc (max_matrix)
!    ! pre allocate variables
!    implicit none
!    ! in/out variables
!    integer,intent(in) :: max_matrix
!    ! process
!    if (max_matrix<1) stop "maximum number of matrix error, stop!"
!    maxfct = max_matrix
!    ! allocate internal solver memory pointer pt and iparm vector
!    allocate ( pt(maxfct,64) )
!    allocate ( iparm(maxfct,64) )
!    ! allocate other solver variables
!    allocate(mtype(maxfct),  n(maxfct),nnz(maxfct))
!    end subroutine pardiso_alloc
!    subroutine pardiso_direct_1 (symm,n_in,nnz_in,ia_in,ja_in,a_in_double)
!    IMPLICIT NONE
!    logical,intent(in) :: symm
!    integer,intent(in) :: n_in,nnz_in
!    INTEGER,ALLOCATABLE,intent(in) :: ia_in(:),ja_in(:)
!    REAL(8),ALLOCATABLE,intent(in) :: a_in_double(:)
!    if (symm) then
!        mtype_par = 2
!    else
!        mtype_par = 11
!    end if
!    n_par = n_in
!    nnz_par = nnz_in
!    ALLOCATE (pt_par(64),iparm_par(64))
!    call pardisoinit (pt_par, mtype_par, iparm_par)
!    iparm_par(2) = 3
!    iparm_par(24) = 1
!    nrhs_par=1; maxfct_par=1; mnum_par=1; error_par=0; msglvl_par=0
!    phase_par = 12
!    CALL pardiso (pt_par,maxfct_par,mnum_par,mtype_par,phase_par,n_par,a_in_double,ia_in,ja_in,idum_par,nrhs_par,iparm_par,msglvl_par,ddum_par,ddum_par,error_par)
!    IF (error_par /= 0) THEN
!        if (symm) then
!            mtype_par = -2
!            print*,"trying indefinite frac..."
!            CALL pardiso (pt_par,maxfct_par,mnum_par,mtype_par,phase_par,n_par,a_in_double,ia_in,ja_in,idum_par,nrhs_par,iparm_par,msglvl_par,ddum_par,ddum_par,error_par)
!            IF (error_par /= 0) then
!                print*, 'ERROR - The following error was detected during analysis step of pardiso solver: ', error_par
!                stop
!            end if
!        ELSE
!            print*, 'ERROR - The following error was detected during analysis step of pardiso solver: ', error_par
!            stop
!        END IF
!    END IF
!    END subroutine pardiso_direct_1
!
!    subroutine pardiso_direct_2 (b_in_double,x_out_double)
!    implicit none
!    real(8),allocatable :: b_in_double(:,:),x_out_double(:,:)
!    phase_par = 33
!    CALL pardiso (pt_par,maxfct_par,mnum_par,mtype_par,phase_par,n_par,ddum_par,idum_par,idum_par,idum_par,size(b_in_double,2),iparm_par,msglvl_par,b_in_double,x_out_double,error_par)
!    IF (error_par/= 0) THEN
!        print*, 'ERROR - The following error was detected during solve step of pardiso solver: ', error_par
!        stop
!    END IF
!    END subroutine pardiso_direct_2
!
!    subroutine pardiso_direct_3
!    implicit none
!    phase_par = -1
!    CALL pardiso (pt_par,maxfct_par,mnum_par,mtype_par,phase_par,n_par,ddum_par,idum_par,idum_par,idum_par,nrhs_par,iparm_par,msglvl_par,ddum_par,ddum_par,error_par)
!    IF(ALLOCATED(pt_par)) DEALLOCATE(pt_par)
!    IF(ALLOCATED(iparm_par)) DEALLOCATE(iparm_par)
!    IF (error_par /= 0) THEN
!        print*, 'ERR - The following error was detected during terminal step of pardiso solver: ',error_par
!        stop
!    END IF
!    END subroutine pardiso_direct_3
!    subroutine pardiso_dealloc
!    ! decallocate all woring arrays of PARDISO after all calls to pardiso
!    implicit none
!    IF(ALLOCATED(pt)) DEALLOCATE(pt)
!    IF(ALLOCATED(iparm)) DEALLOCATE(iparm)
!    IF(ALLOCATED(mtype)) DEALLOCATE(mtype)
!    IF(ALLOCATED(n)) DEALLOCATE(n)
!    IF(ALLOCATED(nnz)) DEALLOCATE(nnz)
!    END subroutine pardiso_dealloc
!    ! pardiso_upper_checker
!    subroutine pardiso_checker (sym,n_in,nnz_in,ia_in,ja_in,checker)
!    USE MKL_SPARSE_HANDLE
!    implicit none
!    ! variables
!    logical,intent(in) :: sym
!    integer,intent(in) :: n_in,nnz_in
!    INTEGER,intent(in) :: ia_in(n_in+1),ja_in(nnz_in)
!    logical,intent(out) :: checker
!    TYPE(SPARSE_STRUCT) :: PT
!    integer :: CHECK_RESULT_CODE
!    ! process
!    CALL SPARSE_MATRIX_CHECKER_INIT(PT)
!    PT % N = n_in
!    PT % CSR_IA = LOC(ia_in)
!    PT % CSR_JA = LOC(ja_in)
!    PT % INDEXING         = MKL_ONE_BASED
!    if (sym) then
!        PT % MATRIX_STRUCTURE = MKL_UPPER_TRIANGULAR
!    else
!        PT % MATRIX_STRUCTURE = MKL_GENERAL_STRUCTURE
!    end if
!    PT % PRINT_STYLE      = MKL_C_STYLE
!    PT % MESSAGE_LEVEL    = MKL_PRINT
!    CHECK_RESULT_CODE = SPARSE_MATRIX_CHECKER(PT)
!    IF (CHECK_RESULT_CODE/=MKL_SPARSE_CHECKER_SUCCESS) THEN
!        WRITE(*,'("Matrix check details: (",i0,", ",i0,", ",i0,")")') PT % CHECK_RESULT
!        WRITE(*,*)'matrix check failed, stop!'
!        checker = .false.
!    else
!        checker = .true.
!    end if
!    end subroutine pardiso_checker
!    end module pardiso_multi_real
    
    module source_mod
    !source: ind  pulse_type (1 fenoglio 2 Gaussian pulse) t0 t1 current
    !   51     1    1.0e-7    0.0    1.0
    !   52     2    5.0e-9    0.0    1.e-11
    !   53     3    1.0e-7    0.0    1.0
    implicit none
    private
    public :: time_const_p
    contains
    subroutine time_const_p (pulse_type,t,t0,cur,tmp)
    implicit none
    integer,intent(in) :: pulse_type
    real(8),intent(in) :: t,t0,cur
    real(8),intent(out) :: tmp
    if (pulse_type==2) then
        tmp = deltafunc (1,t,t0) * cur
    else if (pulse_type==3) then
        tmp = diffdeltafunc (1,t,t0) * cur
    else if (pulse_type==4) then
        tmp = lighting (1,t) * cur
    else
        stop 'no such kind of source'
    end if
      end subroutine time_const_p
      

   
    real(8) function deltafunc (deriv,t,t0)
    implicit none
    integer,intent(in) :: deriv
    real(8),intent(in) :: t,t0
    if (deriv==0) then
        deltafunc = exp(-1.0D0/t0**2*(t*5.436563656918091D0-t0*2.718281828459046D&
            &0)**2)
    else if (deriv==1 ) then
            deltafunc = -1.0D0/t0**2*exp(-1.0D0/t0**2*(t*5.436563656918091D0-t0*2.718&
              &281828459046D0)**2)*(t*5.911244879144521D+1-t0*2.955622439572261D+&
              &1)
    else if (deriv==2) then
        deltafunc = 1.0D0/t0**2*exp(-1.0D0/t0**2*(t*5.436563656918091D0-t0*2.7182&
            &81828459046D0)**2)*(-5.911244879144521D+1)+1.0D0/t0**4*exp(-1.0D0/&
            &t0**2*(t*5.436563656918091D0-t0*2.718281828459046D0)**2)*(t*5.9112&
            &44879144521D+1-t0*2.955622439572261D+1)**2
    else if (deriv==3) then
        deltafunc = 1.0D0/t0**4*exp(-1.0D0/t0**2*(t*5.436563656918091D0-t0*2.7182&
            &81828459046D0)**2)*(t*5.911244879144521D+1-t0*2.955622439572261D+1&
            &)*5.911244879144521D+1+1.0D0/t0**4*exp(-1.0D0/t0**2*(t*5.436563656&
            &918091D0-t0*2.718281828459046D0)**2)*(t*6.988563204242466D+3-t0*3.&
            &494281602121233D+3)-1.0D0/t0**6*exp(-1.0D0/t0**2*(t*5.436563656918&
            &091D0-t0*2.718281828459046D0)**2)*(t*5.911244879144521D+1-t0*2.955&
            &622439572261D+1)**3
    else
        stop "deltafunc error"
    end if
    end function deltafunc
    
    real(8) function diffdeltafunc (deriv,t,t0)
    implicit none
    integer,intent(in) :: deriv
    real(8),intent(in) :: t,t0
    if (deriv==0) then
        diffdeltafunc = (t*exp(t**2*1.0D0/t0**2*3.141592653589793D0*(-4.0D0)))/t0
    else if (deriv==1 ) then
        diffdeltafunc = exp(t**2*1.0D0/t0**2*3.141592653589793D0*(-4.0D0))/t0-t**2*1.&
            &0D0/t0**3*3.141592653589793D0*exp(t**2*1.0D0/t0**2*3.1415926535897&
            &93D0*(-4.0D0))*8.0D0
    else
        stop "deltafunc error"
    end if
    end function diffdeltafunc
    
    real(8) function my (deriv,t,t1,t2)
    implicit none
    integer,intent(in) :: deriv
    real(8),intent(in) :: t,t1,t2
    if (deriv==0) then
        my = 1
    else if (deriv==1 ) then
        if (t<=t1) then
            my = 1.0D0/t1**2*exp(1.0D0/t1**2*(t-t1)**2*(-7.389056098930652D0))&
              &*(t*2.0D0-t1*2.0D0)*(-7.389056098930652D0)/(t1+t2)
        else
            my = 1.0D0/t2**2*exp(1.0D0/t2**2*(t-t1)**2*(-7.389056098930652D0))&
              &*(t*2.0D0-t1*2.0D0)*(-7.389056098930652D0)/(t1+t2)
        end if
    else
        stop "deltafunc error"
    end if
    end function my
    
    real(8) function lighting (deriv,t)
    implicit none
    integer,intent(in) :: deriv
    real(8),intent(in) :: t
    if (deriv==0) then
        lighting = 1
    else if (deriv==1 ) then
        if (t<=5.d-4) then
            lighting = 2000
        else
            lighting =  -exp(-(t-5.0D-4)/0.001)/0.001
        end if
    else
        stop "deltafunc error"
    end if
    end function lighting
    
    end module source_mod
    
    module loop_solution
    use MKL_SPBLAS
    implicit none
    contains
    
    function itoa(i) result(res)
	character(:),allocatable :: res
	integer,intent(in) :: i
	character(range(i)+2) :: tmp
	write(tmp,'(i0)') i
	res = trim(tmp)
    end function
    
    subroutine mkl_add_3 (num,av,ai_mkl,bv,bi_mkl,cv,ci_mkl,r_mkl)
    implicit none
    ! IN/OUT VARIABLES
    type(SPARSE_MATRIX_T),intent(in) :: ai_mkl,bi_mkl,ci_mkl
    integer,intent(in) :: num
    real(8),intent(in) :: av,bv,cv
    type(SPARSE_MATRIX_T),intent(out) :: r_mkl
    ! other variables
    integer :: info
    type(SPARSE_MATRIX_T) :: z_mkl,zt_mkl,tmp_mkl_1,tmp_mkl_2
    integer :: z_I(1),z_j(1)
    real(8) :: z_v(1)
    ! z_mkl
    z_i(1) = 1; z_j(1) = 1; z_v(1) = 0.d0
    info = mkl_sparse_d_create_coo (zt_mkl, SPARSE_INDEX_BASE_one, num, num, 1, z_i, z_j, z_v)
    if (info/=0) stop '2325'
    info = mkl_sparse_convert_csr (zt_mkl, SPARSE_OPERATION_NON_TRANSPOSE, z_mkl)
    if (info/=0) stop '2327'
    ! process
    info = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, Ai_mkl, av , z_mkl, tmp_mkl_1)
    if (info/=0) stop '2331'
    info = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, bi_mkl, bv , tmp_mkl_1, tmp_mkl_2)
    if (info/=0) stop '2333'
    info = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, ci_mkl, cv , tmp_mkl_2, r_mkl)
    if (info/=0) stop '2335'
    ! destroy
    info = mkl_sparse_destroy (z_mkl)
    info = mkl_sparse_destroy (zt_mkl)
    info = mkl_sparse_destroy (tmp_mkl_1)
    info = mkl_sparse_destroy (tmp_mkl_2)
    end subroutine mkl_add_3
    
    subroutine mkl_to_csr (ai_mkl,n_csr,NNZ_CSR,a_csr,i_csr,j_csr)
    USE, INTRINSIC :: ISO_C_BINDING, only : C_PTR,c_f_pointer
    implicit none
    ! IN/OUT VARIABLES
    type(SPARSE_MATRIX_T),intent(in) :: ai_mkl
    integer,intent(out) :: n_csr,NNZ_CSR
    integer,allocatable,intent(out) :: i_csr(:),j_csr(:)
    real(8),allocatable,intent(out) :: a_csr(:)
    ! LOCAL VARIABLES
    integer :: info_1,info_2
    TYPE(C_PTR) :: isys_c_1,isys_c_2,jsys_c_1,vsys_c_1
    integer,pointer :: isys_1(:),isys_2(:),jsys_1(:)
    real(8),pointer :: vsys_1(:)
    ! process
    info_1 = mkl_sparse_d_export_csr (ai_mkl, info_2, n_csr, n_csr, isys_c_1, isys_c_2, jsys_c_1, vsys_c_1)
    if (info_1/=0) stop
    call c_f_pointer(isys_c_1, isys_1,[n_csr])
    call c_f_pointer(isys_c_2, isys_2,[n_csr])
    allocate(i_csr(n_csr+1))
    i_csr (1:n_csr) = isys_1
    i_csr (2:n_csr+1) = isys_2
    NNZ_CSR = i_csr(n_csr+1)-1
    call c_f_pointer(jsys_c_1, jsys_1,[NNZ_CSR])
    call c_f_pointer(vsys_c_1, vsys_1,[NNZ_CSR])
    allocate (j_csr(NNZ_CSR),a_csr(NNZ_CSR))
    j_csr = jsys_1
    a_csr = vsys_1
    end subroutine mkl_to_csr
    

    
    real(8) function two_error (a,b)
    ! output: ||a-b||/||b||
    implicit none
    include "mkl_vml.f90"
    ! in/out variables
    real(8),intent(in),allocatable :: a(:),b(:)
    ! local variablea
    real(8),allocatable :: tmp_arr(:)
    real(8) :: val_1,val_2
    integer :: length
    ! process
    if (size(a)/=size(b)) stop "ERROR - array size do not consist, stop!"
    length = size(a)
    allocate(tmp_arr(length))
    call vdsub (length,a,b,tmp_arr)
    val_1 = two_norm (tmp_arr)
    val_2 = two_norm (b)
    two_error = val_1 / val_2
    end function two_error
    
    real(8) function two_norm (a)
    implicit none
    include "mkl_vml.f90"
    ! in/out variables
    real(8),intent(in),allocatable :: a(:)
    ! local variablea
    real(8),allocatable :: tmp_arr(:)
    integer :: length
    ! process
    length = size(a)
    allocate(tmp_arr(length))
    call vdsqr ( length, a ,tmp_arr )
    two_norm = sqrt(sum(tmp_arr))
    end function two_norm
    
    end module loop_solution