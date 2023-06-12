module msh_struct
implicit none
private
public mesh
!-----------------------------------------------------------------------
! defining kernel data structure
!-----------------------------------------------------------------------
type mesh
    real(8), dimension(:,:), allocatable :: nds
    integer :: totnds, order, dim
    integer :: totsurfs, surfshape, surfnds
    integer :: totvols, volshape, volnds
    integer, dimension(:,:), allocatable :: surfs, vols
    real(8) :: dx
end type

end module msh_struct

module msh_ops
implicit none
contains

subroutine gen_gridpnts3d(p,ln,n,msh)
!-----------------------------------------------------------------------
! Generate 3D grid points
!-----------------------------------------------------------------------
    use msh_struct
!-----------------------------------------------------------------------   
    integer, dimension(3), intent(IN) :: n
    real(8), dimension(3), intent(IN) :: p, ln
!-----------------------------------------------------------------------
    type(mesh), intent(OUT) :: msh
!-----------------------------------------------------------------------
    integer :: i, j, k, l
    real(8) :: dx, dy, dz
!-----------------------------------------------------------------------

    msh%dim = 3
    dx = ln(1)/(n(1)-1)
    dy = ln(2)/(n(2)-1)
    dz = ln(3)/(n(3)-1)
    msh%totnds = n(1)*n(2)*n(3);
    allocate(msh%nds(msh%totnds+1,3))
    msh%nds(1,1) = p(1)
    msh%nds(1,2) = p(2)
    msh%nds(1,3) = p(3)
    l = 1;
    do k=1,n(3)
        do j=1,n(2)
            do i=1,n(1)
                msh%nds(l+1,1) = msh%nds(l,1) + dx;
                msh%nds(l+1,2) = msh%nds(l,2);
                msh%nds(l+1,3) = msh%nds(l,3)
                l = l+1;
            end do
            msh%nds(l,1) = msh%nds(1,1)
            msh%nds(l,2) = msh%nds(l,2) + dy
            msh%nds(l,3) = msh%nds(l,3)
        end do
        msh%nds(l,1) = p(1)
        msh%nds(l,2) = p(2)
        msh%nds(l,3) = msh%nds(l,3) + dz
    end do
    write(*,'(a,i0)') "totnds = ", msh%totnds

end subroutine gen_gridpnts3d

subroutine gen_quadgrid(p, ln, n, msh)
!-----------------------------------------------------------------------
! Generate first-order quadrilateral grid
!-----------------------------------------------------------------------
    use msh_struct
!-----------------------------------------------------------------------
    integer, dimension(3), intent(IN) :: n
    real(8), dimension(3), intent(IN) :: p, ln
!-----------------------------------------------------------------------
    type(mesh), intent(OUT) :: msh
!-----------------------------------------------------------------------
    integer :: elems_x, elems_y, i, j, k, l
    integer, dimension(:), allocatable :: nd
    real(8) :: dx, dy
!-----------------------------------------------------------------------

    msh%dim = 2
    dx = ln(1)/(n(1)-1)
    dy = ln(2)/(n(2)-1)
    msh%totnds = n(1)*n(2);
    elems_x = n(1)-1;
    elems_y = n(2)-1;
    msh%order = 1;
    msh%surfshape = 2;
    msh%surfnds=4;
    msh%totsurfs = elems_x * elems_y
    allocate(msh%nds(msh%totnds+1,3))
    allocate(nd(msh%totnds))
    allocate(msh%surfs(msh%totsurfs+1,msh%surfnds))
    
    msh%nds(1,1) = p(1)
    msh%nds(1,2) = p(2)
    msh%nds(1,3) = 0
    l = 1;
    do j=1,n(2)
        do i=1,n(1)
            msh%nds(l+1,1) = msh%nds(l,1) + dx;
            msh%nds(l+1,2) = msh%nds(l,2);
            msh%nds(l+1,3) = 0
            l = l+1;
        end do
        msh%nds(l,1) = msh%nds(1,1)
        msh%nds(l,2) = msh%nds(l,2) + dy
        msh%nds(l,3) = 0
    end do

    do i=1,msh%totnds
       nd(i) = i
    end do
    k = 1;
    l = 1;
    do j = 1,elems_y
        do i = 1,elems_x
            msh%surfs(l,1) = nd(k)
            msh%surfs(l,2) = nd(k+1)
            msh%surfs(l,3) = nd(k+elems_x+2)
            msh%surfs(l,4) = nd(k+elems_x+1)
            k = k+1
            l = l+1
        end do
        k = k+1
    end do
    write(*,1) "totnds = ", msh%totnds
    write(*,1) "totelems = ", msh%totsurfs
    write(*,1) "elemnds = ", msh%surfnds
    write(*,1) "order = ", msh%order
    write(*,1) "elemshape = ", msh%surfshape
    1 format(a,i0)
end subroutine gen_quadgrid

subroutine gen_hexgrid(p,ln,n,msh)
!-----------------------------------------------------------------------
! Generate first_order hexahedral grid
!-----------------------------------------------------------------------
    use msh_struct
!-----------------------------------------------------------------------
    integer, dimension(3), intent(IN) :: n
    real(8), dimension(3), intent(IN) :: p, ln
!-----------------------------------------------------------------------
    type(mesh), intent(OUT) :: msh
!-----------------------------------------------------------------------
    integer,dimension(:,:,:), allocatable :: ndnum
    integer :: i, j, k, l, e, s
!-----------------------------------------------------------------------

    call gen_gridpnts3d(p,ln,n,msh)
    msh%dim = 3
    msh%order = 1
    msh%volshape = 2
    msh%volnds = 8
    msh%surfshape = 2
    msh%surfnds = 4

    msh%totvols = (n(1)-1)*(n(2)-1)*(n(3)-1)
    msh%totsurfs = msh%totvols*6

    allocate(ndnum(n(1),n(2),n(3)))
    allocate(msh%surfs(msh%totsurfs,msh%surfnds))
    allocate(msh%vols(msh%totvols,msh%volnds))

    l = 1
    do k=1,n(3)
        do j=1,n(2)
            do i=1,n(1)
                ndnum(i,j,k) = l
                l = l + 1
            end do
        end do
    end do

    e = 1
    do k=1,n(3)-1
        do j=1,n(2)-1
            do i=1,n(1)-1
                msh%vols(e,1) = ndnum(i,j,k)
                msh%vols(e,2) = ndnum(i+1,j,k)
                msh%vols(e,3) = ndnum(i+1,j+1,k)
                msh%vols(e,4) = ndnum(i,j+1,k)
                msh%vols(e,5) = ndnum(i,j,k+1)
                msh%vols(e,6) = ndnum(i+1,j,k+1)
                msh%vols(e,7) = ndnum(i+1,j+1,k+1)
                msh%vols(e,8) = ndnum(i,j+1,k+1)
                e = e + 1
            end do
        end do
    end do

    s = 1
    do i = 1,msh%totvols
        msh%surfs(s,1) = msh%vols(i,1)
        msh%surfs(s,2) = msh%vols(i,2)
        msh%surfs(s,3) = msh%vols(i,3)
        msh%surfs(s,4) = msh%vols(i,4)

        msh%surfs(s+1,1) = msh%vols(i,5)
        msh%surfs(s+1,2) = msh%vols(i,6)
        msh%surfs(s+1,3) = msh%vols(i,7)
        msh%surfs(s+1,4) = msh%vols(i,8)

        msh%surfs(s+2,1) = msh%vols(i,2)
        msh%surfs(s+2,2) = msh%vols(i,6)
        msh%surfs(s+2,3) = msh%vols(i,7)
        msh%surfs(s+2,4) = msh%vols(i,3)

        msh%surfs(s+3,1) = msh%vols(i,1)
        msh%surfs(s+3,2) = msh%vols(i,5)
        msh%surfs(s+3,3) = msh%vols(i,8)
        msh%surfs(s+3,4) = msh%vols(i,4)

        msh%surfs(s+4,1) = msh%vols(i,1)
        msh%surfs(s+4,2) = msh%vols(i,2)
        msh%surfs(s+4,3) = msh%vols(i,6)
        msh%surfs(s+4,4) = msh%vols(i,5)

        msh%surfs(s+5,1) = msh%vols(i,4)
        msh%surfs(s+5,2) = msh%vols(i,3)
        msh%surfs(s+5,3) = msh%vols(i,7)
        msh%surfs(s+5,4) = msh%vols(i,8)

        s = s+6
    end do
end subroutine gen_hexgrid

subroutine gen_tetgrid(p,ln,n,msh)
!-----------------------------------------------------------------------
! Generate first_order tetrahedral grid
!-----------------------------------------------------------------------
    use msh_struct
!-----------------------------------------------------------------------
    integer, dimension(3), intent(IN) :: n
    real(8), dimension(3), intent(IN) :: p, ln
!-----------------------------------------------------------------------
    type(mesh), intent(OUT) :: msh
!-----------------------------------------------------------------------
    integer,dimension(:,:,:), allocatable :: ndnum
    integer :: i, j, k, l, e, s
    integer :: tets_per_cube
!-----------------------------------------------------------------------

    call gen_gridpnts3d(p,ln,n,msh)
    msh%dim = 3
    msh%order = 1
    msh%volshape = 1
    msh%volnds = 4
    msh%surfshape = 1
    msh%surfnds = 3
    tets_per_cube = 5
    select case (tets_per_cube)
        case (5)
            msh%totvols =  5*(n(1)-1)*(n(2)-1)*(n(3)-1)
        case (6)
            msh%totvols =  6*(n(1)-1)*(n(2)-1)*(n(3)-1)
    end select
 
    msh%totsurfs = msh%totvols*4;

    allocate(ndnum(n(1),n(2),n(3)))
    allocate(msh%surfs(msh%totsurfs,msh%surfnds))
    allocate(msh%vols(msh%totvols,msh%volnds))

    l = 1
    do k=1,n(3)
        do j=1,n(2)
            do i=1,n(1)
                ndnum(i,j,k) = l
                l = l + 1
            end do
        end do
    end do
    
    e = 1
    select case (tets_per_cube)
        case (5)
            do k = 1, n(3)-1
                do j = 1, n(2)-1
                    do i = 1, n(1)-1
                        msh%vols(e,1) = ndnum(i,j,k)
                        msh%vols(e,2) = ndnum(i+1,j,k)
                        msh%vols(e,3) = ndnum(i,j+1,k)
                        msh%vols(e,4) = ndnum(i,j,k+1)

                        msh%vols(e+1,1) = ndnum(i+1,j,k)
                        msh%vols(e+1,2) = ndnum(i+1,j+1,k)
                        msh%vols(e+1,3) = ndnum(i,j+1,k)
                        msh%vols(e+1,4) = ndnum(i+1,j+1,k+1)

                        msh%vols(e+2,1) = ndnum(i,j,k+1)
                        msh%vols(e+2,2) = ndnum(i,j+1,k)
                        msh%vols(e+2,3) = ndnum(i,j+1,k+1)
                        msh%vols(e+2,4) = ndnum(i+1,j+1,k+1)

                        msh%vols(e+3,1) = ndnum(i,j,k+1)
                        msh%vols(e+3,2) = ndnum(i+1,j,k)
                        msh%vols(e+3,3) = ndnum(i+1,j+1,k+1)
                        msh%vols(e+3,4) = ndnum(i+1,j,k+1)

                        msh%vols(e+4,1) = ndnum(i,j+1,k)
                        msh%vols(e+4,2) = ndnum(i+1,j,k)
                        msh%vols(e+4,3) = ndnum(i,j,k+1)
                        msh%vols(e+4,4) = ndnum(i+1,j+1,k+1)

                        e = e+5;
                    end do
                end do
            end do
        case (6)
            do k = 1, n(3)-1
                do j = 1, n(2)-1
                    do i = 1, n(1)-1
                        msh%vols(e,1) = ndnum(i,j,k);
                        msh%vols(e,2) = ndnum(i+1,j,k);
                        msh%vols(e,3) = ndnum(i+1,j+1,k);
                        msh%vols(e,4) = ndnum(i+1,j,k+1);

                        msh%vols(e+1,1) = ndnum(i,j,k);
                        msh%vols(e+1,2) = ndnum(i+1,j,k+1);
                        msh%vols(e+1,3) = ndnum(i,j,k+1);
                        msh%vols(e+1,4) = ndnum(i,j+1,k+1);

                        msh%vols(e+2,1) = ndnum(i,j+1,k);
                        msh%vols(e+2,2) = ndnum(i+1,j+1,k);
                        msh%vols(e+2,3) = ndnum(i+1,j+1,k+1);
                        msh%vols(e+2,4) = ndnum(i+1,j,k+1);

                        msh%vols(e+3,1) = ndnum(i,j,k);
                        msh%vols(e+3,2) = ndnum(i+1,j+1,k);
                        msh%vols(e+3,3) = ndnum(i,j+1,k);
                        msh%vols(e+3,4) = ndnum(i+1,j,k+1);

                        msh%vols(e+4,1) = ndnum(i,j,k);
                        msh%vols(e+4,2) = ndnum(i,j+1,k);
                        msh%vols(e+4,3) = ndnum(i,j+1,k+1);
                        msh%vols(e+4,4) = ndnum(i+1,j,k+1);

                        msh%vols(e+5,1) = ndnum(i,j+1,k);
                        msh%vols(e+5,2) = ndnum(i+1,j+1,k+1);
                        msh%vols(e+5,3) = ndnum(i,j+1,k+1);
                        msh%vols(e+5,4) = ndnum(i+1,j,k+1);
                        e = e+6;

                    end do
                end do
            end do
    end select

    s = 1
    do i = 1, msh%totvols
        msh%surfs(s,1) = msh%vols(i,1)
        msh%surfs(s,2) = msh%vols(i,2)
        msh%surfs(s,3) = msh%vols(i,3)

        msh%surfs(s+1,1) = msh%vols(i,1)
        msh%surfs(s+1,2) = msh%vols(i,2)
        msh%surfs(s+1,3) = msh%vols(i,4)

        msh%surfs(s+2,1) = msh%vols(i,2)
        msh%surfs(s+2,2) = msh%vols(i,3)
        msh%surfs(s+2,3) = msh%vols(i,4)

        msh%surfs(s+3,1) = msh%vols(i,1)
        msh%surfs(s+3,2) = msh%vols(i,3)
        msh%surfs(s+3,3) = msh%vols(i,4)
        s = s + 4
    end do
end subroutine gen_tetgrid

subroutine gen_trigrid(p,ln,n,msh)
!-----------------------------------------------------------------------
! Generate first_order triangular grid
!-----------------------------------------------------------------------
    use msh_struct
!-----------------------------------------------------------------------
    integer, dimension(3), intent(IN) :: n
    real(8), dimension(3), intent(IN) :: p, ln
!-----------------------------------------------------------------------
    type(mesh), intent(OUT) :: msh
!-----------------------------------------------------------------------
    integer :: elems_x, elems_y, i, j, k, l
    integer, dimension(:), allocatable :: nd
    real(8) :: dx, dy
!-----------------------------------------------------------------------

    msh%dim = 2
    dx = ln(1)/(n(1)-1)
    dy = ln(2)/(n(2)-1)
    elems_x = n(1)-1
    elems_y = n(2)-1
    msh%order = 1
    msh%surfshape = 1
    msh%surfnds = 3
    msh%totsurfs = 2*elems_x*elems_y
    allocate(msh%surfs(msh%totsurfs, msh%surfnds))

    msh%totnds = (elems_x+1)*(elems_y+1)
    allocate(msh%nds(msh%totnds+1,3))
    msh%nds(1,1) = p(1)
    msh%nds(1,2) = p(2)
    msh%nds(1,3) = 0
    l = 1;
    do j=1,n(2)
        do i=1,n(1)
            msh%nds(l+1,1) = msh%nds(l,1) + dx;
            msh%nds(l+1,2) = msh%nds(l,2);
            msh%nds(l+1,3) = 0
            l = l+1;
        end do
        msh%nds(l,1) = msh%nds(1,1)
        msh%nds(l,2) = msh%nds(l,2) + dy
        msh%nds(l,3) = 0
    end do

    allocate(nd(msh%totnds))
    do i = 1, msh%totnds
        nd(i) = i
    end do

    k = 1
    l = 1
    do j = 1, elems_y
        do i = 1, elems_x
            if ((mod(i,2) == 0 .and. mod(j,2) /= 0) .or. (mod(i,2) /= 0 .and. mod(j,2) == 0)) then
                msh%surfs(l,1) = nd(k)
                msh%surfs(l,2) = nd(k+1)
                msh%surfs(l,3) = nd(k+elems_x+1)

                msh%surfs(l+1,1) = nd(k+1)
                msh%surfs(l+1,2) = nd(k+elems_x+2)
                msh%surfs(l+1,3) = nd(k+elems_x+1)
            else  if ((mod(i,2) /= 0 .and. mod(j,2) /= 0) .or. (mod(i,2) == 0 .and. mod(j,2)== 0)) then
                msh%surfs(l,1) = nd(k)
                msh%surfs(l,2) = nd(k+1)
                msh%surfs(l,3) = nd(k+elems_x+2)

                msh%surfs(l+1,1) = nd(k)
                msh%surfs(l+1,2) = nd(k+elems_x+2)
                msh%surfs(l+1,3) = nd(k+elems_x+1)
            end if
            k = k+1
            l = l+2
        end do
        k = k+1
    end do
end subroutine

subroutine gen_trigrid2(p, ln, n, msh)
!-----------------------------------------------------------------------
! Generate first_order triangular grid
!-----------------------------------------------------------------------
    use msh_struct
!-----------------------------------------------------------------------
    integer, dimension(3), intent(IN) :: n
    real(8), dimension(3), intent(IN) :: p, ln
!-----------------------------------------------------------------------
    type(mesh), intent(OUT) :: msh
!-----------------------------------------------------------------------
    integer :: elems_x, elems_y, i, j, k, l
    integer, dimension(:), allocatable :: nd
    real(8) :: dx, dy
!-----------------------------------------------------------------------

    msh%dim = 2
    dx = ln(1)/(n(1)-1)
    dy = ln(2)/(n(2)-1)
    elems_x = n(1)-1
    elems_y = n(2)-1
    msh%order = 1
    msh%surfshape = 1
    msh%surfnds = 3
    msh%totsurfs = 2*elems_x*elems_y
    allocate(msh%surfs(msh%totsurfs, msh%surfnds))

    msh%totnds = (elems_x+1)*(elems_y+1)
    allocate(msh%nds(msh%totnds+1,3))
    msh%nds(1,1) = p(1)
    msh%nds(1,2) = p(2)
    msh%nds(1,3) = 0
    l = 1;
    do j=1,n(2)
        do i=1,n(1)
            msh%nds(l+1,1) = msh%nds(l,1) + dx;
            msh%nds(l+1,2) = msh%nds(l,2);
            msh%nds(l+1,3) = 0
            l = l+1;
        end do
        msh%nds(l,1) = msh%nds(1,1)
        msh%nds(l,2) = msh%nds(l,2) + dy
        msh%nds(l,3) = 0
    end do

    allocate(nd(msh%totnds))
    do i = 1, msh%totnds
        nd(i) = i
    end do

    k = 1
    l = 1
    do j = 1, elems_y
        do i = 1, elems_x
            msh%surfs(l,1) = nd(k)
            msh%surfs(l,2) = nd(k+1)
            msh%surfs(l,3) = nd(k+elems_x+1)

            msh%surfs(l+1,1) = nd(k+1)
            msh%surfs(l+1,2) = nd(k+elems_x+2)
            msh%surfs(l+1,3) = nd(k+elems_x+1)
            l = l+2;
            k = k+1;
        end do
        k = k+1;
    end do
end subroutine gen_trigrid2

subroutine prnt_vtk(u, msh, fname, t)
!-----------------------------------------------------------------------
! Print field variable u and mesh in VTK format
!-----------------------------------------------------------------------
    use msh_struct
!-----------------------------------------------------------------------
    type(mesh), intent(IN) :: msh                  
    integer, intent(IN) :: t
    real(8),dimension(:), intent(IN) :: u
    character(len=50), intent(IN) :: fname
!-----------------------------------------------------------------------
    integer :: i, j, sum
    character(len=50) :: file_name
    character(len=10) :: t_s
    integer :: elemorder
    character*29 :: path='sim/out/'
!-----------------------------------------------------------------------
    call system('mkdir -p ' // adjustl(trim(path)) // trim(adjustl(fname)) // "_vtk/")
    !Construct the filename:
    write(t_s, '(i0)') t
    file_name = trim(adjustl(path)) // trim(adjustl(fname)) // "_vtk/" // trim(adjustl(fname))  // trim(adjustl(t_s)) // ".vtk"
   
    elemorder = 1
    if (msh%surfnds == 6 .or. msh%surfnds == 8) then
        elemorder = 2
    end if
    sum = ((msh%surfnds/elemorder)+1)*msh%totsurfs;

    write(*,'(a)') "Printing results to: " // file_name
    open(UNIT=1, FILE=file_name)
    write(1,'(a)') "# vtk DataFile Version 1.0"
    write(1,'(a)') "Cube example"
    write(1,'(a)') "ASCII"
    write(1,*)
    write(1,'(a)') "DATASET POLYDATA"
    write(1,1) "POINTS ",msh%totnds," double"
    1 format(A,i0,A)
     
    do i = 1,msh%totnds
      do j = 1,3
       write(1,'(2X,F12.6)',advance='no') msh%nds(i,j)      
      end do
       write(1,*)
    end do
    write(1,*) 

    write(1,2) "POLYGONS ",msh%totsurfs," ",sum
    2 format(A,i0,A,i0)
    do i = 1,msh%totsurfs
        write(1,'(i0)',advance='no') msh%surfnds/elemorder
        do j = 1,msh%surfnds/elemorder
            write(1,'(4X,i0)',advance='no') msh%surfs(i,j)-1
        end do
        write(1,*)
    end do
    write(1,*)
    
    write(1,3) "POINT_DATA ",msh%totnds
    3 format(A,i0)
    write(1,'(a)') "SCALARS myscalars double"
    write(1,'(a)') "LOOKUP_TABLE custom_table"
    do i=1,msh%totnds
        write(1,'(F18.6)',advance='yes') u(i)
    end do
    close(unit=1)
end subroutine prnt_vtk

subroutine prnt_vec_vtk(vel, msh, fname, t)
!-----------------------------------------------------------------------
! Print vectors vel and mesh in VTK format
!-----------------------------------------------------------------------
    use msh_struct
!-----------------------------------------------------------------------
    type(mesh), intent(IN) :: msh                  
    integer, intent(IN) :: t
    real(8),dimension(:,:), intent(IN) :: vel
    character(len=50), intent(IN) :: fname
!-----------------------------------------------------------------------
    integer :: i, j, sum
    character(len=50) :: file_name
    character(len=10) :: t_s
    character*29 :: path='sim/out/'
    real(8), dimension(:), allocatable :: v
!-----------------------------------------------------------------------
    call system('mkdir -p ' // adjustl(trim(path)) // trim(adjustl(fname)) // "_vtk/")
    write(t_s, '(i0)') t
    file_name = trim(adjustl(path)) // trim(adjustl(fname)) // "_vtk/" // trim(adjustl(fname))  // trim(adjustl(t_s)) // ".vtk"
   
    sum = (msh%surfnds+1)*msh%totsurfs

    allocate(v(msh%totnds))
    v = sqrt(vel(:,1)*vel(:,1)+vel(:,2)*vel(:,2)+vel(:,3)*vel(:,3))

    write(*,'(a)') "Printing results to: " // file_name
    open(UNIT=1, FILE=file_name)
    write(1,'(a)') "# vtk DataFile Version 1.0"
    write(1,'(a)') "Cube example"
    write(1,'(a)') "ASCII"
    write(1,*)
    write(1,'(a)') "DATASET POLYDATA"
    write(1,1) "POINTS ",msh%totnds," double"
    1 format(A,i0,A)

    do i = 1,msh%totnds
    do j = 1,3
    write(1,'(2X,F12.6)',advance='no') msh%nds(i,j)      
    end do
       write(1,*)
    end do
    write(1,*) 
    
    write(1,2) "POLYGONS ",msh%totsurfs," ",sum
    2 format(A,i0,A,i0)
    do i = 1,msh%totsurfs
        write(1,'(i0)',advance='no') msh%surfnds
        do j = 1,msh%surfnds
            write(1,'(4X,i0)',advance='no') msh%surfs(i,j)-1
        end do
        write(1,*)
    end do
    write(1,*)

    write(1,3) "POINT_DATA ",msh%totnds
    3 format(A,i0)
    write(1,'(a)') "SCALARS myscalars double"
    write(1,'(a)') "LOOKUP_TABLE custom_table"
    do i=1,msh%totnds
        write(1,'(F18.6)',advance='yes') v(i)
    end do
    write(1,*)

    write(1,'(a)') "vec vec double"
    do i = 1, msh%totnds
        write(1,'(F18.6, a, F18.6, a, F18.6)') vel(i,1), " ", vel(i,2), " ", vel(i,3)
    end do

    close(unit=1)
end subroutine prnt_vec_vtk

subroutine prnt_pnts_txt(u, nds, totnds, fname, t)
!-----------------------------------------------------------------------
! Print field variable u and mesh nodes  TXT format
!-----------------------------------------------------------------------
    integer, intent(IN) :: totnds, t
    real(8),dimension(:), intent(IN) :: u
    real(8), dimension(:,:), intent(IN) :: nds
    character(len=50), intent(IN) :: fname
!-----------------------------------------------------------------------
    integer:: i, j
    character(len=50) :: file_name
    character(len=10) :: t_s
    character*29 :: path='sim/out/'
!-----------------------------------------------------------------------

    call system('mkdir -p ' // adjustl(trim(path)) // trim(adjustl(fname)) // "_txt/")
    write(t_s, '(i0)') t
    file_name = trim(adjustl(path)) // trim(adjustl(fname)) // "_txt/" // trim(adjustl(fname))  // trim(adjustl(t_s)) // ".txt"
    
    write(*,'(a)') "Printing results to: " // file_name
    open(UNIT=1, FILE=file_name)
    
    do i = 1,totnds
      do j = 1,3
       write(1,'(2X,F12.6)',advance='no') nds(i,j)      
      end do
       write(1,'(F12.6)',advance='no') u(i)
       write(1,*)
    end do
    write(1,*) 
    close(unit=1)
end subroutine prnt_pnts_txt

subroutine read_pnts(mesh_fname, ps)
!-----------------------------------------------------------------------
! Read mesh nodes
!-----------------------------------------------------------------------
    use msh_struct
!-----------------------------------------------------------------------
    character(len=50), intent(IN) :: mesh_fname
!-----------------------------------------------------------------------
    type(mesh), intent(OUT) :: ps
!-----------------------------------------------------------------------
    real(8), dimension(:,:), allocatable :: read_txt
    integer :: k
    real(8) :: maxy, maxz
    character(len=50) :: file_name
!-----------------------------------------------------------------------
    file_name = trim(adjustl(mesh_fname))

    open(unit=2,file=file_name,access='sequential',status='old',action='read') !read from the existing file
    allocate (read_txt(1,2))
    read(2, *) read_txt(1,:)
    ps%totnds = int(read_txt(1,1))
    allocate(ps%nds(ps%totnds,3))
    maxy = 0.0
    maxz = 0.0
    do k = 1,ps%totnds
        read(2,*) ps%nds(k,:)
        if (abs(ps%nds(k,2)) >= maxy) maxy = ps%nds(k,2)
        if (abs(ps%nds(k,3)) >= maxz) maxz= ps%nds(k,3)
    end do
    write(*,'(a, i0)') "Total Points = ", ps%totnds
    if (maxy == 0.0 .and. maxz == 0.0) then
        ps%dim = 1
    else if (maxy > 0.0 .and. maxz == 0.0) then
        ps%dim = 2
    else
        ps%dim = 3
    end if
    write(*,'(a, i0)') "Spatial dim = ", ps%dim
    close(unit=2)
end subroutine read_pnts

subroutine read_surfs(mesh_fname, ps)
!-----------------------------------------------------------------------
! Read 2D mesh elements (surfaces)
!-----------------------------------------------------------------------
    use msh_struct
!-----------------------------------------------------------------------
    character(len=50), intent(IN) :: mesh_fname
!-----------------------------------------------------------------------
    type(mesh), intent(INOUT) :: ps
!-----------------------------------------------------------------------
    real(8), dimension(:,:), allocatable :: read_txt
    integer :: k
    character(len=50) :: file_name
!-----------------------------------------------------------------------
    file_name = trim(adjustl(mesh_fname))
    
    open(unit=2,file=file_name,access='sequential',status='old',action='read') !read from the existing file
    allocate (read_txt(1,4))
    read(2, *) read_txt(1,:)
    ps%totsurfs = int(read_txt(1,1))
    ps%surfnds = int(read_txt(1,2))
    ps%surfshape = int(read_txt(1,3))
    ps%order = int(read_txt(1,4))

    allocate(ps%surfs(ps%totsurfs,ps%surfnds))
    do k = 1,ps%totsurfs
        read(2,*) ps%surfs(k,:)
    end do
    write(*,'(a, i0)') "Total Elements = ", ps%totsurfs
    write(*,'(a, i0)') "Element Nds = ", ps%surfnds
    write(*,'(a, i0)') "Element Shape = ", ps%surfshape
    write(*,'(a, i0)') "Element Order = ", ps%order
    close(unit=2)
end subroutine read_surfs

subroutine read_vols(mesh_fname, ps)
!-----------------------------------------------------------------------
! Read 3D mesh elements (volumes)
!-----------------------------------------------------------------------
    use msh_struct
!-----------------------------------------------------------------------
    character(len=50), intent(IN) :: mesh_fname
!-----------------------------------------------------------------------
    type(mesh), intent(INOUT) :: ps
!-----------------------------------------------------------------------
    real(8), dimension(:,:), allocatable :: read_txt
    integer :: k
    character(len=50) :: file_name
!-----------------------------------------------------------------------

    file_name = trim(adjustl(mesh_fname))

    open(unit=2,file=file_name,access='sequential',status='old',action='read') !read from the existing file
    allocate (read_txt(1,4))
    read(2, *) read_txt(1,:)
    ps%totvols = int(read_txt(1,1))
    ps%volnds = int(read_txt(1,2))
    ps%volshape = int(read_txt(1,3))
    ps%order = int(read_txt(1,4))

    allocate(ps%vols(ps%totvols,ps%volnds))
    do k = 1,ps%totvols
        read(2,*) ps%vols(k,:)
    end do
    write(*,'(a, i0)') "Total Elements = ", ps%totvols
    write(*,'(a, i0)') "Element Nds = ", ps%volnds
    write(*,'(a, i0)') "Element Shape = ", ps%volshape
    write(*,'(a, i0)') "Element Order = ", ps%order
    close(unit=2)
end subroutine read_vols

subroutine prnt_pltctrl(pnts, nt, prnt_freq)
!-----------------------------------------------------------------------
! Print plotting controls to a text file
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: nt, prnt_freq
!-----------------------------------------------------------------------
    real(8) :: xmin, xmax, ymin, ymax
    character(len=50) :: fname
    character*29 :: path = 'sim/out/'
!-----------------------------------------------------------------------

    !Construct the filename:
    fname = trim(adjustl(path)) // "pltctrl.txt"
    xmin = minval(pnts(:,1))
    ymin = minval(pnts(:,2))
    xmax = maxval(pnts(:,1))
    ymax = maxval(pnts(:,2))
    
    open(UNIT=1, FILE=fname)
    write(1,'(i0, a, i0, a, f10.3, a, f10.3, a, f10.3, a, f10.3)') nt, " ", prnt_freq, " ", xmin, " ", xmax, " ", ymin, " ", ymax
    close(UNIT=1)
end subroutine prnt_pltctrl

subroutine prnt_elems_txt(elems, totelems)
!-----------------------------------------------------------------------
! Print element connectivity to a text file
!-----------------------------------------------------------------------
    integer, dimension(:,:), intent(in) :: elems
    integer, intent(in):: totelems
!-----------------------------------------------------------------------
    integer :: i
    character(len=50) :: fname
    character*29 :: path='sim/out/'
!-----------------------------------------------------------------------

    !Construct the filename:
    fname = trim(adjustl(path)) // "elems.txt"
    
    open(UNIT=1, FILE=fname)
    !nd number is minus 1 to be able to plot it in matplotlib
    do i = 1, totelems
        write(1, '(i0, a, i0, a, i0)') elems(i,1)-1, " ", elems(i,2)-1, " ", elems(i,3)-1
    end do
    close(UNIT=1)
end subroutine prnt_elems_txt

subroutine get_avg_dx(pnts, avg_dx)
!-------------------------------------------------------------------------
    use omp_lib
!-------------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: pnts
!-------------------------------------------------------------------------
    real(8), intent(out) :: avg_dx
!-------------------------------------------------------------------------
    real(8) :: d_min, dminf = 1e6, dx, dy, dz, d
    integer :: i, j
    real(8) :: sum_d = 0
    integer :: totsamples = 50
!-------------------------------------------------------------------------
    do i=1, totsamples
        d_min = 1e6
        do j = 1, TotSamples
            dx = pnts(i, 1) - pnts(j, 1)
            dy = pnts(i, 2) - pnts(j, 2)
            dz = pnts(i, 3) - pnts(j, 3)
            d = sqrt(dx * dx + dy * dy + dz * dz)
            if (d <= d_min .and. d > 1.0e-6) d_min = d
        end do
        if (d_min <= dminf) dminf = d_min
        sum_d = sum_d + d_min
    end do
    !avg_dx = sum_d / TotSamples;
    avg_dx = dminf;
    write(*,'(a,f10.3)') "Average dx = ", avg_dx
end subroutine get_avg_dx
end module msh_ops


