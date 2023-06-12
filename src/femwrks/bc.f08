module bc_struct
implicit none
private
public bc
!-----------------------------------------------------------------------
! defining kernel data structure
!-----------------------------------------------------------------------
type bc
    integer, dimension(:), allocatable :: pnts
    real(8), dimension(:), allocatable :: vals
    integer :: tot
end type

end module bc_struct

module bc_ops
implicit none
contains

subroutine set_dbc(pnts, totpnts, edge, loc, value, bc_out)
!-----------------------------------------------------------------------
! Collecting boundary-condition nodes in a specified location
! and assigning a value to them
!-----------------------------------------------------------------------
    use bc_struct
    implicit none
    real(8), dimension(:,:), intent(IN) :: pnts
    integer,intent(IN) :: totpnts, edge
    real(8), intent(IN) :: loc, value
!-----------------------------------------------------------------------
    type(bc), intent(INOUT) :: bc_out
!-----------------------------------------------------------------------
    integer :: i, m = 0
    logical :: cond1, cond2
    integer, dimension(:), allocatable :: ndnum
!-----------------------------------------------------------------------

    allocate(ndnum(0))
    write(*,'(a)') "Setting grid direchlet boundary conditions..."
    !x axis
    if (edge == 1) then
        m = 0
        do i=1,totpnts
            cond1 = pnts(i,1) >= loc-0.01*abs(loc)
            cond2 = pnts(i,1) <= loc+0.01*abs(loc)
            if (cond1 .and. cond2) then
                m = m+1
                ndnum = [ndnum, i]
            end if
        end do
    end if
    !y axis
    if (edge == 2) then
        m = 0
        do i=1,totpnts
            cond1 = pnts(i,2) >= loc-0.01*abs(loc)
            cond2 = pnts(i,2) <= loc+0.01*abs(loc)
            if (cond1 .and. cond2) then
                m = m+1
                ndnum = [ndnum, i]
            end if
        end do
    end if
    !z axis
    if (edge == 3) then
        m = 0
        do i=1,totpnts
            cond1 = pnts(i,3) >= loc-0.01*abs(loc)
            cond2 = pnts(i,3) <= loc+0.01*abs(loc)
            if (cond1 .and. cond2) then
                m = m+1
                ndnum = [ndnum, i]
            end if
        end do
    end if
    do i = 1,m
        bc_out%tot = bc_out%tot + 1
        bc_out%pnts = [bc_out%pnts, ndnum(i)]
        bc_out%vals = [bc_out%vals, value]
    end do

end subroutine set_dbc

subroutine read_dbc(bc_fname, value, mshpnts, mshtotpnts, bc_out)
!-----------------------------------------------------------------------
! Reading the set of boundary-condition nodes from an input file
!-----------------------------------------------------------------------
    use msh_ops
    use msh_struct
    use bc_struct
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: mshpnts
    integer, intent(IN) :: mshtotpnts
    real(8), intent(IN) :: value
    character(len=50), intent(IN) :: bc_fname
!-----------------------------------------------------------------------
    type(bc), intent(INOUT) :: bc_out
!-----------------------------------------------------------------------
    integer :: i, j
    real(8) :: dx, dy, dz, d, tol = 1e-4
    type(mesh) :: bcs
!-----------------------------------------------------------------------
    call read_pnts(bc_fname, bcs)
    do i = 1, bcs%totnds
        do j = 1, mshtotpnts
            dx = bcs%nds(i,1)-mshpnts(j,1)
            dy = bcs%nds(i,2)-mshpnts(j,2)
            dz = bcs%nds(i,3)-mshpnts(j,3)
            d = sqrt(dx*dx+dy*dy+dz*dz)
            if (d <= tol) then
                bc_out%tot = bc_out%tot + 1
                bc_out%pnts = [bc_out%pnts, j]
                bc_out%vals = [bc_out%vals, value]
            end if
        end do
    end do
end subroutine

subroutine set_var(u, bcs)
!-----------------------------------------------------------------------
! Sets the value of an array to the specified Direchlet boundary 
! condition values
!-----------------------------------------------------------------------
    use bc_struct
!-----------------------------------------------------------------------
    type(bc), intent(in) :: bcs
!-----------------------------------------------------------------------
    real(8), dimension(:), intent(INOUT) :: u
!-----------------------------------------------------------------------
    if (bcs%tot > 0) u(bcs%pnts(:)) = bcs%vals(:)
end subroutine set_var

end module bc_ops

