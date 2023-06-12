module krnl_struct
implicit none
private
public kernel
!-----------------------------------------------------------------------
! defining kernel data structure
!-----------------------------------------------------------------------
type kernel
    integer :: totnbrs, elem_num
    real(8) :: quad_wt, pntx, pnty, pntz
    integer, dimension(:), allocatable :: nbrs
    real(8), dimension(:), allocatable :: n, nx, ny, nz
    real(8), dimension(:,:), allocatable :: nx_nx, ny_ny, nz_nz, n_nx, n_ny, n_nz, n_n
end type
end module krnl_struct

module krnl_ops
implicit none
contains

subroutine get_intrps(dim0, nds, elems, totelems, elemnds, eps, eta, zta, w, elemquadpnts, krnls, totquadpnts)
!-----------------------------------------------------------------------
! Calculate and store Lagrange interpolants for each element
!-----------------------------------------------------------------------
    use krnl_struct
    use eq_slvrs
!-----------------------------------------------------------------------
    integer, intent(IN) :: dim0, totelems, elemnds, elemquadpnts
    real(8), dimension(:,:), intent(IN) :: nds
    integer, dimension(:,:), intent(IN) :: elems
    real(8), dimension(:), intent(IN) :: eps, eta, zta, w
!-----------------------------------------------------------------------
    type(kernel), dimension(:), allocatable, intent(OUT)    :: krnls
    integer, intent(OUT) :: totquadpnts
!-----------------------------------------------------------------------
    integer :: m = 0, i, j, k
    real(8) :: det_jac
!-----------------------------------------------------------------------

    totquadpnts = elemquadpnts * totelems
    write(*,'(a,i0)') "Total quad pnts generated = ", totquadpnts
    allocate(krnls(totquadpnts))

    do k = 1, totelems
        do i = 1, elemquadpnts
            m = m + 1
            call get_lagrange(dim0, nds, elems(k,:), elemnds, eps(i), eta(i), zta(i), krnls(m)%n, &
                krnls(m)%nx, krnls(m)%ny, krnls(m)%nz, det_jac)
            
            krnls(m)%quad_wt = w(i)*abs(det_jac)
            krnls(m)%elem_num = k

            krnls(m)%pntx = 0
            krnls(m)%pnty = 0
            krnls(m)%pntz = 0
            do j = 1, elemnds
                krnls(m)%pntx = krnls(m)%pntx + krnls(m)%n(j)*nds(elems(k,j),1)
                krnls(m)%pnty = krnls(m)%pnty + krnls(m)%n(j)*nds(elems(k,j),2)
                krnls(m)%pntz = krnls(m)%pntz + krnls(m)%n(j)*nds(elems(k,j),3)
            end do
            call mul_trnsp(krnls(m)%nx, krnls(m)%nx, krnls(m)%nx_nx, elemnds)
            call mul_trnsp(krnls(m)%ny, krnls(m)%ny, krnls(m)%ny_ny, elemnds)
            call mul_trnsp(krnls(m)%nz, krnls(m)%nz, krnls(m)%nz_nz, elemnds)
            call mul_trnsp(krnls(m)%n, krnls(m)%nx, krnls(m)%n_nx, elemnds)
            call mul_trnsp(krnls(m)%n, krnls(m)%ny, krnls(m)%n_ny, elemnds)
            call mul_trnsp(krnls(m)%n, krnls(m)%nz, krnls(m)%n_nz, elemnds)
            call mul_trnsp(krnls(m)%n, krnls(m)%n, krnls(m)%n_n, elemnds)
        end do
    end do

end subroutine get_intrps

subroutine get_lagrange(dim, nds, nbrnds, totnbrnds, eps, eta, zta, n, nx, ny, nz, det_jac)
!-----------------------------------------------------------------------
! Calculate Lagrange inteprolants for a given set of nodes
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: nds
    integer, dimension(:), intent(IN) :: nbrnds
    integer, intent(IN) :: dim, totnbrnds
    real(8), intent(IN) :: eps, eta, zta
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: n, nx, ny, nz
    real(8), intent(OUT) :: det_jac
!-----------------------------------------------------------------------

    select case (dim)
        case (2)
            select case (totnbrnds)
                case(4)
                    call Q4(nds, nbrnds, totnbrnds, eps, eta, n, nx, ny, det_jac)
                case(8)
                    call Q8(nds, nbrnds, totnbrnds, eps, eta, n, nx, ny, det_jac)
                case (3)
                    call T3(nds, nbrnds, totnbrnds, eps, eta, n, nx, ny, det_jac)
                case (6)
                    call T6(nds, nbrnds, totnbrnds, eps, eta, n, nx, ny, det_jac)
            end select
            allocate(nz(totnbrnds))
            nz(:) = 0
        case (3)
            select case (totnbrnds)
                case (4)
                    call T4(nds, nbrnds, totnbrnds, eps, eta, zta, n, nx, ny, nz, det_jac)
                case (10)
                    call T10(nds, nbrnds, totnbrnds, eps, eta, zta, n, nx, ny, nz, det_jac)
                case (8)
                    call H8(nds, nbrnds, totnbrnds, eps, eta, zta, n, nx, ny, nz, det_jac)
                case (20)
                    call H20(nds, nbrnds, totnbrnds, eps, eta, zta, n, nx, ny, nz, det_jac)
            end select
    end select
end subroutine get_lagrange

subroutine H8(nds, nbrnds, totnbrnds, eps, eta, zta, n, nx, ny, nz, det_jac)
!-----------------------------------------------------------------------
! Calculate Lagrange inteprolants for a 3D first order Hexahedral
! element
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: nds
    integer, dimension(:), intent(IN) :: nbrnds
    integer, intent(IN) :: totnbrnds
    real(8), intent(IN) :: eps, eta, zta
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: n, nx, ny, nz
    real(8), intent(OUT) :: det_jac
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable :: dNdeps, dNdeta, dNdzta
!-----------------------------------------------------------------------
    
    allocate (n(totnbrnds))
    allocate(dNdeps(totnbrnds))
    allocate(dNdeta(totnbrnds))
    allocate(dNdzta(totnbrnds))

    n(1) = 0.125*(1-eps)*(1-eta)*(1+zta)
    n(2) = 0.125*(1+eps)*(1-eta)*(1+zta)
    n(3) = 0.125*(1+eps)*(1+eta)*(1+zta)
    n(4) = 0.125*(1-eps)*(1+eta)*(1+zta)
    n(5) = 0.125*(1-eps)*(1-eta)*(1-zta)
    n(6) = 0.125*(1+eps)*(1-eta)*(1-zta)
    n(7) = 0.125*(1+eps)*(1+eta)*(1-zta)
    n(8) = 0.125*(1-eps)*(1+eta)*(1-zta)

    dNdeps(1) = -0.125*(1-eta)*(1+zta)
    dNdeta(1) = -0.125*(1-eps)*(1+zta)
    dNdzta(1) = 0.125*(1-eps)*(1-eta)

    dNdeps(2) = 0.125*(1-eta)*(1+zta)
    dNdeta(2) = -0.125*(1+eps)*(1+zta)
    dNdzta(2) = 0.125*(1+eps)*(1-eta)

    dNdeps(3) = 0.125*(1+eta)*(1+zta)
    dNdeta(3) = 0.125*(1+eps)*(1+zta)
    dNdzta(3) = 0.125*(1+eps)*(1+eta)

    dNdeps(4) = -0.125*(1+eta)*(1+zta)
    dNdeta(4) = 0.125*(1-eps)*(1+zta)
    dNdzta(4) = 0.125*(1-eps)*(1+eta)

    dNdeps(5) = -0.125*(1-eta)*(1-zta)
    dNdeta(5) = -0.125*(1-eps)*(1-zta)
    dNdzta(5) = -0.125*(1-eps)*(1-eta)

    dNdeps(6) = 0.125*(1-eta)*(1-zta)
    dNdeta(6) = -0.125*(1+eps)*(1-zta)
    dNdzta(6) = -0.125*(1+eps)*(1-eta)

    dNdeps(7) = 0.125*(1+eta)*(1-zta)
    dNdeta(7) = 0.125*(1+eps)*(1-zta)
    dNdzta(7) = -0.125*(1+eps)*(1+eta)

    dNdeps(8) = -0.125*(1+eta)*(1-zta)
    dNdeta(8) = 0.125*(1-eps)*(1-zta)
    dNdzta(8) = -0.125*(1-eps)*(1+eta)

    call jacobian_volume(nds, nbrnds, totnbrnds, dNdeps, dNdeta, dNdzta, n, nx, ny, nz, det_jac)

end subroutine H8

subroutine H20(nds, nbrnds, totnbrnds, eps, eta, zta, n, nx, ny, nz, det_jac)
!-----------------------------------------------------------------------
! Calculate Lagrange inteprolants for a 3D second order Hexahedral
! element
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: nds
    integer, dimension(:), intent(IN) :: nbrnds
    integer, intent(IN) :: totnbrnds
    real(8), intent(IN) :: eps, eta, zta
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: n, nx, ny, nz
    real(8), intent(OUT) :: det_jac
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable :: dNdeps, dNdeta, dNdzta
!----------------------------------------------------------------------- 
    
    allocate (n(totnbrnds))
    allocate(dNdeps(totnbrnds))
    allocate(dNdeta(totnbrnds))
    allocate(dNdzta(totnbrnds))

    n(1) = 0.125*(1-eps)*(1-eta)*(1-zta)*(-eps-eta-zta-2)
    n(2) = 0.125*(1-eps)*(1-eta)*(1-zta)*(eps-eta-zta-2)
    n(3) = 0.125*(1+eps)*(1+eta)*(1-zta)*(eps+eta-zta-2)
    n(4) = 0.125*(1-eps)*(1+eta)*(1-zta)*(-eps+eta-zta-2)
    n(5) = 0.125*(1-eps)*(1-eta)*(1+zta)*(-eps-eta+zta-2)
    n(6) = 0.125*(1+eps)*(1-eta)*(1+zta)*(eps-eta+zta-2)
    n(7) = 0.125*(1+eps)*(1+eta)*(1+zta)*(eps+eta+zta-2)
    n(8) = 0.125*(1-eps)*(1+eta)*(1+zta)*(-eps+eta+zta-2)
    n(9) = 0.125*2*(1-eps*eps)*(1-eta)*(1-zta)
    n(10) = 0.125*2*(1+eps)*(1-eta*eta)*(1-zta)
    n(11) = 0.125*2*(1-eps*eps)*(1+eta)*(1-zta)
    n(12) = 0.125*2*(1-eps)*(1-eta*eta)*(1-zta)
    n(13) = 0.125*2*(1-eps)*(1-eta)*(1-zta*zta)
    n(14) = 0.125*2*(1+eps)*(1-eta)*(1-zta*zta)
    n(15) = 0.125*2*(1+eps)*(1+eta)*(1-zta*zta)
    n(16) = 0.125*2*(1-eps)*(1+eta)*(1-zta*zta)
    n(17) = 0.125*2*(1-eps*eps)*(1-eta)*(1+zta)
    n(18) = 0.125*2*(1+eps)*(1-eta*eta)*(1+zta)
    n(19) = 0.125*2*(1-eps*eps)*(1+eta)*(1+zta)
    n(20) = 0.125*2*(1-eps)*(1-eta*eta)*(1+zta)

    dNdeps(1) = -0.125*(1-eta)*(1-zta)*(-eps-eta-zta-2)-0.125*(1-eps)*(1-eta)*(1-zta)
    dNdeta(1) = -0.125*(1-eps)*(1-zta)*(-eps-eta-zta-2)-0.125*(1-eps)*(1-eta)*(1-zta)
    dNdzta(1) = -0.125*(1-eps)*(1-eta)*(-eps-eta-zta-2)-0.125*(1-eps)*(1-eta)*(1-zta)

    dNdeps(2) = -0.125*(1-eta)*(1-zta)*(eps-eta-zta-2)+0.125*(1-eps)*(1-eta)*(1-zta)
    dNdeta(2) = -0.125*(1-eps)*(1-zta)*(eps-eta-zta-2)-0.125*(1-eps)*(1-eta)*(1-zta)
    dNdzta(2) = -0.125*(1-eps)*(1-eta)*(eps-eta-zta-2)-0.125*(1-eps)*(1-eta)*(1-zta)

    dNdeps(3) = 0.125*(1+eta)*(1-zta)*(eps+eta-zta-2)+0.125*(1+eps)*(1+eta)*(1-zta)
    dNdeta(3) = 0.125*(1+eps)*(1-zta)*(eps+eta-zta-2)+0.125*(1+eps)*(1+eta)*(1-zta)
    dNdzta(3) = -0.125*(1+eps)*(1+eta)*(eps+eta-zta-2)-0.125*(1+eps)*(1+eta)*(1-zta)

    dNdeps(4) = -0.125*(1+eta)*(1-zta)*(-eps+eta-zta-2)-0.125*(1-eps)*(1+eta)*(1-zta)
    dNdeta(4) =  0.125*(1-eps)*(1-zta)*(-eps+eta-zta-2)+ 0.125*(1-eps)*(1+eta)*(1-zta)
    dNdzta(4) =  -0.125*(1-eps)*(1+eta)*(-eps+eta-zta-2)-0.125*(1-eps)*(1+eta)*(1-zta)

    dNdeps(5) = -0.125*(1-eta)*(1+zta)*(-eps-eta+zta-2)-0.125*(1-eps)*(1-eta)*(1+zta)
    dNdeta(5) = -0.125*(1-eps)*(1+zta)*(-eps-eta+zta-2)-0.125*(1-eps)*(1-eta)*(1+zta)
    dNdzta(5) = 0.125*(1-eps)*(1-eta)*(-eps-eta+zta-2)+0.125*(1-eps)*(1-eta)*(1+zta)

    dNdeps(6) = 0.125*(1-eta)*(1+zta)*(eps-eta+zta-2)+0.125*(1+eps)*(1-eta)*(1+zta)
    dNdeta(6) = -0.125*(1+eps)*(1+zta)*(eps-eta+zta-2)-0.125*(1+eps)*(1-eta)*(1+zta)
    dNdzta(6) = 0.125*(1+eps)*(1-eta)*(eps-eta+zta-2)+0.125*(1+eps)*(1-eta)*(1+zta)

    dNdeps(7) = 0.125*(1+eta)*(1+zta)*(eps+eta+zta-2)+0.125*(1+eps)*(1+eta)*(1+zta)
    dNdeta(7) = 0.125*(1+eps)*(1+zta)*(eps+eta+zta-2)+0.125*(1+eps)*(1+eta)*(1+zta)
    dNdzta(7) = 0.125*(1+eps)*(1+eta)*(eps+eta+zta-2)+0.125*(1+eps)*(1+eta)*(1+zta)

    dNdeps(8) = -0.125*(1+eta)*(1+zta)*(-eps+eta+zta-2)-0.125*(1-eps)*(1+eta)*(1+zta)
    dNdeta(8) = 0.125*(1-eps)*(1+zta)*(-eps+eta+zta-2)+0.125*(1-eps)*(1+eta)*(1+zta)
    dNdzta(8) = 0.125*(1-eps)*(1+eta)*(-eps+eta+zta-2)+0.125*(1-eps)*(1+eta)*(1+zta)

    dNdeps(9) = 0.125*2*(-2*eps)*(1-eta)*(1-zta)
    dNdeta(9) = -0.125*2*(1-eps*eps)*(1-zta)
    dNdzta(9) = -0.125*2*(1-eps*eps)*(1-eta)

    dNdeps(10) = 0.125*2*(1-eta*eta)*(1-zta)
    dNdeta(10) = 0.125*2*(1+eps)*(-2*eta)*(1-zta)
    dNdzta(10) = -0.125*2*(1+eps)*(1-eta*eta)

    dNdeps(11) = 0.125*2*(-2*eps)*(1+eta)*(1-zta)
    dNdeta(11) = 0.125*2*(1-eps*eps)*(1-zta)
    dNdzta(11) = -0.125*2*(1-eps*eps)*(1+eta)

    dNdeps(12) =  -0.125*2*(1-eta*eta)*(1-zta)
    dNdeta(12) =  0.125*2*(1-eps)*(-2*eta)*(1-zta)
    dNdzta(12) =  -0.125*2*(1-eps)*(1-eta*eta)

    dNdeps(13) = -0.125*2*(1-eta)*(1-zta*zta)
    dNdeta(13) = -0.125*2*(1-eps)*(1-zta*zta)
    dNdzta(13) = 0.125*2*(1-eps)*(1-eta)*(-2*zta)

    dNdeps(14) = 0.125*2*(1-eta)*(1-zta*zta)
    dNdeta(14) =  -0.125*2*(1+eps)*(1-zta*zta)
    dNdzta(14) =  0.125*2*(1+eps)*(1-eta)*(-2*zta)

    dNdeps(15) = 0.125*2*(1+eta)*(1-zta*zta)
    dNdeta(15) = 0.125*2*(1+eps)*(1-zta*zta)
    dNdzta(15) = 0.125*2*(1+eps)*(1+eta)*(-2*zta)

    dNdeps(16) = -0.125*2*(1+eta)*(1-zta*zta)
    dNdeta(16) = 0.125*2*(1-eps)*(1-zta*zta)
    dNdzta(16) = 0.125*2*(1-eps)*(1+eta)*(-2*zta)

    dNdeps(17) = 0.125*2*(-2*eps)*(1-eta)*(1+zta)
    dNdeta(17) = -0.125*2*(1-eps*eps)*(1+zta)
    dNdzta(17) = 0.125*2*(1-eps*eps)*(1-eta)

    dNdeps(18) = 0.125*2*(1-eta*eta)*(1+zta)
    dNdeta(18) = 0.125*2*(1+eps)*(-2*eta)*(1+zta)
    dNdzta(18) = 0.125*2*(1+eps)*(1-eta*eta)

    dNdeps(19) = 0.125*2*(-2*eps)*(1+eta)*(1+zta)
    dNdeta(19) = 0.125*2*(1-eps*eps)*(1+zta)
    dNdzta(19) = 0.125*2*(1-eps*eps)*(1+eta)

    dNdeps(20) = -0.125*2*(1-eta*eta)*(1+zta)
    dNdeta(20) = 0.125*2*(1-eps)*(-2*eta)*(1+zta)
    dNdzta(20) = 0.125*2*(1-eps)*(1-eta*eta)

    call jacobian_volume(nds, nbrnds, totnbrnds, dNdeps, dNdeta, dNdzta, n, nx, ny, nz, det_jac)
end subroutine H20

subroutine T4(nds, nbrnds, totnbrnds, eps, eta, zta, n, nx, ny, nz, det_jac)
!-----------------------------------------------------------------------
! Calculate Lagrange inteprolants for a 3D first order Tetrahedral 
! element
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: nds
    integer, dimension(:), intent(IN) :: nbrnds
    integer, intent(IN) :: totnbrnds
    real(8), intent(IN) :: eps, eta, zta
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: n, nx, ny, nz
    real(8), intent(OUT) :: det_jac
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable :: dNdeps, dNdeta, dNdzta
!----------------------------------------------------------------------- 
    
    allocate (n(totnbrnds))
    allocate(dNdeps(totnbrnds))
    allocate(dNdeta(totnbrnds))
    allocate(dNdzta(totnbrnds))

    n(1) = 1-eps-eta-zta
    n(2) = eps
    n(3) = eta
    n(4) = zta

    dNdeps(1) = -1.0
    dNdeta(1) = -1.0
    dNdzta(1) = -1.0

    dNdeps(2) = 1.0
    dNdeta(2) = 0.0
    dNdzta(2) = 0.0

    dNdeps(3) = 0.0
    dNdeta(3) = 1.0
    dNdzta(3) = 0.0

    dNdeps(4) = 0.0
    dNdeta(4) = 0.0
    dNdzta(4) = 1.0

     call jacobian_volume(nds, nbrnds, totnbrnds, dNdeps, dNdeta, dNdzta, n, nx, ny, nz, det_jac)

end subroutine T4

subroutine T10(nds, nbrnds, totnbrnds, eps, eta, zta, n, nx, ny, nz, det_jac)
!-----------------------------------------------------------------------
! Calculate Lagrange inteprolants for a 3D second order Tetrahedral 
! element
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: nds
    integer, dimension(:), intent(IN) :: nbrnds
    integer, intent(IN) :: totnbrnds
    real(8), intent(IN) :: eps, eta, zta
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: n, nx, ny, nz
    real(8), intent(OUT) :: det_jac
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable :: dNdeps, dNdeta, dNdzta
!----------------------------------------------------------------------- 
    
    allocate (n(totnbrnds))
    allocate(dNdeps(totnbrnds))
    allocate(dNdeta(totnbrnds))
    allocate(dNdzta(totnbrnds))

    n(1) = (1 - eps - eta - zta)*(2*(1 - eps - eta - zta)-1)
    n(2) = eps*(2*eps-1)
    n(3) = eta*(2*eta-1)
    n(4) = zta*(2*zta-1)
    n(5) = 4*(1 - eps - eta - zta) * eps
    n(6) = 4*eps*eta
    n(7) = 4*(1 - eps - eta - zta) * eta
    n(8) = 4*(1 - eps - eta - zta) * zta
    n(9) = 4*eps*zta
    n(10) = 4*eta*zta

    dNdeps(1) = -1*(2*(1 - eps - eta - zta)-1)+(1 - eps - eta - zta)*(-2)
    dNdeta(1) = -1*(2*(1 - eps - eta - zta)-1)+(1 - eps - eta - zta)*(-2)
    dNdzta(1) = -1*(2*(1 - eps - eta - zta)-1)+(1 - eps - eta - zta)*(-2)

    dNdeps(2) = 4*eps-1
    dNdeta(2) = 0
    dNdzta(2) = 0

    dNdeps(3) = 0
    dNdeta(3) = 4*eta-1
    dNdzta(3) = 0

    dNdeps(4) = 0
    dNdeta(4) = 0
    dNdzta(4) = 4*zta-1

    dNdeps(5) = -4*eps+4*(1 - eps - eta - zta)
    dNdeta(5) = -4*eps
    dNdzta(5) = -4*eps

    dNdeps(6) = 4*eta
    dNdeta(6) = 4*eps
    dNdzta(6) = 0

    dNdeps(7) = -4*eta
    dNdeta(7) = -4*eta+4*(1 - eps - eta - zta)
    dNdzta(7) = -4*eta

    dNdeps(8) = -4*zta
    dNdeta(8) = -4*zta
    dNdzta(8) = -4*zta+4*(1 - eps - eta - zta)

    dNdeps(9) = 4*zta
    dNdeta(9) = 0
    dNdzta(9) = 4*eps

    dNdeps(10) = 0
    dNdeta(10) = 4*zta
    dNdzta(10) = 4*eta

    call jacobian_volume(nds, nbrnds, totnbrnds, dNdeps, dNdeta, dNdzta, n, nx, ny, nz, det_jac)

end subroutine T10

subroutine Q4(nds, nbrnds, totnbrnds, eps, eta, n, nx, ny, det_jac)
!-----------------------------------------------------------------------
! Calculate Lagrange inteprolants for a 2D first-order Quadrilateral
! element
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: nds
    integer, dimension(:), intent(IN) :: nbrnds
    integer, intent(IN) :: totnbrnds
    real(8), intent(IN) :: eps, eta
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: n, nx, ny
    real(8), intent(OUT) :: det_jac
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable :: dNdeps, dNdeta
!-----------------------------------------------------------------------

    allocate(n(totnbrnds))
    allocate(dNdeps(totnbrnds))
    allocate(dNdeta(totnbrnds))

    n(1) = 0.25*(1-eps)*(1-eta)
    n(2) = 0.25*(1+eps)*(1-eta)
    n(3) = 0.25*(1+eps)*(1+eta)
    n(4) = 0.25*(1-eps)*(1+eta)

    dNdeps(1) = -0.25*(1-eta)
    dNdeta(1) = -0.25*(1-eps)
    
    dNdeps(2) = 0.25*(1-eta)
    dNdeta(2) = -0.25*(1+eps)
    
    dNdeps(3) = 0.25*(1+eta)
    dNdeta(3) = 0.25*(1+eps)
    
    dNdeps(4) = -0.25*(1+eta)
    dNdeta(4) = 0.25*(1-eps)

    call jacobian_surface(nds, nbrnds, totnbrnds, dNdeps, dNdeta, n, nx, ny, det_jac)

end subroutine Q4

subroutine Q8(nds, nbrnds, totnbrnds, eps, eta, n, nx, ny, det_jac)
!-----------------------------------------------------------------------
! Calculate Lagrange inteprolants for a 2D second-order Quadrilateral
! element
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: nds
    integer, dimension(:), intent(IN) :: nbrnds
    integer, intent(IN) :: totnbrnds
    real(8), intent(IN) :: eps, eta
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: n, nx, ny
    real(8), intent(OUT) :: det_jac
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable :: dNdeps, dNdeta
!-----------------------------------------------------------------------

    allocate(n(totnbrnds))
    allocate(dNdeps(totnbrnds))
    allocate(dNdeta(totnbrnds))

    n(1) = 0.25*(eps-1)*(1-eta)*(eps+eta+1)
    n(2) = 0.25*(eps+1)*(1-eta)*(eps-eta-1)
    n(3) = 0.25*(1+eps)*(1+eta)*(eps+eta-1)
    n(4) = 0.25*(1-eps)*(1+eta)*(-eps+eta-1)
    n(5) = 0.5*(1-eps*eps)*(1-eta)
    n(6) = 0.5*(1+eps)*(1-eta*eta)
    n(7) = 0.5*(1-eps*eps)*(1+eta)
    n(8) = 0.5*(1-eps)*(1-eta*eta)

    dNdeps(1) = 0.25*(eps-1)*(1-eta) + 0.25*(1-eta)*(eps+eta+1)
    dNdeta(1) = -0.25*(eps-1)*(eps+eta+1) + 0.25*(eps-1)*(1-eta)

    dNdeps(2) = 0.25*(1-eta)*(eps-eta-1) + 0.25*(eps+1)*(1-eta)
    dNdeta(2) = -0.25*(eps+1)*(eps-eta-1) - 0.25*(eps+1)*(1-eta)

    dNdeps(3) = 0.25*(1+eta)*(eps+eta-1) + 0.25*(1+eps)*(1+eta)
    dNdeta(3) = 0.25*(1+eps)*(eps+eta-1) + 0.25*(1+eps)*(1+eta)

    dNdeps(4) = -0.25*(1+eta)*(-eps+eta-1) - 0.25*(1-eps)*(1+eta)
    dNdeta(4) = 0.25*(1-eps)*(-eps+eta-1) + 0.25*(1-eps)*(1+eta)

    dNdeps(5) = -0.5*(1-eta)*(2*eps)
    dNdeta(5) = -0.5*(1-eps*eps)

    dNdeps(6) = 0.5*(1-eta*eta)
    dNdeta(6) = -0.5*(1+eps)*(2*eta)

    dNdeps(7) = -0.5*(1+eta)*(2*eps)
    dNdeta(7) = 0.5*(1-eps*eps)

    dNdeps(8) = -0.5*(1-eta*eta)
    dNdeta(8) = -0.5*(1-eps)*(2*eta)

    call jacobian_surface(nds, nbrnds, totnbrnds, dNdeps, dNdeta, n, nx, ny, det_jac)
end subroutine Q8

subroutine T3(nds, nbrnds, totnbrnds, eps, eta, n, nx, ny, det_jac)
!-----------------------------------------------------------------------
! Calculate Lagrange inteprolants for a 2D first-order Triangular
! element
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: nds
    integer, dimension(:), intent(IN) :: nbrnds
    integer, intent(IN) :: totnbrnds
    real(8), intent(IN) :: eps, eta
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: n, nx, ny
    real(8), intent(OUT) :: det_jac
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable :: dNdeps, dNdeta
!-----------------------------------------------------------------------

    allocate(n(totnbrnds))
    allocate(dNdeps(totnbrnds))
    allocate(dNdeta(totnbrnds))

    n(1) = (1-eps-eta)
    n(2) = eps
    n(3) = eta

    dNdeps(1) = -1
    dNdeta(1) = -1
    
    dNdeps(2) = 1
    dNdeta(2) = 0
    
    dNdeps(3) = 0
    dNdeta(3) = 1

    call jacobian_surface(nds, nbrnds, totnbrnds, dNdeps, dNdeta, n, nx, ny, det_jac)
end subroutine T3

subroutine T6(nds, nbrnds, totnbrnds, eps, eta, n, nx, ny, det_jac)
!-----------------------------------------------------------------------
! Calculate Lagrange inteprolants for a 2D second order Triangular
! element
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: nds
    integer, dimension(:), intent(IN) :: nbrnds
    integer, intent(IN) :: totnbrnds
    real(8), intent(IN) :: eps, eta
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: n, nx, ny
    real(8), intent(OUT) :: det_jac
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable :: dNdeps, dNdeta
!-----------------------------------------------------------------------

    allocate(n(totnbrnds))
    allocate(dNdeps(totnbrnds))
    allocate(dNdeta(totnbrnds))

    n(1) = (2*(1-eps-eta)-1)*(1-eps-eta)
    n(2) = eps*(2*eps - 1)
    n(3) = eta*(2*eta - 1)
    n(4) = 4*eps*(1-eps-eta)
    n(5) = 4*eps*eta
    n(6) = 4*eta*(1-eps-eta)

    dNdeps(1) = (-(2*(1-eps-eta)-1))-2*(1-eps-eta)
    dNdeta(1) = (-(2*(1-eps-eta)-1))-2*(1-eps-eta)
    
    dNdeps(2) = 4*eps - 1
    dNdeta(2) = 0
    
    dNdeps(3) = 0
    dNdeta(3) = 4*eta - 1
    
    dNdeps(4) = -4*eps + 4*(1-eps-eta)
    dNdeta(4) = -4*eps
    
    dNdeps(5) = 4*eta
    dNdeta(5) = 4*eps
    
    dNdeps(6) = -4*eta
    dNdeta(6) = -4*eta + 4*(1-eps-eta)

    call jacobian_surface(nds, nbrnds, totnbrnds, dNdeps, dNdeta, n, nx, ny, det_jac)
end subroutine T6

subroutine jacobian_surface(nds, nbrnds, totnbrnds, dNdeps, dNdeta, n, nx, ny, det_jac)
!-----------------------------------------------------------------------
! Calculate the jacobian matrix for a 2D element
!-----------------------------------------------------------------------
    use eq_slvrs
    real(8), dimension(:,:), intent(IN) :: nds
    integer, dimension(:), intent(IN) :: nbrnds
    integer, intent(IN) :: totnbrnds
    real(8), dimension(:), intent(IN) :: n, dNdeps, dNdeta
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: nx, ny
    real(8), intent(OUT) :: det_jac
!-----------------------------------------------------------------------
    real(8), dimension(2,2) :: jac, inv_jac
    real(8) :: dxdeps, dxdeta, dydeps, dydeta, sumn, pntx, pnty, pntz
!-----------------------------------------------------------------------
    allocate(nx(totnbrnds))
    allocate(ny(totnbrnds))

    pntx = sum(n(:) * nds(nbrnds(:),1))
    pnty = sum(n(:) * nds(nbrnds(:),2))
    pntz = sum(n(:) * nds(nbrnds(:),3))

    dxdeps = sum(dNdeps(:) * nds(nbrnds(:),1))
    dydeps = sum(dNdeps(:) * nds(nbrnds(:),2))

    dxdeta = sum(dNdeta(:) * nds(nbrnds(:),1))
    dydeta = sum(dNdeta(:) * nds(nbrnds(:),2))

    jac(1,1) = dxdeps
    jac(2,1) = dydeps

    jac(1,2) = dxdeta
    jac(2,2) = dydeta

    call det_2x2(jac, det_jac)
    call inv_2x2(jac, inv_jac)

    nx = dNdeps*inv_jac(1,1) + dNdeta*inv_jac(2,1)
    ny = dNdeps*inv_jac(1,2) + dNdeta*inv_jac(2,2)

    sumn = sum(n)
    if (sumn <= 0.9) then
        write(*,'(a, F7.5)') "ERROR in constructed Lagrange interpolants. Partition of Unity is Not satisfied. sum = ", sumn
    end if

end subroutine jacobian_surface

subroutine jacobian_volume(nds, nbrnds, totnbrnds, dNdeps, dNdeta, dNdzta, n, nx, ny, nz, det_jac)
!-----------------------------------------------------------------------
! Calculate the jacobian matrix for a 3D element
!-----------------------------------------------------------------------
    use eq_slvrs
    real(8), dimension(:,:), intent(IN) :: nds
    integer, dimension(:), intent(IN) :: nbrnds
    integer, intent(IN) :: totnbrnds
    real(8), dimension(:), intent(IN) :: dNdeps, dNdeta, dNdzta, N
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: nx, ny, nz
    real(8), intent(OUT) :: det_jac
!-----------------------------------------------------------------------
    real(8), dimension(3,3) :: jac, inv_jac
    real(8) :: dxdeps, dxdeta, dxdzta, dydeps, dydeta, dydzta
    real(8) :: dzdeps, dzdeta, dzdzta, sumn, pntx, pnty, pntz
!-----------------------------------------------------------------------
    allocate(nx(totnbrnds))
    allocate(ny(totnbrnds))
    allocate(nz(totNbrNds))

    pntx = sum(n(:)*nds(nbrnds(:),1))
    pnty = sum(n(:)*nds(nbrnds(:),2))
    pntz = sum(n(:)*nds(nbrnds(:),3))

    dxdeps = sum(dNdeps(:)*nds(nbrnds(:),1))
    dydeps = sum(dNdeps(:)*nds(nbrnds(:),2))
    dzdeps = sum(dNdeps(:)*nds(nbrnds(:),3))

    dxdeta = sum(dNdeta(:)*nds(nbrnds(:),1))
    dydeta = sum(dNdeta(:)*nds(nbrnds(:),2))
    dzdeta = sum(dNdeta(:)*nds(nbrnds(:),3))

    dxdzta = sum(dNdzta(:)*nds(nbrnds(:),1))
    dydzta = sum(dNdzta(:)*nds(nbrnds(:),2))
    dzdzta = sum(dNdzta(:)*nds(nbrnds(:),3))

    jac(1,1) = dxdeps
    jac(1,2) = dydeps
    jac(1,3) = dzdeps

    jac(2,1) = dxdeta
    jac(2,2) = dydeta
    jac(2,3) = dzdeta

    jac(3,1) = dxdzta
    jac(3,2) = dydzta
    jac(3,3) = dzdzta

    call det_3x3(jac, det_jac)
    call inv_3x3(jac, inv_jac)

    nx = dNdeps*inv_jac(1,1) + dNdeta*inv_jac(1,2) + dNdzta*inv_jac(1,3)
    ny = dNdeps*inv_jac(2,1) + dNdeta*inv_jac(2,2) + dNdzta*inv_jac(2,3)
    nz = dNdeps*inv_jac(3,1) + dNdeta*inv_jac(3,2) + dNdzta*inv_jac(3,3)
    
    sumn = sum(n)
    if (sumn <= 0.9) then
        write(*,'(a, F7.5)') "ERROR in constructed Lagrange interpolants. Partition of Unity is Not satisfied. sum = ", sumn
    end if
end subroutine jacobian_volume

end module krnl_ops
