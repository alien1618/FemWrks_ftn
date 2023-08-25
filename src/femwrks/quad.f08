module quad_lib
implicit none
contains

subroutine gen_quadpnts(dim, elemtype, level, eps, eta, zta, w, elemquadpnts)
!-----------------------------------------------------------------------
! Generate quadrature points depending on element dimension, type and 
! the level or order of quadrature points to be generated per element
!-----------------------------------------------------------------------
    integer, intent(IN) :: dim, elemtype, level
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: eps, eta, zta, w
    integer, intent(OUT) :: elemquadpnts
!-----------------------------------------------------------------------
    select case (dim)
        case (3)
            select case (elemtype)
                case (1)
                    select case (level)
                        case (1)
                            call tetra(1, eps, eta, zta, w, elemquadpnts)
                        case (2)
                            call tetra(4, eps, eta, zta, w, elemquadpnts)
                    end select
                case (2)
                    select case (level)
                        case (1)
                            call hexa(2, eps, eta, zta, w, elemquadpnts)
                        case (2)
                            call hexa(3, eps, eta, zta, w, elemquadpnts)
                    end select
            end select
        case (2)
            select case (elemtype)
                case (1)
                    select case (level)
                        case (1)
                            call tri(1, eps, eta, zta, w, elemquadpnts)
                        case (2)
                            call tri(3, eps, eta, zta, w, elemquadpnts)
                    end select
                case (2)
                    select case (level)
                        case (1)
                            call quad(2, eps, eta, zta, w, elemquadpnts)
                        case (2)
                            call quad(3, eps, eta, zta, w, elemquadpnts)
                    end select
            end select
    end select

end subroutine gen_quadpnts

subroutine hexa(num_pts, eps, eta, zta, w, integ_pts)
!-----------------------------------------------------------------------
! Generate quadrature points for a hexahedral element
!-----------------------------------------------------------------------
    integer, intent(IN) :: num_pts
    real(8), dimension(:), allocatable, intent(OUT) :: eps, eta, zta, w
    integer, intent(OUT) :: integ_pts
!-----------------------------------------------------------------------
    real(8) :: a1, a2, a3, w1, w2, w3
    real(8), dimension(:), allocatable :: eps_pnt, w_pnt
    integer :: i, j, k, m
!-----------------------------------------------------------------------
    allocate(eps_pnt(num_pts))
    allocate(w_pnt(num_pts))
    select case (num_pts)
        case (1)
            a1 = 0.000000000000000
            w1 = 2.000000000000000
            eps_pnt(1) = a1
            w_pnt(1) = w1
        case (2)
            a1 = 0.577350269189626
            w1 = 1.000000000000000
            eps_pnt(1) = a1
            eps_pnt(2) = -a1
            w_pnt(1) = w1
            w_pnt(2) = w1
        case (3)
            a1 = 0.774596669241483
            a2 = 0.000000000000000
            w1 = 0.555555555555556
            w2 = 0.888888888888889
            eps_pnt(1) = a1
            eps_pnt(2) = a2
            eps_pnt(3) = -a1
            w_pnt(1) = w1
            w_pnt(2) = w2
            w_pnt(3) = w1
        case (4)
            a1 = 0.861136311594053
            a2 = 0.339981043584856
            w1 = 0.347854845137454
            w2 = 0.652145154862546
            eps_pnt(1) = a1
            eps_pnt(2) = a2
            eps_pnt(3) = -a2
            eps_pnt(4) = -a1
            w_pnt(1) = w1
            w_pnt(2) = w2
            w_pnt(3) = w2
            w_pnt(4) = w1
        case (5)
            a1 = 0.906179845938664
            a2 = 0.538469310105683
            a3 = 0.000000000000000
            w1 = 0.236926885056189
            w2 = 0.478628670499366
            w3 = 0.568888888888889
            eps_pnt(1) = a1
            eps_pnt(2) = a2
            eps_pnt(3) = a3
            eps_pnt(4) = -a2
            eps_pnt(5) = -a1
            w_pnt(1) = w1
            w_pnt(2) = w2
            w_pnt(3) = w3
            w_pnt(4) = w2
            w_pnt(5) = w1
        case (6)
            a1 = 0.932469514203152
            a2 = 0.661209386466265
            a3 = 0.238619186083197
            w1 = 0.171324492379170
            w2 = 0.360761573048139
            w3 = 0.467913934572691
            eps_pnt(1) = a1
            eps_pnt(2) = a2
            eps_pnt(3) = a3
            eps_pnt(4) = -a3
            eps_pnt(5) = -a2
            eps_pnt(6) = -a1
            w_pnt(1) = w1;
            w_pnt(2) = w2;
            w_pnt(3) = w3;
            w_pnt(4) = w3;
            w_pnt(5) = w2;
            w_pnt(6) = w1;
    end select

    integ_pts = num_pts*num_pts*num_pts;
    allocate(eps(integ_pts))
    allocate(eta(integ_pts))
    allocate(zta(integ_pts))
    allocate(w(integ_pts))

    m = 0
    do k = 1,num_pts
        do j = 1, num_pts
            do i = 1, num_pts
                m = m + 1
                eps(m) = eps_pnt(i)
                eta(m) = eps_pnt(j)
                zta(m) = eps_pnt(k)
                w(m) = w_pnt(i)*w_pnt(j)*w_pnt(k)
            end do
       end do
    end do
end subroutine hexa

subroutine quad(num_pts, eps, eta, zta, w, integ_pts)
!-----------------------------------------------------------------------
! Generate quadrature points for a quadrilateral element
!-----------------------------------------------------------------------
    integer, intent(IN) :: num_pts
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: eps, eta, zta, w
    integer, intent(OUT) :: integ_pts
!-----------------------------------------------------------------------
    real(8) :: a1, a2, a3, w1, w2, w3
    real(8), dimension(:), allocatable :: eps_pnt, w_pnt
    integer :: i, j, m
!-----------------------------------------------------------------------
    
    allocate(eps_pnt(num_pts))
    allocate(w_pnt(num_pts))
    select case (num_pts)
        case (1)
            a1 = 0.000000000000000
            w1 = 2.000000000000000
            eps_pnt(1) = a1
            w_pnt(1) = w1
        case (2)
            a1 = 0.577350269189626
            w1 = 1.000000000000000
            eps_pnt(1) = a1
            eps_pnt(2) = -a1
            w_pnt(1) = w1
            w_pnt(2) = w1
        case (3)
            a1 = 0.774596669241483
            a2 = 0.000000000000000
            w1 = 0.555555555555556
            w2 = 0.888888888888889
            eps_pnt(1) = a1
            eps_pnt(2) = a2
            eps_pnt(3) = -a1
            w_pnt(1) = w1
            w_pnt(2) = w2
            w_pnt(3) = w1
        case (4)
            a1 = 0.861136311594053
            a2 = 0.339981043584856
            w1 = 0.347854845137454
            w2 = 0.652145154862546
            eps_pnt(1) = a1
            eps_pnt(2) = a2
            eps_pnt(3) = -a2
            eps_pnt(4) = -a1
            w_pnt(1) = w1
            w_pnt(2) = w2
            w_pnt(3) = w2
            w_pnt(4) = w1
        case (5)
            a1 = 0.906179845938664
            a2 = 0.538469310105683
            a3 = 0.000000000000000
            w1 = 0.236926885056189
            w2 = 0.478628670499366
            w3 = 0.568888888888889
            eps_pnt(1) = a1
            eps_pnt(2) = a2
            eps_pnt(3) = a3
            eps_pnt(4) = -a2
            eps_pnt(5) = -a1
            w_pnt(1) = w1
            w_pnt(2) = w2
            w_pnt(3) = w3
            w_pnt(4) = w2
            w_pnt(5) = w1
        case (6)
            a1 = 0.932469514203152
            a2 = 0.661209386466265
            a3 = 0.238619186083197
            w1 = 0.171324492379170
            w2 = 0.360761573048139
            w3 = 0.467913934572691
            eps_pnt(1) = a1
            eps_pnt(2) = a2
            eps_pnt(3) = a3
            eps_pnt(4) = -a3
            eps_pnt(5) = -a2
            eps_pnt(6) = -a1
            w_pnt(1) = w1
            w_pnt(2) = w2
            w_pnt(3) = w3
            w_pnt(4) = w3
            w_pnt(5) = w2
            w_pnt(6) = w1
    end select

    integ_pts = num_pts*num_pts
    allocate(eps(integ_pts))
    allocate(eta(integ_pts))
    allocate(zta(integ_pts))
    allocate(w(integ_pts))
    
    m = 0
    do j = 1, num_pts
        do i = 1, num_pts
            m = m + 1
            eps(m) = eps_pnt(i)
            eta(m) = eps_pnt(j)
            zta(m) = 1
            w(m) = w_pnt(i)*w_pnt(j)
        end do
    end do

end subroutine quad

subroutine tetra(num_pts, eps, eta, zta, w, integ_pts)
!-----------------------------------------------------------------------
! Generate quadrature points for a tetrahedral element
!-----------------------------------------------------------------------
    integer, intent(IN) :: num_pts
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: eps, eta, zta, w
    integer, intent(OUT) :: integ_pts
!-----------------------------------------------------------------------
    real(8) :: a1, a2, w1
!-----------------------------------------------------------------------
                            
    integ_pts = num_pts
    allocate(eps(num_pts))
    allocate(eta(num_pts))
    allocate(zta(num_pts))
    allocate(w(num_pts))
    select case(num_pts)
        case (1)
            eps(1) = 0.250
            eta(1) =  0.250
            zta(1) =  0.250
            w(1) = 1
        case (4)
            a1 = 0.58541020
            a2 = 0.13819660
            w1 = 0.250

            eps(1) = a1
            eta(1) = a2
            zta(1) = a2
            w(1) = w1

            eps(2) = a2
            eta(2) = a1
            zta(2) = a2
            w(2) = w1

            eps(3) = a2
            eta(3) = a2
            zta(3) = a1
            w(3) = w1

            eps(4) = a2
            eta(4) = a2
            zta(4) = a2
            w(4) = w1
    end select
end subroutine tetra

subroutine tri(num_pts, eps, eta, zta, w, integ_pts)
!-----------------------------------------------------------------------
! Generate quadrature points for a triangular element
!-----------------------------------------------------------------------
    integer, intent(IN) :: num_pts
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: eps, eta, zta, w
    integer, intent(OUT) :: integ_pts
!-----------------------------------------------------------------------
    real(8) :: a1, a2, a3, a4, a5, w1, w2, w3
!-----------------------------------------------------------------------

    integ_pts = num_pts
    allocate(eps(integ_pts))
    allocate(eta(integ_pts))
    allocate(zta(integ_pts))
    allocate(w(integ_pts))
    select case(num_pts)
        case (1)
            eps(1) = 0.333333333333
            eta(1) =  0.333333333333
            zta(1) =  1
            w(1) = 1
        case (3)
            a1 = 0.1666666666667
            a2 = 0.6666666666667
            w1 = 0.1666666666667
            eps(1) = a1
            eps(2) = a2
            eps(3) = a1
            eta(1) = a1
            eta(2) = a1
            eta(3) = a2
            zta(1) = 1
            zta(2) = 1
            zta(3) = 1
            w(1) = w1
            w(2) = w1
            w(3) = w1
        case (7)
            a1 = 0.1012865073235
            a2 = 0.7974269853531
            a3 = 0.4701420641051
            a4 = 0.0597158717898
            a5 = 0.3333333333333
            w1 = 0.0629695902724
            w2 = 0.0661970763942
            w3 = 0.1125000000000

            eps(1) = a1
            eps(2) = a2
            eps(3) = a1
            eps(4) = a3
            eps(5) = a3
            eps(6) = a4
            eps(7) = a5
            eta(1) = a1
            eta(2) = a1
            eta(3) = a2
            eta(4) = a4
            eta(5) = a3
            eta(6) = a3
            eta(7) = a5
            zta(1) = 1
            zta(2) = 1
            zta(3) = 1
            zta(4) = 1
            zta(5) = 1
            zta(6) = 1
            zta(7) = 1
            w(1) = w1
            w(2) = w1
            w(3) = w1
            w(4) = w2
            w(5) = w2
            w(6) = w2
            w(7) = w3
    end select
end subroutine tri
end module quad_lib

