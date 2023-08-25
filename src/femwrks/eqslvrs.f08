module eqslvrs_lib
implicit none
contains

pure subroutine det_2X2(mat, det)
!-----------------------------------------------------------------------
! Calculate the determinant of a 2X2 matrix
!-----------------------------------------------------------------------
    real(8), dimension(2,2), intent(IN) :: mat
!-----------------------------------------------------------------------
    real(8), intent(OUT) :: det
!-----------------------------------------------------------------------
    det = (mat(1,1) * mat(2,2)) - (mat(2,1) * mat(1,2))
end subroutine det_2X2

pure subroutine det_3X3(mat, det)
!-----------------------------------------------------------------------
! Calculate the determinant of a 3X3 matrix
!-----------------------------------------------------------------------
    real(8), dimension(3,3), intent(IN) :: mat
!-----------------------------------------------------------------------
    real(8), intent(OUT) :: det
!-----------------------------------------------------------------------
    det = mat(1,1) * ((mat(2,2)*mat(3,3))-(mat(3,2)*mat(2,3)))
    det = det - (mat(2,1) * ((mat(1,2)*mat(3,3))-(mat(3,2)*mat(1,3))))
    det = det + (mat(3,1) * ((mat(1,2)*mat(2,3))-(mat(2,2)*mat(1,3))))
end subroutine

pure subroutine inv_2x2(A, A_inv)
!-----------------------------------------------------------------------
! Calculate the inverse of a 2X2 matrix
!-----------------------------------------------------------------------
    real(8), dimension(2,2), intent(IN) :: A
!-----------------------------------------------------------------------
    real(8), dimension(2,2), intent(OUT) :: A_inv
!-----------------------------------------------------------------------
    real(8) :: det_A, num
!-----------------------------------------------------------------------

    call det_2X2(A, det_A)
    num = 1/det_A
    A_inv(1,1) = num*A(2,2)
    A_inv(2,1) = -num*A(2,1)
    A_inv(1,2) = -num*A(1,2)
    A_inv(2,2) = num*A(1,1)
end subroutine inv_2x2

pure subroutine inv_3x3(A, A_inv)
!-----------------------------------------------------------------------
! Calculate the determinant of a 3X3 matrix
!-----------------------------------------------------------------------
    real(8), dimension(3,3), intent(IN) :: A
!-----------------------------------------------------------------------
    real(8), dimension(3,3), intent(OUT) :: A_inv
!-----------------------------------------------------------------------
    real(8), dimension(3,3) :: T
    real(8) :: T11, T21, T31, T12, T22, T32, T13, T23, T33
    real(8) :: det_A, num
    integer :: i, j
!-----------------------------------------------------------------------

    call det_3X3(A, det_A)
    num = 1/det_A

    do j = 1, 3
        do i = 1, 3
            T(i,j) = A(j,i)
        end do
    end do

    T11 = T(2,2)*T(3,3) - T(3,2)*T(2,3)
    T21 = -(T(1,2)*T(3,3) - T(3,2)*T(1,3))
    T31 = T(1,2)*T(2,3) - T(2,2)*T(1,3)
    T12 = -(T(2,1)*T(3,3) - T(3,1)*T(2,3))
    T22 = T(1,1)*T(3,3) - T(3,1)*T(1,3)
    T32 = -(T(1,1)*T(2,3) - T(2,1)*T(1,3))
    T13 = T(2,1)*T(3,2) - T(3,1)*T(2,2)
    T23 = -(T(1,1)*T(3,2) - T(3,1)*T(1,2))
    T33 = T(1,1)*T(2,2) - T(2,1)*T(1,2)

    A_inv(1,1) = num*T11
    A_inv(2,1) = num*T21
    A_inv(3,1) = num*T31
    A_inv(1,2) = num*T12
    A_inv(2,2) = num*T22
    A_inv(3,2) = num*T32
    A_inv(1,3) = num*T13
    A_inv(2,3) = num*T23
    A_inv(3,3) = num*T33
end subroutine inv_3x3

subroutine slv_syseq(transient, u, gk, gm, gf, totnds, dt, cgitr, tol)
!-----------------------------------------------------------------------
! Solve the system of equations [gk] * {u} = {gf} to give us a solution
! for {u}
!-----------------------------------------------------------------------
    integer, intent(IN) :: transient, totnds, cgitr
    real(8), dimension(:,:), intent(IN) :: gk, gm
    real(8), dimension(:), intent(IN) :: gf
    real(8), intent(IN) :: dt, tol
!-----------------------------------------------------------------------
    real(8), dimension(:), intent(INOUT)  :: U
!-----------------------------------------------------------------------

    select case(transient)
        case (0)
            call slv_pcg(gk, gf, totnds, cgitr, tol, u)
        case (1)
            call slv_trns_imp(gk, gm, gf, u, totnds, dt, cgitr, tol)
        case (2)
            call slv_trns_exp(gk, gm, gf, u, totnds, dt, cgitr, tol)
    end select
end subroutine slv_syseq

subroutine slv_trns_imp(K, N, F, u, totnds, dt, cgitr, tol)
!-----------------------------------------------------------------------
! Calculates the implicit transient solution for {u_new} by solving
! ([N] + [K] * dt) u_new = {F} * dt + [N]{u}
! using the pre-conditioned conjugate gradient method
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: K, N
    real(8), dimension(:), intent(IN) :: F
    integer, intent(IN) :: totnds, cgitr
    real(8), intent(IN) :: dt, tol
!-----------------------------------------------------------------------
    real(8), dimension(:), intent(INOUT) :: u
!-----------------------------------------------------------------------
    real(8), dimension(:,:), allocatable :: Term1
    real(8), dimension(:), allocatable :: Term2
!-----------------------------------------------------------------------

    allocate(term1(totnds,totnds))
    allocate(term2(totnds))
	
    term1 = N + K*dt
    term2 = F*dt + matmul(N, u)

    call slv_pcg(term1, term2, totnds, cgitr, tol, u)
    
end subroutine slv_trns_imp

subroutine slv_trns_exp(K, N, F, u, totnds, dt, cgitr, tol)
!-----------------------------------------------------------------------
! Calculates the explicit transient solution for {u_new} by solving
! [N] u_new = ([N] + [K] * dt) {u} + {F} * dt 
! using the pre-conditioned conjugate gradient method
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: K, N
    real(8), dimension(:), intent(IN) :: F
    integer, intent(IN) :: totnds, cgitr
    real(8), intent(IN) :: dt, tol
!-----------------------------------------------------------------------
    real(8), dimension(:), intent(INOUT) :: u
    real(8), dimension(:), allocatable :: term
!-----------------------------------------------------------------------

    allocate(term(totnds))
    term =  matmul((N + K*dt), u) + F*dt

    call slv_pcg(N, term, totnds, cgitr, tol, u)

end subroutine slv_trns_exp

subroutine slv_pcg(A, b, mat_size, cgitr, tol, x)
!-----------------------------------------------------------------------
! Calculates the solution of [A] {x} = {b} using the pre-conditioned
! conjugate gradient method where [A] is a square and sparse matrix 
! and {x} is an unknown.
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: A
    real(8), dimension(:), intent(IN) :: b
    integer, intent(IN) :: mat_size, cgitr
    real(8), intent(IN) :: tol
!-----------------------------------------------------------------------
    real(8), dimension(:), intent(INOUT) :: x
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable :: x_old, A_p, r, r_new, p, precond
    integer :: i, t
    real(8) :: alpha, betta, p_A_p, r_new_r_new, r_r, res
!-----------------------------------------------------------------------
    ! Serial Code with no OpenMP
    allocate(x_old(mat_size))
    allocate(A_p(mat_size))
    allocate(r(mat_size))
    allocate(r_new(mat_size))
    allocate(p(mat_size))
    allocate(precond(mat_size))

    ! initial guess of x
    x = b

    ! calculate r = b-A*x;
    r = b - matmul(A,x)
    
    do i = 1, mat_size
        precond(i) = 1.0/A(i,i)
    end do
    
    p = precond*r

    do t = 1, cgitr

        x_old = x;

        !calculate r'*r
        r_r = sum(r*r*precond)

        ! calculate p'*A*p
        A_p = matmul(A,p)
        p_A_p = sum(p*A_p)

        alpha=0
        if (p_A_p == 0) then
            p_A_p = 1;
            write(*,'(a)') "WARNING: MATRIX SINGULARITY ISSUE"
        else
            alpha = r_r/p_A_p
        end if

        ! calculate x = x + alpha*p;
        x = x + alpha * p
        r_new = r - alpha * A_p
        !if (r_new[1] <= 0.00001)
        !{ break;}

        ! calculate r_new'*r_new
        r_new_r_new = sum(r_new*r_new*precond)

        betta = 0
        if (r_r <= 0.0) then
            r_r = 1
            write(*,'(a)') "ERROR: MATRIX SINGULARITY ISSUE"
        else
            betta = r_new_r_new/r_r
        end if
        
        p = precond * r_new + betta * p
        r = r_new
        res = sum(abs(x_old-x))
        if (res <= tol) then
            exit
        end if
    end do
    if (t < cgitr) then
        write(*, '(a, i0, a)') "Converged matrix solution obtained after ", t, " itrations..."
    else
        write(*,'(a)') "ERROR: No convergence within the specified tolerance was obtained for the PCG slvr"
        call exit(0)
    end if
end subroutine slv_pcg

subroutine mul_trnsp(A, B, C, s)
!-----------------------------------------------------------------------
! Calculates the solution of {A} * transpose{B} = [C]
! where {A} and {B} are arrays of rank 1 and [C] is a square matrix
!-----------------------------------------------------------------------
    real(8), dimension(:), intent(IN) :: A, B
    integer, intent(IN) :: s
!-----------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(OUT)  :: C
!-----------------------------------------------------------------------
    integer :: i, j
!-----------------------------------------------------------------------

    allocate(C(s,s))
    
    do j = 1, s
        do i = 1, s
            C(i,j) = A(j) * B(i)
        end do
    end do
end subroutine mul_trnsp

end module eqslvrs_lib
