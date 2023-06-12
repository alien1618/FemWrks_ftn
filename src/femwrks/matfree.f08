module lmat_struct
!--------------------------------------------------------------------------------------------
! Defining local matrix data structure
!--------------------------------------------------------------------------------------------
implicit none
private
public local

type local
    integer::tot
    real(8), dimension(:,:), allocatable :: mat
end type
end module lmat_struct

module nbr_struct
!--------------------------------------------------------------------------------------------
! Defining neighbours to an entity (i.e. a node or element) data structure
!--------------------------------------------------------------------------------------------
implicit none
private
public nbr

type nbr
    integer::totnbrs
    integer, dimension(:), allocatable :: nbrs
end type
end module nbr_struct

module matfree_ops
implicit none
contains

subroutine get_ndkrnls(totnds, elems, totelems, elemnds, ndkrnls)
!--------------------------------------------------------------------------------------------
! This function collects all the elements that contribute to a node
!--------------------------------------------------------------------------------------------
    use nbr_struct    
!--------------------------------------------------------------------------------------------
    integer, intent(IN) :: totnds, totelems, elemnds
    integer, dimension(:,:), intent(IN) :: elems
!--------------------------------------------------------------------------------------------
    type(nbr), dimension(:), allocatable, intent(OUT) :: ndkrnls
!--------------------------------------------------------------------------------------------
    integer :: i, j, k, ndnum
!--------------------------------------------------------------------------------------------

    allocate(ndkrnls(totnds))
    do i = 1, totnds
        ndkrnls(i)%totnbrs = 0
        allocate(ndkrnls(i)%nbrs(0))
    end do
    do j = 1, totelems
        do k = 1, elemnds
            ndnum = elems(j,k)
            ndkrnls(ndnum)%totnbrs = ndkrnls(ndnum)%totnbrs + 1
            ndkrnls(ndnum)%nbrs = [ndkrnls(ndnum)%nbrs,j]
        end do
    end do
end subroutine get_ndkrnls

subroutine get_lmat(k, v, krnls, totquadpnts, elems, totelems, elemnds, &
                    dt, transient, supg, stab_prmtr, lmat)
!--------------------------------------------------------------------------------------------
! This function calculates a compact and dense local matrices at each element without
! the need for calculating global stiffness or mass matrices.
! only the force matrix is calculated globaly since it's just a 1D array
!--------------------------------------------------------------------------------------------
    use krnl_struct
    use lmat_struct
    use nbr_struct
    use eq_slvrs
!--------------------------------------------------------------------------------------------
    real(8), dimension(:,:), intent(IN) :: v                
    real(8), intent(IN) :: k, dt, stab_prmtr
    integer, intent(IN) :: totquadpnts, totelems, elemnds
    integer, intent(IN) :: transient
    integer, dimension(:,:), intent(IN) :: elems
    type(kernel), dimension(:), intent(IN) :: krnls
!--------------------------------------------------------------------------------------------
    type(local), dimension(:), allocatable, intent(out)  :: lmat !local matrices at element level
!--------------------------------------------------------------------------------------------
    logical :: supg
    real(8), dimension(:), allocatable :: dN_v
    real(8), dimension(:,:), allocatable :: lk_supg, lm_supg
    integer :: i, g, cj, cm, e, mm, n
    real(8) :: vx, vy, vz, norm_V, tau_supg
!--------------------------------------------------------------------------------------------
     
    allocate(lmat(totelems))
    
    do i = 1, totelems
        allocate(lmat(i)%mat(elemnds, elemnds))
        lmat(i)%mat(:,:) = 0.0
    end do
    
    if (supg .eqv. .true.) then
        allocate (dn_v(elemnds))
        allocate (lk_supg(elemnds, elemnds))
        allocate (lm_supg(elemnds, elemnds))
    end if
    
    do g = 1, totquadpnts
        !----------------------------------------------------------
        !interpolate material properties at quad point
        !----------------------------------------------------------
        vx = 0;
        vy = 0;
        vz = 0;
        e = krnls(g)%elem_num
        do mm = 1, elemnds
            cm = elems(e,mm)
            vx = vx + krnls(g)%N(mm)*v(cm,1)
            vy = vy + krnls(g)%N(mm)*v(cm,2)
            vz = vz + krnls(g)%N(mm)*v(cm,3)
        end do

        !----------------------------------------------------------
        ! cout << "get quad STIFFNESS AND FORCE MATRICES" << endl;
        !----------------------------------------------------------
        if (transient == 0) then
            !CALCULATES IMPLICIT FORM: 0 = (K*U + Q)
            
            !First: calculate K
            lmat(e)%mat = lmat(e)%mat + k*(krnls(g)%nx_nx+krnls(g)%ny_ny+krnls(g)%nz_nz)           
            lmat(e)%mat = lmat(e)%mat + (krnls(g)%n_nx*vx + krnls(g)%n_ny*vy+ krnls(g)%n_nz*vz)
            lmat(e)%mat = lmat(e)%mat * krnls(g)%quad_wt
        else
            !CALCULATES IMPLICIT FORM: U_new(M + K*dt) = ((M*U_now) + (Q*dt))
            
            !First: calculate (M + K*dt)
            lmat(e)%mat = lmat(e)%mat + dt*k*(krnls(g)%nx_nx + krnls(g)%ny_ny + krnls(g)%nz_nz)           
            lmat(e)%mat = lmat(e)%mat + dt*(krnls(g)%n_nx*vx + krnls(g)%n_ny*vy + krnls(g)%n_nz*vz)
            lmat(e)%mat = lmat(e)%mat + krnls(g)%n_n
            lmat(e)%mat = lmat(e)%mat * krnls(g)%quad_wt

            if (supg .eqv. .true.) then
                norm_V = ((vx*vx + vy*vy + vz*vz)**0.5)+1e-5
                            tau_supg = 0.0
            if (norm_v == 0.0) then
                tau_supg = 0.5*stab_prmtr
            else
                tau_supg = stab_prmtr/norm_V
            end if
                lk_supg(:,:) = 0.0
                do n = 1, elemnds
                    cj = elems(e,n)
                    dN_V(n) = v(cj,1)*krnls(g)%Nx(n) + v(cj,2)*krnls(g)%Ny(n) + v(cj,3)*krnls(g)%Nz(n)
                end do
                call mul_trnsp(dn_v, dn_v, lk_supg, elemnds)
                lmat(e)%mat = lmat(e)%mat + krnls(g)%quad_wt*dt*tau_supg*lk_supg
            end if
        end if
    end do
end subroutine get_lmat

subroutine asmbl_gf_ebe(u, q, totnds, krnls, totquadpnts, elems, elemnds, dt, transient, gf)
!--------------------------------------------------------------------------------------------
! This function calculates a compact and dense local matrices at each element without
! the need for calculating global stiffness or mass matrices.
! only the force matrix is calculated globaly since it's just a 1D array
!--------------------------------------------------------------------------------------------
    use krnl_struct
    use lmat_struct
    use nbr_struct
    use eq_slvrs
    use omp_lib
!--------------------------------------------------------------------------------------------
    real(8), dimension(:), intent(IN) :: u
    real(8), dimension(:), intent(IN) :: q
    real(8), intent(IN) :: dt
    integer, intent(IN) :: totnds, totquadpnts, elemnds
    integer, intent(IN) :: transient
    integer, dimension(:,:), intent(IN) :: elems
    type(kernel), dimension(:), intent(IN) :: krnls
!--------------------------------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(out)   :: gf !global force matrix
!--------------------------------------------------------------------------------------------    
    integer :: i, j, g, ci, cj, cm, e, mm
    real(8) :: lq, row_sum
!--------------------------------------------------------------------------------------------

    allocate(gf(totnds))
    gf(:) = 0
    
  !  !$omp parallel do private(g)
    do g = 1, totquadpnts
        lq = 0;
        e = krnls(g)%elem_num
        do mm = 1, elemnds
            cm = elems(e,mm)
            lq = lq + krnls(g)%N(mm)*q(cm)
        end do

        !----------------------------------------------------------
        ! cout << "get quad STIFFNESS AND FORCE MATRICES" << endl;
        !----------------------------------------------------------
        if (transient == 0) then
            !Second: calculate (Q)
            do j = 1, elemnds
                cj = elems(e,j)
                gf(cj) = gf(cj) + krnls(g)%quad_wt*(lq)
            end do
        else
            !for (int n = 1; n <= c_totnds; n++)
            !{F_quad_SUPG[n] = dN_V[n]*Q;}
            !Second: calculate ((M * U_now) + (Q*dt))
            do j = 1,elemnds
                cj = elems(e,j)
                row_sum = 0;
                do i = 1, elemnds
                    ci = elems(e,i)
                    row_sum = row_sum + krnls(g)%n_n(i,j)*u(ci)
                end do
                gf(cj) = gf(cj) + krnls(g)%quad_wt*(row_sum + dt*lq)
            end do
        end if
    end do
  !  !$omp end parallel do
end subroutine asmbl_gf_ebe

subroutine get_lmat_elst(dim, load_type, moe, nu, thick, krnls, totquadpnts, elems, totelems, elemnds, lmat)
!--------------------------------------------------------------------------------------------
! Calculates the local matrices for elasticity
!--------------------------------------------------------------------------------------------
    use krnl_struct
    use nbr_struct
    use eq_slvrs
    use gm_ops
    use lmat_struct
!--------------------------------------------------------------------------------------------
    integer, intent(IN) :: dim, load_type
    real(8), intent(IN) :: moe, nu, thick
    integer, intent(IN) :: totquadpnts, totelems, elemnds
    integer, dimension(:,:), intent(IN) :: elems
    type(kernel), dimension(:), intent(IN) :: krnls
!--------------------------------------------------------------------------------------------
    type(local), dimension(:), allocatable, intent(OUT) :: lmat
!--------------------------------------------------------------------------------------------
    integer :: dof, i, j, e
    integer, dimension(:), allocatable :: dof_mat
    real(8), dimension(:,:), allocatable :: dmat, b
    integer :: itr, kk, s
!--------------------------------------------------------------------------------------------
    !this function slvs the steady state or transient solution of
    !the general trnsprt equation using the finite element method
    !without generating global matrices.
    !it slvs it element-by-element itratively using the pre-conditioned
    !conjugate gradient method. Function assumes there is only one
    !degree of freedom per node. Used for heat/fluid flow simulations.

    dof = dim
    allocate(dof_mat(elemnds*dim))

    allocate(lmat(totelems))
    do i = 1, totelems
        allocate(lmat(i)%mat(dof*elemnds, dof*elemnds))
        lmat(i)%mat(:,:) = 0.0
    end do

    select case (dim)
        case (2)
            allocate (dmat(3, 3))
            allocate (b(dof*elemnds, 3))
        case (3)
            allocate (dmat(6, 6))
            allocate (b(dof*elemnds, 6))
    end select    

    do s = 1, totquadpnts
        e = krnls(s)%elem_num
        kk = 1
        
        do j = 1, elemnds
            do itr = dof-1, 0, -1
                dof_mat(kk) =  dof*elems(e,j)-itr
                kk = kk + 1
            end do
        end do
        
        select case (dim)
            case (2)
                call get_d2d(load_type,moe,nu,thick, dmat)
                call get_b2d(elemnds, krnls(s)%nx, krnls(s)%ny, b)
            case (3)
                call get_d3d(moe,nu, dmat)
                call get_b3d(elemnds, krnls(s)%nx, krnls(s)%ny, krnls(s)%nz, b)
        end select
        lmat(e)%mat = lmat(e)%mat + krnls(s)%quad_wt*matmul(matmul(b, dmat), transpose(b))
    end do

end subroutine get_lmat_elst

subroutine slv_pcg_ebe(dof, F0, totnds, elems, totelems, elemnds, lmat, &
                bcs, cgitr, conv_tol, x)
!--------------------------------------------------------------------------------------------
! Function solves for [K] {x} = {F} to obtain a solution for x using the pre-conditioned 
! gradient method and without constructing global matrices for [k]
!--------------------------------------------------------------------------------------------
    use lmat_struct
    use bc_struct
!--------------------------------------------------------------------------------------------
    integer, intent(IN) :: dof, totnds, elemnds, totelems, cgitr
    real(8), dimension(:), intent(IN) :: F0
    type(local), dimension(:), intent(IN) :: lmat
    integer, dimension(:,:), intent(IN) ::  elems
    type(bc), intent(IN) :: bcs
    real(8), intent(IN) :: conv_tol
!--------------------------------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(out) :: x
!--------------------------------------------------------------------------------------------
    real(8)                                    :: alpha, betta, p_A_p, product_quad, r_new_r_new, r_r, res
    integer                                 :: e, k, i, m, j, itr, h, ff, kk, t
    integer, dimension(:), allocatable      :: dof_mat
    real(8), dimension(:), allocatable         :: b, x_old, A_p, r, r_new, p, precond
!--------------------------------------------------------------------------------------------

    allocate(b(dof*totnds))
    allocate(precond(dof*totnds))
    allocate(x(dof*totnds))
    allocate(x_old(dof*totnds))
    allocate(A_p(dof*totnds))
    allocate(r(dof*totnds))
    allocate(r_new(dof*totnds))
    allocate(p(dof*totnds))
    allocate(dof_mat(dof*elemnds))
    
    A_p(:) = 0.0
    b = F0
    precond(:) = 1.0
    
    do h = 1, bcs%tot
        do m = 1, dof
            itr = dof-m
            b(dof*bcs%pnts(h)-itr) = bcs%vals(h)
        end do
    end do

    !Conjugate Gradient slvr(Serial Code with no OpenMP)
    x = b; !initial guess of x
    do e = 1, totelems
        kk = 1
        do ff = 1, elemnds
            do itr = dof-1,0,-1
                dof_mat(kk) =  dof*elems(e,ff)-itr;
                kk = kk + 1;
            end do
        end do
        do j = 1, dof*elemnds    
            product_quad = 0; 
            do i = 1, dof*elemnds
                product_quad = product_quad + lmat(e)%mat(i,j)*x(dof_mat(i))
            end do
            A_p(dof_mat(j)) = A_p(dof_mat(j)) + product_quad
        end do
    end do

    do k = 1, bcs%tot
        do m = 1, dof
            itr = dof-m
            A_p(dof*bcs%pnts(k)-itr) = x(dof*bcs%pnts(k)-itr)
        end do
    end do
    
    !calculate r = b-A*x;
    r = b - A_p
    p = precond * r 
    do t = 0, cgitr
        x_old = x

        !calculate r'*r
        r_r = sum(precond*r*r)

        !calculate p'*A*p
        A_p(:) = 0
        do e = 1, totelems
            kk = 1
            do ff = 1, elemnds
                do itr = dof-1, 0, -1
                    dof_mat(kk) =  dof*elems(e,ff)-itr
                    kk = kk + 1
                end do
            end do

            do j = 1, dof*elemnds        
                product_quad = 0
                do i = 1,dof*elemnds
                    product_quad = product_quad + lmat(e)%mat(i,j)*p(dof_mat(i))
                end do
                A_p(dof_mat(j)) = A_p(dof_mat(j)) + product_quad
            end do
        end do

        do m = 1, dof
            itr = dof-m
            do k = 1, bcs%tot
                A_p(dof*bcs%pnts(k)-itr) = p(dof*bcs%pnts(k)-itr)
            end do
        end do

        p_A_p = sum(p*A_p)

        alpha=1
        if (p_A_p == 0.0) then
            p_A_p = 1.0
            write(*,*) "WARNING: MATRIX SINGULARITY ISSUE"
        else
            alpha = r_r/p_A_p
        end if

        !calculate x = x + alpha*p;
        x = x + alpha * p
        r_new = r - alpha * A_p
        !if (r_new[1] <= 0.00001)
        !{ break;}

        !calculate r_new'*r_new
        r_new_r_new = sum(precond*r_new*r_new)

        betta = 0
        if (r_r == 0.0) then
            r_r = 1.0
            !write(*,*) "ERROR: MATRIX SINGULARITY ISSUE"
            !exit(0)
        else
            betta = r_new_r_new/r_r
        end if
        
        p = precond*r_new + betta * p
        r = r_new
        res = sum(abs(x_old - x))

        !write(*,*) "RES = ", res
        if (res <= conv_tol) then
            exit
        end if
    end do

    if (t < cgitr) then
        write(*, '(a, i0, a)') "Converged matrix solution obtained after ", t, " itrations..."
    else
        write(*,'(a)') "ERROR: No convergence within the specified tolerance was obtained for the PCG slvr"
        call exit(0)
    end if

end subroutine slv_pcg_ebe

subroutine get_dfrm(u, dof, totnds, nds, deformation)
!--------------------------------------------------------------------------------------------
! calculates 2d and 3d deformation using the calculated elasticity solution
!--------------------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: dof, totnds
    real(8), dimension(:), intent(in) :: u
!--------------------------------------------------------------------------------------------
    real(8), dimension(:,:), intent(inout) :: nds
    real(8), dimension(:), allocatable, intent(out) :: deformation
!--------------------------------------------------------------------------------------------
    integer :: i
!--------------------------------------------------------------------------------------------

    allocate(deformation(totnds))
    select case(dof)
        case (2)
            do i = 1, totnds
                nds(i,1) = nds(i,1) + u(dof*i-1)
                nds(i,2) = nds(i,2) + u(dof*i)
                deformation(i) = sqrt(u(dof*i-1)*u(dof*i-1) + u(dof*i)*u(dof*i))
            end do
        case (3)
            do i = 1, totnds
                nds(i,1) = nds(i,1) + u(dof*i-2)
                nds(i,2) = nds(i,2) + u(dof*i-1)
                nds(i,3) = nds(i,3) + u(dof*i)
                deformation(i) = sqrt(u(dof*i-2)*u(dof*i-2) + u(dof*i-1)*u(dof*i-1) + u(dof*i)*u(dof*i))
            end do
    end select
end subroutine get_dfrm

subroutine get_vms(dof, mat, msh, krnls, totquadpnts, u, vms)
!--------------------------------------------------------------------------------------------
! calculates 2d and 3d von-misses stress the calculated elasticity solution
!--------------------------------------------------------------------------------------------
    use msh_struct
    use bc_struct 
    use krnl_struct
    use nbr_struct
    use mat_struct
    use gm_ops
!--------------------------------------------------------------------------------------------    
    implicit none
    type(mesh), intent(in) :: msh
    integer, intent(in) :: totquadpnts, dof
    type(kernel), dimension(:), intent(in) :: krnls
    type(material), intent(in) :: mat
    real(8), dimension(:), intent(in) :: u
!--------------------------------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(out) :: vms
!--------------------------------------------------------------------------------------------
    real(8), dimension(:), allocatable :: vms_local, u_supp
    integer, dimension(:), allocatable :: dof_mat
    real(8) :: strain_xx, strain_yy, strain_zz, gamma_xy, gamma_xz, gamma_yz
    real(8) :: stress_xx, stress_yy, stress_zz, stress_xy, stress_xz, stress_yz
    real(8) :: dx, dy, dz, dist, ww, sum_w, sum_vms
    real(8), dimension(6) :: strain
    real(8), dimension(:,:), allocatable :: d_mat, lb
    type(nbr), dimension(:), allocatable :: ndkrnls, elems
    integer :: e, i, j, k, g = 0, q, itr
!--------------------------------------------------------------------------------------------

    select case (msh%dim)
        case (3)
            call get_ndkrnls(msh%totnds, msh%vols, msh%totvols, msh%volnds, ndkrnls)
            g = msh%totvols
        case (2)

            call get_ndkrnls(msh%totnds, msh%surfs, msh%totsurfs, msh%surfnds, ndkrnls)
            g = msh%totsurfs
    end select

    allocate(vms_local(totquadpnts))
    select case(msh%dim)
        case (3)
            allocate(d_mat(6,6))
            call get_d3d(mat%moe, mat%nu, d_mat)
            allocate (u_supp(dof*msh%volnds))
            allocate (dof_mat(dof*msh%volnds))
            do g = 1, totquadpnts
                e = krnls(g)%elem_num;
                call get_b3d(msh%volnds, krnls(g)%nx, krnls(g)%ny, krnls(g)%nz, lb);
                k = 1
                do j = 1, msh%volnds
                    do  itr = dof-1, 0, -1
                        dof_mat(k) =  dof*msh%vols(e,j)-itr
                        k = k + 1
                    end do
                end do

                strain_xx = sum(lb(:,1)*u(dof_mat(:)))
                strain_yy = sum(lb(:,2)*u(dof_mat(:)))
                strain_zz = sum(lb(:,3)*u(dof_mat(:)))
                gamma_xy = sum(lb(:,4)*u(dof_mat(:)))
                gamma_xz = sum(lb(:,5)*u(dof_mat(:)))
                gamma_yz = sum(lb(:,6)*u(dof_mat(:)))
                
                strain = [strain_xx, strain_yy, strain_zz, gamma_xy, gamma_xz, gamma_yz]
                
                stress_xx = sum(d_mat(:,1)*strain)
                stress_yy = sum(d_mat(:,2)*strain)
                stress_zz = sum(d_mat(:,3)*strain)
                stress_xy = gamma_xy*d_mat(4,4)
                stress_xz = gamma_xz*d_mat(5,5)
                stress_yz = gamma_yz*d_mat(6,6)

                vms_local(g) = sqrt(((stress_xx-stress_yy)**2)+((stress_yy-stress_zz)**2)+((stress_xx-stress_zz)**2) &
                + 6*(gamma_xy*gamma_xy+gamma_yz*gamma_yz+gamma_xz*gamma_xz))/sqrt(2.0)
            end do
        case (2)
            allocate(d_mat(3,3))
            call get_d2d(mat%plane_type, mat%moe, mat%nu, mat%thickness, d_mat)
            allocate (u_supp(dof*msh%surfnds))
            allocate (dof_mat(dof*msh%surfnds))
            do g = 1, TotQuadPnts
                e = krnls(g)%elem_num
                call get_b2d(msh%surfnds, krnls(g)%nx, krnls(g)%ny, lb)
                k = 1
                do j = 1, msh%surfnds
                    do itr = dof-1, 0, -1
                        dof_mat(k) =  dof*msh%surfs(e,j)-itr
                        k = k + 1
                    end do
                end do
            
                strain_xx = sum(lb(:,1)*u(dof_mat(:)))
                strain_yy = sum(lb(:,2)*u(dof_mat(:)))
                gamma_xy = sum(lb(:,3)*u(dof_mat(:)))
                
                stress_xx = d_mat(1,1)*strain_xx + d_mat(2,1)*strain_yy + d_mat(3,1)*gamma_xy
                stress_yy = d_mat(1,2)*strain_xx + d_mat(2,2)*strain_yy + d_mat(3,2)*gamma_xy
                vms_local(g) = sqrt(((stress_xx-stress_yy)**2)+ stress_xx*stress_xx &
                    + stress_yy*stress_yy + 6*gamma_xy*gamma_xy)/sqrt(2.0)
            end do
    end select

    allocate(elems(g))
    do i = 1, g
        elems(i)%totnbrs = 0
        allocate(elems(i)%nbrs(0))
    end do
    do i = 1, totquadpnts
        elems(krnls(i)%elem_num)%totnbrs = elems(krnls(i)%elem_num)%totnbrs + 1
        elems(krnls(i)%elem_num)%nbrs = [elems(krnls(i)%elem_num)%nbrs, i]
    end do

    allocate(vms(msh%totnds))
    do i = 1, msh%totnds 
        sum_w = 0.0
        sum_vms = 0.0
        do j = 1, ndkrnls(i)%totnbrs
            e = ndkrnls(i)%nbrs(j)
            do g = 1, elems(e)%totnbrs
                q = elems(e)%nbrs(g)
                dx = msh%nds(i,1) - krnls(q)%pntx
                dy = msh%nds(i,2) - krnls(q)%pnty
                dz = msh%nds(i,3) - krnls(q)%pntz
                dist  = sqrt(dx*dx + dy*dy + dz*dz)
                ww = 1.0/dist
                sum_w = sum_w + ww
                sum_vms = sum_vms + ww*vms_local(q)
            end do
        end do
        vms(i) = sum_vms/sum_w
    end do
end subroutine

end module matfree_ops
