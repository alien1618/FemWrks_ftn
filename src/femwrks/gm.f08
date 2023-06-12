module gm_ops
implicit none
contains

subroutine asmbl_diff_mat(gk, supg, totnds, &
            elems, elemnds, intrps, totquadpnts, d, v, stab_prmtr)
!----------------------------------------------------------------------------------
! Assemble the global advection-diffusion matrix
!----------------------------------------------------------------------------------
    use krnl_struct
    use eq_slvrs
!----------------------------------------------------------------------------------
    logical, intent(IN) :: supg
    integer, intent(IN) :: totnds, elemnds, totquadpnts
    real(8), intent(IN) :: d, stab_prmtr
    real(8), dimension(:,:), intent(IN) :: v
    integer, dimension(:,:), intent(IN) :: elems
    type(kernel), dimension(:), intent(IN) :: intrps
!----------------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(OUT)  :: gk !global stiffness matrix
!----------------------------------------------------------------------------------
    integer :: g, e, mm, n, cm, cn, cj
    real(8) :: vx, vy, vz, tau_supg, norm_v
    real(8), dimension(:,:), allocatable :: lk_diff, lk_adv, lk_supg !quadrature stiffness matrices
    real(8), dimension(:), allocatable :: dn_v
!----------------------------------------------------------------------------------
   
    allocate(gk(totnds, totnds))
    allocate(lk_diff(elemnds, elemnds))
    allocate(lk_adv(elemnds, elemnds))
    allocate(lk_supg(elemnds, elemnds))
    allocate(dn_v(elemnds))
            
    gk(:,:) = 0

    do g = 1, totquadpnts
        e = intrps(g)%elem_num
        vx = 0.0
        vy = 0.0
        vz = 0.0
        do mm = 1, elemnds
            cm = elems(e,mm)
            vx = vx + intrps(g)%n(mm)*v(cm,1)
            vy = vy + intrps(g)%n(mm)*v(cm,2)
            vz = vz + intrps(g)%n(mm)*v(cm,3)
        end do

        !get local quad matrices
        lk_diff = d*(intrps(g)%nx_nx + intrps(g)%ny_ny + intrps(g)%nz_nz)
        lk_adv = vx*intrps(g)%n_nx + vy * intrps(g)%n_ny + vz * intrps(g)%n_nz

        !asmbl local quad matrices into global matrices
        do n = 1, elemnds
            do mm = 1, elemnds
                cm = elems(e,mm)
                cn = elems(e,n)
                gk(cm,cn) = gk(cm,cn) + intrps(g)%quad_wt * (lk_diff(mm,n) - lk_adv(mm,n))
            end do
        end do
        if (supg .eqv. .true.) then
            
            !calculate local stabilization matrices
            norm_v = sqrt(vx*vx + vy*vy + vz*vz)
            tau_supg = 0.0
            if (norm_v == 0.0) then
                tau_supg = 0.5*stab_prmtr
            else
                tau_supg = stab_prmtr/norm_V
            end if
            
            lk_supg(:,:) = 0.0

            do n = 1, elemnds
                cj = elems(e,n)
                dn_v(n) = v(cj,1)*intrps(g)%nx(n) + v(cj,2)*intrps(g)%ny(n) + v(cj,3)*intrps(g)%nz(n)
            end do
            call mul_trnsp(dn_v, dn_v, lk_supg, elemnds)
            !M_quad_SUPG = mul_trnsp(dN_V, intrps[g].N, elem_nds);
            !for (int n = 1; n <= elem_nds; n++)
            !{F_quad_SUPG[n] = dN_V[n]*Q;}

            !asmbl local stabilization matrices into global matrices
            do n = 1, elemnds
                do mm = 1, elemnds
                    cm = elems(e,mm)
                    cn = elems(e,n)
                    gk(cm,cn) = gk(cm,cn) + intrps(g)%quad_wt * (tau_supg * lk_supg(mm,n))
                end do
            end do
        end if
    end do
end subroutine asmbl_diff_mat

subroutine asmbl_stiff_mat(gk, dim, load_type, moe, nu, t, totnds, krnls, totquadpnts, elems, elemnds)
!----------------------------------------------------------------------------------
! Assemble the global stiffness matrix
!----------------------------------------------------------------------------------
    use eq_slvrs
    use krnl_struct
!----------------------------------------------------------------------------------
    integer, intent(IN) :: dim, load_type, totnds, elemnds, totquadpnts
    real(8), intent(IN) :: moe, nu, t
    integer, dimension(:,:), intent(IN) :: elems
    type(kernel), dimension(:), intent(IN) :: krnls
!----------------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(OUT)    :: gk !global stiffness matrix
!----------------------------------------------------------------------------------    
    integer :: dof, j, kk, l, itr, s, e
    integer, dimension(:), allocatable :: dof_mat
    real(8), dimension(:,:), allocatable :: dmat, b, lk_diff
!----------------------------------------------------------------------------------
    dof = dim
    allocate(gk(dof*totnds, dof*totnds))
    allocate(dof_mat(elemnds*dim))
    gk(:,:) = 0

    select case (dim)
        case (2)
            allocate (dmat(3, 3))
            allocate (b(dof*elemnds, 3))
        case (3)
            allocate (dmat(6, 6))
            allocate (b(dof*elemnds, 6))
    end select    
    allocate (lk_diff(dof*elemnds, dof*elemnds))

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
                call get_d2d(load_type,moe,nu,t, dmat)
                call get_b2d(elemnds, krnls(s)%nx, krnls(s)%ny, b)
            case (3)
                call get_d3d(moe,nu, dmat)
                call get_b3d(elemnds, krnls(s)%nx, krnls(s)%ny, krnls(s)%nz, b)
        end select
        
        lk_diff = matmul(matmul(b, dmat), transpose(b))
        do l = 1, dof*elemnds
            do j = 1, dof*elemnds
                gk(dof_mat(j),dof_mat(l)) = gk(dof_mat(j),dof_mat(l)) + lk_diff(l,j)*krnls(s)%quad_wt
            end do
        end do
    end do
end subroutine asmbl_stiff_mat

subroutine asmbl_mass_mat(gm, supg, totnds, &
            elems, elemnds, intrps, totquadpnts, v, stab_prmtr)
!----------------------------------------------------------------------------------
! Assemble the global mass matrix
!----------------------------------------------------------------------------------
    use krnl_struct
    use eq_slvrs
!----------------------------------------------------------------------------------
    logical, intent(IN) :: supg
    integer, intent(IN) :: totnds, elemnds, totquadpnts
    real(8), intent(IN) :: stab_prmtr
    real(8), dimension(:,:), intent(IN) :: v
    integer, dimension(:,:), intent(IN) :: elems
    type(kernel), dimension(:), intent(IN) :: intrps
!----------------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(OUT) :: gm !global mass matrix
!----------------------------------------------------------------------------------
    integer :: g, e, mm, n, cm, cn, cj
    real(8) :: vx, vy, vz, tau_supg, norm_v
    real(8), dimension(:,:), allocatable :: lm, lm_supg !quadrature mass matrices
    real(8), dimension(:), allocatable :: dn_v
!----------------------------------------------------------------------------------

    allocate(gm(totnds,totnds))
    allocate(lm(elemnds, elemnds))
    allocate(dn_v(elemnds))
    allocate(lm_supg(elemnds, elemnds))
             
    gm(:,:) = 0

    do g = 1, totquadpnts
        e = intrps(g)%elem_num
        vx = 0.0
        vy = 0.0
        vz = 0.0
        do mm = 1, elemnds
            cm = elems(e,mm)
            vx = vx + intrps(g)%N(mm)*v(cm,1)
            vy = vy + intrps(g)%N(mm)*v(cm,2)
            vz = vz + intrps(g)%N(mm)*v(cm,3)
        end do

        !get global mass matrix
        do n = 1, elemnds
            do mm = 1, elemnds
                cm = elems(e,mm)
                cn = elems(e,n)
                gm(cm,cn) = gm(cm,cn) + intrps(g)%quad_wt * intrps(g)%n_n(mm,n)
            end do
        end do

        if (supg .eqv. .true.) then
            !calculate local stabilization matrices
            norm_v = sqrt(vx*vx + vy*vy + vz*vz)
            tau_supg = 0.0
            if (norm_v == 0.0) then
                tau_supg = 0.5*stab_prmtr
            else
                tau_supg = stab_prmtr/norm_V
            end if
            
            lm_supg(:,:) = 0.0

            do n = 1, elemnds
                cj = elems(e,n)
                dn_v(n) = v(cj,1)*intrps(g)%nx(n) + v(cj,2)*intrps(g)%ny(n) + v(cj,3)*intrps(g)%nz(n)
            end do
            !lm_supg = mul_trnsp(dn_v, intrps[g].n, elem_nds);
            !for (int n = 1; n <= elem_nds; n++)
            !{lf_supg[n] = dn_v[n]*q;}

            !asmbl local stabilization matrices into global matrices
            do n = 1, elemnds
                do mm = 1, elemnds
                    cm = elems(e,mm)
                    cn = elems(e,n)
                    gm(cm,cn) = gm(cm,cn) + intrps(g)%quad_wt * tau_supg * lm_supg(mm,n)
                end do
            end do
        end if
    end do
end subroutine asmbl_mass_mat

subroutine asmbl_force_mat(gf, q, totnds, krnls, totquadpnts, elems, elemnds)
!----------------------------------------------------------------------------------
! Assemble the global force matrix
!----------------------------------------------------------------------------------
    use krnl_struct
!----------------------------------------------------------------------------------
    integer, intent(IN) :: elemnds, totquadpnts, totnds
    integer, dimension(:,:), intent(IN) :: elems
    real(8), dimension(:), intent(IN) :: q
    type(kernel), dimension(:), intent(IN) :: krnls
!----------------------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(OUT) :: gf   !global force matrix
!----------------------------------------------------------------------------------
    real(8), dimension(:), allocatable :: lf   !quadrature force matrix
    real(8) :: lq
    integer :: mm, cm, n, cn, s, e
!----------------------------------------------------------------------------------
    
    allocate(gf(totnds))
    allocate(lf(elemnds))

    gf(:) = 0

    do s = 1, totquadpnts
        e = krnls(s)%elem_num
        lq = 0.0;
        do mm = 1, elemnds
            cm = elems(e,mm)
            lq = lq + krnls(s)%n(mm)*q(cm)
        end do

        lf = krnls(s)%n*lq
        do n = 1, elemnds
            cn = elems(e,n)
            gf(cn) = gf(cn) + krnls(s)%quad_wt * (lf(n))
        end do
    end do
end subroutine asmbl_force_mat

subroutine set_dbc_trans(totnds, bcs, gk, gm, gf)
!----------------------------------------------------------------------------------
! Assign boundary conditions to the global stiffness and force matrices
!----------------------------------------------------------------------------------
    use bc_struct
!----------------------------------------------------------------------------------
    integer, intent(IN) :: totnds
    type(bc), intent(IN) :: bcs
!----------------------------------------------------------------------------------
    real(8), dimension(:,:), intent(INOUT) :: gk, gm
    real(8), dimension(:), intent(INOUT) :: gf
!----------------------------------------------------------------------------------
    integer :: i, j, h
!----------------------------------------------------------------------------------

    do j = 1, bcs%tot
        do i = 1, totnds
            gf(i) = gf(i) - (gk(bcs%pnts(j),i)*bcs%vals(j))
        end do
    end do
    do h = 1, bcs%tot
        gf(bcs%pnts(h)) = bcs%vals(h)
    end do
    do h = 1, bcs%tot
        do j = 1, totnds
            gk(bcs%pnts(h),j) = 0.0
            gk(j,bcs%pnts(h)) = 0.0
            gm(bcs%pnts(h),j) = 0.0
            gm(j,bcs%pnts(h)) = 0.0
        end do
        gk(bcs%pnts(h),bcs%pnts(h)) = 1.0
        gm(bcs%pnts(h),bcs%pnts(h)) = 1.0
    end do
end subroutine set_dbc_trans

subroutine set_dbc_ss(dof_mat, totnds, bcs, gk, gf)
!----------------------------------------------------------------------------------
! Assign boundary conditions to the global stiffness, mass and 
! force matrices
!----------------------------------------------------------------------------------
    use bc_struct
!----------------------------------------------------------------------------------
    integer, intent(IN) :: dof_mat, totnds
    type(bc), intent(IN) :: bcs
!----------------------------------------------------------------------------------
    real(8), dimension(:,:), intent(INOUT) :: gk
    real(8), dimension(:), intent(INOUT) :: gf
!----------------------------------------------------------------------------------
    integer :: i, j, s, h, itr
!----------------------------------------------------------------------------------
    do j = 1, bcs%tot
        do s = 1, dof_mat
            itr = dof_mat-s
            do i = 1, dof_mat*totnds
                gf(i) = gf(i)-(gk(dof_mat*bcs%pnts(j)-itr,i)*bcs%vals(j))
            end do
        end do
    end do
        
    do h = 1, bcs%tot
        do s = 1, dof_mat
            itr = dof_mat-s
            gf(dof_mat*bcs%pnts(h)-itr) = bcs%vals(h)
        end do
    end do
        
    do h = 1, bcs%tot
        do s = 1, dof_mat
            itr = dof_mat-s
            do i = 1, dof_mat*totnds
                gk(dof_mat*bcs%pnts(h)-itr,i) = 0.0
                gk(i,dof_mat*bcs%pnts(h)-itr) = 0.0
            end do
            gk(dof_mat*bcs%pnts(h)-itr,dof_mat*bcs%pnts(h)-itr) = 1.0
        end do
    end do
end subroutine set_dbc_ss

subroutine get_d2d(load_type, moe, nu, t, dmat)
!----------------------------------------------------------------------------------
! Constructs the elasticity matrix in 2D
!----------------------------------------------------------------------------------
    real(8), intent(IN) :: moe, nu, t
    integer, intent(IN) :: load_type
!----------------------------------------------------------------------------------
    real(8), dimension(3,3), intent(OUT) :: dmat
!----------------------------------------------------------------------------------
    real(8) :: d1, d2, d3
!----------------------------------------------------------------------------------

    d1 = 1
    d2 = 1
    d3 = 1
    select case (load_type)
        case (1) !plane stress
            d1 = moe*t/(1-(nu*nu))
            d2 = nu*D1
            d3 = D1*(1-nu)/2
        case (2) !plane strain
            d1 = (1-nu)*moe*t/((1+nu)*(1-2*nu))
            d2 = nu*moe*t/((1+nu)*(1-2*nu))
            d3 = (1-2*nu)*moe*t/2*((1+nu)*(1-2*nu))
        case default
            write(*,'(a)') "Error: load type is invalid"
    end select
    dmat(1,1) = d1; dmat(2,1) = d2; dmat(3,1) = 0;
    dmat(1,2) = d2; dmat(2,2) = d1; dmat(3,2) = 0;
    dmat(1,3) = 0; dmat(2,3) = 0; dmat(3,3) = d3;
end subroutine get_d2d

subroutine get_d3d(moe, nu, dmat)
!----------------------------------------------------------------------------------
! Constructs the elasticity matrix in 3D
!----------------------------------------------------------------------------------
    real(8), intent(IN) :: moe, nu
!----------------------------------------------------------------------------------
    real(8), dimension(6,6), intent(OUT) :: dmat
!----------------------------------------------------------------------------------
    real(8) :: t1, d1, d2, d3
!----------------------------------------------------------------------------------

    t1 = (moe/((1+nu)*(1-2*nu)))
    D1 = (1-nu)*t1
    D2 = nu*t1
    D3 = 0.5*(1-2*nu)*t1
    dmat(1,1) = D1; dmat(2,1) = D2; dmat(3,1) = D2;
    dmat(1,2) = D2; dmat(2,2) = D1; dmat(3,2) = D2;
    dmat(1,3) = D2; dmat(2,3) = D2; dmat(3,3) = D1;
    dmat(4,4) = D3; dmat(5,5) = D3; dmat(6,6) = D3;
end subroutine get_d3d

subroutine get_b3d(volnds, dNdx, dNdy, dNdz, B)
!----------------------------------------------------------------------------------
! Constructs the elasticity constitutive matrix in 3D
!----------------------------------------------------------------------------------
    integer, intent(IN) :: volnds
    real(8), dimension(:), intent(IN) :: dNdx, dNdy, dNdz
!----------------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(OUT) :: B
!----------------------------------------------------------------------------------
    integer :: k, dof, i
!----------------------------------------------------------------------------------

    k = 1;
    dof = 3;
    allocate(b(dof*volnds, 6))

    do i = 1, volnds
        b(k,1) = dNdx(i)
        b(dof*i-1,2) = dNdy(i)
        b(dof*i,3) = dNdz(i)
        b(k,4) = dNdy(i)
        b(k+1,4) = dNdx(i)
        b(k+1,5) = dNdz(i)
        b(k+2,5) = dNdy(i)
        b(k,6) = dNdz(i)
        b(k+2,6) = dNdx(i)
        k = k + 3;
    end do
end subroutine get_b3d

subroutine get_b2d(surfnds, dNdx, dNdy, B)
!----------------------------------------------------------------------------------
! Constructs the elasticity constitutive matrix in 3D
!----------------------------------------------------------------------------------
    integer, intent(IN) :: surfnds
    real(8), dimension(:), intent(IN) :: dNdx, dNdy
!----------------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(OUT) :: B
!----------------------------------------------------------------------------------
    integer :: k, dof, i
!----------------------------------------------------------------------------------

    k = 1;
    dof = 2;
    allocate(b(dof*surfnds, 3))

    do i = 1, surfnds
        b(k,1) = dNdx(i)
        b(dof*i,2) = dNdy(i)
        b(dof*i,3) = dNdx(i)
        b(k,3) = dNdy(i)
        k = k + 2
    end do
end subroutine get_b2d

end module gm_ops
