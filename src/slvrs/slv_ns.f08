module ns
implicit none
contains

subroutine run_ns()
!---------------------------------------------------------------------------------------------
! Subroutine solves special cases of the general transport equation using
! the finite element method.
!---------------------------------------------------------------------------------------------
    use msh_lib
    use bc_lib
    use prmtrs_lib
    implicit none
!---------------------------------------------------------------------------------------------
    type(mesh)                              :: msh
    type(bc), dimension(6)                  :: bcs      !boundary conditions data structure
    real(8), dimension(:,:), allocatable    :: vel, f
    real(8)                                 :: ro, nu
    type(slvr_prmtrs)                       :: sp 
!---------------------------------------------------------------------------------------------
    !assigning mesh
    call set_msh(msh)

    !assigning material properties
    call set_mat('nu', nu)
    call set_mat('ro', ro)
    
    !assigning boundary conditions
    call set_bc(msh%nds, msh%totnds, bcs)
    
    !initialize field variables
    call set_v0(msh%nds, msh%totnds, vel)
    vel = vel + 1e-6 !adding a small value to avoid singular matrices
    allocate(f(msh%totnds,3))      
    f(:,:) = 0.0
    
    !assigning solver parameters
    call set_ctrls(sp)
    
    !running solver
    call slv_ns(msh, bcs, sp, nu, ro, vel, f)
end subroutine run_ns

subroutine slv_ns(msh, bcs, sp, nu, ro, vel, f)
!---------------------------------------------------------------------------------------------
! Subroutine solves special cases of the general transport equation using
! the finite element method.
!---------------------------------------------------------------------------------------------
    use msh_lib
    use bc_lib
    use krnl_lib
    use quad_lib
    use gm_lib
    use matfree_lib
    use eqslvrs_lib
    use prmtrs_lib
    implicit none
!---------------------------------------------------------------------------------------------
    type(mesh), intent(in)                    :: msh
    type(bc), dimension(6), intent(in)        :: bcs      !boundary conditions data structure
    type(slvr_prmtrs), intent(in)             :: sp
    real(8), dimension(:,:), intent(in)       :: f
    real(8), intent(in)                       :: ro, nu
    real(8), dimension(:,:), intent(inout)    :: vel
!---------------------------------------------------------------------------------------------
    integer                                   :: t, cntr = 1
    integer                                   :: elemquadpnts, totquadpnts
    real(8), dimension(:), allocatable        :: eps, eta, zta, w, gf, v
    real(8), dimension(:,:), allocatable      :: gk, gm
    type(kernel), dimension(:), allocatable   :: krnls
    real(8), dimension(:), allocatable        :: p, dpdx, dpdy, dpdz
    character(len=50)                         :: v_fname = 'v', p_fname = 'p', vec_fname = 'vec'
!---------------------------------------------------------------------------------------------
    write(*,'(a)') "Initializing field variables..."
    allocate(v(msh%totnds))   
    allocate(p(msh%totnds))
    allocate(dpdx(msh%totnds))
    allocate(dpdy(msh%totnds))
    allocate(dpdz(msh%totnds))
    v(:) = sqrt(vel(:,1)**2 + vel(:,2)**2 + vel(:,3)**3)
    p(:) = 0.0
    dpdx(:) = 0.0
    dpdy(:) = 0.0
    dpdz(:) = 0.0
!---------------------------------------------------------------------------------------------
    !quadrature
    write(*,'(a)') "Generating quadrature points..."
    call gen_quadpnts(msh%dim, msh%surfshape, msh%order, eps, eta, zta, w, elemquadpnts)   
!---------------------------------------------------------------------------------------------
    !interpolation
    write(*,'(a)') "Computing interpolants..."
    select case (msh%dim)
        case (2)
            call get_intrps(msh%dim, msh%nds, msh%surfs, msh%totsurfs, &
                msh%surfnds, eps, eta, zta, w, elemquadpnts, krnls, totquadpnts)
        case (3)
            call get_intrps(msh%dim, msh%nds, msh%vols, msh%totvols, &
                msh%volnds, eps, eta, zta, w, elemquadpnts, krnls, totquadpnts)
    end select

    if (sp%vtk .eqv. .true.) then
        call prnt_vtk(v, msh, v_fname, 0) 
        call prnt_vtk(p, msh, p_fname, 0)
        call prnt_vec_vtk(vel, msh, vec_fname, 0)   
    else
        call prnt_pltctrl(msh%nds, sp%nt, sp%prnt_frq)
        call prnt_elems_txt(msh%surfs, msh%totsurfs)
        call prnt_pnts_txt(v, msh%nds, msh%totnds, v_fname, 0)
        call prnt_pnts_txt(p, msh%nds, msh%totnds, p_fname, 0)
    end if
    !-----------------------------------------------------------------------------------------
    !constructing and solving global matrices
    allocate(gk(msh%totnds, msh%totnds))
    allocate(gm(msh%totnds, msh%totnds))
    allocate(gf(msh%totnds))

    do t = 1, sp%nt
        write(*,'(a,i0,a,i0)')"Processing Time Step ", t ," of ", sp%nt
        select case(msh%dim)
            case (2)               
                !write(*,'(a)') "Assembling global matrices..."
                call asmbl_diff_mat(gk, sp%supg, msh%totnds, &
                    msh%surfs, msh%surfnds, krnls, totquadpnts, nu, vel, sp%stab_prmtr)
                call asmbl_mass_mat(gm, sp%supg, msh%totnds, &
                    msh%surfs, msh%surfnds, krnls, totquadpnts, vel, sp%stab_prmtr)
                    
                !write(*,'(a)') "Calculating Intermediate Vx..."
                gf = f(:,1)
                call set_dbc_trans(msh%totnds, bcs(1), gk, gm, gf)
                call slv_syseq(1, vel(:,1), gk, gm, gf, msh%totnds, sp%dt, sp%cgitrs, sp%tol)

                !write(*,'(a)') "Calculating Intermediate Vy..."
                gf = f(:,2)
                call set_dbc_trans(msh%totnds, bcs(2), gk, gm, gf)
                call slv_syseq(1, vel(:,2), gk, gm, gf, msh%totnds, sp%dt, sp%cgitrs, sp%tol)

                !write(*,'(a)') "Solving the Poisson Equation for Pressure..."
                call slv_p(P, dpdx, dpdy, dpdz, bcs(4), sp%dt, vel, &
                    msh%totnds, msh%surfs, msh%surfnds, krnls, totquadpnts, ro, sp%cgitrs, sp%tol)
            case (3)
                !write(*,'(a)') "Assembling global matrices..."
                call asmbl_diff_mat(gk, sp%supg, msh%totnds, &
                    msh%vols, msh%volnds, krnls, totquadpnts, nu, vel, sp%stab_prmtr)
                call asmbl_mass_mat(gm, sp%supg, msh%totnds, &
                    msh%vols, msh%volnds, krnls, totquadpnts, vel, sp%stab_prmtr)
                                
                !write(*,'(a)') "Calculating Intermediate Vx..."
                gf = f(:,1)
                call set_dbc_trans(msh%totnds, bcs(1), gk, gm, gf)
                call slv_syseq(1, vel(:,1), gk, gm, gf, msh%totnds, sp%dt, sp%cgitrs, sp%tol)

                !write(*,'(a)') "Calculating Intermediate Vy..."
                gf = f(:,2)
                call set_dbc_trans(msh%totnds, bcs(2), gk, gm, gf)
                call slv_syseq(1, vel(:,2), gk, gm, gf, msh%totnds, sp%dt, sp%cgitrs, sp%tol)

                !write(*,'(a)') "Calculating Intermediate Vz..."
                gf = f(:,3)
                call set_dbc_trans(msh%totnds, bcs(3), gk, gm, gf)
                call slv_syseq(1, vel(:,3), gk, gm, gf, msh%totnds, sp%dt, sp%cgitrs, sp%tol)

                !write(*,'(a)') "Solving the Poisson Equation for Pressure..."
                call slv_p(P, dpdx, dpdy, dpdz, bcs(4), sp%dt, vel, &
                msh%totnds, msh%vols, msh%volnds, krnls, totquadpnts, ro, sp%cgitrs, sp%tol)
        end select

        !write(*,'(a)') "Computing Corrected Flow Velocities Vx, Vy, Vz..."
        vel(:,1) = vel(:,1) - dpdx * (sp%dt / ro)
        vel(:,2) = vel(:,2) - dpdy * (sp%dt / ro)
        vel(:,3) = vel(:,3) - dpdz * (sp%dt / ro)

        call set_var(vel(:,1), bcs(1))
        call set_var(vel(:,2), bcs(2))
        call set_var(vel(:,3), bcs(3))

        v = sqrt(vel(:,1)*vel(:,1) + vel(:,2)*vel(:,2) + vel(:,3)*vel(:,3))
            
        if (t/sp%prnt_frq == cntr) then
            cntr = cntr + 1
            write(*,'(a)') "Printing Solutions to File..."
            if (sp%vtk .eqv. .true.) then
                call prnt_vtk(v, msh, v_fname, t) 
                call prnt_vtk(p, msh, p_fname, t)
                call prnt_vec_vtk(vel, msh, vec_fname, t)
            else
                call prnt_pnts_txt(v, msh%nds, msh%totnds, v_fname, t)
                call prnt_pnts_txt(p, msh%nds, msh%totnds, p_fname, t)
            end if
        end if
    end do
end subroutine slv_ns

subroutine slv_p(p, dpdx, dpdy, dpdz, bc_p, dt, &
    vel, totnds, elems, elemnds, intrps, totquadpnts, ro, cgitr, tol)
!---------------------------------------------------------------------------------------------
! CONSTRUCTION AND SOLUTION OF GLOBAL MATRICES FOR PRESSURE
!---------------------------------------------------------------------------------------------  
    use krnl_lib
    use bc_lib
    use eqslvrs_lib
    use gm_lib

    real(8), dimension(:,:), intent(IN)                :: vel
    integer, intent(IN)                             :: totnds, elemnds, totquadpnts, cgitr
    real(8), intent(IN)                                :: dt, tol, ro
    integer, dimension(:,:), intent(IN)             :: elems
    type(kernel), dimension(:), intent(IN)          :: intrps
    type(bc), intent(IN)                            :: bc_p   
    real(8), dimension(:), allocatable, intent(INOUT)  :: P           
    real(8), dimension(:), allocatable, intent(OUT)    :: dpdx, dpdy, dpdz
    integer                                         :: i, j, g, e, m, n, cm, cn, nbr
    real(8)                                            :: dvxdx, dvydy, dvzdz, px, py, pz, stab_prmtr
    real(8), dimension(:,:), allocatable               :: gk, gm, lk
    real(8), dimension(:), allocatable                 :: gf, lf
    logical                                         :: acm, supg
    
    allocate(dpdx(totnds))
    allocate(dpdy(totnds))
    allocate(dpdz(totnds))
    allocate(gk(totnds, totnds))
    allocate(gf(totnds))
    allocate(lk(elemnds,elemnds))
    allocate(lf(elemnds))

    acm = .false.
    supg = .false.
    stab_prmtr = 0.0

    if (acm .eqv. .true.) then
        allocate(gm(totnds, totnds))
        gm(:,:) = 0
    end if
    gk(:,:) = 0.0
    gf(:) = 0.0

    do g = 1, totquadpnts
        e = intrps(g)%elem_num

        !donot change...pressure will not be calculated correctly if you modify.
        dvxdx = 0.0
        dvydy = 0.0
        dvzdz = 0.0
        do j = 1, elemnds
            dvxdx = dvxdx+intrps(g)%nx(j)*vel(elems(e,j),1)
            dvydy = dvydy+intrps(g)%ny(j)*vel(elems(e,j),2)
            dvzdz = dvzdz+intrps(g)%nz(j)*vel(elems(e,j),3)
        end do

        do j = 1, elemnds
            do i = 1, elemnds
                lk(i,j) = (1.0/ro)*(intrps(g)%nx_nx(i,j) + intrps(g)%ny_ny(i,j) + intrps(g)%nz_nz(i,j))
            end do
            lf(j) = (-1.0/dt) * (dvxdx + dvydy + dvzdz)
        end do

        do n = 1, elemnds
            cn = elems(e,n)
            do m = 1, elemnds
                cm = elems(e,m)
                gk(cm,cn) = gk(cm,cn) + intrps(g)%quad_wt * lk(m,n)
            end do
            gf(cn) = gf(cn) + intrps(g)%quad_wt * lf(n)
        end do
    end do

    call set_var(P, bc_p)
    if (acm .eqv. .true.) then
        call asmbl_mass_mat(GM, supg, totnds, &
            elems, elemnds, intrps, totquadpnts, vel, stab_prmtr)
        call set_dbc_trans(totnds, bc_p, gk, gm, gf)
        call slv_syseq(1, p, gk, gm, gf, totnds, dt, cgitr, tol)
    else
        call set_dbc_ss(1, totnds, bc_p, gk, gf)
        call slv_pcg(gk, gf, totnds, cgitr, tol, p);
    end if

    dpdx(:) = 0.0
    dpdy(:) = 0.0
    dpdz(:) = 0.0
    do i = 1, totquadpnts
        e = intrps(i)%elem_num
        !donot change...pressure will not be calculated correctly if you modify.
        px = 0.0
        py = 0.0
        pz = 0.0
        do j = 1, elemnds
            nbr = elems(e,j)
            px = px + intrps(i)%nx(j) * p(nbr)
            py = py + intrps(i)%ny(j) * p(nbr)
            pz = pz + intrps(i)%nz(j) * p(nbr)
        end do

        do n = 1, elemnds
            cn = elems(e,n)
            dpdx(cn) = dpdx(cn) + intrps(i)%quad_wt * px
            dpdy(cn) = dpdy(cn) + intrps(i)%quad_wt * py
            dpdz(cn) = dpdz(cn) + intrps(i)%quad_wt * pz
        end do
    end do
end subroutine slv_p

end module ns

