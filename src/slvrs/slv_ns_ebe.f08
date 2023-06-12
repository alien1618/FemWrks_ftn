module ns_ebe
implicit none
contains

subroutine run_ns_ebe()
!---------------------------------------------------------------------------------------------
! Subroutine solves special cases of the general transport equation using
! the element-by-element finite element method without the need to form global matrices
!---------------------------------------------------------------------------------------------
    use msh_struct
    use bc_struct
    use slvr_prmtrs_struct
    use prmtrs
    implicit none
!---------------------------------------------------------------------------------------------
    type(mesh)                             :: msh
    type(bc), dimension(6)                 :: bcs      !boundary conditions data structure
    real(8), dimension(:,:), allocatable   :: vel, f
    real(8)                                :: ro, nu
    type(slvr_prmtrs)                      :: sp 
!---------------------------------------------------------------------------------------------
    !assigning mesh
    call set_msh(msh)

    !assigning material properties
    call set_mat('nu', nu)
    call set_mat('ro', ro)
    
    !assigning boundary conditions
    call set_bc(msh%nds, msh%totnds, bcs)
    
    !initialization of field variables
    write(*,'(a)') "Initializing field variables..."
    call set_v0(msh%nds, msh%totnds, vel)
    vel = vel + 1e-6 !adding a small value to avoid singular matrices
    allocate(f(msh%totnds,3))      
    f(:,:) = 0.0

    !assigning solver parameters
    call set_ctrls(sp)
    
    !running solver
    call slv_ns_ebe(msh, bcs, sp, nu, ro, vel, f)
end subroutine run_ns_ebe

subroutine slv_ns_ebe(msh, bcs, sp, nu, ro, vel, f)
!---------------------------------------------------------------------------------------------
! Subroutine solves special cases of the general transport equation using
! the element-by-element finite element method without the need to form global matrices
!---------------------------------------------------------------------------------------------
    use msh_struct
    use bc_struct
    use krnl_struct
    use nbr_struct
    use lmat_struct
    use slvr_prmtrs_struct
    use msh_ops
    use quad_ops
    use krnl_ops
    use gm_ops
    use bc_ops
    use matfree_ops
    use eq_slvrs
    implicit none
!---------------------------------------------------------------------------------------------
    type(mesh), intent(in)                       :: msh
    type(bc), dimension(6), intent(in)           :: bcs      !boundary conditions data structure
    type(slvr_prmtrs), intent(in)                :: sp 
    real(8), intent(in)                          :: ro, nu
    real(8), dimension(:,:), intent(in)          :: f
    real(8), dimension(:,:), intent(inout)       :: vel
!---------------------------------------------------------------------------------------------
    integer                                      :: t, cntr = 1
    integer                                      :: elemquadpnts, totquadpnts
    real(8), dimension(:), allocatable           :: eps, eta, zta, w, gf, v
    type(kernel), dimension(:), allocatable      :: krnls
    character(len=50)                            :: v_fname = 'v', p_fname = 'p', vec_fname = 'vec'
    real(8), dimension(:), allocatable           :: p, dpdx, dpdy, dpdz
    type(local), dimension(:), allocatable       :: lmat !local matrices at element level 
    real(8), dimension(:), allocatable           :: vx, vy, vz
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
    !solving element-by-element without creating global matrices    
        do t = 1, sp%nt
            write(*,'(a,i0,a,i0)')"Processing ebe time step ", t ," of ", sp%nt
            select case (msh%dim)
                case (2)
                    call get_lmat(nu, vel, krnls, totquadpnts, msh%surfs,&
                        msh%totsurfs, msh%surfnds, sp%dt, 1, sp%supg, sp%stab_prmtr, lmat)

                    call asmbl_gf_ebe(vel(:,1), f(:,1), msh%totnds, krnls, totquadpnts, msh%surfs, msh%surfnds, sp%dt, 1, gf)

                    call slv_pcg_ebe(1, gf, msh%totnds, msh%surfs, msh%totsurfs, msh%surfnds, lmat, &
                        bcs(1), sp%cgitrs, sp%tol, vx)
                    
                    call asmbl_gf_ebe(vel(:,2), f(:,2), msh%totnds, krnls, totquadpnts, msh%surfs, msh%surfnds, sp%dt, 1, gf)
                    
                    call slv_pcg_ebe(1, gf, msh%totnds, msh%surfs, msh%totsurfs, msh%surfnds, lmat, &
                        bcs(2), sp%cgitrs, sp%tol, vy)
                    
                    vel(:,1) = vx
                    vel(:,2) = vy     
                    call slv_p_ebe(p, dpdx, dpdy, dpdz, vel, msh%totnds, krnls, totquadpnts, &
                        msh%surfs, msh%totsurfs, msh%surfnds, ro, sp%dt, sp%cgitrs, sp%tol, bcs(4))
                    
                case (3)
                    call get_lmat(nu, vel, krnls, totquadpnts, msh%vols,&
                        msh%totvols, msh%volnds, sp%dt, 1, sp%supg, sp%stab_prmtr, lmat)

                    call asmbl_gf_ebe(vel(:,1), f(:,1), msh%totnds, krnls, totquadpnts, msh%vols, msh%volnds, sp%dt, 1, gf)
                    
                    call slv_pcg_ebe(1, gf, msh%totnds, msh%vols, msh%totvols, msh%volnds, lmat, &
                        bcs(1), sp%cgitrs, sp%tol, vx)                   
                  
                    call asmbl_gf_ebe(vel(:,2), f(:,2), msh%totnds, krnls, totquadpnts, msh%vols, msh%volnds, sp%dt, 1, gf)
                    
                    call slv_pcg_ebe(1, gf, msh%totnds, msh%vols, msh%totvols, msh%volnds, lmat, &
                        bcs(2), sp%cgitrs, sp%tol, vy) 
                    
                    call asmbl_gf_ebe(vel(:,3), f(:,3), msh%totnds, krnls, totquadpnts, msh%vols, msh%volnds, sp%dt, 1, gf)
                    
                    call slv_pcg_ebe(1, gf, msh%totnds, msh%vols, msh%totvols, msh%volnds, lmat, &
                        bcs(3), sp%cgitrs, sp%tol, vz) 
                    
                    vel(:,1) = vx
                    vel(:,2) = vy     
                    vel(:,3) = vz

                    call slv_p_ebe(p, dpdx, dpdy, dpdz, vel, msh%totnds, krnls, totquadpnts, &
                         msh%vols, msh%totvols, msh%volnds, ro, sp%dt, sp%cgitrs, sp%tol, bcs(4))
                    
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
end subroutine slv_ns_ebe

subroutine slv_p_ebe(p, dpdx, dpdy, dpdz, v, totnds, krnls, totquadpnts, elems, totelems, elemnds, &
                    ro, dt, cgitr, conv_tol, bcs)
!--------------------------------------------------------------------------------------------
! Function solves for fluid flow pressure using the element-by-element finite elements
! without the need to construct global matrices
!--------------------------------------------------------------------------------------------
    use krnl_struct
    use lmat_struct
    use nbr_struct
    use bc_struct
    use eq_slvrs
    use bc_ops
    use matfree_ops

    real(8), dimension(:), allocatable, intent(INOUT) :: p
    real(8), dimension(:,:), intent(IN)               :: v                
    type(bc), intent(IN)                              :: bcs
    real(8), intent(IN)                               :: ro, dt, conv_tol
    integer, intent(IN)                               :: totnds, totquadpnts, totelems, elemnds
    integer, intent(IN)                               :: cgitr
    integer, dimension(:,:), intent(IN)               :: elems
    type(kernel), dimension(:), intent(IN)            :: krnls
    integer                                           :: i, j, g, cj, cm, e, mm, n
    real(8), dimension(:), allocatable                :: gf !global force matrix
    type(local), dimension(:), allocatable            :: lmat
    real(8), dimension(:), allocatable, intent(OUT)   :: dpdx, dpdy, dpdz
    real(8)                                           :: px, py, pz, dvxdx, dvydy, dvzdz

    allocate(dpdx(totnds))
    allocate(dpdy(totnds))
    allocate(dpdz(totnds))
    allocate(gf(totnds))
    gf(:) = 0

    allocate(lmat(totelems))
    do i = 1, totelems
        allocate(lmat(i)%mat(elemnds, elemnds))
        lmat(i)%mat(:,:) = 0.0
    end do

    do g = 1, totquadpnts
        !----------------------------------------------------------
        !interpolate material properties at quad point
        !----------------------------------------------------------
        dvxdx = 0;
        dvydy = 0;
        dvzdz = 0;
        e = krnls(g)%elem_num
        do mm = 1, elemnds
            cm = elems(e,mm)
            dvxdx = dvxdx + krnls(g)%nx(mm)*v(cm,1)
            dvydy = dvydy + krnls(g)%ny(mm)*v(cm,2)
            dvzdz = dvzdz + krnls(g)%nz(mm)*v(cm,3)
        end do

        !----------------------------------------------------------
        ! cout << "get quad STIFFNESS AND FORCE MATRICES" << endl;
        !----------------------------------------------------------
        !CALCULATES IMPLICIT FORM: 0 = (K*U + Q)
            
        !First: calculate K
        lmat(e)%mat = lmat(e)%mat + krnls(g)%quad_wt * (krnls(g)%nx_nx+krnls(g)%ny_ny+krnls(g)%nz_nz)           

        !Second: calculate (Q)
        do j = 1, elemnds
            cj = elems(e,j)
            gf(cj) = gf(cj) + krnls(g)%quad_wt*((-ro/dt)*(dvxdx+dvydy+dvzdz))
        end do
    end do

    call slv_pcg_ebe(1, gf, totnds, elems, totelems, elemnds, lmat, bcs, cgitr, conv_tol, p)
    call set_var(p, bcs)
    dpdx(:) = 0.0
    dpdy(:) = 0.0
    dpdz(:) = 0.0
    do i = 1, totquadpnts
        e = krnls(i)%elem_num
        !donot change...pressure will not be calculated correctly if you modify.
        px = 0.0
        py = 0.0
        pz = 0.0
        do j = 1, elemnds
            cm = elems(e,j)
            px = px + krnls(i)%nx(j) * p(cm)
            py = py + krnls(i)%ny(j) * p(cm)
            pz = pz + krnls(i)%nz(j) * p(cm)
        end do

        do n = 1, elemnds
            cm = elems(e,n)
            dpdx(cm) = dpdx(cm) + krnls(i)%quad_wt * px
            dpdy(cm) = dpdy(cm) + krnls(i)%quad_wt * py
            dpdz(cm) = dpdz(cm) + krnls(i)%quad_wt * pz
        end do
    end do

end subroutine slv_p_ebe

end module ns_ebe

