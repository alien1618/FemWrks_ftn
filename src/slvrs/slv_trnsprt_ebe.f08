module trnsprt_ebe
implicit none
contains

subroutine run_trnsprt_ebe()
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
    type(mesh)                              :: msh      !mesh data structure
    real(8), dimension(:), allocatable      :: u, q, phi     !temperature and generation term
    real(8), dimension(:,:), allocatable    :: v        !advective velocity
    type(bc), dimension(6)                  :: bcs      !boundary conditions data structure
    real(8)                                 :: k
    type(slvr_prmtrs)                       :: sp 
    real(8)                                 :: start, finish
!---------------------------------------------------------------------------------------------
    !assigning mesh
    call set_msh(msh)

    !assigning thermal conductivity
    call set_mat('k', k)
    
    !assigning boundary conditions
    call set_bc(msh%nds, msh%totnds, bcs)
    
    !initialize field variables
    call set_phi0(msh%nds, msh%totnds, msh%dx, phi)
    call set_u0(phi, msh%totnds, u)
    call set_v0(msh%nds, msh%totnds, v)
    allocate(q(msh%totnds))
    q(:) = 0
    
    !set solver parameters
    call set_ctrls(sp)
    
    !run solver
    call cpu_time(start)
    call slv_trnsprt_ebe(msh, bcs, sp, k, u, v, q)
    call cpu_time(finish)
    print '("Time = ",f8.1," seconds.")',finish-start

end subroutine run_trnsprt_ebe

subroutine slv_trnsprt_ebe(msh, bcs, sp, k, u, v, q)
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
    type(mesh), intent(in)                              :: msh
    type(bc), dimension(6), intent(in)                  :: bcs
    type(slvr_prmtrs), intent(in)                       :: sp 
    real(8), intent(in)                                 :: k
    real(8), dimension(:), allocatable, intent(inout)   :: u
    real(8), dimension(:), intent(in)                   :: q
    real(8), dimension(:,:), intent(in)                 :: v
!---------------------------------------------------------------------------------------------
    integer                                 :: t, cntr = 1  
    integer                                 :: elemquadpnts, totquadpnts
    real(8), dimension(:), allocatable      :: eps, eta, zta, w
    real(8), dimension(:), allocatable      :: gf       !global force matrix
    type(kernel), dimension(:), allocatable :: krnls    !interpolation krnls data structure
    character(len=50)                       :: fname = 'u'
    type(local), dimension(:), allocatable  :: lmat !local matrices at element level 
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
!-----------------------------------------------------------------------------------------
    !computing local matrices
        select case (msh%dim)
            case (2)
                call get_lmat(k, v, krnls, totquadpnts, msh%surfs, msh%totsurfs, msh%surfnds, &
                    sp%dt, sp%transient, sp%supg, sp%stab_prmtr, lmat)
            case (3)
                call get_lmat(k, v, krnls, totquadpnts, msh%vols, msh%totvols, msh%volnds, &
                    sp%dt, sp%transient, sp%supg, sp%stab_prmtr, lmat)
        end select       
!-----------------------------------------------------------------------------------------
    !solving the element-by-element system of equations  
        if (sp%transient == 0) then
            write(*,'(a)') "Computing steady-state solution..."

            select case (msh%dim)
            case (2)
                call asmbl_gf_ebe(u, q, msh%totnds, krnls, totquadpnts, msh%surfs, msh%surfnds, sp%dt, sp%transient, gf)
                
                call slv_pcg_ebe(1, gf, msh%totnds, msh%surfs, msh%totsurfs, msh%surfnds, lmat, bcs(1), sp%cgitrs, sp%tol, u)
            case (3)
                call asmbl_gf_ebe(u, q, msh%totnds, krnls, totquadpnts, msh%vols, msh%volnds, sp%dt, sp%transient, gf)
                
                call slv_pcg_ebe(1, gf, msh%totnds, msh%vols, msh%totvols, msh%volnds, lmat, bcs(1), sp%cgitrs, sp%tol, u)
            end select

            if (sp%vtk .eqv. .true.) then
                call prnt_vtk(u, msh, fname, 1)
            else
                call prnt_pltctrl(msh%nds, 1, 1)
                call prnt_elems_txt(msh%surfs, msh%totsurfs)
                call prnt_pnts_txt(u, msh%nds, msh%totnds, fname, 1)
            end if
        else
            cntr = 1;
            if (sp%vtk .eqv. .true.) then
                call prnt_vtk(u, msh, fname, 0)
            else
                call prnt_pltctrl(msh%nds, sp%nt, sp%prnt_frq)
                call prnt_elems_txt(msh%surfs, msh%totsurfs)
                call prnt_pnts_txt(u, msh%nds, msh%totnds, fname, 0)
            end if
            
            do t = 1, sp%nt
                write(*,'(a, i0, a, i0)') "Processing time step ", t, " of ", sp%nt

                write(*,'(a)') "Solving the system of equations..."
                select case (msh%dim)
                case (2)
                    call asmbl_gf_ebe(u, q, msh%totnds, krnls, totquadpnts, msh%surfs, msh%surfnds, sp%dt, sp%transient, gf)
                    
                    call slv_pcg_ebe(1, gf, msh%totnds, msh%surfs, msh%totsurfs, msh%surfnds, lmat, bcs(1), sp%cgitrs, sp%tol, u)
                case (3)
                    call asmbl_gf_ebe(u, q, msh%totnds, krnls, totquadpnts, msh%vols, msh%volnds, sp%dt, sp%transient, gf)
                    
                    call slv_pcg_ebe(1, gf, msh%totnds, msh%vols, msh%totvols, msh%volnds, lmat, bcs(1), sp%cgitrs, sp%tol, u)
                end select
                
                if (t/sp%prnt_frq == cntr) then
                    cntr = cntr + 1
                    if (sp%vtk .eqv. .true.) then
                        call prnt_vtk(u, msh, fname, t)
                    else
                        call prnt_pnts_txt(u, msh%nds, msh%totnds, fname, t)
                    end if
                end if
            end do
        end if
end subroutine slv_trnsprt_ebe

end module trnsprt_ebe

