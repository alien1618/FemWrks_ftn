module elst_ebe
implicit none
contains

subroutine run_elst_ebe()
!---------------------------------------------------------------------------------------------
! Subroutine solves linear elasticity using the element-by-element finite element method
! without the need to form global matrices
!---------------------------------------------------------------------------------------------
    use msh_struct
    use bc_struct
    use mat_struct
    use slvr_prmtrs_struct
    use prmtrs
!---------------------------------------------------------------------------------------------
    implicit none
    type(mesh)               :: msh
    type(material)           :: mat
    type(bc), dimension(6)   :: bcs
    type(slvr_prmtrs)        :: sp 
    real(8)                  :: start, finish
!---------------------------------------------------------------------------------------------
    !assigning mesh
    call set_msh(msh)

    !assigning material properties
    call set_elst(mat)
    
    !assigning boundary conditions
    call set_bc(msh%nds, msh%totnds, bcs)
    
    !assigning solver parameters
    call set_ctrls(sp)
    
    !running solver
    call cpu_time(start)
    call slv_elst_ebe(msh, mat, bcs, sp)
    call cpu_time(finish)
    print '("Time = ",f8.1," seconds.")',finish-start
end subroutine run_elst_ebe

subroutine slv_elst_ebe(msh, mat, bcs, sp)
!---------------------------------------------------------------------------------------------
! Subroutine solves linear elasticity using the element-by-element finite element method
! without the need to form global matrices
!---------------------------------------------------------------------------------------------
    use msh_struct
    use bc_struct
    use krnl_struct
    use nbr_struct
    use lmat_struct
    use mat_struct
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
    type(mesh), intent(inout)               :: msh
    type(material), intent(in)              :: mat
    type(bc), dimension(6), intent(in)      :: bcs
    type(slvr_prmtrs), intent(in)           :: sp
!---------------------------------------------------------------------------------------------
    integer                                 :: itr, j
    integer                                 :: ndnum
    integer                                 :: elemquadpnts, totquadpnts, dof
    real(8), dimension(:), allocatable      :: eps, eta, zta, w, u, gf, deformation, vms
    type(kernel), dimension(:), allocatable :: krnls
    character(len=50)                       :: u_fname = 'u', vms_fname = 'vms'
    type(local), dimension(:), allocatable  :: lmat !local matrices at element level 
!---------------------------------------------------------------------------------------------
    !generating quadrature points
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
!---------------------------------------------------------------------------------------------
    !calculating the force matrix
    dof = msh%dim
    allocate(gf(dof*msh%totnds))
    gf(:) = 0.0
    select case(msh%dim)
        case (2)
            write(*,'(a)') "Applying boundary conditions..."
            itr = 1;
            do j = 1, bcs(2)%tot
                ndnum = bcs(2)%pnts(j)
                gf(dof*ndnum-itr) = bcs(2)%vals(j)
            end do
            do j = 1, bcs(3)%tot
                ndnum = bcs(3)%pnts(j)
                gf(dof*ndnum) = bcs(3)%vals(j)
            end do

        case (3)
            write(*,'(a)') "Applying boundary conditions..."
            itr = 2
            do j = 1, bcs(2)%tot
                ndnum = bcs(2)%pnts(j)
                gf(dof*ndnum-itr) = bcs(2)%vals(j)
            end do
            itr = 1
            do j = 1, bcs(3)%tot
                ndnum = bcs(3)%pnts(j)
                gf(dof*ndnum-itr) = bcs(3)%vals(j)
            end do
            do j = 1, bcs(4)%tot
                ndnum = bcs(4)%pnts(j)
                gf(dof*ndnum) = bcs(4)%vals(j)
            end do
    end select
!---------------------------------------------------------------------------------------------
    !element-by-element solver without forming global matrices  
    write(*,'(a)') "Computing steady-state solution..."
    select case (msh%dim)
    case (2)
        call get_lmat_elst(msh%dim, mat%plane_type, mat%moe, mat%nu, mat%thickness, krnls, &
                totquadpnts, msh%surfs, msh%totsurfs, msh%surfnds, lmat)
         
        call slv_pcg_ebe(dof, gf, msh%totnds, msh%surfs, msh%totsurfs, msh%surfnds, lmat, &
                bcs(1), sp%cgitrs, sp%tol, u)
    case (3)
        call get_lmat_elst(msh%dim, mat%plane_type, mat%moe, mat%nu, mat%thickness, krnls, &
                totquadpnts, msh%vols, msh%totvols, msh%volnds, lmat)
         
        call slv_pcg_ebe(dof, gf, msh%totnds, msh%vols, msh%totvols, msh%volnds, lmat, &
                bcs(1), sp%cgitrs, sp%tol, u)
    end select
!---------------------------------------------------------------------------------------------
    !calculating deformation
    write(*,'(a)') "Computing deformation..."
    call get_dfrm(u, dof, msh%totnds, msh%nds, deformation)
!---------------------------------------------------------------------------------------------
    !calculating strains and stresses
    write(*,'(a)') "Computing strains and stresses..."
    call get_vms(dof, mat, msh, krnls, totquadpnts, u, vms)
!---------------------------------------------------------------------------------------------
    !printing data to file
    if (sp%vtk .eqv. .true.) then
        call prnt_vtk(deformation, msh, u_fname, 1)  
        call prnt_vtk(vms, msh, vms_fname, 1)  
    else
        call prnt_pltctrl(msh%nds, 1, 1)
        call prnt_elems_txt(msh%surfs, msh%totsurfs)
        call prnt_pnts_txt(deformation, msh%nds, msh%totnds, u_fname, 1)
        call prnt_pnts_txt(vms, msh%nds, msh%totnds, vms_fname, 1)
    end if
end subroutine slv_elst_ebe

end module elst_ebe

