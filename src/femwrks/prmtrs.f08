module slvr_prmtrs_struct
!-----------------------------------------------------------------------
! defining solver parameters data structure
!-----------------------------------------------------------------------
implicit none
private
public slvr_prmtrs

type slvr_prmtrs
    integer :: slvr, transient, nt, cgitrs, prnt_frq
    real(8) :: dt, tol, stab_prmtr
    logical :: vtk, supg
end type

end module slvr_prmtrs_struct

module mat_struct
!-----------------------------------------------------------------------
! defining material properties data structure
!-----------------------------------------------------------------------
implicit none
private
public material

type material
    real(8) :: moe, nu, thickness
    integer :: plane_type
end type

end module mat_struct

module prmtrs
implicit none
contains

subroutine read_array(fname, data, file_exists)
    use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none
!-----------------------------------------------------------------------
    character*50, intent(in) :: fname
!-----------------------------------------------------------------------
    logical, intent(out) :: file_exists
    real(8), dimension(:), allocatable, intent(out) :: data
!-----------------------------------------------------------------------
    integer :: error
    real(8) :: x
!-----------------------------------------------------------------------

    allocate(data(0))
    inquire(file = fname, exist = file_exists)
    if (file_exists .eqv. .true.) then
        open (1, file = fname, status = 'old')
        do
            read(1, *, iostat = error) x
            select case(error)
            case (0)
                data = [data , x]
            case(iostat_end)
                exit
            case Default
                write(*,*) 'Error in reading file'
                stop
            end select
        end Do
        close(1)
    end if
end subroutine read_array

subroutine set_msh(msh)
    use msh_struct
    use msh_ops
    implicit none
!-----------------------------------------------------------------------
    type(mesh), intent(out) :: msh
!-----------------------------------------------------------------------
    integer  :: i
    real(8), dimension(3) :: trans, b, l
    integer, dimension(3) :: n
    real(8) :: scale
    character(len=50) :: pnts_fname, surf_fname, vol_fname
    character(len=50) :: fname1, fname2, fname3
    real(8), dimension(:), allocatable :: data
    logical :: msh_exists, grd_exists, trns_exists
!-----------------------------------------------------------------------
    fname1 = 'sim/in/1_msh/nds.txt'
    fname2 = 'sim/in/1_msh/grd.txt'
    fname3 = 'sim/in/1_msh/trns.txt'

    inquire(file = fname1, exist = msh_exists)
    call read_array(fname2, data, grd_exists)
    if (msh_exists .eqv. .true.) then
        pnts_fname = 'sim/in/1_msh/nds.txt'
        surf_fname = 'sim/in/1_msh/surfs.txt'   
        vol_fname = 'sim/in/1_msh/vols.txt'           
        call read_pnts(pnts_fname, msh)  
        call read_surfs(surf_fname, msh)
        if (msh%dim == 3) then
            call read_vols(vol_fname, msh)
        else
            msh%totvols = 0
            msh%volnds = 0
            msh%volshape = 0
        end if
    else if (grd_exists .eqv. .true.) then
        b = data(1:3)          !base point
        l = data(4:6)          !domain length
        n = int(data(7:9))    !domain discritization
        if (n(3) <= 1) then
            call gen_quadgrid(b, l, n, msh)
        else
            call gen_hexgrid(b, l, n, msh)
        end if
    else
        write(*,*) "ERROR: no msh.txt or grd.txt exist in sim/in/1_msh/"
        write(*,*) "Solver can not proceed without either file. Exiting..."
        call exit()
    end if

    call read_array(fname3, data, trns_exists)
    if (trns_exists .eqv. .true.) then
        scale = data(1)       !geometry scale factor
        trans = data(2:4)    !geometry translation
         msh%nds =  msh%nds * scale
        do i = 1, 3
             msh%nds(:, i) =  msh%nds(:, i) + trans(i)
        end do
    end if
    call get_avg_dx(msh%nds, msh%dx)
end subroutine set_msh

subroutine set_slvr(slvr)
!-----------------------------------------------------------------------
    use slvr_prmtrs_struct
    implicit none
!-----------------------------------------------------------------------
    integer, intent(out) :: slvr
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable      :: data
    character*50 :: fname = 'sim/in/5_ctrl/slvr.txt'
    logical :: file_exists
!-----------------------------------------------------------------------
    call read_array(fname, data, file_exists)

    if (file_exists .eqv. .true.) then
        slvr = int(data(1))
    else
        write(*,'(a)') "ERROR: sim/in/5_ctrl/slvr.txt does not exist. Can not proceed without this file. Exiting..."
        call exit()
    end if
end subroutine set_slvr

subroutine set_ctrls(sp)
!-----------------------------------------------------------------------
    use slvr_prmtrs_struct
    implicit none
!-----------------------------------------------------------------------
    type(slvr_prmtrs), intent(inout)        :: sp   !solver parameters data structure
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable      :: data
    character*50 :: fname = 'sim/in/5_ctrl/ctrls.txt'
    logical :: file_exists
!-----------------------------------------------------------------------
    call read_array(fname, data, file_exists)

    if (file_exists .eqv. .true.) then
        sp%transient = int(data(1))
        sp%dt = data(2)          !time step size
        sp%nt = int(data(3))           !total time steps
        sp%cgitrs = int(data(4))             !iterations
        sp%tol = data(5)                    !convergence tolerance
        sp%prnt_frq = int(data(6))       !print frequency
        sp%vtk = .false.         !output data format vtk or txt
        if (int(data(7)) == 1) sp%vtk = .true.
        sp%supg = .false.
        if (int(data(8)) == 1) sp%supg = .true.
        sp%stab_prmtr = data(9)
    else
        write(*,'(a)') "ERROR: sim/in/5_ctrl/ctrls.txt does not exist. Can not proceed without this file. Exiting..."
        call exit()
    end if
end subroutine set_ctrls

subroutine set_elst(mat)

    use slvr_prmtrs_struct
    use mat_struct
    implicit none
!-----------------------------------------------------------------------
    type(material), intent(out) :: mat   !solver parameters data structure
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable      :: data
    character*50 :: fname = 'sim/in/2_mat/elst.txt'
    logical :: file_exists
!-----------------------------------------------------------------------
    call read_array(fname, data, file_exists)

    if (file_exists .eqv. .true.) then
        mat%moe = data(1)    !modulous of elasticity
        mat%nu = data(2)      !poisson ratio
        mat%thickness = data(3)           !material thickness
        mat%plane_type = int(data(4))             !plane stress or plane strain
    else
        write(*,'(a)') "ERROR: sim/in/2_mat/elst.txt does not exist. Can not proceed without this file. Exiting..."
        call exit()
    end if
end subroutine set_elst

subroutine set_bc(pnts, totpnts, bcs)
    use bc_struct
    use bc_ops
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: totpnts
!-----------------------------------------------------------------------
    type(bc), dimension(6), intent(out) :: bcs
!-----------------------------------------------------------------------
    real(8), dimension(:,:), allocatable :: data
    real(8) :: bcval, bcloc
    integer :: i, totbcs, bcaxis, bctype
    logical :: check
    character*100, dimension(:,:), allocatable :: fname
    character*100 :: file_name
    character*50 :: path='sim/in/4_bc/'
!-----------------------------------------------------------------------
    !bc types: 1 = u, 2 = phi, 3 = vx, 4 = vy, 5 = vz, 6 = p
    do i = 1, 6
        bcs(i)%tot = 0
    end do
    inquire( file = 'sim/in/4_bc/bc.txt', exist = check)
    if (check .eqv. .true.) then
        open(unit = 1, file = 'sim/in/4_bc/bc.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(1, *) totbcs
        allocate(data(totbcs, 4))
        do i = 1, totbcs
            read(1,*) data(i,:)
        end do
        close (1)
        
        !maximum number of bc types possible is 6 (vx, vy, vz, p, u, phi) or (ux, uy, uz, fx, fy, fz)
        do i = 1, 6
            bcs(i)%tot = 0
            allocate(bcs(i)%pnts(0))
            allocate(bcs(i)%vals(0))
        end do

        do i = 1, totbcs
            bctype = int(data(i, 1))
            bcaxis = int(data(i, 2))
            bcloc = data(i, 3)
            bcval = data(i, 4)
            call set_dbc(pnts, totpnts, bcaxis, bcloc, bcval, bcs(bctype))
        end do
        close (1)
    end if
    inquire(file = 'sim/in/4_bc/bc_f.txt', exist = check)
    if (check .eqv. .true.) then
        open(unit = 2, file = 'sim/in/4_bc/bc_f.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(2, *) totbcs                   
        allocate(fname(totbcs, 3))
        do i = 1, totbcs
            read(2,*) fname(i,:)
        end do

        do i = 1, 6 !total possible bc types is 6 (vx, vy, vz, p, u, phi) or (ux, uy, uz, fx, fy, fz)
            bcs(i)%tot = 0
            allocate(bcs(i)%pnts(0))
            allocate(bcs(i)%vals(0))
        end do

        do i = 1, totbcs
            read(fname(i,1),*)  bctype
            read(fname(i,3), '(f10.0)' )  bcval 
            file_name = trim(adjustl(path)) // trim(adjustl(fname(i,2)))
            call read_dbc(file_name, bcval, pnts, totpnts, bcs(bctype))
        end do
    end if

end subroutine set_bc

subroutine set_u0(phi, totpnts, u)

    use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none
!-----------------------------------------------------------------------
    real(8), dimension(:), intent(in) :: phi
    integer, intent(in) :: totpnts
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(out) :: u
!-----------------------------------------------------------------------
    integer                                 :: i
    real(8)                                 :: u1, u2
    character*50 :: fname = 'sim/in/3_init/u0.txt'
    real(8), dimension(:), allocatable      :: data
    logical :: file_exists
!-----------------------------------------------------------------------
    call read_array(fname, data, file_exists)

    allocate(u(totpnts))
    if (file_exists .eqv. .true.) then
        if (size(data) > 1)  then
            u1 = data(1)
            u2 = data(2)
            u(:) = 0.0
            do i = 1, totpnts
                if (phi(i) <= 0.5) then
                    u(i) = u1
                else 
                    u(i) = u2
                end if
            end do
        else
            u(:) = data(1)
        end if
    else
        write(*,'(a)') "sim/in/3_init/u0.txt does not exist. Setting u(:) = 0.0"
        u(:) = 0.0
    end if
end subroutine set_u0

subroutine set_v0(pnts, totpnts, v)
    implicit none
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: totpnts
!-----------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: v
!-----------------------------------------------------------------------
    integer :: i, adv_v
    real(8) :: vx1, vy1, vz1, damp
    character*50 :: fname = 'sim/in/3_init/v0.txt'
    real(8), dimension(:), allocatable      :: data
    logical :: file_exists
!-----------------------------------------------------------------------
    call read_array(fname, data, file_exists)

    allocate(v(totpnts, 3))
    v(:,:) = 0.0   

    if (file_exists .eqv. .true.) then
        adv_v = int(data(1))
        damp = data(2)
        vx1 = data(3)
        vy1 = data(4)
        vz1 = data(5)

        select case (adv_v)
        case (1)
            do i = 1, totpnts
                v(i,1) = vx1
                v(i,2) = vy1
                v(i,3) = vz1
            end do
        case (2)
            call set_vrtx_vel(v, pnts, totpnts)
        case (3)
            call set_shearing_vel(v, pnts, totpnts)
            v = damp * v
        case (4)
            call set_rot_vel(v, pnts, totpnts, vx1, vy1)
        end select
    else
        write(*,'(a)') "sim/in/3_init/v0.txt does not exist. Setting v(:,:) = 0.0"
    end if
end subroutine set_v0

subroutine set_mat(matname, k)
    use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none
!-----------------------------------------------------------------------
    character(len=*),  intent(in) :: matname
!-----------------------------------------------------------------------
    real(8), intent(out) :: k
!-----------------------------------------------------------------------
    character(len=50) :: fname
    real(8) :: x       
    logical :: file_exists  
    integer :: error
!-----------------------------------------------------------------------
    fname = 'sim/in/2_mat/' // trim(adjustl(matname)) // '.txt'

    inquire(file = fname, exist = file_exists)
    if (file_exists .eqv. .true.) then
        open (1, file = fname, status = 'old')
        do
            read(1, *, iostat = error) x
            select case(error)
            case (0)
                k = x
            case(iostat_end)
                exit
            case Default
                write(*,*) 'Error in reading file'
                stop
            end select
        end Do
        close(1)
    end if
end subroutine set_mat

subroutine set_phi0(pnts, totpnts, dx, phi)
    implicit none
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: totpnts
    real(8), intent(in) :: dx
!-----------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(out) :: phi
!-----------------------------------------------------------------------
    real(8), dimension(:,:), allocatable :: qb, sphr, cyl, cosf
    real(8), dimension(:), allocatable :: rand
    real(8) :: r, a0, b0, c0, d0, e0, noise, num !random noise parameters
    real(8), dimension(3) :: c, sz
    real(8) :: sgn, axis, val
    integer :: k, totseeds
    logical :: c1, c2, c3, c4, c5
    character(len=50) :: fname
!-----------------------------------------------------------------------

    allocate(phi(totpnts))
    phi(:) = 1.0

    inquire( file = 'sim/in/3_init/phi0_pln.txt', exist = c1)
    if (c1 .eqv. .true.) then
        open(unit = 1, file = 'sim/in/3_init/phi0_pln.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(1, *) sgn
        read(1, *) axis
        read(1, *) val
        close (1)
        phi = sgn*(pnts(:,int(axis)) + val)
    end if

    inquire( file = 'sim/in/3_init/phi0_sphr.txt', exist = c1)
    if (c1 .eqv. .true.) then
        open(unit = 1, file = 'sim/in/3_init/phi0_sphr.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(1, *) totseeds
        allocate(sphr(totseeds, 5))
        do k = 1, totseeds
            read(1,*) sphr(k,:)
            c = sphr(k, 2:4)
            r = sphr(k, 5)
            if (int(sphr(k,1)) == 1) then
                call add_sphr(c, r, pnts, totpnts, phi)
            else
                call sub_sphr(c, r, pnts, totpnts, phi)
            end if
        end do
        close (1)
    end if

    inquire( file = 'sim/in/3_init/phi0_cyl.txt', exist = c2)
    if (c2 .eqv. .true.) then
        open(unit = 1, file = 'sim/in/3_init/phi0_cyl.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(1, *) totseeds
        allocate(cyl(totseeds, 5))
        do k = 1, totseeds
            read(1,*) cyl(k,:)
            c = cyl(k, 2:4)
            r = cyl(k, 5)
            if (int(cyl(k,1)) == 1) then
                call add_cyl_z(c, r, pnts, totpnts, phi)
            else
                call sub_cyl_z(c, r, pnts, totpnts, phi)
            end if
        end do
        close (1)
    end if

    inquire( file = 'sim/in/3_init/phi0_qb.txt', exist = c3)
    if (c3 .eqv. .true.) then
        open(unit = 2, file = 'sim/in/3_init/phi0_qb.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(2, *) totseeds
        allocate(qb(totseeds, 7))
        do k = 1, totseeds
            read(2,*) qb(k,:)
            c = qb(k, 2:4)
            sz = qb(k, 5:7)
            if (int(qb(k,1)) == 1) then
                call add_cube(c, sz, pnts, totpnts, phi)
            else
                call sub_cube(c, sz, pnts, totpnts, phi)
            end if
        end do
        close (2)
    end if

    inquire( file = 'sim/in/3_init/phi0_cos.txt', exist = c4)
    if (c4 .eqv. .true.) then
        open(unit = 2, file = 'sim/in/3_init/phi0_cos.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(2, *) totseeds
        allocate(cosf(totseeds, 5))
        do k = 1, totseeds
            read(2,*) cosf(k,:)
            a0 = cosf(k, 1)
            b0 = cosf(k, 2)
            c0 = cosf(k, 3)
            d0 = cosf(k, 4)
            e0 = cosf(k, 5)
            if (k == 1) then
                !phi = ((ps%pnts(:,2) - (0.006 * cos(25.0 * ps%pnts(:,1) + 0.1) + 0.0)) - 0.5)
                phi = ((pnts(:,2) - (a0 * cos(b0 * pnts(:,1) + c0) + d0)) - e0)
            else
                !phi = phi - (0.006 * cos(25.0 * ps%pnts(:,3) + 0.1) + 0.0)
                phi = phi - (a0 * cos(b0 * pnts(:,3) + c0) + d0) - e0
            end if
        end do
        close (2)
    end if
    !if (c1 .eqv. .true. .or. c2 .eqv. .true. .or. c3 .eqv. .true. .or. c4 .eqv. .true.) then
        call get_vof(phi, totpnts, dx)
    !end if

    fname = 'sim/in/3_init/phi0_rand.txt'
    call read_array(fname, rand, c5)
    if (c5 .eqv. .true.) then
        c0 = rand(1)   
        noise = rand(2)      
        do k = 1, totpnts
            call random_number(num)
            phi(k) = c0 + noise * (0.5 - num)
        end do
    end if

end subroutine set_phi0

pure subroutine get_vof(phi, totpnts, w)
!-----------------------------------------------------------
! convert a signed distance function to a heaveside function
! from 0 to 1
!-----------------------------------------------------------
    real(8), intent(in)                    ::  w
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i
!-----------------------------------------------------------

    do i = 1, totpnts
        phi(i) = 0.5 * (1 - tanh(phi(i) / w))
        if (phi(i) > 1) phi(i) = 1
        if (phi(i) < 0) phi(i) = 0
    end do
end subroutine get_vof

pure subroutine add_sphr(c, a, pnts, totpnts, phi)
!-----------------------------------------------------------
!  add an implicit sphere to phi
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)    ::  pnts
    real(8), dimension(3), intent(in)      ::  c
    real(8), intent(in)                    ::  a
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i
!-----------------------------------------------------------
    do i = 1, totpnts
        phi(i) = min(phi(i), ((pnts(i, 1) - c(1))**2) + ((pnts(i, 2) - c(2))**2) &
                     + ((pnts(i, 3) - c(3))**2) - (a**2));
    end do
end subroutine add_sphr

pure subroutine sub_sphr(c, a, pnts, totpnts, phi)
!-----------------------------------------------------------
!  add an implicit sphere to phi
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)    ::  pnts
    real(8), dimension(3), intent(in)      ::  c
    real(8), intent(in)                    ::  a
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i
!-----------------------------------------------------------
    do i = 1, totpnts
        phi(i) = max(phi(i), ((pnts(i, 1) - c(1))**2) + ((pnts(i, 2) - c(2))**2) &
                     + ((pnts(i, 3) - c(3))**2) - (a**2));
    end do
end subroutine sub_sphr

pure subroutine add_cyl_z(c, a, pnts, totpnts, phi)
!-----------------------------------------------------------
!  add an implicit cylinder to phi along the z axis
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)    ::  pnts
    real(8), dimension(3), intent(in)      ::  c
    real(8), intent(in)                    ::  a
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i   
!-----------------------------------------------------------
    do i = 1, totpnts
        phi(i) = min(phi(i), ((pnts(i, 1) - c(1))**2) + ((pnts(i, 2) - c(2))**2) - (a**2))
    end do
end subroutine add_cyl_z

pure subroutine sub_cyl_z(c, a, pnts, totpnts, phi)
!-----------------------------------------------------------
!  add an implicit cylinder to phi along the z axis
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)    ::  pnts
    real(8), dimension(3), intent(in)      ::  c
    real(8), intent(in)                    ::  a
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i   
!-----------------------------------------------------------
    do i = 1, totpnts
        phi(i) = max(phi(i), ((pnts(i, 1) - c(1))**2) + ((pnts(i, 2) - c(2))**2) - (a**2))
    end do
end subroutine sub_cyl_z

pure subroutine add_cube(c, l, pnts, totpnts, phi)
!-----------------------------------------------------------
!  add an implicit cube to phi
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)    ::  pnts
    real(8), dimension(3), intent(in)      ::  c, l
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i  
    real(8)                                ::  x, y, z, width, height, depth 
!-----------------------------------------------------------

    x = c(1) + l(1) / 2.0
    y = c(2) + l(2) / 2.0
    z = c(3) + l(3) / 2.0
    width = l(1) / 2.0
    height = l(2) / 2.0
    depth = 1.0
    if (l(3) > 0.0) depth = l(3) / 2.0

    do i = 1, totpnts
        phi(i) = min(phi(i), (((pnts(i, 1) - x) / width)**100) + (((pnts(i, 2) - y) / height)**100) &
                    + (((pnts(i,3) - z) / depth)**6) - 1.0)
    end do
end subroutine add_cube

pure subroutine sub_cube(c, l, pnts, totpnts, phi)
!-----------------------------------------------------------
!  subtract an implicit cube from phi
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)    ::  pnts
    real(8), dimension(3), intent(in)      ::  c, l
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i  
    real(8)                                ::  x, y, z, width, height, depth
!-----------------------------------------------------------

    x = c(1) + l(1) / 2.0
    y = c(2) + l(2) / 2.0
    z = c(3) + l(3) / 2.0

    width = l(1) / 2.0
    height = l(2) / 2.0
    depth = 1.0
    if (l(3) > 0.0) depth = l(3) / 2.0

    do i = 1, totpnts
        phi(i) = max(phi(i), -((((pnts(i, 1) - x) / width)**100) + (((pnts(i, 2)-y) / height)**100) & 
                    + (((pnts(i, 3) - z) / depth)**6) - 1.0))
    end do
end subroutine sub_cube


pure subroutine set_vrtx_vel(vel, pnts, totpnts)
!--------------------------------------------------------------------
! assign vrtx flow field distribution
!--------------------------------------------------------------------
    real(8), dimension(:,:), intent(in)         ::  pnts
    integer, intent(in)                         ::  totpnts
!--------------------------------------------------------------------
    real(8), dimension(totpnts,3), intent(out)  ::  vel
!--------------------------------------------------------------------
    integer                                     ::  i
!--------------------------------------------------------------------

    do i = 1, totpnts
        vel(i, 1) = -(sin(3.14 * pnts(i, 1))**2) * sin(2 * 3.14 * pnts(i, 2))
        vel(i, 2) =  (sin(3.14 * pnts(i, 2))**2) * sin(2 * 3.14 * pnts(i, 1))
    end do
end subroutine set_vrtx_vel

pure subroutine set_shearing_vel(vel, pnts, totpnts)
!-----------------------------------------------------------
! assign shearing flow field distribution
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)                 ::  pnts
    integer, intent(in)                                 ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(totpnts,3), intent(out)   ::  vel
!-----------------------------------------------------------
    integer                                             ::  i
!-----------------------------------------------------------

    do i = 1, totpnts
        vel(i, 1) = sin(4 * 3.14 * (pnts(i, 1) + 0.5)) * sin(4 * 3.14 * (pnts(i, 2) + 0.5))
        vel(i, 2) = cos(4 * 3.14 * (pnts(i, 1) + 0.5)) * cos(4 * 3.14 * (pnts(i, 2) + 0.5))
    end do
end subroutine set_shearing_vel

pure subroutine set_rot_vel(vel, pnts, totpnts, u, v)
!-----------------------------------------------------------
! assign rotating flow field distribution
!-----------------------------------------------------------
    real(8), intent(in)                         :: u, v
    real(8), dimension(:,:), intent(in)         ::  pnts
    integer, intent(in)                         ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(totpnts,3), intent(out)  ::  vel
!-----------------------------------------------------------
    integer                                     ::  i
!-----------------------------------------------------------
    do i = 1, totpnts
        vel(i, 1) = u * pnts(i, 2)
        vel(i, 2) = -v * pnts(i, 1)
    end do
end subroutine set_rot_vel

subroutine get_avg_dx(pnts, avg_dx)
!-------------------------------------------------------------------------
    use omp_lib
!-------------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: pnts
!-------------------------------------------------------------------------
    real(8), intent(out) :: avg_dx
!-------------------------------------------------------------------------
    real(8) :: d_min, dminf = 1e6, dx, dy, dz, d
    integer :: i, j
    real(8) :: sum_d = 0
    integer :: totsamples = 50
!-------------------------------------------------------------------------
    do i=1, totsamples
        d_min = 1e6
        do j = 1, TotSamples
            dx = pnts(i, 1) - pnts(j, 1)
            dy = pnts(i, 2) - pnts(j, 2)
            dz = pnts(i, 3) - pnts(j, 3)
            d = sqrt(dx * dx + dy * dy + dz * dz)
            if (d <= d_min .and. d > 1.0e-6) d_min = d
        end do
        if (d_min <= dminf) dminf = d_min
        sum_d = sum_d + d_min
    end do
    !avg_dx = sum_d / TotSamples;
    avg_dx = dminf;
    write(*,'(a,f10.3)') "Average dx = ", avg_dx
end subroutine get_avg_dx

end module prmtrs



