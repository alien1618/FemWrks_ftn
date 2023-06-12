program FemWrks
    use trnsprt
    use elst
    use ns
    use trnsprt_ebe
    use elst_ebe
    use ns_ebe
    use prmtrs
!---------------------------------------------------------------------------------------------
    integer :: n
!---------------------------------------------------------------------------------------------
    
    call set_slvr(n) 
    select case (n)
        case (1)
            call run_trnsprt_ebe()
        case (2)
            call run_trnsprt()
        case (3)
            call run_elst_ebe()
        case (4)
            call run_elst()
        case (5)
            call run_ns_ebe()
        case (6)
            call run_ns()
    end select
end program FemWrks
