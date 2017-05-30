module module_snapshots

    USE class_inputdata
    USE class_grid
    USE class_junction
    USE class_density
    USE class_system
    USE class_coupling
    USE class_ltensor
    USE class_observables
    USE class_rhox
    USE class_bosonic_bath

    implicit none

    type(grid)::x_grid
    type(densityvector)::density
    type(system):: mol_system
    type(leads)::theleads
    type(coupling)::couplings
    type(ltensor)::tensor
    type(observables)::observable_set
    type(rhox)::rhox_representation
    type(bath)::hbath

contains

    subroutine initialize_entire_system(input)

        CLASS(inputdata), intent(in)::input
        INTEGER                     :: i !< Loop Variable
        !Initialize ...
        !Compute the Position Grid in Hermite Base for the DVR Calculations
        CALL x_grid%init_grid(input)
        !----the density vector
        CALL density%init_densityvector(input, "S")
        !----the internal system
        CALL mol_system%init_system(input, x_grid, density)
        !----the Leads
        CALL theleads%init_leads(input, density, mol_system)
        !----the Couplings
        CALL couplings%init_coupling(mol_system, x_grid, density, input)
        !--- Harmonic Heath Bath
        CALL hbath%init_bath(input, density, mol_system)
        !--- tensor object
        CALL tensor%init_tensor(input, density)
        !--- Observable Object
        CALL observable_set%init_observables(input, density, mol_system)
        !--- RHOX Object for the rho(x) represention
        CALL rhox_representation%init_rhox(input, density, x_grid, mol_system)
        !---Franck Condon
        CALL couplings%write_couplings(1, "stat", 25)


    end subroutine

!
!Routine for taking a simple snaptshot of the data before calculations begin. It is generated
!by a call of the WRITE_SUMMARY_TOFILE routine in the OBSERVABLE object.
!Optionally it generates a summary file after the proposed switching
!

    subroutine take_snapshot(input, switched)

        CLASS(inputdata), intent(in)::input
        INTEGER, intent(in)         :: switched

        CALL initialize_entire_system(input)

        CALL observable_set%write_summary_tofile(mol_system, couplings, 20, "snaptshot_status_1")

        IF(switched.eq.1) THEN
            WRITE(*,*) "Switching for second Snaptshot"
            CALL mol_system%change_gate(input%end_gate_voltage)
            CALL observable_set%write_summary_tofile(mol_system, couplings, 20, "snaptshot_status_2")
        END IF


    end subroutine take_snapshot

!Subroutine drives the ltensor routines.

    SUBROUTINE theoretical_validity_hbath_test(input)

        CLASS(inputdata), intent(in)::input
        INTEGER                     ::test_pass

        CALL initialize_entire_system(input)

            CALL tensor%get_stationary_rho(density, mol_system,  theleads, couplings, hbath, input)

            CALL tensor%check_bath_ltensor(density)
            CALL tensor%analyse_bath_ltensor(mol_system, test_pass)
            CALL tensor%write_checked_bath_ltensor(mol_system)

            CALL observable_set%write_summary_tofile(mol_system, couplings, 20 , "hbath")

    END SUBROUTINE

end module module_snapshots
