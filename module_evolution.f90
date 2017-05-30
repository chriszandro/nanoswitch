module module_evolution

    USE class_inputdata
    USE class_grid
    USE class_junction
    USE class_density
    USE class_system
    USE class_coupling
    USE class_ltensor
    USE class_observables
    USE class_rhox
    USE module_base_transformation
    USE class_bosonic_bath

    implicit none

    DOUBLE COMPLEX, parameter ::    RPAR=DCMPLX(0.0D0,1.0D0)
    DOUBLE COMPLEX    WTRU
    DOUBLE COMPLEX    ERR

     !--- Arrays
    DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: ZWORK
    REAL(8), DIMENSION(:), ALLOCATABLE ::  RWORK
    INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK

    REAL(8) :: ABERR
    REAL(8) :: AEMAX=0.0d0
    INTEGER :: zdim

    !--- Tolerance
    INTEGER  :: ITOL = 1
    REAL(8) :: ATOL = 1.d-8
    REAL(8) :: RTOL = 1.d-8

    INTEGER :: IPAR
    INTEGER :: ITASK=1
    INTEGER :: ISTATE=1
    INTEGER :: IOPT=0
    INTEGER :: LZW
    INTEGER :: LRW
    INTEGER :: LWI
    INTEGER             :: MF=22!< 10: Non-stiff System 
                                !< 22:stiff method, internally generated full Jacobian
    


    !--- Right side of the differential equation d rho(t)/dt = rhs_v_d_p
    EXTERNAL rhs_v_d_p
    !--- Jacobi Matrix
    EXTERNAL jacmat_v_d_p

    !!****************************************************************************
    !Expokit
    INTEGER                              :: ldh !hm....dont know
    INTEGER                              :: lwsp !workspace dimension
    INTEGER                              :: ideg

    INTEGER                              :: ns !Number of scaling-squaring used
    INTEGER                              :: iflag !iexit flag. ! 0 - no problem  ! <0 - problem
    INTEGER, DIMENSION(:), ALLOCATABLE   :: ipiv

    INTEGER, DIMENSION(:), ALLOCATABLE :: iwsp
    COMPLEX(8), DIMENSION(:), ALLOCATABLE :: wps
    COMPLEX(8), DIMENSION(:), ALLOCATABLE :: wsp
    !!****************************************************************************

    type(grid)::x_grid

    type(densityvector)::density_status_1
    type(system):: mol_system_status_1
    type(leads)::theleads_status_1
    type(coupling)::couplings_status_1
    type(ltensor)::tensor_status_1
    type(observables)::observables_status_1
    type(bath)::hbath_status_1

    type(densityvector)::density_status_2
    type(system):: mol_system_status_2
    type(leads)::theleads_status_2
    type(coupling)::couplings_status_2
    type(ltensor)::tensor_status_2
    type(observables)::observables_status_2
    type(bath)::hbath_status_2

    type(densityvector)::density_status_3
    type(observables)::observables_status_3

    type(densityvector)::density_t
    type(observables)::observables_t

    type(observables)::observable_set
    type(rhox)::rhox_representation

    COMPLEX(8), ALLOCATABLE, DIMENSION(:,:)     :: Ltensor_public

    REAL(8), DIMENSION(:,:)        :: preparation(1:5, 1:5)!A matrix to store the observables for the stationary states in  status 1 and status 2
    REAL(8), ALLOCATABLE, DIMENSION(:)           :: zeit
    
    !Output Handling
    CHARACTER(len=100)  ::output_string
    CHARACTER(len=100)  ::output_gate
CONTAINS

    subroutine start_initial_state_preparation(input)

        CLASS(inputdata)::input

        !Compute the Position Grid in Hermite Base for the DVR Calculations
        CALL x_grid%init_grid(input)

        !>Preparate the states in status 1 and status 2
        CALL initialize_system_in_status_1(input, x_grid)
        CALL initialize_system_in_status_2(input, x_grid)
        !<

        SELECT CASE(input%initial_state)

        CASE(1)
        WRITE(*,*) '**********************************************************'
        WRITE(*,*) 'Initial State is stationary state with the Reservoirs'
        WRITE(*,*) '**********************************************************'
        !>Init rho(t)
        CALL density_t%init_densityvector(input, "S")
        !> Perform the base transformation for rho in status 1 into the the system in status 2
        CALL transform_to_newbase_eq(density_status_1,density_t, mol_system_status_1, mol_system_status_2)
        
        !Pure eigenstate
        CASE(2)
        WRITE(*,*) '**********************************************************'
        WRITE(*,*) 'Pure Initial State prepared'
        WRITE(*,*) '**********************************************************'
        CALL density_t%init_densityvector(input, "S")
        CALL density_t%set_rho_pure(input%initial_state_number, input%initial_occupation)

        !Symmetric Case 
        CASE(3)
        WRITE(*,*) '**********************************************************'
        WRITE(*,*) 'Localized initial state prepared'
        WRITE(*,*) '**********************************************************'
        !This Case may become obsolet in future
        !CALL density_t%set_rho_superposition_pure(input%initial_state_number, input%initial_state_number2, input%initial_occupation)
        CALL density_t%init_densityvector(input, "S")
        CALL density_t%set_rho_superposition_symmetric(0, input%initial_state_number, input%initial_state_number2, input%initial_occupation)
        
        !Antisymmtric Case 
        CASE(4) 
        WRITE(*,*) '**********************************************************'
        WRITE(*,*) 'Localized initial state prepared'
        WRITE(*,*) '**********************************************************'
        CALL density_t%init_densityvector(input, "S")
        CALL density_t%set_rho_superposition_symmetric(1, input%initial_state_number, input%initial_state_number2, input%initial_occupation)
        
        !Get stationary state with fermionic bath and then switch the bosonic
        !bath on
        CASE(5) 
        !>Init rho(t)
        CALL density_t%init_densityvector(input, "S")
        WRITE(*,*) '**********************************************************'
        WRITE(*,*) 'Switch Enviroment on with state in fermionic equilibrium'
        WRITE(*,*) '**********************************************************'
        !Uncouple bosonic bath and calculate everything again 
        hbath_status_1%eta=0d0
        CALL tensor_status_1%get_stationary_rho(density_status_1, mol_system_status_1,  theleads_status_1, couplings_status_1, hbath_status_1, input)
        !> Perform the base transformation for rho in status 1 into the the system in status 2
        CALL transform_to_newbase_eq(density_status_1, density_t, mol_system_status_1, mol_system_status_2)
        
        END SELECT

        !> Initialize the observalbe(t) object
        CALL observables_t%init_observables(input, density_t, mol_system_status_2)
        !>
        CALL observables_t%get_and_putin_array(preparation, 3, mol_system_status_2, theleads_status_2, density_t, couplings_status_2)
        !> Check the accuracy of the base transformation +
        CALL evaluate_preparation_intoconsole(preparation)

        CALL rhox_representation%init_rhox(input, density_t, x_grid, mol_system_status_2)


    end subroutine start_initial_state_preparation

    subroutine run_evolution_zvode(input)

        CLASS(inputdata), intent(in) :: input

        REAL(8)         :: t0 !< Initial point in time (copy)
        INTEGER         :: i,j !< Loop Variable

        CALL start_initial_state_preparation(input)


        WRITE(*,*) '**************************************************'
        WRITE(*,*) 'Evolution with ZVODE Mode'
        WRITE(*,*) '**************************************************'

        CALL set_working_arrays_zvode(density_t%medim)

        t0 = input%time_start

        ASSOCIATE(gate_volt => input%sec_parameter_grid, &
                  gate_grid_end => input%sec_parameter_grid_length,&
                  time => input%parameter_grid, &
                  grid_end => input%parameter_grid_length)

            SELECT CASE(input%evo_method)

                CASE(1)

                    DO i = 1, grid_end

                        !CALL observables_t%get_performance_time(i,"start")

                        CALL ZVODE(rhs_v_d_p, zdim , density_t%rho, t0, time(i)  , ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, ZWORK, LZW,RWORK,LRW,IWORK,LWI,jacmat_v_d_p,MF,RPAR,IPAR)

                        !CALL density_t%check_trace

                        CALL observables_t%get_and_store_observables(i, density_t, mol_system_status_2, theleads_status_2, couplings_status_2)
                        CALL observables_t%get_state_population(i, density_t)
                        
                        !SubMatrix 
                        !CALL observables_t%get_sub_matrix(i,  density_t)

                        !CALL observables_t%get_performance_time(i,"end")

                        !CALL observables_t%put_console_obs(i)

                    END DO

                    !Writing Data Output
                    !---Standard
                    !CALL observables_t%write_performance_tofile("ztime")
                    CALL observables_t%write_population_tofile("ztime", 1)
                    CALL observables_t%write_population_tofile("ztime", 0)
                    !CALL observables_t%write_summary_tofile(mol_system_status_2, couplings_status_2, 25 , "ztime")
                    CALL observables_t%write_results_tofile("ztime")
                    !CALL observables_t%write_sub_matrix("ztime")
                
                !Only Rhox
                CASE(2)
                    DO i = 1, grid_end
!                       CALL observables_t%get_performance_time(i,"start")

                        CALL ZVODE(rhs_v_d_p, zdim , density_t%rho, t0, time(i)  , ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, ZWORK, LZW,RWORK,LRW,IWORK,LWI,jacmat_v_d_p,MF,RPAR,IPAR)
                        CALL rhox_representation%get_and_store_rhox(i, density_t, x_grid, mol_system_status_2)

!                       CALL observables_t%get_performance_time(i,"end")
                    END DO

                    !Writing Data Output
                    !---Standard
                    CALL observables_t%write_performance_tofile("ztime")
                    CALL rhox_representation%write_rhox_tofile_xyz("ztime", x_grid)
                    CALL rhox_representation%write_rhox_tofile_mformat("ztime", x_grid)
                    CALL rhox_representation%write_rhox_negative_elements("ztime")


                CASE(3)

                    DO i = 1, grid_end

                        CALL observables_t%get_performance_time(i,"start")

                        CALL ZVODE(rhs_v_d_p, zdim , density_t%rho, t0, time(i)  , ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, ZWORK, LZW,RWORK,LRW,IWORK,LWI,jacmat_v_d_p,MF,RPAR,IPAR)

                        !CALL density_t%check_trace

                        CALL observables_t%get_and_store_observables(i, density_t, mol_system_status_2, theleads_status_2, couplings_status_2)
                        CALL observables_t%get_state_population(i, density_t)
                        
                        !SubMatrix 
                        !CALL observables_t%get_sub_matrix(i,  density_t)

                        CALL observables_t%get_performance_time(i,"end")

                        !CALL observables_t%put_console_obs(i)
                        !Rhox
                        CALL rhox_representation%get_and_store_rhox(i, density_t, x_grid, mol_system_status_2)
                    END DO

                    !Writing Data Output
                    !---Standard
                    CALL observables_t%write_performance_tofile("ztime")
                    CALL observables_t%write_population_tofile("ztime", 1)
                    CALL observables_t%write_population_tofile("ztime", 0)
                    CALL observables_t%write_summary_tofile(mol_system_status_2, couplings_status_2, 25 , "ztime")
                    CALL observables_t%write_results_tofile("ztime")
                    !CALL observables_t%write_sub_matrix("ztime")

                    !Rhox
                    CALL rhox_representation%write_rhox_tofile_xyz("ztime", x_grid)
                    CALL rhox_representation%write_rhox_tofile_mformat("ztime", x_grid)
                    CALL rhox_representation%write_rhox_negative_elements("ztime")

                CASE(4)
               
                WRITE(*,*) "*******************************"
                WRITE(*,*) "Time evolution Heatmap Mode"
                WRITE(*,*) "*******************************"

                DO j = 1, gate_grid_end

                t0 = input%time_start
                
                !UPDATE Ltensor
                CALL update_system_in_status_2(input, x_grid, gate_volt(j))
                CALL transform_to_newbase_eq(density_status_1,density_t, mol_system_status_1, mol_system_status_2)
                CALL observables_t%set_X_auxil_vector(mol_system_status_2, density_t)
                !Reset ZVODE Routine
                ISTATE=1
                !Generate string what gate is on
                output_string = "ztime_"//trim(adjustl(output_gate))               
                WRITE(*,*) output_string 
                
                DO i = 1, grid_end
                                CALL observables_t%get_performance_time(i,"start")

                                CALL ZVODE(rhs_v_d_p, zdim , density_t%rho, t0, time(i)  , ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, ZWORK, LZW,RWORK,LRW,IWORK,LWI,jacmat_v_d_p,MF,RPAR,IPAR) 
                                !CALL density_t%check_trace

                                CALL observables_t%get_and_store_observables(i, density_t, mol_system_status_2, theleads_status_2, couplings_status_2)
                                CALL observables_t%get_state_population(i, density_t)
                                
                                !SubMatrix 
                                !CALL observables_t%get_sub_matrix(i,  density_t)

                                CALL observables_t%get_performance_time(i,"end")

                                !CALL observables_t%put_console_obs(i)
                                !Rhox
                                CALL rhox_representation%get_and_store_rhox(i, density_t, x_grid, mol_system_status_2)

                            END DO
             
                
                           !Writing Data Output
                            !---Standard
                            CALL observables_t%write_performance_tofile(output_string)
                            CALL observables_t%write_population_tofile(output_string, 1)
                            CALL observables_t%write_population_tofile(output_string, 0)
                            CALL observables_t%write_summary_tofile(mol_system_status_2, couplings_status_2, 25 , output_string) 
                            CALL observables_t%write_results_tofile(output_string)
                            !CALL observables_t%write_sub_matrix(output_string)

                            !Rhox
                            CALL rhox_representation%write_rhox_tofile_xyz(output_string, x_grid)
                            CALL rhox_representation%write_rhox_tofile_mformat(output_string, x_grid)
                            CALL rhox_representation%write_rhox_negative_elements(output_string)
                            CALL rhox_representation%write_prob_integral(output_string)

                            !Write Tunnel Probability
                            CALL rhox_representation%write_prob_integral(output_string)

                            !---Observables and Tunnel probability in Matrix Form
                            CALL observable_set%put_in_observable_matrix(j)
                            CALL rhox_representation%put_tunnel_probability_matrix(j)
       
                        END DO
                        
                        CALL observables_t%write_parameter_grid
                        CALL mol_system_status_2%potential_surface("potential", 1)

                        !Write Heatmaps Data 
                        !--- XYZ Format
                        CALL observables_t%write_observable_energy_matrix_tofile_xyz("obs_matrix_time")
                        CALL observables_t%write_observable_position_matrix_tofile_xyz("obs_matrix_time")
                        CALL observables_t%write_observable_current_matrix_tofile_xyz("obs_matrix_time")
                        CALL observables_t%write_observable_population_matrix_tofile_xyz("obs_matrix_time", 1)
                        CALL observables_t%write_observable_population_matrix_tofile_xyz("obs_matrix_time", 0)

                        !--- M-Format
                        CALL observables_t%write_observable_current_matrix_heatmap("heatmap_time")
                        CALL observables_t%write_observable_position_matrix_heatmap("heatmap_time")
                        CALL observables_t%write_observable_energy_matrix_heatmap("heatmap_time")
                        CALL observables_t%write_observable_population_matrix_heatmap("heatmap_time",1)
                        CALL observables_t%write_observable_population_matrix_heatmap("heatmap_time",0)

                        !State_populations 
                        CALL observables_t%write_observable_state_population_matrix_heatmap("heatmap_time", 1)
                        CALL observables_t%write_observable_state_population_matrix_heatmap("heatmap_time", 0)

                        !Tunnel Probability Heatmap
                        CALL rhox_representation%write_tunnel_probability_matrix("heatmap_time")

                CASE DEFAULT

                    WRITE(*,*)'No proper evo method for zvode chosen'

            END SELECT


        END ASSOCIATE
    end subroutine run_evolution_zvode

    subroutine run_evolution_expokit(input)

    CLASS(inputdata), intent(in):: input
    INTEGER                     :: ideg
    INTEGER                     :: i

    COMPLEX(8), ALLOCATABLE, DIMENSION(:) :: rho0

      CALL start_initial_state_preparation(input)

        WRITE(*,*) '**************************************************'
        WRITE(*,*) 'Evolution with EXPOKIT Mode'
        WRITE(*,*) '**************************************************'

        ideg = input%ideg

        CALL set_working_arrays_expokit(density_t%medim)

        ALLOCATE(rho0(density_t%medim))
        rho0 = density_t%rho

        ASSOCIATE(time => input%parameter_grid, &
            grid_end => input%parameter_grid_length)


        SELECT CASE(input%evo_method)

        CASE(1)

        WRITE(*,'(A)') "    ++++++++++++++++++++++++++++++++++++++++++++++++++"
        WRITE(*,'(A)') "    Chebyshev Method Full Mode"
        WRITE(*,'(A)') "    ++++++++++++++++++++++++++++++++++++++++++++++++++"

                DO i = 1, grid_end
                        CALL observables_t%get_performance_time(i,"start")

                        !The routine ZGCHBV overrides the input density in every Call. Therefore
                        !a copy of the initial rho is copied.
                        CALL ZCOPY(density_t%medim, rho0, 1, density_t%rho, 1)

                        !Calculate exp(A*t)
                        CALL ZGCHBV(density_t%medim, time(i), LTensor_public, ldh, density_t%rho, wps, iwsp, iflag)

                        !CALL density_t%check_trace

                        CALL observables_t%get_and_store_observables(i, density_t, mol_system_status_2, theleads_status_2, couplings_status_2)
                        !CALL observables_t%get_state_population(i, density_t)
                        !CALL rhox_representation%get_and_store_rhox(i, density_t, x_grid, mol_system_status_2)

                        CALL observables_t%get_performance_time(i,"end")

                        !CALL observables_t%put_console_obs(i)
                END DO

                    !Writing Data Output
                    !---Standard
                    CALL observables_t%write_performance_tofile("exptime")
                    CALL observables_t%write_population_tofile("exptime", 1)
                    CALL observables_t%write_summary_tofile(mol_system_status_2, couplings_status_2, 25 , "exptime")
                    CALL observables_t%write_results_tofile("exptime")
                    CALL rhox_representation%write_rhox_tofile_xyz("exptime", x_grid)
                    CALL rhox_representation%write_rhox_negative_elements("exptime")

        CASE(2)
        WRITE(*,'(A)') "    ++++++++++++++++++++++++++++++++++++++++++++++++++"
        WRITE(*,'(A)') "    Chebyshev Method Reduced Mode 2"
        WRITE(*,'(A)') "    ++++++++++++++++++++++++++++++++++++++++++++++++++"

                DO i = 1, grid_end
                        CALL observables_t%get_performance_time(i,"start")

                        !The routine ZGCHBV overrides the input density in every Call. Therefore
                        !a copy of the initial rho copied.
                        CALL ZCOPY(density_t%medim, rho0, 1, density_t%rho, 1)

                        !Calculate exp(A*t)
                        CALL ZGCHBV(density_t%medim, time(i), LTensor_public, ldh, density_t%rho, wps, iwsp, iflag)

                        CALL rhox_representation%get_and_store_rhox(i, density_t, x_grid, mol_system_status_2)

                        CALL observables_t%get_performance_time(i,"end")

                END DO

                    !Writing Data Output
                    !---Standard
                    CALL observables_t%write_performance_tofile("exptime")
                    CALL rhox_representation%write_rhox_tofile_xyz("exptime", x_grid)
                    CAlL rhox_representation%write_rhox_tofile_mformat("exptime", x_grid)
                    CALL rhox_representation%write_rhox_negative_elements("exptime")

        CASE DEFAULT

                    WRITE(*,*)'No proper evo method for expokit chosen'

        END SELECT

        END ASSOCIATE

    end subroutine run_evolution_expokit

    subroutine evaluate_preparation_intoconsole(preparation)
        ! Preparation Matrix
        ! (1,1) (1,2) (1,3) (1,4) ! pos1 cu1 ener1 occu1  <- stationary state in status 1
        ! (2,1) (2,2) (2,3) (2, 4) ! pos2 cu2 ener2 occu2 <-stationary state in status 2
        ! (3,1) (3,2) (3,3) (3,4) ! pos3 cu3 ener3 occu3 <- stationary rho in status 1 pluged into system 2
        ! (4,1) (4,2) (4,3) (4,4) ! pos4f cu4f ener4f occu4f


        REAL(8), intent(inout)         :: preparation(1:, 1:)

        INTEGER i

        DO i = 1, 4
            preparation(4, i)= ABS(preparation(1,i) - preparation(3,i))
        END DO

        !Ausgabe
        WRITE(*,'(A)') "Base Transformation Check"
        WRITE(*,'(A)') "----------------------------------------------"
        WRITE(*,'(T20, A, T40, A, T60, A, T80, A)') "Old Base", "New Base", "Error"
        WRITE(*, '(A, T22, E, T50, E, T80, E, T110, F10.3)') "Position", preparation(1,1), preparation(3,1), preparation(4,1)
        WRITE(*, '(A, T22, E, T50, E, T80, E, T110, F10.3)') "Current", preparation(1,2), preparation(3,2)
        WRITE(*, '(A, T22, E, T50, E, T80, E, T110, F10.3)') "Energy",  preparation(1,3), preparation(3,3)
        WRITE(*, '(A, T22, E, T50, E, T80, E, T110, F10.3)') "Occupation", preparation(1,4), preparation(3,4), preparation(4,4)

        WRITE(*,'(A)') "Stationary State the system is evaluating"
        WRITE(*,'(A)') "----------------------------------------------"
        WRITE(*,'(T20, A, T40, A, T60, A, T80, A)') "Initial Value", "Stationary Value"
        WRITE(*, '(A, T22, E, A, T50, E, T80, E, T110, F10.3)') "Position", preparation(1,1), "-->", preparation(2,1)
        WRITE(*, '(A, T22, E, A, T50, E, T80, E, T110, F10.3)') "Current",preparation(1,2), "-->", preparation(2,2)
        WRITE(*, '(A, T22, E, A, T50, E, T80, E, T110, F10.3)') "Energy", preparation(1,3), "-->",preparation(2,3)
        WRITE(*, '(A, T22, E, A, T50, E, T80, E, T110, F10.3)') "Occupation", preparation(1,4), "-->",preparation(2,4)

    end subroutine

    subroutine initialize_system_in_status_1(input, x_grid)

        CLASS(inputdata), intent(in) :: input
        CLASS(grid), intent(in) :: x_grid

        WRITE(*,*) 'Preparing density matrix for initial state'
        !Initialize ...
        !----the densityvector
        CALL density_status_1%init_densityvector(input, "S")
        !--- set system with gate voltage 1
        CALL mol_system_status_1%init_system(input, x_grid, density_status_1, gate_in=1)
        !----the Leads
        CALL theleads_status_1%init_leads(input, density_status_1, mol_system_status_1, voltage_in=1)
        !----the Couplings
        CALL couplings_status_1%init_coupling(mol_system_status_1, x_grid, density_status_1, input)
        !--- Harmonic Heath Bath
        CALL hbath_status_1%init_bath(input, density_status_1, mol_system_status_1)
        !---Observables
        CALL observables_status_1%init_observables(input, density_status_1, mol_system_status_1)
        !--- tensor object
        CALL tensor_status_1%init_tensor(input, density_status_1)

        !stationary rho
        CALL tensor_status_1%get_stationary_rho(density_status_1, mol_system_status_1,  theleads_status_1, couplings_status_1, hbath_status_1, input)

        !Stationary rho, observables in the potential before switching
        CALL observables_status_1%get_and_putin_array(preparation, 1, mol_system_status_1, theleads_status_1, density_status_1, couplings_status_1)

        CALL observables_status_1%write_summary_tofile(mol_system_status_1, couplings_status_1 , 20 , "switch_off")

        !Franck-Condon I
        CALL couplings_status_1%write_couplings(1, "evolution_status1", 25)

    end subroutine initialize_system_in_status_1

    subroutine initialize_system_in_status_2(input, x_grid)

        CLASS(inputdata), intent(in) :: input
        CLASS(grid), intent(in) :: x_grid

        WRITE(*,*) 'Preparing Ltensor for time evolution'
        !Initialize ...
        !----the densityvector
        CALL density_status_2%init_densityvector(input, "S")
        !--- set system with gate voltage 1
        CALL mol_system_status_2%init_system(input, x_grid, density_status_2, gate_in=2)
        !----the Leads
        CALL theleads_status_2%init_leads(input, density_status_2, mol_system_status_2, voltage_in=2)
        !----the Couplings
        CALL couplings_status_2%init_coupling(mol_system_status_2, x_grid, density_status_2, input)
        !--- Harmonic Heath Bath
        CALL hbath_status_2%init_bath(input, density_status_2, mol_system_status_2)
        !---Observables
        CALL observables_status_2%init_observables(input, density_status_2, mol_system_status_2)
        !--- tensor object
        CALL tensor_status_2%init_tensor(input, density_status_2)

        !stationary rho
        CALL tensor_status_2%get_stationary_rho(density_status_2, mol_system_status_2,  theleads_status_2, couplings_status_2, hbath_status_2, input )

        !Stationary rho, to have a idea where the observables are going to
        CALL observables_status_2%get_and_putin_array(preparation, 2, mol_system_status_2, theleads_status_2, density_status_2, couplings_status_2)

        !CALL observables_status_2%write_summary_tofile(mol_system_status_2, couplings_status_2 , 20 , "switch_on")

        ALLOCATE(Ltensor_public(tensor_status_2%dimension, tensor_status_2%dimension))
        Ltensor_public = tensor_status_2%Ltensor

        !Franck Condon Matrix Elements status 2
        CALL couplings_status_2%write_couplings(1, "evolution_status2", 25)

    end subroutine initialize_system_in_status_2

    subroutine update_system_in_status_2(input, x_grid, gate)

        CLASS(inputdata), intent(in)    ::input
        CLASS(grid), intent(in)         ::x_grid
        REAL(8), intent(in)             ::gate

        CHARACTER(len=100)  ::gate_string
        CHARACTER(len=100)  ::output_string
        
        !RUN Compoment UPDATES
        CALL mol_system_status_2%change_gate(gate)
        CALL theleads_status_2%update(mol_system_status_2)
        CALL couplings_status_2%update(mol_system_status_2)
        CALL hbath_status_2%update_bath(mol_system_status_2)

        !Run Update Ltensor
        !CALL tensor_status_2%get_stationary_rho(density_status_2, mol_system_status_2, theleads_status_2, couplings_status_2, hbath_status_2, input )
        CALL tensor_status_2%set_tensor(density_status_2, mol_system_status_2, theleads_status_2, couplings_status_2, hbath_status_2)
        Ltensor_public = tensor_status_2%Ltensor

        WRITE(gate_string, '(F10.3)' ) gate
        output_gate = "switch_on_gate_at_"//trim(adjustl(gate_string))
 
        !Output
        !CALL observables_status_2%write_summary_tofile(mol_system_status_2, couplings_status_2 , 20 , output_string)

        !Franck Condon Matrix Elements status 2
        !CALL couplings_status_2%write_couplings(1, "evolution_status2", 25)

    end subroutine update_system_in_status_2

    subroutine set_working_arrays_zvode(m)

        INTEGER, intent(in)         :: m !< Dimension of the Probblem necessary to calculate the array dimensions

        zdim = m
        !---------------Routine Variables
        !---- LRW
        !LRW    = Declared length of RWORK (in user's DIMENSION statement).
        LRW = 2*m + 2
        !---- LIW
        !LIW    = Declared length of IWORK (in user's DIMENSION statement).
        LWI = 2*m + 30
        !--- ZWORK
        SELECT CASE(MF)
            CASE(10)
                !---- LZW
                !LZW    = Declared length of ZWORK (in user's DIMENSION statement).
                LZW = 15*m
            CASE(21)
                LZW = 8*m + 2*m**2
            CASE(22)
                LZW = 8*m + 2*m**2
            CASE(24)
                WRITE(*,*) "User defined jacobian band matrix"
                !Not clear what MU is
                !ALLOCATE(ZWORK(10*m + (3*ML + 2*MU)*m))
            CASE(25)
                WRITE(*,*) "User defined jacobian band matrix"
                !Not clear what MU is
                !ALLOCATE(ZWORK(10*m + (3*ML + 2*MU)*m))
        END SELECT

        !--- ZWORK
        ALLOCATE(ZWORK(LZW))
        !--- RWORK
        ALLOCATE(RWORK(LRW))
        !--- IWORK
        ALLOCATE(IWORK(LWI))

        WRITE(*,*) 'Working Arrays set to'
        WRITE(*,'(A,I8,A,I8,A,I8)') 'ZWORK: ',  SIZE(ZWORK, 1), "RWORK: ", SIZE(RWORK,1), "IWORK: ",  SIZE(IWORK,1)

    end subroutine set_working_arrays_zvode


    subroutine set_working_arrays_expokit(m)

        INTEGER, intent(in)         ::m

        lwsp = 8*(4*m*m+ideg+1)
        ldh = m

        ALLOCATE(wps(lwsp))
        ALLOCATE(wsp(lwsp))
        ALLOCATE(iwsp(m))
        ALLOCATE(ipiv(m))

    end subroutine set_working_arrays_expokit

end module module_evolution

!****************************************************************************************************
!Subroutine get_function_rhs
!Right hand side for time evolution
!
!
!****************************************************************************************************

SUBROUTINE rhs_v_d_p( n , time_in, w, w_t, RPAR_in , IPAR_in)

    USE module_evolution

    IMPLICIT NONE

    INTEGER, intent(in) :: n
    INTEGER, intent(in) :: IPAR_in
    DOUBLE COMPLEX, intent(in) :: RPAR_in

    DOUBLE COMPLEX, intent(in) :: w(n)
    DOUBLE COMPLEX, intent(out) :: w_t(n)

    REAL(8) :: time_in


    !Multiplikation
    !Fortran inbuilt function
    !w_t =  matmul(Ltensor_public, w)

    !With MKL
    CALL  zgemv('N',density_t%medim, density_t%medim, (1.0d0,0.0d0), Ltensor_public, density_t%medim, w , 1, (0.0d0,0.0d0),w_t, 1)



    RETURN
END SUBROUTINE rhs_v_d_p

!****************************************************************************************************
!Subroutine get_jacobi
!Analytical Jacoby matrix for time evolution
!
!
!************* ***************************************************************************************
SUBROUTINE jacmat_v_d_p(m, time, Y, ML, MU, PD, NRPD, RPAR, IPAR)

    IMPLICIT NONE

    REAL(8) :: time
    INTEGER m

    DOUBLE COMPLEX :: Y(m)
    DOUBLE COMPLEX :: PD(NRPD, m)
    DOUBLE COMPLEX, intent(in) :: RPAR

    INTEGER, intent(in) :: ML
    INTEGER, intent(in) :: MU

    INTEGER, intent(in) :: NRPD
    INTEGER, intent(in) :: IPAR

    ! PD = LTensor_public

    RETURN
END SUBROUTINE  jacmat_v_d_p
