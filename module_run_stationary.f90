module module_run_stationary

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

    subroutine run_current_voltage_characteristics(input)

        CLASS(inputdata), intent(in)::input
        INTEGER                     :: i !< Loop Variable
        INTEGER                     :: test_pass 

        CALL initialize_entire_system(input)

        !Loop over the voltages and calculation of the stationary state
        ASSOCIATE(volt => input%parameter_grid, &
            grid_end => input%parameter_grid_length)

            DO i = 1, grid_end

                !Performance
                CALL observable_set%get_performance_time(i,"start")

                !Variation and Computation of rho
                CALL theleads%change_voltage( volt(i))
                CALL tensor%get_stationary_rho(density, mol_system,  theleads, couplings, hbath, input)

                !Processing rho
                !---Observables
                CALL observable_set%get_and_store_observables(i, density, mol_system, theleads, couplings)
                CALL observable_set%get_state_population(i, density)
                !---Rho(x)
                CALL rhox_representation%get_and_store_rhox(i, density, x_grid, mol_system)

                !Performance
                CALL observable_set%get_performance_time(i,"end")

                WRITE(*,*) volt(i), i
            END DO

            CALL tensor%check_bath_ltensor(density)
            CALL tensor%analyse_bath_ltensor(mol_system, test_pass)

            !Writing Data Output
            !---Standard
            CALL observable_set%write_performance_tofile("stat")
            CALL observable_set%write_results_tofile("stat")
            CALL observable_set%write_population_tofile("stat",1)
            CALL observable_set%write_population_tofile("stat",0)
            CALL observable_set%write_summary_tofile(mol_system, couplings, 20 , "stat")
            CALL observable_set%write_population_histo_tofile("stat",0)
            CALL observable_set%write_population_histo_tofile("stat",1)
            !---rhox
            CALL rhox_representation%write_rhox_tofile_xyz("stat", x_grid)
            CALL rhox_representation%write_rhox_tofile_mformat("stat", x_grid)
            CALL rhox_representation%write_prob_integral("stat")
            CALL mol_system%write_system_to_file(25, "stat")          
            CALL mol_system%write_x_energy(1, "stat", 25)
        END ASSOCIATE

    end subroutine  run_current_voltage_characteristics

    subroutine run_gate_voltage(input)

        CLASS(inputdata), intent(in)::input

        INTEGER                     :: i,j!< Loop Variables
        CHARACTER(len=20)           :: bias_string
        CHARACTER(len=20)           :: gate_string
        CHARACTER(len=100)          :: output_string

        INTEGER                     :: test_pass

        CALL initialize_entire_system(input)
       
        !This optional, only the object properties are set. 
        CALL density%set_rho_superposition_symmetric(0, input%initial_state_number, input%initial_state_number2, input%initial_occupation)
        
        ASSOCIATE(gate_volt => input%parameter_grid, &
            gate_grid_end => input%parameter_grid_length,&
            bias_volt => input%sec_parameter_grid, &
            bias_grid_end => input%sec_parameter_grid_length)

            DO i = 1, bias_grid_end
                    !Change voltage of the leads
                    CALL theleads%change_voltage(bias_volt(i))
                
                    DO j = 1, gate_grid_end

                    WRITE(*,*)
                    WRITE(*,'(A, E10.3, I)') "Gate Voltage: ", gate_volt(j), j
                    WRITE(*,'(A, E10.3, I)') "Bias Voltage: ", bias_volt(i), i
                    WRITE(*,*)

                    !Performance
                    CALL observable_set%get_performance_time(j,"start")

                    !Variation of the gate and new Computation of the system, leads and the couplings
                    !Eigensystem and Eigenenergies 
                    CALL mol_system%change_gate(gate_volt(j))
                    !Leads
                    CALL couplings%update(mol_system)
                    CALL theleads%update(mol_system)
                    !Bath | Remark: Elements <n|x|m> is calculated in mol_system%change_gate
                    CALL hbath%update_bath(mol_system)

                    !Frank Matrix, <n | s(x) | m>^2 
                    CALL observable_set%put_in_frank_matrix(couplings, j)

                    !Position in energy representation <n | x | m> equivalent of the frank
                    !matrix for the harmonich bat  
                    CALL observable_set%put_in_x_energy_matrix(mol_system, j)

                    !''''''''''''''''''''''''''''''''''''''' 
                    !Calculate Rho
                    CALL tensor%get_stationary_rho(density, mol_system,  theleads, couplings, hbath, input)
                    !''''''''''''''''''''''''''''''''''''''' 
                    
                    !''''''''''''''''''''''''''''''''''''''
                    !Check if Rho is calculated within the Theorie
                    !''''''''''''''''''''''''''''''''''''''
                    CALL tensor%check_bath_ltensor(density)
                    CALL tensor%analyse_bath_ltensor(mol_system, test_pass)
                    WRITE(*,*) '**********Testpass****************' 
                    WRITE(*,*) test_pass
                    WRITE(*,*) '**********************************' 
                    CALL observable_set%get_hbath_validity(test_pass, j)
                    !''''''''''''''''''''''''''''''''''''''
                    
                    !Processing rho
                    !---Observables
                    CALL observable_set%get_and_store_observables_mod_system(j, density, mol_system, theleads, couplings)
                    CALL observable_set%get_state_population(j, density)
                    !---Rho(x)
                    CALL rhox_representation%update_rhox(density, mol_system, x_grid)
                    CALL rhox_representation%get_and_store_rhox(j, density, x_grid, mol_system)

                    !Optional
                    CALL rhox_representation%get_wave_integral(x_grid, mol_system,j)

                    !Energys
                    CALL observable_set%get_state_energy(j, mol_system, density)

                    !Performance
                    CALL observable_set%get_performance_time(j,"end")
                    !State Eigenenergie 
                    CALL observable_set%get_state_energy(j, mol_system, density)
                    
                    WRITE(gate_string, '(F10.3)' ) gate_volt(j)
                    output_string= "gate_at_"//trim(adjustl(gate_string))
                    !CALL mol_system%write_system_to_file(15, output_string)          

                    !Optonal
                    CALL observable_set%get_sub_matrix(j,  density)
                
                END DO

                !Generate suffix for the outputfiles for better specification
                WRITE(bias_string, '(F10.3)' ) bias_volt(i)
                output_string= "switchmod_bias_at"//trim(adjustl(bias_string))

                CALL mol_system%potential_surface(output_string, 1)
                !Writing Data Output
                !---Standard
                CALL observable_set%write_performance_tofile(output_string)
                CALL observable_set%write_results_tofile(output_string)
                CALL observable_set%write_state_energy_tofile(output_string)
                CALL observable_set%write_population_tofile(output_string, 1)
                CALL observable_set%write_population_tofile(output_string, 0)
                !---rhox
                CALL rhox_representation%write_rhox_tofile_xyz(output_string,  x_grid)
                !--- Tunnelprobability 
                CALL rhox_representation%write_prob_integral(output_string)
                CALL rhox_representation%write_wave_integral(output_string)
               
                !Write the transition elements <n|s(x)|m>^2 for every gate voltage
                CALL observable_set%write_frank_matrix(output_string, 12)

                !Write the transition elements <n|x|m>^2 for every gate voltage
                CALL observable_set%write_x_energy_matrix(output_string, 12)

                !Peformance for runtime estimation
                CALL observable_set%get_performance_sum(gate_grid_end)
                
                !Optional
                CALL observable_set%write_sub_matrix(output_string)
                
                !Validity Matrix
                CALL observable_set%write_hbath_validity(output_string)

            END DO

        END ASSOCIATE
    end subroutine run_gate_voltage

subroutine run_current_voltage_characteristics_3d(input)
    CLASS(inputdata), intent(in)::input

    INTEGER                     :: i,j!< Loop Variables
    CHARACTER(len=20)           :: gate_string
    CHARACTER(len=100)          :: output_string

    CALL initialize_entire_system(input)

    !Loop over the voltages and calculation of the stationary state
    ASSOCIATE(gate_volt => input%sec_parameter_grid, &
        gate_grid_end => input%sec_parameter_grid_length,&
        bias_volt => input%parameter_grid, &
        bias_grid_end => input%parameter_grid_length)

        DO i = 1, gate_grid_end

            !Variation of the gate and new Computation of the system, leads and the couplings
            CALL mol_system%change_gate(gate_volt(i))
            CALL theleads%update(mol_system)
            CALL couplings%update(mol_system)
            CALL hbath%update_bath(mol_system)

            DO j = 1, bias_grid_end
                WRITE(*,*)
                WRITE(*,'(A, E10.3, I)') "Bias Voltage: ", bias_volt(j), j
                WRITE(*,'(A, E10.3, I)') "Gate Voltage: ", gate_volt(i), i
                WRITE(*,*)

                !Performance
                !CALL observable_set%get_performance_time(j,"start")
                !Change voltage of the leads
                CALL theleads%change_voltage(bias_volt(j))
                CALL tensor%get_stationary_rho(density, mol_system,  theleads, couplings, hbath, input)
                !Processing rho
                !---Observables
                CALL observable_set%get_and_store_observables_mod_system(j, density, mol_system, theleads, couplings)
                CALL observable_set%get_state_population(j, density)
                !---Rho(x)
                CALL rhox_representation%update_rhox(density, mol_system, x_grid)
                CALL rhox_representation%get_and_store_rhox(j, density, x_grid, mol_system)
                !Performance
                !CALL observable_set%get_performance_time(j,"end")

            END DO

            !Generate suffix for the outputfiles for better specification
            WRITE(gate_string, '(F10.3)' ) gate_volt(i)
            output_string= "gate_at_"//trim(adjustl(gate_string))

            !Writing Data Output
            !---Standard
            CALL observable_set%write_performance_tofile(output_string)
            CALL observable_set%write_results_tofile(output_string)
            CALL observable_set%write_population_tofile(output_string, 1)
            CALL observable_set%write_population_tofile(output_string, 0)
            CALL observable_set%write_summary_tofile(mol_system, couplings, 20 , output_string)
            !---rhox
            CALL rhox_representation%write_rhox_tofile_xyz(output_string,  x_grid)
            !--- Tunnelprobability 
            CALL rhox_representation%write_prob_integral(output_string)
            !---Observables and Tunnel probability in Matrix Form
            CALL observable_set%put_in_observable_matrix(i)
            CALL rhox_representation%put_tunnel_probability_matrix(i)
            !Peformance for runtime estimation
            CALL observable_set%get_performance_sum(gate_grid_end)
        END DO
        
        !Writing the grid
        CALL observable_set%write_parameter_grid
        CALL mol_system%potential_surface("potential", 1)
        
        !In and Output of matrix data
        CALL observable_set%write_observable_energy_matrix_tofile_xyz("obs_matrix")
        CALL observable_set%write_observable_position_matrix_tofile_xyz("obs_matrix")
        CALL observable_set%write_observable_current_matrix_tofile_xyz("obs_matrix")
        CALL observable_set%write_observable_population_matrix_tofile_xyz("obs_matrix", 1)
        CALL observable_set%write_observable_population_matrix_tofile_xyz("obs_matrix", 0)

        !Output of heatmpas
        CALL observable_set%write_observable_current_matrix_heatmap("heatmap")
        CALL observable_set%write_observable_position_matrix_heatmap("heatmap")
        CALL observable_set%write_observable_energy_matrix_heatmap("heatmap")
        CALL observable_set%write_observable_population_matrix_heatmap("heatmap",1)
        CALL observable_set%write_observable_population_matrix_heatmap("heatmap",0)
        !State_populations 
        CALL observable_set%write_observable_state_population_matrix_heatmap("heatmap", 1)
        CALL observable_set%write_observable_state_population_matrix_heatmap("heatmap", 0)
        !Tunnel Probabilitys
        CALL rhox_representation%write_tunnel_probability_matrix("heatmap")
    END ASSOCIATE

end subroutine run_current_voltage_characteristics_3d

    subroutine run_cvc_with_bias_field(input)

        CLASS(inputdata), intent(in)::input

        INTEGER                     :: i,j!< Loop Variables
        CHARACTER(len=20)           :: gate_string
        CHARACTER(len=20)           :: bias_string
        CHARACTER(len=100)          :: output_string

        CALL initialize_entire_system(input)
        
        !Corrent potential is used
       
        !Loop over the voltages and calculation of the stationary state
        ASSOCIATE(gate_volt => input%sec_parameter_grid, &
            gate_grid_end => input%sec_parameter_grid_length,&
            bias_volt => input%parameter_grid, &
            bias_grid_end => input%parameter_grid_length)

                !Writing the grid
                CALL observable_set%write_parameter_grid

                CALL mol_system%write_system_to_file(25, "diode_start")          

            DO i = 1, gate_grid_end

                DO j = 1, bias_grid_end
                    WRITE(*,*)
                    WRITE(*,'(A, E10.3, I)') "Bias Voltage: ", bias_volt(j), j
                    WRITE(*,'(A, E10.3, I)') "Gate Voltage: ", gate_volt(i), i
                    WRITE(*,*)

                    !Performance
                    CALL observable_set%get_performance_time(j,"start")
                    !Variation of the gate and new Computation of the system, leads and the couplings
                    CALL mol_system%change_gate_and_bias_field(gate_volt(i), bias_volt(j))
                    CALL theleads%update(mol_system)
                    CALL couplings%update(mol_system)
                    CALL hbath%update_bath(mol_system)

                    !Change voltage of the leads
                    CALL theleads%change_voltage(bias_volt(j))
                    CALL tensor%get_stationary_rho(density, mol_system,  theleads, couplings, hbath, input)

                    !Frank Matrix, <n | s(x) | m>^2 
                    CALL observable_set%put_in_frank_matrix(couplings, j)

                    !Processing rho
                    !---Observables
                    CALL observable_set%get_and_store_observables_mod_system(j, density, mol_system, theleads, couplings)
                    CALL observable_set%get_state_population(j, density)
                    CALL observable_set%get_state_energy(j, mol_system, density)
                    !---Rho(x)
                    CALL rhox_representation%update_rhox(density, mol_system, x_grid)
                    CALL rhox_representation%get_and_store_rhox(j, density, x_grid, mol_system)
                    !Performance
                    CALL observable_set%get_performance_time(j,"end")
                
                    !Generate suffix for the outputfiles for better specification
                    WRITE(bias_string, '(F10.3)' ) bias_volt(i)
                    output_string= "bias_field__"//trim(adjustl(bias_string))
                   
                    !BEWARE THAT THIS PRODUCES A HUGE AMOUNT OF DATA IF
                    !summary_bool=1
!                    CALL observable_set%write_summary_tofile(mol_system, couplings, 20 , output_string)

                END DO

                !Generate suffix for the outputfiles for better specification
                WRITE(gate_string, '(F10.3)' ) gate_volt(i)
                output_string= "gate_w_bias_"//trim(adjustl(gate_string))

                CALL mol_system%potential_surface(output_string, 1)
                !Writing Data Output
                !---Standard
                CALL observable_set%write_performance_tofile(output_string)
                CALL observable_set%write_results_tofile(output_string)
                CALL observable_set%write_population_tofile(output_string, 1)
                CALL observable_set%write_population_tofile(output_string, 0)
                CALL observable_set%write_state_energy_tofile(output_string)
                CALL observable_set%write_summary_tofile(mol_system, couplings, 20 , output_string)
                !---rhox
                CALL rhox_representation%write_rhox_tofile_xyz(output_string,  x_grid)
                CALL rhox_representation%write_rhox_tofile_mformat(output_string, x_grid)
                !---Observables in Matrix Form
                CALL observable_set%put_in_observable_matrix(i)
                CALL rhox_representation%put_tunnel_probability_matrix(j)
                !Peformance for runtime estimation
                CALL observable_set%get_performance_sum(gate_grid_end)
                !Write the transition elements <n|s(x)|m>^2 for every gate voltage
                CALL observable_set%write_frank_matrix(output_string, 12)

            END DO

            !In and Output of matrix data
            CALL observable_set%write_observable_energy_matrix_tofile_xyz("obs_matrix")
            CALL observable_set%write_observable_position_matrix_tofile_xyz("obs_matrix")
            CALL observable_set%write_observable_current_matrix_tofile_xyz("obs_matrix")
            CALL observable_set%write_observable_population_matrix_tofile_xyz("obs_matrix", 1)

            !Output of heatmpas
            CALL observable_set%write_observable_current_matrix_heatmap("heatmap")
            CALL observable_set%write_observable_position_matrix_heatmap("heatmap")
            CALL observable_set%write_observable_energy_matrix_heatmap("heatmap")
            CALL observable_set%write_observable_population_matrix_heatmap("heatmap",1)
            !Tunnel Probability
            CALL rhox_representation%write_tunnel_probability_matrix("heatmap")
        END ASSOCIATE

    END SUBROUTINE run_cvc_with_bias_field

    subroutine run_gate_temperature_characteristics_3d(input)
        CLASS(inputdata), intent(in)::input

        INTEGER                     :: i,j!< Loop Variables
        CHARACTER(len=20)           :: gate_string
        CHARACTER(len=20)           :: bias_string
        CHARACTER(len=100)          :: output_string

        CALL initialize_entire_system(input)

        !Loop over the voltages and calculation of the stationary state
        ASSOCIATE(gate_volt => input%sec_parameter_grid, &
            gate_grid_end => input%sec_parameter_grid_length,&
            temperature => input%parameter_grid, &
            temperature_end => input%parameter_grid_length)

            DO i = 1, gate_grid_end

                !Variation of the gate and new Computation of the system, leads and the couplings
                CALL mol_system%change_gate(gate_volt(i))
                CALL theleads%update(mol_system)
                CALL couplings%update(mol_system)
                CALL hbath%update_bath(mol_system)
               
                DO j = 1, temperature_end
                    WRITE(*,*)
                    WRITE(*,'(A, E10.3, I)') "Temperature: ", temperature(j), j
                    WRITE(*,'(A, E10.3, I)') "Gate Voltage: ", gate_volt(i), i
                    WRITE(*,*)

                    !Performance
                    CALL observable_set%get_performance_time(j,"start")
                    !Change Temperature of the enviroment 
                    !---Leads
                    CALL theleads%change_temperature(temperature(j))
                    !---Harmonic Bath
                    CALL  hbath%update_and_change_temperature(mol_system,temperature(j))   
                    
                    CALL tensor%get_stationary_rho(density, mol_system,  theleads, couplings, hbath, input)
                    !Processing rho
                    !---Observables
                    CALL observable_set%get_and_store_observables_mod_system(j, density, mol_system, theleads, couplings)
                    CALL observable_set%get_state_population(j, density)
                    !---Rho(x)
                    CALL rhox_representation%update_rhox(density, mol_system, x_grid)
                    CALL rhox_representation%get_and_store_rhox(j, density, x_grid, mol_system)
                    !Performance
                    CALL observable_set%get_performance_time(j,"end")

                END DO

                !Generate suffix for the outputfiles for better specification
                WRITE(gate_string, '(F10.3)' ) gate_volt(i)
                output_string= "temp3d_gate_at_"//trim(adjustl(gate_string))

                !Writing Data Output
                !---Standard
                CALL observable_set%write_performance_tofile(output_string)
                CALL observable_set%write_results_tofile(output_string)
                CALL observable_set%write_population_tofile(output_string, 1)
                CALL observable_set%write_population_tofile(output_string, 0)
                CALL observable_set%write_summary_tofile(mol_system, couplings, 20 , output_string)
                !---rhox
                CALL rhox_representation%write_rhox_tofile_xyz(output_string,  x_grid)
                !--- Tunnelprobability 
                CALL rhox_representation%write_prob_integral(output_string)
                !---Observables and Tunnel probability in Matrix Form
                CALL observable_set%put_in_observable_matrix(i)
                CALL rhox_representation%put_tunnel_probability_matrix(i)
                !Peformance for runtime estimation
                CALL observable_set%get_performance_sum(gate_grid_end)
            END DO
            
            !Writing the grid
            CALL observable_set%write_parameter_grid
            CALL mol_system%potential_surface("potential", 1)
            
            !In and Output of matrix data
            CALL observable_set%write_observable_energy_matrix_tofile_xyz("obs_matrix")
            CALL observable_set%write_observable_position_matrix_tofile_xyz("obs_matrix")
            CALL observable_set%write_observable_current_matrix_tofile_xyz("obs_matrix")
            CALL observable_set%write_observable_population_matrix_tofile_xyz("obs_matrix", 1)
            CALL observable_set%write_observable_population_matrix_tofile_xyz("obs_matrix", 0)

            !Output of heatmpas
            CALL observable_set%write_observable_current_matrix_heatmap("heatmap")
            CALL observable_set%write_observable_position_matrix_heatmap("heatmap")
            CALL observable_set%write_observable_energy_matrix_heatmap("heatmap")
            CALL observable_set%write_observable_population_matrix_heatmap("heatmap",1)
            CALL observable_set%write_observable_population_matrix_heatmap("heatmap",0)
            !State_populations 
            CALL observable_set%write_observable_state_population_matrix_heatmap("heatmap", 1)
            CALL observable_set%write_observable_state_population_matrix_heatmap("heatmap", 0)
            !Tunnel Probabilitys
            CALL rhox_representation%write_tunnel_probability_matrix("heatmap")
        END ASSOCIATE

    end subroutine run_gate_temperature_characteristics_3d

subroutine run_current_frank_characteristics_3d(input)
    CLASS(inputdata), intent(in)::input

    INTEGER                     :: i,j!< Loop Variables
    CHARACTER(len=20)           :: gate_string
    CHARACTER(len=100)          :: output_string

    CALL initialize_entire_system(input)

    !Loop over the voltages and calculation of the stationary state
    ASSOCIATE(gate_volt => input%sec_parameter_grid, &
        gate_grid_end => input%sec_parameter_grid_length,&
        bias_volt => input%parameter_grid, &
        bias_grid_end => input%parameter_grid_length)

        DO i = 1, gate_grid_end

            !Variation of the gate and new Computation of the system, leads and the couplings
            CALL mol_system%change_gate(gate_volt(i))
            CALL theleads%update(mol_system)
            CALL couplings%update(mol_system)
            CALL hbath%update_bath(mol_system)
           
            DO j = 1, bias_grid_end
                WRITE(*,*)
                WRITE(*,'(A, E10.3, I)') "Bias Voltage: ", bias_volt(j), j
                WRITE(*,'(A, E10.3, I)') "Gate Voltage: ", gate_volt(i), i
                WRITE(*,*)

                !Performance
                CALL observable_set%get_performance_time(j,"start")
                !Change voltage of the leads
                CALL theleads%change_voltage(bias_volt(j))
                CALL tensor%get_stationary_rho(density, mol_system,  theleads, couplings, hbath, input)
                !Processing rho
                !---Observables
                CALL observable_set%get_and_store_observables_mod_system(j, density, mol_system, theleads, couplings)
                CALL observable_set%get_state_population(j, density)
                !---Rho(x)
                CALL rhox_representation%update_rhox(density, mol_system, x_grid)
                CALL rhox_representation%get_and_store_rhox(j, density, x_grid, mol_system)
                !Performance
                CALL observable_set%get_performance_time(j,"end")

            END DO

            !Generate suffix for the outputfiles for better specification
            WRITE(gate_string, '(F10.3)' ) gate_volt(i)
            output_string= "gate_at_"//trim(adjustl(gate_string))

            !Writing Data Output
            !---Standard
            CALL observable_set%write_performance_tofile(output_string)
            CALL observable_set%write_results_tofile(output_string)
            CALL observable_set%write_population_tofile(output_string, 1)
            CALL observable_set%write_population_tofile(output_string, 0)
            CALL observable_set%write_summary_tofile(mol_system, couplings, 20 , output_string)
            !---rhox
            CALL rhox_representation%write_rhox_tofile_xyz(output_string,  x_grid)
            !--- Tunnelprobability 
            CALL rhox_representation%write_prob_integral(output_string)
            !---Observables and Tunnel probability in Matrix Form
            CALL observable_set%put_in_observable_matrix(i)
            CALL rhox_representation%put_tunnel_probability_matrix(i)
            !Peformance for runtime estimation
            CALL observable_set%get_performance_sum(gate_grid_end)
        END DO
        
        !Writing the grid
        CALL observable_set%write_parameter_grid
        CALL mol_system%potential_surface("potential", 1)
        
        !In and Output of matrix data
        CALL observable_set%write_observable_energy_matrix_tofile_xyz("obs_matrix")
        CALL observable_set%write_observable_position_matrix_tofile_xyz("obs_matrix")
        CALL observable_set%write_observable_current_matrix_tofile_xyz("obs_matrix")
        CALL observable_set%write_observable_population_matrix_tofile_xyz("obs_matrix", 1)
        CALL observable_set%write_observable_population_matrix_tofile_xyz("obs_matrix", 0)

        !Output of heatmpas
        CALL observable_set%write_observable_current_matrix_heatmap("heatmap")
        CALL observable_set%write_observable_position_matrix_heatmap("heatmap")
        CALL observable_set%write_observable_energy_matrix_heatmap("heatmap")
        CALL observable_set%write_observable_population_matrix_heatmap("heatmap",1)
        CALL observable_set%write_observable_population_matrix_heatmap("heatmap",0)
        !State_populations 
        CALL observable_set%write_observable_state_population_matrix_heatmap("heatmap", 1)
        CALL observable_set%write_observable_state_population_matrix_heatmap("heatmap", 0)
        !Tunnel Probabilitys
        CALL rhox_representation%write_tunnel_probability_matrix("heatmap")
    END ASSOCIATE

end subroutine run_current_frank_characteristics_3d

    SUBROUTINE test(input)

        CLASS(inputdata), intent(in)::input

        CALL hbath%init_bath(input, density, mol_system)
        CALL hbath%write_bath_functions
        CALL hbath%write_bath

    END SUBROUTINE test

end module module_run_stationary
