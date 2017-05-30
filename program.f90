program program

    USE class_inputdata
    USE module_run_stationary
    USE module_evolution

    USE module_snapshots
    IMPLICIT NONE

    REAL(8)               :: clock_program_start
    REAL(8)               :: clock_program_end
    REAL(8)               :: program_time

    type(inputdata)::input

    INTEGER :: mode !< Program mode received from inputdata object

    WRITE(*,*) '***********************************************************'
    WRITE(*,*) 'Proton Transfer in a Molecular Junction as switch mechanism'
    WRITE(*,*) 'Author: Chriszandro Hofmeister                             '
    WRITE(*,*) 'Last update: 2015/02/25                                    '
    WRITE(*,*) '***********************************************************'

    !Get the data from the files and process it
    CALL input%init_inputdata
    !Get the initial time for the calucation
    CALL CPU_time(clock_program_start)

    CALL input%get_mode(mode)
    SELECT CASE(mode)

        CASE(1) !Stationary Mode
            CALL run_current_voltage_characteristics(input)

        CASE(2) !Gate Variation Mode or "Slow switching mode"
            WRITE(*,*) "***********************"
            WRITE(*,*) "Mode 2 Current-Gate-Characteristics (Switch)"
            WRITE(*,*) "***********************"

            CALL run_gate_voltage(input)

        CASE(3)
            WRITE(*,*) "***********************"
            WRITE(*,*) "Mode 3 Current voltage characteristics in 3d"
            WRITE(*,*) "***********************"

            CALL run_current_voltage_characteristics_3d(input)

        CASE(4)

            WRITE(*,*) "***********************"
            WRITE(*,*) "Mode 4 CVC with bias dependent electric field"
            WRITE(*,*) "***********************"

            CALL run_cvc_with_bias_field(input)

        CASE(5)

            WRITE(*,*) "***************************"
            WRITE(*,*) "Mode 5 gate-temperature 3D"
            WRITE(*,*) "***************************"
           
            CALL run_gate_temperature_characteristics_3d(input)

        CASE(20) !ZVODE Evolution

            WRITE(*,*) "***********************"
            WRITE(*,*) "Mode 20 Time Evolution with ZVODE"
            WRITE(*,*) "***********************"

            CALL run_evolution_zvode(input)

        CASE(40, 30) !Expokit Evolution

            WRITE(*,*) "***********************"
            WRITE(*,*) "Mode 30 and 40 Time Evolution with Expokit"
            WRITE(*,*) "Mode 40 has a magnitude grid insteas of an equidistant" 
            WRITE(*,*) "***********************"

            CALL run_evolution_expokit(input)

        !CASE(40)

            !CALL test(input)

        CASE(50)

            WRITE(*,*) "***********************"
            WRITE(*,*) "Mode 50 Taking System Snapshot"
            WRITE(*,*) "***********************"

            CALL take_snapshot(input,0)

        CASE(42) !Secret Mode

            WRITE(*,*) "All right you figured it out. You are done with everything since &
                        you pondered the meaning of life"            

        CASE(51)

            WRITE(*,*) "***********************"
            WRITE(*,*) "Mode 51 Theoretical Validity Harmonic Bath"
            WRITE(*,*) "***********************"

        CALL theoretical_validity_hbath_test(input)

        CASE DEFAULT
    END SELECT

    CALL CPU_time(clock_program_end)
    program_time = ABS(clock_program_end - clock_program_start)
    WRITE(*,*) 'Program End'

    !Performance
    WRITE(*,'(A, E)') 'Time  in s: ',  program_time
    WRITE(*, '(A, E, E)') 'Time in min:', program_time/60.d0


end program program
