module class_inputdata
    implicit none
    private

    TYPE, public :: parameter_grid_structure
         !Parameter Grid
        INTEGER                 :: length
        REAL(8)                 :: start
        REAL(8)                 :: end
        REAL(8)                 :: step

    END TYPE

    type, public :: inputdata

        CHARACTER(len=2000) :: sInputpath                  !< Path location for the inputfile and the computation file

        CHARACTER(len=1000) :: sInputparameters             !Name of the inputfile file for physical parameters
        CHARACTER(len=1000) :: sInputparameters_copy       !A Copy: Name of the inputfile file for physical parameters
        CHARACTER(len=1000) :: sInputfile


        CHARACTER(len=1000) :: sComputationparameters       !Name of inputfile handling computational parameters
        CHARACTER(len=1000) :: sComputationparameters_copy !Copy:
        CHARACTER(len=1000) :: sComputationfile !Copy:

        INTEGER                        ::inputNumber   !Number of inputfile data

        CHARACTER(len=40):: timestamp

        !Input File
        !----------------------------------------------------------------------------------------
        !Variable Set for file management
        CHARACTER(len=1000)                           :: sInputfiles               !Path of the directory containing the input-files
        CHARACTER(len=2000)                           :: sOutputpath               !output Path

        !Potential Parameters
        REAL(8)                                    :: l              !Value V(-l)=0 and V(l)=0 for Vb = 0
        REAL(8)                                    :: delta   !Shift of the right minima
        REAL(8)                                    :: Vb                     !Value V(0)=Vb
        REAL(8)                                    :: shiftEnergy
        REAL(8)                                    :: x_shift

        !Second Set of Potential Parameters
        REAL(8)                                    :: l_2              !Value V(-l)=0 and V(l)=0 for Vb = 0
        REAL(8)                                    :: delta_2   !Shift of the right minima
        REAL(8)                                    :: Vb_2                     !Value V(0)=Vb

        !Switch Function
        REAL(8)                                    :: A                         !upper step
        REAL(8)                                    :: B                         !low step
        REAL(8)                                    :: F                         !Smoothing Factor
        REAL(8)                                    :: switchshift
        !Harmonic Potential
        REAL(8)                                    :: Omega
        REAL(8)                                    :: lambda
        REAL(8)                                    :: parabola_shift

        !---GRIDS
        !Voltage Grid
        INTEGER                                    :: grid_bias_voltage
        REAL(8)                                    :: start_bias_voltage
        REAL(8)                                    :: end_bias_voltage
        !Field Grid
        INTEGER                                    :: grid_gate_voltage
        REAL(8)                                    :: start_gate_voltage
        REAL(8)                                    :: end_gate_voltage

        !Gate
        REAL(8)                                    :: d
        !Junction
        REAL(8)                                    :: Lj
        REAL(8)                                    :: angle
        REAL(8)                                    :: trans_shift

        !Junctions
        !---Left
        REAL(8)                                    ::beta1L
        REAL(8)                                    ::beta2L
        !---Right
        REAL(8)                                    ::beta1R
        REAL(8)                                    ::beta2R
        !---Both
        REAL(8)                                    ::T
        REAL(8)                                    ::fermi_level
        !---Time Evolution
        INTEGER                                    ::time_grid
        REAL(8)                                    ::time_end
        REAL(8)                                    ::time_start
        !--- Various Parameter Grid
        INTEGER                                    :: rhox_step
        REAL(8)                                    :: parameter1
        REAL(8)                                    :: parameter2

        !---Harmonic Bath
        REAL(8)                                    :: wcut !!Cutoff Frequency
        REAL(8)                                    :: eta !!Coupling Parameter
        REAL(8)                                    :: tharmonic !!Temperatur
        !--- Swap Variable
        REAL(8)                                    :: swap
        !Computation File
        !----------------------------------------------------------------------------------------
        !Variable set for program flow control
        !---Program Mode
        INTEGER                                     :: mode !< Computation Mode.
        INTEGER                                     :: plot_bool !< Pseudo Boolean for outputfiles. 0 - no output 1 - output.
        INTEGER                                     :: rhox_bool !< Rhox Boolean if rhox output is generated. 0- no output 1 - output
        INTEGER                                     :: surf_bool !< Bolean for surf plots
        INTEGER                                     :: pop_bool !< Boolean for the populations
        INTEGER                                     :: performance_bool!< Boolean for the performace
        INTEGER                                     :: summary_bool!< Boolean for the summary
        INTEGER                                     :: coupling_bool!< Boolean for the coupling
        INTEGER                                     :: pop_number!< Boolean for the coupling

        !--- DVR
        INTEGER                                     :: N !< Number of Gridpoints in the DVR
        INTEGER                                     :: id !< id Number of potential to use in the DVR Calculations
        INTEGER                                     :: dvr_method
        REAL(8)                                     :: fourier !< genauigkeit

        !---Dimension
        INTEGER                :: medim1 !< dimension of the subspace 11
        INTEGER                :: medim0 !< dimension of the subspace 00
        INTEGER                :: meBND !< band width for the density matrix
        INTEGER                :: meBnd_small !< band width for the large density matrix for the time evolution calculations

        CHARACTER(len=20)               :: temp1 !<temporal Integer 1 for input from console
        CHARACTER(len=20)                :: temp2 !<temporal Integer 2 for input from console
        CHARACTER(len=20)                :: temp3 !<temporal Integer 3 for input from console
        CHARACTER(len=20)                :: temp4 !<temporal Integer 4 for input from console

        !---Paths
        CHARACTER(len=1000)       ::system_output_filepath!< the path for outputfiles
        CHARACTER(len=1000)       ::summary_output_filepath!< summary outpufile path
        CHARACTER(len=1000)       ::result_output_filepath

        !Bath String
        CHARACTER(len=1000)       ::bath_string

        !---Logicals
        INTEGER                 ::evo_method!< this Integer selects the Method for the time evolution.
        INTEGER                 ::solve_method !< Integer for the diagonalization method in
        REAL(8)                 ::abstol !< genauigkeit
        INTEGER                 ::ideg

        INTEGER                 ::timebool
        REAL(8)                 ::xrange!< specifies the range from -xrange to +xrange for the density matrix in position space
        INTEGER                 ::performancebool !< Pseudo Bool for an append of the performaceg
        REAL(8)                 ::rhox_tolerance !< tolerance value for rhox computations. Negative values of rho(x) whose absolute value is smaller than this value will be set to 0.
        !--- Ensembles
        INTEGER                 ::fermion_bath
        INTEGER                 ::boson_bath
        !--- Future
        INTEGER                 :: computation_modell
        INTEGER                 :: computation_modell_2

        !--- Initial State Preparation
        INTEGER                            ::initial_state !< 1 - Normal State. 2 - Pure State
        INTEGER                            ::initial_occupation
        INTEGER                            ::initial_state_number
        INTEGER                            ::initial_state_number2
        COMPLEX(8)                         :: c1_coeff
        COMPLEX(8)                         :: c2_coeff

        !Parameter Grid
        INTEGER                 :: parameter_grid_length
        REAL(8)                 :: paramter_grid_start
        REAL(8)                 :: parameter_grid_end
        REAL(8)                 :: paramter_grid_step
        CHARACTER(len=1)        :: mode_specifier

        REAL(8), ALLOCATABLE, DIMENSION(:):: parameter_grid

        !Parameter Grid
        INTEGER                 :: sec_parameter_grid_length
        REAL(8)                 :: sec_paramter_grid_start
        REAL(8)                 :: sec_parameter_grid_end
        REAL(8)                 :: sec_paramter_grid_step
        CHARACTER(len=1)        :: sec_mode_specifier

        REAL(8), ALLOCATABLE, DIMENSION(:):: sec_parameter_grid

    contains

        procedure, private :: read_inputfile
        procedure, private :: read_computation
        procedure, private :: convert_to_au
        procedure, private :: get_timestamp
        procedure, private :: get_bath_string
        procedure, public :: print_inputdata
        procedure, public :: init_inputdata
        procedure, private :: get_filepath_strings
        procedure, private :: set_grid
        procedure, private :: compute_grid
        procedure, private :: compute_sec_grid
        procedure, private :: compute_magnitude_grid
        procedure, public::get_mode

    end type inputdata

contains


     !****************************************************************************************************
     !Subroutine new_input
     !This is the constructor of the "inputdata" class.
     !Subroutine reads input from console
     !****************************************************************************************************
    subroutine init_inputdata(this)
        class(inputdata), intent(inout) :: this


        this%inputNumber = iargc()             !iarg() reading the number of inputfiles
        !Additional Parameters for control flow

        !Get the paths to the Inputfile and the Computationsfile
        CALL getarg(1,this%sInputparameters)
        CALL getarg(2, this%sComputationparameters)
        CALL getarg(3, this%sInputpath)
        CALL getarg(4, this%temp1)
        CALL getarg(5, this%temp2)
        CALL getarg(6, this%temp3)
        CALL getarg(7, this%temp4)

        !Locate where the Inputfile and the Computation File is
        this%sComputationfile = trim(adjustl(this%sInputpath))//trim(adjustl(this%sComputationparameters))
        this%sInputfile = trim(adjustl(this%sInputpath))//trim(adjustl(this%sInputparameters))

        CALL read_computation(this)
        CALL read_inputfile(this)
        CALL convert_to_au(this)
        CALL get_timestamp(this)
        CALL get_bath_string(this)
        CALL get_filepath_strings(this)
        CALL set_grid(this)

    end subroutine init_inputdata

    subroutine read_inputfile(this)
        class(inputdata), intent(inout) :: this

        OPEN(unit=100,file=this%sInputfile, status='old', action='read')  !INPUT-FILME

        READ(100, *)!Douple Well Paramters
        READ(100,'(1D20.10)') this%l                                    !
        READ(100,'(1D20.10)') this%delta                                !
        READ(100,'(1D20.10)') this%Vb                                   !
        READ(100,'(1D20.10)') this%shiftEnergy
        READ(100,'(1D20.10)') this%x_shift
        READ(100,*)!Parameters for the harmonic system
        READ(100,'(1D20.10)') this%Omega                                 !
        READ(100,'(1D20.10)') this%lambda                                 !
        READ(100,'(1D20.10)') this%parabola_shift
        READ(100,*)!Switch function parameters
        READ(100,'(1D20.10)') this%A                                    !
        READ(100,'(1D20.10)') this%B                                 !
        READ(100,'(1D20.10)') this%F                                 !
        READ(100,'(1D20.10)') this%switchShift
        READ(100,*)!Bias Voltage Grid
        READ(100,'(1I20)') this%grid_bias_voltage
        READ(100,'(1D20.10)') this%start_bias_voltage                                   !
        READ(100,'(1D20.10)') this%end_bias_voltage                                    !
        READ(100,*)!Gate Voltage Grid
        READ(100,'(1I20)') this%grid_gate_voltage
        READ(100,'(1D20.10)') this%start_gate_voltage
        READ(100,'(1D20.10)') this%end_gate_voltage
        READ(100, '(1D20.10)') this%d
        READ(100,*)!Dependance of E with the bias voltage
        READ(100,'(1D20.10)') this%Lj
        READ(100,'(1D20.10)') this%angle
        READ(100,'(1D20.10)') this%trans_shift
        READ(100,*)
        READ(100,*)!Coupling Paramter left lead
        READ(100,'(1D20.10)') this%beta1L                                 !
        READ(100,'(1D20.10)') this%beta2L
        READ(100,*)
        READ(100,*)!Coupling Parameters right lead
        READ(100,'(1D20.10)') this%beta1R                                 !
        READ(100,'(1D20.10)') this%beta2R
        READ(100,*)
        READ(100,*)!Temperature and fermi level
        READ(100,'(1D20.10)') this%T
        READ(100,'(1D20.10)') this%fermi_level
        READ(100,*)!Time Grid
        READ(100,'(1I20)') this%time_grid
        READ(100,'(1D20.10)') this%time_start
        READ(100,'(1D20.10)') this%time_end
        READ(100,*)!Rhox parameters
        READ(100,'(1I20)') this%rhox_step
        READ(100,'(1D20.10)') this%parameter1
        READ(100,'(1D20.10)') this%parameter2
        READ(100,*)!Harmonic Bath Operators
        READ(100,'(1D20.10)') this%wcut
        READ(100,'(1D20.10)') this%eta
        READ(100,'(1D20.10)') this%tharmonic
        READ(100,*)!Pure State Management
        READ(100,'(10I20)') this%initial_occupation
        READ(100,'(10I20)') this%initial_state_number
        READ(100,'(10I20)') this%initial_state_number2
        READ(100,*)!Second Parameter Set
        READ(100,'(1D20.10)') this%l_2                                    !
        READ(100,'(1D20.10)') this%delta_2                                !
        READ(100,'(1D20.10)') this%Vb_2                                   !
        READ(100,*)!Complex Coefficients for the initial pure state
!       READ(100,*) this%c1_coeff
!       READ(100,*) this%c2_coeff
        CLOSE(100)
    end subroutine read_inputfile


    SUBROUTINE convert_to_au(this)
        USE UnitConv

        class(inputdata), intent(inout) :: this

        !Unit Conversion
        !--1st Set: Double Well Potential
        this%l = this%l*ang2au
        this%delta = this%delta*ev2au
        this%Vb = this%Vb*ev2au

        !--2nd Set: Double Well Potential
        this%l_2 = this%l_2*ang2au
        this%delta_2 = this%delta_2*ev2au
        this%Vb_2 = this%Vb_2*ev2au

        this%shiftEnergy = this%shiftEnergy*ev2au
        !---Parabolas
        this%Omega = this%Omega*ev2au
        this%lambda = this%lambda*ev2au
        this%parabola_shift = this%parabola_shift*ev2au

        !--- Electric Field of the Gate Voltage
        !-- conversion in nano meters
        this%d = this%d * 1d-9

        !--- The Leads
        !---Left
        this%beta1L = this%beta1L * ev2au
        this%beta2L = this%beta2L * ev2au
        !---Right
        this%beta1R = this%beta1R * ev2au
        this%beta2R = this%beta2R * ev2au
        !---Left + Right
        this%T = this%T * kb
        this%fermi_level = this%fermi_level * ev2au

        !--- Electric Field of the leads
        !--- for bias dependent calculations
        this%angle = this%angle*D_to_R
        this%Lj = this%Lj*1d-9
        this%trans_shift = this%trans_shift*ang2au

        !--- Harmonic Bath
        this%wcut = this%wcut*ev2au
        this%Tharmonic = this%Tharmonic * kb

    END SUBROUTINE convert_to_au


    subroutine read_computation(this)
        class(inputdata), intent(inout) :: this

        OPEN(unit=101,file=this%sComputationfile, status='old', action='read')  !INPUT-FILME

        !System
        READ(101,*)
        READ(101,'(4I20)') this%mode
        READ(101,'(4I20)') this%plot_bool
        READ(101,'(4I20)') this%rhox_bool
        READ(101,'(4I20)') this%surf_bool
        READ(101,'(4I20)') this%pop_bool
        READ(101,'(4I20)') this%summary_bool
        READ(101,'(4I20)') this%performance_bool
        READ(101,'(4I20)') this%coupling_bool
        READ(101,'(4I20)') this%pop_number
        READ(101,*)                                                   !program-mode
        READ(101,'(10I20)') this%N
        READ(101,'(10I20)') this%dvr_method
        READ(101,'(1D20.10)') this%fourier
        READ(101,'(4I20)') this%id
        READ(101,*)
        READ(101,'(10I20)') this%medim0
        READ(101,'(10I20)') this%medim1
        READ(101,'(10I20)') this%meBnd
        READ(101,'(10I20)') this%meBnd_small
        READ(101, *)
        READ(101, *)
        READ(101,'(A)') this%system_output_filepath
        READ(101, *)
        READ(101,'(A)') this%summary_output_filepath
        READ(101, *)
        READ(101,'(A)') this%result_output_filepath
        READ(101, *)
        READ(101,'(10I20)') this%evo_method
        READ(101,'(10I20)') this%ideg !!< For Pade Approximation.
        READ(101, *)
        READ(101,'(4I20)') this%solve_method !!< Method switch for Diagonalization of the Schroedinger Eqauation in systemdvr.f90 module
        READ(101,'(1D20.10)') this%abstol
        READ(101, *)
        READ(101,'(4I20)') this%timebool !!<Boolean for timestamp.
        READ(101,'(1D20.10)') this%xrange !!<Range parameter for the position
        READ(101,'(4I20)') this%performancebool!!Boolean for timestamp.
        READ(101,'(1D20.10)') this%rhox_tolerance!!Boolean for timestamp.
        READ(101, *)
        READ(101,'(10I20)') this%fermion_bath
        READ(101,'(10I20)') this%boson_bath
        READ(101, *)
        READ(101,'(10I20)') this%initial_state
        READ(101, *)
        READ(101,'(10I20)') this%computation_modell
        READ(101,'(10I20)') this%computation_modell_2
        CLOSE(101)

    end subroutine read_computation

    subroutine get_timestamp(this)

        USE timestamper
        class(inputdata), intent(inout) :: this

        CALL timestring(this%timestamp)


    end subroutine get_timestamp

    subroutine get_bath_string(self)

        class(inputdata), intent(inout) ::self

        CHARACTER(len=10) bathstring
        CHARACTER(len=10) bosonstring
        CHARACTER(len=10) fermionstring

        bosonstring = ""
        fermionstring = ""

        IF(self%fermion_bath.eq.1) THEN
            fermionstring = "Fe"
        END IF

        IF(self%boson_bath.eq.1) THEN
            bosonstring = "Bo"
        END IF

        self%bath_string = "_"//adjustl(trim(fermionstring))//adjustl(trim(bosonstring))//"_"

    end subroutine



    subroutine get_filepath_strings(self)

        class(inputdata), intent(inout) ::self

        CHARACTER(len=2000) swap_system_output_filepath
        CHARACTER(len=2000) swap_summary_output_filepath
        CHARACTER(len=2000) swap_result_output_filepath

        CHARACTER(len=2000) swaptime_system_output_filepath
        CHARACTER(len=2000) swaptime_summary_output_filepath
        CHARACTER(len=2000) swaptime_result_output_filepath

        !Generate the the file paths + prefixes for the different paths
        swap_system_output_filepath = trim(adjustl(self%system_output_filepath))//"/"//trim(adjustl(self%sInputparameters))//trim(adjustl(self%bath_string))
        swap_summary_output_filepath = trim(adjustl(self%summary_output_filepath))//"/"//trim(adjustl(self%sInputparameters))//trim(adjustl(self%bath_string))
        swap_result_output_filepath = trim(adjustl(self%result_output_filepath))//"/"//trim(adjustl(self%sInputparameters))//trim(adjustl(self%bath_string))

        !Generate the filepaths + prefixes (inputfilename + timestamp)
        swaptime_system_output_filepath = trim(adjustl(swap_system_output_filepath))//"_"//trim(adjustl(self%timestamp))
        swaptime_summary_output_filepath = trim(adjustl(swap_summary_output_filepath))//"_"//trim(adjustl(self%timestamp))
        swaptime_result_output_filepath = trim(adjustl(swap_result_output_filepath))//"_"//trim(adjustl(self%timestamp))

        IF(self%timebool.eq.0) THEN

            self%system_output_filepath = swap_system_output_filepath
            self%summary_output_filepath = swap_summary_output_filepath
            self%result_output_filepath = swap_result_output_filepath

        ELSE

            self%system_output_filepath = swaptime_system_output_filepath
            self%summary_output_filepath = swaptime_summary_output_filepath
            self%result_output_filepath = swaptime_result_output_filepath

        END IF

    end subroutine get_filepath_strings

    subroutine check_validity_of_data(self)

        CLASS(inputdata)::self
        INTEGER         :: i

        IF(self%mode.eq.4) THEN
            self%id = 10
            WRITE(*,*) "Selected Mode 4. Potential ID is set to id=10"
        END IF


    end subroutine check_validity_of_data

    SUBROUTINE set_grid(self)
        USE UnitConv
        CLASS(inputdata)::self

        !Parameters
        REAL(8)                     :: diff

        !Loop
        INTEGER                     :: i

        SELECT CASE(self%mode)

            
            CASE(20, 30, 21) !Time evolution Mode with normal time-grid

                CALL self%compute_grid(self%time_start, self%time_end, self%time_grid)
                CALL self%compute_sec_grid(self%start_gate_voltage, self%end_gate_voltage, self%grid_gate_voltage)
                self%mode_specifier = "T"
                !For Heatmaps    
                IF(self%evo_method.eq.4) THEN 
                        self%mode_specifier = "K"
                END IF

            CASE(40) !Time evolution Mode with a magnitude grid calculated with expokit 
                CALL self%compute_magnitude_grid(self%time_grid)
                self%mode_specifier = "U"

            CASE(1) !Current-Voltage-Characteristics
                CALL self%compute_grid(self%start_bias_voltage, self%end_bias_voltage, self%grid_bias_voltage)
                self%mode_specifier = "V"

            CASE(2) !Current-Gate-Characteristics with heatmap
                self%mode_specifier = "W"
                CALL self%compute_sec_grid(self%start_bias_voltage, self%end_bias_voltage, self%grid_bias_voltage)
                CALL self%compute_grid(self%start_gate_voltage, self%end_gate_voltage, self%grid_gate_voltage)

            CASE(3) !Current Gate Characteristics with heatmap in 3D
                CALL self%compute_grid(self%start_bias_voltage, self%end_bias_voltage, self%grid_bias_voltage)
                CALL self%compute_sec_grid(self%start_gate_voltage, self%end_gate_voltage, self%grid_gate_voltage)
                self%mode_specifier = "G"

            CASE(4) !CVC with biased Field with heatmap

                CALL self%compute_grid(self%start_bias_voltage, self%end_bias_voltage, self%grid_bias_voltage)
                CALL self%compute_sec_grid(self%start_gate_voltage, self%end_gate_voltage, self%grid_gate_voltage)
                self%mode_specifier = "B"

                WRITE(*,*) "Computation Mode 4 is chosen. The potential is set to id = 10, ignoring the id given in the computation file to prevent wrong calculations!"
                self%id=10

            CASE(5) !Gate-Temperature Heatmap 
                !Temperature and bias voltage swap 
                self%swap = self%T
                self%Tharmonic = self%start_bias_voltage* kb
                self%T = self%start_bias_voltage* kb
                self%start_bias_voltage = self%swap

                CALL self%compute_grid(self%start_bias_voltage, self%end_bias_voltage, self%grid_bias_voltage)
                CALL self%compute_sec_grid(self%start_gate_voltage, self%end_gate_voltage, self%grid_gate_voltage)
                self%mode_specifier = "P"
                WRITE(*,*) "Computation Mode 5 is chosen! Gate-Temperature-Heatmap"
        
            END SELECT

    end subroutine set_grid


    subroutine compute_magnitude_grid(self, length)

        CLASS(inputdata)::self
        !Calculating the magnitude Grid
        REAL(8) :: mantisse
        REAL(8) :: division
        REAL(8) :: expo

        INTEGER, intent(in) :: length

        INTEGER :: sub_length
        INTEGER :: start_sub
        INTEGER :: allocation_counter
        REAL(8) :: refinement
        REAL(8) :: grid_point
        REAL(8) :: end_first_grid

        ! Loop 
        INTEGER             :: i, j


        !DO i = 1, self%parameter_grid_length

            !division = REAL(i,8)/10.d0
            !expo = REAL(int(division))
            !mantisse = mod( REAL(i,8), 10.d0)

            !self%parameter_grid(i) = mantisse*10**expo

            !IF (mantisse.eq.0.0d0) THEN 
                !self%parameter_grid(i) = (0.95d0)*10.d0**expo
            !END IF

            !!Debugging Output
            !!WRITE(*,*) i, self%parameter_grid(i)
            !!WRITE(*,*) division, expo, mantisse
            !!WRITE(*,*) 10**expo

        !END DO
        
        !self%parameter_grid_length = length
        refinement = 0.10d0
        start_sub = int(1.0d0 / refinement)
        sub_length = int(10.0d0/refinement) - 1 
        end_first_grid = 80

        !Calculate Dimension for Allocation | Quick and Dirty
        allocation_counter = 0
        WRITE(*,*) allocation_counter

        DO i = 1, end_first_grid 
                allocation_counter = allocation_counter + 1
        END DO

        DO i = 8, 14 
                DO j = start_sub, sub_length 
                        allocation_counter = allocation_counter + 1
                END DO
        END DO

        WRITE(*,*) allocation_counter
        
        self%parameter_grid_length = allocation_counter 
        ALLOCATE(self%parameter_grid(self%parameter_grid_length))

        WRITE(*,*) 'First Grid'
        allocation_counter = 0
        !1-80 Loop
        DO i = 1, end_first_grid 
                
            allocation_counter = allocation_counter + 1

            division = REAL(i,8)/10.d0
            expo = REAL(int(division))
            mantisse = mod( REAL(i,8), 10.d0)

            self%parameter_grid(i) = mantisse*10**expo

            IF (mantisse.eq.0.0d0) THEN 
                self%parameter_grid(i) = (0.95d0)*10.d0**expo
            END IF

            !Debugging Output
            WRITE(*,'(I,E)') i, self%parameter_grid(i)
            !WRITE(*,*) division, expo, mantisse
            !WRITE(*,*) 10**expo

        END DO
         
        WRITE(*,*) 'Second Grid' 
        DO i = 8, 14 
            
                expo = REAL(i,8)

                DO j = start_sub, sub_length 
                    
                    allocation_counter = allocation_counter + 1

                    grid_point = (REAL(j,8)*refinement)*10**expo

                    self%parameter_grid(allocation_counter) = grid_point

                    WRITE(*,'(I,E)') allocation_counter, self%parameter_grid(allocation_counter)
                END DO
        END DO


    end subroutine compute_magnitude_grid


    subroutine compute_grid(self, start, end, length)

        CLASS(inputdata)::self
        REAL(8), intent(in) :: start
        REAL(8), intent(in) :: end
        INTEGER :: length

        REAL(8)             :: diff
        INTEGER             :: i
        
        !This +1 is maybe to much since in the allocation the arraylength is exteded by 1 as well
        self%parameter_grid_length = length + 1
        self%paramter_grid_start = start
        self%parameter_grid_end = end

        diff = ABS(start - end)
        self%paramter_grid_step = diff/ REAL(length,8)

        ALLOCATE(self%parameter_grid(self%parameter_grid_length + 1 ))

        DO i = 1, self%parameter_grid_length

            self%parameter_grid(i) =   self%paramter_grid_start + (i-1)*self%paramter_grid_step

        END DO

    end subroutine compute_grid

    subroutine compute_sec_grid(self, start, end, length)

        CLASS(inputdata)::self
        REAL(8), intent(in) :: start
        REAL(8), intent(in) :: end
        INTEGER :: length

        REAL(8)             :: diff
        INTEGER             :: i

        self%sec_parameter_grid_length = length + 1
        self%sec_paramter_grid_start = start
        self%sec_parameter_grid_end = end

        diff = ABS(start - end)
        self%sec_paramter_grid_step = diff/ REAL(length,8)

        ALLOCATE(self%sec_parameter_grid(self%sec_parameter_grid_length + 1 ))

        DO i = 1, self%sec_parameter_grid_length

            self%sec_parameter_grid(i) =   self%sec_paramter_grid_start + (i-1)*self%sec_paramter_grid_step

        END DO

    end subroutine compute_sec_grid


    subroutine get_mode(self, mode)

        CLASS(inputdata), intent(in):: self
        INTEGER, intent(out)::mode

        mode = self%mode

    end subroutine

    subroutine print_inputdata(this)
        class(inputdata), intent(inout) :: this

        WRITE(*,*) "ComputatioN: ", this%sComputationfile
        WRITE(*,*) "Mode: ", this%mode
        WRITE(*,*) "Timestamp: ", this%timestamp
    end subroutine

end module class_inputdata
