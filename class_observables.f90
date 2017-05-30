module class_observables

    USE class_inputdata
    USE class_grid
    USE class_junction
    USE class_density
    USE class_system
    USE class_coupling

    implicit none
    private

    type, public :: observables
        !Copy of the primary paramter grid
        REAL(8), ALLOCATABLE, DIMENSION(:)        :: parameter_grid !< this is the paramter the calculating loop varies
        INTEGER                                   :: parameter_grid_length

        !Copy of the primary paramter grid
        REAL(8), ALLOCATABLE, DIMENSION(:)        :: sec_parameter_grid !< this is the paramter the calculating loop varies
        INTEGER                                   :: sec_parameter_grid_length

        !Observables
        REAL(8), ALLOCATABLE, DIMENSION(:)        :: energy
        REAL(8), ALLOCATABLE, DIMENSION(:)        :: current
        REAL(8), ALLOCATABLE, DIMENSION(:,:)      :: bridge_population
        REAL(8), ALLOCATABLE, DIMENSION(:)        :: position
        REAL(8), ALLOCATABLE, DIMENSION(:,:, :)   :: state_population
        REAL(8), ALLOCATABLE, DIMENSION(:,:)      :: sub_matrix
        REAL(8), ALLOCATABLE, DIMENSION(:,:, :)   :: state_energys
        
        !Validity Hbath Matrix  
        REAL(8), ALLOCATABLE, DIMENSION(:)        :: validity_matrix


        !Heatmaps
        !---Observable Matrices
        REAL(8), ALLOCATABLE, DIMENSION(:,:)        :: energy_matrix
        REAL(8), ALLOCATABLE, DIMENSION(:,:)        :: current_matrix
        REAL(8), ALLOCATABLE, DIMENSION(:,:,:)      :: bridge_population_matrix
        REAL(8), ALLOCATABLE, DIMENSION(:,:)        :: position_matrix
        !---Populations Matrix
        REAL(8), ALLOCATABLE, DIMENSION(:,:,:)      :: state_occupation_pop_matrix
        REAL(8), ALLOCATABLE, DIMENSION(:,:,:)      :: state_unoccupation_pop_matrix
        INTEGER                                     :: state_population_number=10
        !---Frank Matrix
        REAL(8), ALLOCATABLE, DIMENSION(:,:,:)      :: frank_matrix
        INTEGER                                     :: frank_number=25
        INTEGER                                     :: frank_matrix_elements
        !---<n|x|m> Matrix
        REAL(8), ALLOCATABLE, DIMENSION(:,:,:)      :: x_energy_matrix
        INTEGER                                     :: x_energy_number=25
        INTEGER                                     :: x_energy_matrix_elements

        
        !Auxilarys
        COMPLEX(8), ALLOCATABLE, DIMENSION(:)     :: x_auxil_vector
        COMPLEX(8), ALLOCATABLE, DIMENSION(:)     :: current_auxil_matrix

        REAL(4), ALLOCATABLE, DIMENSION(:)        :: performance
        INTEGER                                   :: system_clock_start
        INTEGER                                   :: system_clock_end
        INTEGER                                   :: system_clock_rate

        CHARACTER(len=1000)                       :: sOutput
        CHARACTER(len=1000)                       :: sOutput_summary
        LOGICAL                                   :: init_flag=.false.
        INTEGER                                     :: plot_bool !< Pseudo Boolean for outputfiles. 0 - no output 1 - output.
        INTEGER                                     :: pop_bool !< Pseudo Boolean for population outputfiles. 0 - no output 1 - output.
        INTEGER                                     :: pop_number !< Number of States

        INTEGER                                     :: performance_bool !< Pseudo Boolean for performance outputfiles. 0 - no output 1 - output.
        INTEGER                                     :: summary_bool !< Pseudo Boolean for outputfiles. 0 - no output 1 - output.



        CHARACTER                                 :: mode_specifier
    contains
        procedure, public:: get_hbath_validity
        procedure, private:: get_energy
        procedure, private:: get_bridge_occupation
        procedure, private:: get_position
        procedure, public:: get_and_store_observables
        procedure, public:: init_observables
        procedure, public:: get_state_population
        procedure, public:: get_state_energy
        procedure, public:: get_sub_matrix
        procedure, private:: get_current
        procedure, public:: write_results_tofile
        procedure, public:: write_state_energy_tofile
        procedure, public::write_population_tofile
        procedure, public::get_performance_time
        procedure, public::write_summary_tofile
        procedure, public:: get_and_store_observables_mod_system
        procedure, public:: set_X_auxil_vector
        procedure, private:: set_current_auxil_matrix
        procedure, public:: write_performance_tofile
        procedure, public:: get_and_putin_array
        procedure, public:: put_console_obs
        procedure, public:: put_in_frank_matrix
        procedure, public:: put_in_observable_matrix
        procedure, public:: get_performance_sum
        procedure, public:: write_observable_energy_matrix_tofile_xyz
        procedure, public:: write_observable_position_matrix_tofile_xyz
        procedure, public:: write_observable_current_matrix_tofile_xyz
        procedure, public:: write_observable_population_matrix_tofile_xyz
        procedure, public:: write_observable_state_population_matrix_heatmap
        procedure, public:: write_observable_current_matrix_heatmap
        procedure, public:: write_observable_position_matrix_heatmap
        procedure, public:: write_observable_energy_matrix_heatmap
        procedure, public:: write_observable_population_matrix_heatmap
        procedure, public:: write_parameter_grid
        procedure, public:: write_index_to_console
        procedure, public:: write_sub_matrix 
        procedure, public:: write_frank_matrix
        procedure, public:: write_population_histo_tofile   
        procedure, public:: write_x_energy_matrix 
        procedure, public:: write_hbath_validity
        procedure, public:: put_in_x_energy_matrix 
        end type observables

contains

    subroutine put_console_obs(self, indexx)
        CLASS(observables), intent(in)::self
        INTEGER, intent(in) :: indexx

        WRITE(*,'(A, T15, I3, T35, A, T50, E10.3, T70, A, E10.3, A, E10.3, A)'), "Time (a.u) : ", indexx , "Time (SI): ", indexx*au2second
        WRITE(*,'(A, T15, E10.3, F10.3, A, T45, F10.3, A, T60, F10.6, A, T80, F10.6, T110, F10.6, A, F10.6)') "Position", self%position(indexx)
        WRITE(*,'(A, T15, E10.3, F10.3, A, T45, F10.3, A, T60, F10.6, A, T80, F10.6, T110, F10.6, A, F10.6)') "Energy", self%energy(indexx)
        WRITE(*,'(A, T15, E10.3, F10.3, A, T45, F10.3, A, T60, F10.6, A, T80, F10.6, T110, F10.6, A, F10.6)') "Current", self%current(indexx)
        WRITE(*,'(A, T15, E10.3, F10.3, A, T45, F10.3, A, T60, F10.6, A, T80, F10.6, T110, F10.6, A, F10.6)') "Occupation", self%bridge_population(indexx,0)

    end subroutine put_console_obs


    subroutine init_observables(self, input, density, sys)

        class(observables), intent(inout):: self
        class(inputdata), intent(in):: input
        class(densityvector), intent(in)::density
        class(system), intent(in)::sys

        WRITE(*,*), "Observables Arrays Initialization"

        IF(self%init_flag.eq..false.) THEN

            !Population Number. 
            self%pop_number = input%pop_number

            !Outputfile Bools
            self%plot_bool = input%plot_bool
            self%pop_bool = input%pop_bool
            self%summary_bool = input%summary_bool
            self%performance_bool = input%performance_bool

            WRITE(*,*), "self%paramter_grid"
            ALLOCATE(self%parameter_grid(input%parameter_grid_length)) !< this is the paramter the calculating loop varies
            WRITE(*,*), "self%energy"
            ALLOCATE(self%energy(input%parameter_grid_length))
            WRITE(*,*), "self%current"
            ALLOCATE(self%current(input%parameter_grid_length))
            WRITE(*,*), "self%position"
            ALLOCATE(self%position(input%parameter_grid_length))
            WRITE(*,*), "self%bridge"
            ALLOCATE(self%bridge_population(input%parameter_grid_length, 0:1))
            WRITE(*,*), "self_sub%matrix"
            ALLOCATE(self%sub_matrix(input%parameter_grid_length, 0:3))
            WRITE(*,*), "self_state%population"
            ALLOCATE(self%state_population(1:self%pop_number, input%parameter_grid_length, 0:1))
            WRITE(*,*), "state_energys"
            ALLOCATE(self%state_energys(1:self%pop_number, input%parameter_grid_length, 0:1))
            WRITE(*,*), "performance"
            ALLOCATE(self%performance(1:input%parameter_grid_length))
            WRITE(*,*), "auxil_vector"
            ALLOCATE(self%x_auxil_vector(1:density%medim))

            !Get primary Grid Copy
            self%parameter_grid_length = input%parameter_grid_length
            self%parameter_grid = input%parameter_grid

            !Get the folder and the prefix for all possible outputfiles
            self%sOutput = input%result_output_filepath
            self%sOutput_summary = input%summary_output_filepath

            !Get the mode of the computation for further processing
            self%mode_specifier = input%mode_specifier
            WRITE(*,'(A, A)') "Specifier is ", self%mode_specifier 

            !Auxilarys
            !---Compute the auxilary vector for the calculation of <x>
            CALL self%set_X_auxil_vector(sys, density)


            IF((self%mode_specifier.eq."W").or.(self%mode_specifier.eq."B").or.(self%mode_specifier.eq."G").or.(self%mode_specifier.eq."K")) THEN

                !Get second Grid
                self%sec_parameter_grid_length = input%sec_parameter_grid_length

                ALLOCATE(self%sec_parameter_grid(self%sec_parameter_grid_length)) !< this is the paramter the calculating loop varies
                self%sec_parameter_grid = input%sec_parameter_grid

                !Allocte the observable matrices
                ALLOCATE(self%energy_matrix(self%parameter_grid_length, self%sec_parameter_grid_length))
                ALLOCATE(self%position_matrix(self%parameter_grid_length, self%sec_parameter_grid_length))
                ALLOCATE(self%current_matrix(self%parameter_grid_length, self%sec_parameter_grid_length))
                ALLOCATE(self%bridge_population_matrix(self%parameter_grid_length, self%sec_parameter_grid_length, 0:1))
                !State Populations 
                ALLOCATE(self%state_occupation_pop_matrix(self%parameter_grid_length, self%sec_parameter_grid_length, 0:self%state_population_number))
                ALLOCATE(self%state_unoccupation_pop_matrix(self%parameter_grid_length, self%sec_parameter_grid_length, 0:self%state_population_number))
                !Frank Matrix
                ALLOCATE(self%frank_matrix(self%frank_number, self%frank_number, self%parameter_grid_length))
                !<n|x|m> Matrix
                ALLOCATE(self%x_energy_matrix(self%x_energy_number, self%x_energy_number, self%parameter_grid_length))
                !Validity Matrix
                ALLOCATE(self%validity_matrix(self%parameter_grid_length))
            END IF

            self%init_flag = .true.

        ELSE

            WRITE(*,'(A)') "Observables already initialized yet"

        END IF

    end subroutine init_observables


    subroutine get_and_store_observables_mod_system(self, indexx, density, sys, junction, coup)
        class(observables), intent(inout):: self
        class(system), intent(in)::sys
        class(leads), intent(in)::junction
        class(densityvector), intent(inout)::density
        class(coupling), intent(in)::coup

        INTEGER, intent(in) :: indexx

        !Auxilarys
        CALL self%set_X_auxil_vector(sys, density)

        CALL self%get_and_store_observables(indexx, density, sys, junction, coup)

    end subroutine get_and_store_observables_mod_system


    subroutine get_and_putin_array(self, prep, status, sys, junction, density, coup)

        class(observables), intent(inout):: self
        class(system), intent(in)::sys
        class(leads), intent(in)::junction
        class(densityvector), intent(inout)::density
        class(coupling), intent(in)::coup


        REAL(8), intent(out)       :: prep(1:, 1: )
        INTEGER, intent(in)        :: status

        !Auxilarys
        CALL self%set_X_auxil_vector(sys, density)

        prep(status, 1) = self%get_position(sys, density)
        prep(status, 2)= self%get_current(density, junction, sys, coup)
        prep(status, 3) = self%get_energy(sys, density)
        prep(status, 4) = self%get_bridge_occupation(density, 1)

    end subroutine

    subroutine get_and_store_observables(self, indexx, density, sys, junction, coup)

        class(observables), intent(inout):: self
        class(system), intent(in)::sys
        class(leads), intent(in)::junction
        class(densityvector), intent(inout)::density
        class(coupling), intent(in)::coup

        INTEGER, intent(in) :: indexx

        self%energy(indexx) = self%get_energy(sys, density)
        self%position(indexx) = self%get_position(sys, density)
        self%current(indexx) = self%get_current(density, junction, sys, coup)

        self%bridge_population(indexx, 0) = self%get_bridge_occupation(density, 0)
        self%bridge_population(indexx, 1) = self%get_bridge_occupation(density, 1)

    end subroutine get_and_store_observables

    REAL(8) FUNCTION get_energy(self, sys, density)

        IMPLICIT NONE
        CLASS(observables)::self
        CLASS(system):: sys
        CLASS(densityvector)::density

        !Summation
        COMPLEX(8)                                      :: summation
        !Loop
        INTEGER             ::i,j

        ASSOCIATE( all_m => density%NumDim00,              &
            all_v => density%NumDim11,              &
            diag_m => density%medim00,          &
            diag_v => density%medim11)

            summation = (0.0d0,0.0d0)

            !--Diagonaldensityments in 00
            !           DO j = 1, diag_m
            !               summation = summation + density%rho(j)*sys%hsEN(j,0)
            !           END DO

            summation = dot_product(density%rho(1:diag_m), sys%hsEN(1:diag_m,0))
            !--Diagonalelements in 11
            !    DO j = 1, diag_v
            !        summation = summation + density%rho(j+all_m)*sys%hsEN(j,1)
            !    END DO

            summation = summation + dot_product( density%rho((all_m + 1):(all_m + diag_v)), &
                sys%hsEN(1:diag_v, 1))

            get_energy = REAL(summation, 8)*au2ev

        END ASSOCIATE
    END FUNCTION get_energy

    REAL(8) FUNCTION get_bridge_occupation(self, density, occupation_state) RESULT(occupation)

        IMPLICIT NONE
        CLASS(densityvector), intent(in)              :: density
        CLASS(observables), intent(in)                :: self
        INTEGER, intent(in)                             :: occupation_state

        !Summation
        COMPLEX(8)                                      :: summation
        !Loop
        INTEGER             ::m,v

        ASSOCIATE(all_m => density%NumDim00,          &
            all_v => density%NumDim11,          &
            diag_m => density%medim00,          &
            diag_v => density%medim11)

            SELECT CASE (occupation_state)

                CASE(0)

                    !--Diagonalelements in 00
                    !               DO m = 1, diag_m
                    !                   summation = summation + density%rho(m)
                    !               END DO

                    summation = SUM(density%rho(1:diag_m))

                    occupation = REAL(summation, 8)

                    RETURN

                CASE(1)

                    !--Diagonalelements in 11
                    !                   DO v = 1, diag_v
                    !                       summation = summation + density%rho(v+all_m)
                    !                   END DO

                    summation = SUM(density%rho((all_m  + 1 ):(all_m + diag_v)))
                    occupation = REAL(summation, 8)

                    RETURN

                !Failure. Occupation is not 0 or 1
                CASE DEFAULT

                    WRITE(*,'(A)') "Wrong occupation in occupation Function. 0 and 1 are allowed!"
                    WRITE(*,*) occupation
                    RETURN

            END SELECT

        END ASSOCIATE

    END FUNCTION get_bridge_occupation

    REAL(8) FUNCTION get_position(self, sys, density)


        CLASS(observables)::self
        CLASS(densityvector)::density
        CLASS(system)::sys

        !Summation
        COMPLEX(8)                                           :: summation
        !Loop
        INTEGER             ::i,j


        COMPLEX(8), external  ::zdotu

        ASSOCIATE(  all_m => density%NumDim00,              &
            all_v => density%NumDim11,              &
            diag_m => density%medim00,          &
            diag_v => density%medim11,          &

            !Indize
            m1 => density%tulpVector(:,1,0),        &
            m2 => density%tulpVector(:,2,0),        &
            v1 => density%tulpVector(:,1,1),        &
            v2 => density%tulpVector(:,2,1))

            summation = 0.0d0

            summation = zdotu(density%medim ,density%rho,1,self%x_auxil_vector, 1)

            get_position = REAL(summation, 8)*au2ang

        END ASSOCIATE

    END FUNCTION get_position

    SUBROUTINE get_state_population(self, indexx,  density)

        CLASS(densityvector), intent(in) :: density
        CLASS(observables), intent(inout) :: self

        INTEGER indexx

        !Grid Point Variable
        INTEGER i,j

        ASSOCIATE( all_m => density%NumDim00,              &
            all_v => density%NumDim11,              &
            diag_m => density%medim00,          &
            diag_v => density%medim11)

            !Old Version
            !self%state_population(1:diag_m, indexx, 0 ) = density%rho(1:diag_m)
            !self%state_population(1:diag_v, indexx, 1) =density%rho((all_m + 1 ): (all_m + diag_v))
            
            self%state_population(1:self%pop_number, indexx, 0 ) = density%rho(1:self%pop_number)
            self%state_population(1:self%pop_number, indexx, 1) =density%rho((all_m + 1 ): (all_m + self%pop_number))

        END ASSOCIATE

    END SUBROUTINE get_state_population
    
    SUBROUTINE get_sub_matrix(self, indexx,  density)

        CLASS(densityvector), intent(in) :: density
        CLASS(observables), intent(inout) :: self

        INTEGER indexx

        !Grid Point Variable
        INTEGER i,j

        ASSOCIATE( all_m => density%NumDim00,              &
            all_v => density%NumDim11,              &
            diag_m => density%medim00,          &
            diag_v => density%medim11)
            
            self%sub_matrix(indexx, 0) = density%rho(density%state1)   
            self%sub_matrix(indexx, 1) = density%rho(density%state2)   
            self%sub_matrix(indexx, 2) = REAL(density%rho(density%index_of_coherences), 8)
            self%sub_matrix(indexx, 3) = DIMAG(density%rho(density%index_of_coherences))

        END ASSOCIATE

    END SUBROUTINE get_sub_matrix

    SUBROUTINE get_state_energy(self, indexx, sys, density)

        CLASS(observables), intent(inout) :: self
        CLASS(densityvector), intent(in) :: density
        CLASS(system), intent(in)::sys
        INTEGER, intent(in) :: indexx

        ASSOCIATE( all_m => density%NumDim00,              &
            all_v => density%NumDim11,              &
            diag_m => density%medim00,          &
            diag_v => density%medim11)

            self%state_energys(1:diag_m, indexx, 0) = sys%hsEN(1:diag_m,0)*au2ev 
            self%state_energys(1:diag_v, indexx, 1) =  sys%hsEN(1:diag_v,1)*au2ev

        END ASSOCIATE

    END SUBROUTINE get_state_energy

    REAL(8) FUNCTION get_current(self, density, junction, sys, coup)

        IMPLICIT NONE
        class(observables), intent(inout):: self
        class(system), intent(in)::sys
        class(leads), intent(in)::junction
        class(densityvector), intent(inout)::density
        class(coupling), intent(in)::coup

        COMPLEX(8)                                      :: summation

        !Loop
        INTEGER             ::i,j
        !Loop
        INTEGER             ::m,v

        ASSOCIATE(  all_m => density%NumDim00,              &
            all_v => density%NumDim11,              &

            diag_m => density%medim00,          &
            diag_v => density%medim11,          &

            !Indize
            m1 => density%tulpVector(:,1,0),        &
            m2 => density%tulpVector(:,2,0),        &
            v1 => density%tulpVector(:,1,1),        &
            v2 => density%tulpVector(:,2,1), &
            jnSEM => junction%jnSEM, &
            jnVDM => coup%jnVDM)


            summation = (0.0d0, 0.0d0)
            !Calculation Current <I>
            !--Diagonalelements in 00

            DO i = 1, all_m
                DO v = 1, diag_v
                    summation = summation + jnSEM(m1(i), v, 1, 1)*dcmplx(jnVdM(m1(i),v,1)*jnVdM(m2(i),v, 1), 0.0d0) * density%rho(i)
                END DO
            END DO


            !--Diagonalelements in 11

            DO i = 1, all_v
                DO m = 1, diag_m
                    summation = summation + dconjg(jnSEM(m,v1(i), 1, 2))*dcmplx(jnVdM(m,v1(i),1)*jnVdM(m ,v2(i), 1), 0.0d0) * density%rho(i + all_m)
                END DO
            END DO

            get_current = 2.0d0 * dimag(summation)* au2mua

        END ASSOCIATE
    END FUNCTION get_current

    subroutine get_performance_time(self, indexx, start_stop)

        USE ifport

        CLASS(observables), intent(inout)   ::self
        CHARACTER(len=*), intent(in)        ::start_stop
        INTEGER, intent(in)                 ::indexx

        REAL(4)                             ::elapsed_time

        SELECT CASE(start_stop)
            CASE("start")

                CALL system_clock(self%system_clock_start, self%system_clock_rate)

            CASE("end")

                CALL system_clock(self%system_clock_end)
                elapsed_time = REAL((self%system_clock_end-self%system_clock_start), 4)/REAL(self%system_clock_rate,4)

                self%performance(indexx) = elapsed_time
                WRITE(*,'(A, E10.3, A, I)') "Performance time: ",  elapsed_time, " in point ", indexx
            CASE DEFAULT


        END SELECT

    end subroutine get_performance_time


    subroutine get_hbath_validity(self, test_pass, indexx)

        CLASS(observables), intent(inout)   ::self
        INTEGER, intent(in)                 ::test_pass
        INTEGER, intent(in)                 ::indexx

        self%validity_matrix(indexx) = REAL(test_pass, 8)

    end subroutine get_hbath_validity 


    subroutine write_hbath_validity(self, string)

        CLASS(observables), intent(inout)   ::self
        CHARACTER(len=*), intent(in)                 :: string
        !Outputfile
        CHARACTER(len=1000)                          ::outputFile
        !DO LOOP
        INTEGER                                     ::i,j


            outputFile =trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//".val"

            WRITE(*,'(A,A)') "Validity File is written to: ", outputFile

            OPEN(UNIT=9, FILE=outputfile, ACTION="write", STATUS="replace")

                DO i = 1, self%parameter_grid_length
                    WRITE(9, '(10000E)') self%parameter_grid(i), self%validity_matrix(i)
                END DO

            CLOSE(9)

    end subroutine write_hbath_validity 


    subroutine get_performance_sum(self, grid_length)

        CLASS(observables), intent(inout)   ::self
        REAL(8)                             ::performance_sum
        INTEGER, intent(in)                 :: grid_length

        performance_sum = SUM(self%performance(1:self%parameter_grid_length))

        WRITE(*,'(A, E10.3, A, E10.3)') "Performance sum:", performance_sum, " run time: ", grid_length*performance_sum

    end subroutine get_performance_sum


    subroutine put_in_x_energy_matrix(self, sys, indexx)

        CLASS(system), intent(in)::sys
        CLASS(observables), intent(inout) :: self

        INTEGER, intent(in)               :: indexx

        self%x_energy_matrix(1:self%x_energy_number, 1:self%x_energy_number, indexx) = sys%X_energybase(1:self%x_energy_number, 1:self%x_energy_number, 1) 

    end subroutine put_in_x_energy_matrix


    subroutine write_x_energy_matrix(self, string, diagonal_number)

        CLASS(observables), intent(inout) :: self
        INTEGER, intent(in)               :: diagonal_number
        CHARACTER(len=*), intent(in)      :: string

        CHARACTER(len=1000)                         :: outputFile
        CHARACTER(len=2)                            :: diagonal_number_string 
        !Loop
        INTEGER                                     :: n, i, j 

        IF(diagonal_number.le.self%x_energy_number) THEN

        DO n = 0, diagonal_number 
                
            !Conversoin 
            WRITE(diagonal_number_string, '(I2)') n

              outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_diagonal_"//trim(adjustl(diagonal_number_string))//"_x_energyspace_matrix.xmatrix"

                OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

                    ASSOCIATE(sec_end => self%sec_parameter_grid_length, &
                        sec_param => self%sec_parameter_grid,&

                        prim_end => self%parameter_grid_length, &
                        prim_param => self%parameter_grid)
                          
                        DO  i = 1, prim_end 
                            !WRITE(16, '(10000E)') prim_param(i), (Abs(self%x_energy_matrix(j,j+n,i))* Abs(self%x_energy_matrix(j,j+n,i)) , j=1, self%x_energy_number-n)
                            WRITE(16, '(10000E)') prim_param(i), (self%x_energy_matrix(j,j+n,i) , j=1, self%x_energy_number-n)
                        END DO

                    END ASSOCIATE
                   
                CLOSE(16)

                WRITE(*,*) "Position X in energy matrix written to: ", outputFile

         END DO 

         END IF 
                

    end subroutine write_x_energy_matrix
       
        
    subroutine put_in_frank_matrix(self, couplings, indexx)

        CLASS(observables), intent(inout) :: self
        CLASS(coupling), intent(inout) :: couplings
        INTEGER, intent(in)               :: indexx

        INTEGER                           :: i,j
        
        self%frank_matrix(1:self%frank_number,1:self%frank_number, indexx) = couplings%jnVdM(1:self%frank_number,1:self%frank_number,1) 
    
        !   DO i = 1, self%frank_number
        !       DO j = 1, self%frank_number
        !           self%frank_matrix(i,j,indexx) = couplings%jnVdM(i,j,0) 
        !       END DO
        !   END DO
    
    end subroutine put_in_frank_matrix

    subroutine write_frank_matrix(self, string, diagonal_number)

        CLASS(observables), intent(inout)           :: self
        INTEGER, intent(in)                         :: diagonal_number
        CHARACTER(len=*), intent(in)                :: string
        
        
        CHARACTER(len=1000)                         :: outputFile
        CHARACTER(len=2)                            :: diagonal_number_string 
        !Loop
        INTEGER                                     :: n, i, j 
       
        IF(diagonal_number.le.self%frank_number) THEN

        DO n = 0, diagonal_number 
                
            !Conversoin 
            WRITE(diagonal_number_string, '(I2)') n 

              outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_diagonal_"//trim(adjustl(diagonal_number_string))//"_frank_matrix.frank"

                OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

                    ASSOCIATE(sec_end => self%sec_parameter_grid_length, &
                        sec_param => self%sec_parameter_grid,&

                        prim_end => self%parameter_grid_length, &
                        prim_param => self%parameter_grid)
                          
                        DO  i = 1, prim_end 
                            !WRITE(16, '(10000E)') prim_param(i), (Abs(self%frank_matrix(j,j-n,i))* Abs(self%frank_matrix(j,j-n,i)) , j=1, self%frank_number-n)
                            WRITE(16, '(10000E)') prim_param(i), (Abs(self%frank_matrix(j,j+n,i))* Abs(self%frank_matrix(j,j+n,i)) , j=1, self%frank_number-n)
                        END DO

                    END ASSOCIATE
                   
                CLOSE(16)

                WRITE(*,*) "Frank Matrix written to: ", outputFile

         END DO 

         END IF 
                

    end subroutine write_frank_matrix


    subroutine put_in_observable_matrix(self,indexx)

        CLASS(observables), intent(inout) :: self
        INTEGER, intent(in)               :: indexx

        INTEGER                           :: n

        self%energy_matrix(:,indexx) = self%energy
        self%position_matrix(:,indexx) = self%position
        self%current_matrix(:,indexx) = self%current
        self%bridge_population_matrix(:,indexx,0)= self%bridge_population(:,0)
        self%bridge_population_matrix(:,indexx,1)= self%bridge_population(:,1)

        DO n = 1, self%state_population_number 

            self%state_occupation_pop_matrix(:,indexx, n) = self%state_population(n,:,1)
            self%state_unoccupation_pop_matrix(:,indexx, n) = self%state_population(n,:,0)

        END DO

    end subroutine put_in_observable_matrix

    subroutine write_parameter_grid(self)

        CLASS(observables), intent(in)               :: self

        INTEGER                                      :: i, j

        CHARACTER(len=1000)                         :: outputFile_prim
        CHARACTER(len=1000)                         :: outputFile_sec

        outputFile_prim = trim(adjustl(self%sOutput))//"primary.grid"
        outputFile_sec = trim(adjustl(self%sOutput))//"secondary.grid"

            ASSOCIATE(sec_end => self%sec_parameter_grid_length, &
                sec_param => self%sec_parameter_grid,&
                prim_end => self%parameter_grid_length, &
                prim_param => self%parameter_grid)

        !Primary Grid
            OPEN(16, file = outputFile_prim, ACTION="write", STATUS="replace")

                DO  i = 1, prim_end
                    WRITE(16, '(E)') (prim_param(i))
                END DO

            CLOSE(16)

        !Secondary Grid
            OPEN(17, file = outputFile_sec, ACTION="write", STATUS="replace")

                DO  j = 1, sec_end
                    WRITE(17, '(E)') (sec_param(j))
                END DO

            CLOSE(17)
            END ASSOCIATE

            WRITE(*,*) "Parameter grids for heatmaps are written to: "
            WRITE(*,*) outputFile_prim
            WRITE(*,*) outputFile_sec

    end subroutine write_parameter_grid


    subroutine write_observable_state_population_matrix_heatmap(self, string, population)

        CLASS(observables), intent(in)               :: self
        INTEGER, intent(in)                          :: population

        CHARACTER(len=3)                              :: population_string
        CHARACTER(len=*), intent(in)                 :: string
        CHARACTER(len=1000)                         :: outputFile

        INTEGER                                     :: i, j !< Loop Variable
        !Populationsnumber
        INTEGER                                     :: n !< Loop Variable
        CHARACTER(len=2)                            :: n_string 

        IF((population.eq.0).or.(population.eq.1)) THEN

            IF(population.eq.0) THEN
                population_string = "uno"
            ELSE
                population_string = "occ"
            END IF

            DO n = 1, self%state_population_number
                WRITE(n_string, '(I2)') n

              outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_"//population_string//"_"//trim(adjustl(n_string))//"_heatmap_state.pop"

                OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

                ASSOCIATE(sec_end => self%sec_parameter_grid_length, &
                    sec_param => self%sec_parameter_grid,&

                    prim_end => self%parameter_grid_length, &
                    prim_param => self%parameter_grid)

                        IF(population.eq.0) THEN
                            
                            DO  i = 1, prim_end
                                WRITE(16, '(10000E)') (self%state_unoccupation_pop_matrix(i,j,n), j=1,sec_end)
                            END DO
                       
                        ELSE
                      
                            DO  i = 1, prim_end
                                WRITE(16, '(10000E)') (self%state_occupation_pop_matrix(i,j,n), j=1,sec_end)
                            END DO
                        
                        END IF

                END ASSOCIATE

            CLOSE(16)
            
            END DO 

            WRITE(*,*) "observable heatmap  in computation mode 3 is written in Format 2 (for gnuplot): ", outputFile

        ELSE

            WRITE(*,*) "no valid value for population variable in population matrix method"

        END IF

    end subroutine write_observable_state_population_matrix_heatmap

    subroutine write_observable_population_matrix_heatmap(self, string, population)

        CLASS(observables), intent(in)               :: self
        INTEGER, intent(in)                          :: population

        CHARACTER(len=3)                              :: population_string
        CHARACTER(len=*), intent(in)                 :: string
        CHARACTER(len=1000)                         :: outputFile

        INTEGER                                     :: i, j !< Loop Variable

        IF((population.eq.0).or.(population.eq.1)) THEN

            IF(population.eq.0) THEN
                population_string = "uno"
            ELSE
                population_string = "occ"
            END IF

            outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_"//population_string//"_heatmap.popum"

            OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

            ASSOCIATE(sec_end => self%sec_parameter_grid_length, &
                sec_param => self%sec_parameter_grid,&

                prim_end => self%parameter_grid_length, &
                prim_param => self%parameter_grid)

                DO  i = 1, prim_end
                    WRITE(16, '(10000E)') (self%current_matrix(i,j), j=1,sec_end)
                END DO

            END ASSOCIATE

            CLOSE(16)

            WRITE(*,*) "observable heatmap  in computation mode 3 is written in Format 2 (for gnuplot): ", outputFile

        ELSE

            WRITE(*,*) "no valid value for population variable in population matrix method"

        END IF

    end subroutine write_observable_population_matrix_heatmap


    subroutine write_observable_population_matrix_tofile_xyz(self, string, population)

        CLASS(observables), intent(in)               :: self
        INTEGER, intent(in)                          :: population

        CHARACTER(len=3)                              :: population_string
        CHARACTER(len=*), intent(in)                 :: string
        CHARACTER(len=1000)                         :: outputFile

        INTEGER                                     :: i, j !< Loop Variable

        IF((population.eq.0).or.(population.eq.1)) THEN

            IF(population.eq.0) THEN
                population_string = "uno"
            ELSE
                population_string = "occ"
            END IF

            outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_"//population_string//"_xyz.popum"

            OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

            ASSOCIATE(sec_end => self%sec_parameter_grid_length, &
                sec_param => self%sec_parameter_grid,&

                prim_end => self%parameter_grid_length, &
                prim_param => self%parameter_grid)

                !Time Loop
                DO j = 1, sec_end
                    !Important Blank for gnuplot
                    WRITE(16,'(A)') ' '
                    !Loop over the position on the grid
                    DO  i = 1, prim_end
                        WRITE(16, '(E, E, E)') prim_param(i), sec_param(j), self%current_matrix(i,j)
                    END DO
                END DO

            END ASSOCIATE

            CLOSE(16)

            WRITE(*,*) "observable matrix in computation mode 2 is written in Format 2 (for gnuplot): ", outputFile

        ELSE

            WRITE(*,*) "no valid value for population variable in population matrix method"

        END IF

    end subroutine write_observable_population_matrix_tofile_xyz

    subroutine write_observable_current_matrix_heatmap(self, string)

        CLASS(observables), intent(in)               :: self
        CHARACTER(len=*), intent(in)                 :: string
        CHARACTER(len=1000)                         :: outputFile

        INTEGER                                     :: i, j !< Loop Variable

        outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_heatmap.cum"

        OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

        ASSOCIATE(sec_end => self%sec_parameter_grid_length, &
            sec_param => self%sec_parameter_grid,&

            prim_end => self%parameter_grid_length, &
            prim_param => self%parameter_grid)

            DO  i = 1, prim_end
                WRITE(16, '(10000E)') (self%current_matrix(i,j), j=1,sec_end)
            END DO

        END ASSOCIATE

        CLOSE(16)

        WRITE(*,*) "current heatmap in computation mode 3 is written in Format 2 (for gnuplot): ", outputFile

    end subroutine write_observable_current_matrix_heatmap

    subroutine write_observable_current_matrix_tofile_xyz(self, string)

        CLASS(observables), intent(in)               :: self
        CHARACTER(len=*), intent(in)                 :: string
        CHARACTER(len=1000)                         :: outputFile

        INTEGER                                     :: i, j !< Loop Variable

        outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_xyz.cum"

        OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

        ASSOCIATE(sec_end => self%sec_parameter_grid_length, &
            sec_param => self%sec_parameter_grid,&

            prim_end => self%parameter_grid_length, &
            prim_param => self%parameter_grid)

            !Time Loop
            DO j = 1, sec_end
                !Important Blank for gnuplot
                WRITE(16,'(A)') ' '
                !Loop over the position on the grid
                DO  i = 1, prim_end
                    WRITE(16, '(E, E, E)') prim_param(i), sec_param(j), self%current_matrix(i,j)
                END DO
            END DO

        END ASSOCIATE

        CLOSE(16)

        WRITE(*,*) "current matrix in computation mode 2 is written in Format 2 (for gnuplot): ", outputFile

    end subroutine write_observable_current_matrix_tofile_xyz

    subroutine write_observable_position_matrix_heatmap(self, string)

        CLASS(observables), intent(in)               :: self
        CHARACTER(len=*), intent(in)                 :: string
        CHARACTER(len=1000)                         :: outputFile

        INTEGER                                     :: i, j !< Loop Variable

        outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_heatmap.pom"

        OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

        ASSOCIATE(sec_end => self%sec_parameter_grid_length, &
            sec_param => self%sec_parameter_grid,&

            prim_end => self%parameter_grid_length, &
            prim_param => self%parameter_grid)

            !Loop over the position on the grid
            DO  i = 1,prim_end
                WRITE(16, '(10000E)') (self%position_matrix(i,j), j=1,sec_end)
            END DO

        END ASSOCIATE

        CLOSE(16)

        WRITE(*,*) "position heatmap in computation mode 2 is written in Format 2 (for gnuplot): ", outputFile

    end subroutine write_observable_position_matrix_heatmap

    subroutine write_observable_position_matrix_tofile_xyz(self, string)

        CLASS(observables), intent(in)               :: self
        CHARACTER(len=*), intent(in)                 :: string
        CHARACTER(len=1000)                         :: outputFile

        INTEGER                                     :: i, j !< Loop Variable

        outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_xyz.pom"

        OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

        ASSOCIATE(sec_end => self%sec_parameter_grid_length, &
            sec_param => self%sec_parameter_grid,&

            prim_end => self%parameter_grid_length, &
            prim_param => self%parameter_grid)

            !Time Loop
            DO j = 1, sec_end
                !Important Blank for gnuplot
                WRITE(16,'(A)') ' '
                !Loop over the position on the grid
                DO  i = 1, prim_end
                    WRITE(16, '(E, E, E)') prim_param(i), sec_param(j), self%position_matrix(i,j)
                END DO
            END DO

        END ASSOCIATE

        CLOSE(16)

        WRITE(*,*) "position matrix in computation mode 2 is written in Format 2 (for gnuplot): ", outputFile

    end subroutine write_observable_position_matrix_tofile_xyz

    subroutine write_observable_energy_matrix_tofile_xyz(self, string)

        CLASS(observables), intent(in)               :: self
        CHARACTER(len=*), intent(in)                 :: string
        CHARACTER(len=1000)                         :: outputFile

        INTEGER                                     :: i, j !< Loop Variable

        outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_xyz.enm"

        OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

        ASSOCIATE(sec_end => self%sec_parameter_grid_length, &
            sec_param => self%sec_parameter_grid,&

            prim_end => self%parameter_grid_length, &
            prim_param => self%parameter_grid)

            !Time Loop
            DO j = 1, sec_end
                !Important Blank for gnuplot
                WRITE(16,'(A)') ' '
                !Loop over the position on the grid
                DO  i = 1, prim_end
                    WRITE(16, '(E, E, E)') prim_param(i), sec_param(j), self%energy_matrix(i,j)
                END DO
            END DO

            CLOSE(16)

        END ASSOCIATE


        WRITE(*,*) "energy matrix in computation mode 2 is written in Format 2 (for gnuplot): ", outputFile

    end subroutine write_observable_energy_matrix_tofile_xyz

    subroutine write_observable_energy_matrix_heatmap(self, string)

        CLASS(observables), intent(in)               :: self
        CHARACTER(len=*), intent(in)                 :: string
        CHARACTER(len=1000)                         :: outputFile

        INTEGER                                     :: i, j !< Loop Variable

        outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"heatmap.enm"

        OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

        ASSOCIATE(sec_end => self%sec_parameter_grid_length, &
            sec_param => self%sec_parameter_grid,&

            prim_end => self%parameter_grid_length, &
            prim_param => self%parameter_grid)

            DO  i = 1, prim_end
                WRITE(16, '(1000E)') (self%energy_matrix(i,j), j=1,sec_end)
            END DO

        END ASSOCIATE

        CLOSE(16)
        WRITE(*,*) "energy matrix in computation mode 2 is written in Format 2 (for gnuplot): ", outputFile

    end subroutine write_observable_energy_matrix_heatmap


    subroutine write_index_to_console(self, indexx)

        CLASS(observables), intent(in)  ::self
        INTEGER, intent(in)             :: indexx

        WRITE(*,'(A,I)') "Index: ", indexx
        WRITE(*,'(A,E)') "Position: ", self%position(indexx)
        WRITE(*,'(A,E)') "Current: ", self%current(indexx)
        WRITE(*,'(A,E)') "Energy: ", self%energy(indexx)
        WRITE(*,'(A,E)') "Population on Bridge: ", self%bridge_population(indexx,1)

    end subroutine write_index_to_console

    subroutine write_results_tofile(self, string)

        CLASS(observables)::self

        !Outputfile
        CHARACTER(len=1000)                          ::outputFile
        CHARACTER(len=12)                            ::extension
        !DO LOOP
        INTEGER                                     ::i
        INTEGER                                     ::err
        REAL(8)                                     ::parameter_conversion

        CHARACTER(len=*), intent(in)                             :: string


        !No Conversion for non-time-evolution mode
        parameter_conversion = 1.0d0

        IF(self%plot_bool.eq.1) THEN

            SELECT CASE(self%mode_specifier)
                
                CASE("P")

                    extension = "temp_.dat"

                CASE("T", "K", "U")

                    extension = ".evo"
                    parameter_conversion = au2femto 

                CASE("V")

                    extension = ".dat"

                CASE("W")

                    extension = ".dat"

                CASE("G")

                    extension = "_cvg.dat"

                CASE("B")

                    extension = "_wbias.dat"

                CASE DEFAULT

                    extension = ".unkwnown"
                    WRITE(*,*) "No valid specifier in object found"

            END SELECT

            outputFile =trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//trim(adjustl(extension))

            OPEN(UNIT=9, file= outputFile  , ACTION="write", STATUS="replace", IOSTAT = err)

            DO i = 1, self%parameter_grid_length
                WRITE(9, '(10000E)') self%parameter_grid(i)*parameter_conversion, self%position(i), self%current(i), &
                    self%energy(i), self%bridge_population(i,0), self%bridge_population(i,1)
            END DO

            CLOSE(9)

            WRITE(*,'(A,A)') 'Result has been written to file: ', outputFile

        ELSE

            WRITE(*,'(A)') 'Plotbool is 0, thus no results file is written'

        END IF

    end subroutine write_results_tofile


    !Not Used yet
    subroutine set_current_auxil_matrix(self, sys, junction, density, coup)

        class(observables), intent(inout):: self
        class(system), intent(in)::sys
        class(leads), intent(in)::junction
        class(densityvector), intent(inout)::density
        class(coupling), intent(in)::coup


    end subroutine set_current_auxil_matrix


    subroutine set_X_auxil_vector(self, sys, density)

        CLASS(observables)::self
        CLASS(system)::sys
        CLASS(densityvector)::density

        REAL(8), ALLOCATABLE, DIMENSION(:, :, :)  :: X_energybase
        REAL(8), ALLOCATABLE, DIMENSION(:,:,:) ::hsOV_copy

        !Summation
        REAL(8)                                         :: summation
        INTEGER                                         :: occupation
        !Loop
        INTEGER             ::i,j
        INTEGER             ::k

        ALLOCATE(hsOV_copy(1:sys%N, 1:sys%N, 0:1))
        ALLOCATE(X_energybase(1:density%maxdiagvalue, 1:density%maxdiagvalue, 0:1))

        !Create a working copy of systems hsOV for faster processing
        hsOV_copy = sys%hsOV

        ASSOCIATE(  all_m => density%NumDim00,              &
            all_v => density%NumDim11,              &
            diag_m => density%medim00,          &
            diag_v => density%medim11,          &

            !Indize
            m1 => density%tulpVector(:,1,0),        &
            m2 => density%tulpVector(:,2,0),        &
            v1 => density%tulpVector(:,1,1),        &
            v2 => density%tulpVector(:,2,1))


            !> Computation of X in energybase
            DO occupation = 0,1
                DO i = 1, density%maxdiagvalue
                    DO j = 1, density%maxdiagvalue
                        summation = 0.0d0
                        DO k = 1, sys%N
                            summation = summation + sys%grid(k)*sys%hsOV(k,i,occupation)*hsOV_copy(k,j,occupation)
                        END DO
                        X_energybase(i,j, occupation) = summation
                    END DO
                END DO
            END DO

            !< Preparation for the auxilary vector
            !--Diagonalelements in 00
            DO j = 1, all_m
                self%x_auxil_vector(j) = sys%X_energybase(m1(j), m2(j),0)
            END DO

            !--Diagonalelements in 11
            DO j = 1, all_v
                self%x_auxil_vector(j + all_m) = sys%X_energybase(v1(j),v2(j),1)
            END DO
        END ASSOCIATE

    end subroutine set_X_auxil_vector

    subroutine write_sub_matrix(self, string )

        CLASS(observables), intent(in)                          ::self !< self object
        CHARACTER(len=*), intent(in)                            :: string!< String to specifiy the output file exactly

        !Outputfile

        CHARACTER(len=1000)                          ::outputFile

        !DO LOOP
        INTEGER                                     ::i,j

        outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//".sub"
        OPEN(UNIT=9, FILE= outputFile, ACTION="write", STATUS="replace")

            DO i = 1,self%parameter_grid_length, 1
                WRITE (9,'(E, 10000E)')  self%parameter_grid(i), ( self%sub_matrix(i, j), j = 0, 3 ), 2d0*self%sub_matrix(i, 2), &
                    -2d0*self%sub_matrix(i, 3), (self%sub_matrix(i, 0) - self%sub_matrix(i, 1))
            END DO

        CLOSE(9)

        WRITE(*,'(A, A)') "Submatrix written in", outputFile

    end subroutine write_sub_matrix

    subroutine write_state_energy_tofile(self, string)

        CLASS(observables), intent(in)                          ::self !< self object
        CHARACTER(len=*), intent(in)                            :: string!< String to specifiy the output file exactly

        INTEGER                                                 :: state_number

        !Outputfile
        CHARACTER(len=1000)                          ::outputFile

        !DO LOOP
        INTEGER                                     ::i,j

            state_number = self%pop_number
            outputFile =trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//".dia"

            OPEN(UNIT=10, FILE= outputFile, ACTION="write", STATUS="replace")

                DO i = 1,self%parameter_grid_length, 1
                    WRITE (10,'(E, 10000E)')  self%parameter_grid(i), (self%state_energys(j, i, 1), j = 1, state_number)
                END DO

            CLOSE(10)

            WRITE(*,'(A, A)') "States energys written in", outputFile

    end subroutine write_state_energy_tofile
    
    subroutine write_population_histo_tofile(self, string, occupation)

        CLASS(observables), intent(in)                          ::self !< self object
        CHARACTER(len=*), intent(in)                            :: string!< String to specifiy the output file exactly
        INTEGER, intent(in)                                     :: occupation

        INTEGER                                                 :: state_number

        !Outputfile

        CHARACTER(len=1000)                          ::outputFile

        !DO LOOP
        INTEGER                                     ::i,j

        state_number = self%pop_number

        IF(self%pop_bool.eq.1) THEN
            SELECT CASE(occupation)

                CASE(0)

                    outputFile =trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_histo.p0"
                    OPEN(UNIT=9, FILE= outputFile, ACTION="write", STATUS="replace")

                    !First Line 
                    ! Values of the primary parameter
                    WRITE (9,'(10000E)') (self%parameter_grid(i), i = 1, self%parameter_grid_length) 
                    !Population Matrix 
                    DO i = 1, state_number 
                        WRITE (9,'(10000E)') (self%state_population(i, j, 0), j = 1, self%parameter_grid_length)
                    END DO

                    CLOSE(9)

                    WRITE(*,'(A, I1, A, A)') "States Population in histogram format for occupation ", occupation, " written in", outputFile

                CASE(1)

                    outputFile =trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_histo.p1"

                    !Occupied
                    OPEN(UNIT=10, FILE= outputFile, ACTION="write", STATUS="replace")
                    
                    !First Line 
                    ! Values of the primary parameter
                    WRITE (10,'(10000E)') (self%parameter_grid(i), i = 1, self%parameter_grid_length) 
                    !Population Matrix 
                    DO i = 1, state_number 
                        WRITE (10,'(10000E)') (self%state_population(i, j, 1), j = 1, self%parameter_grid_length)
                    END DO
                    CLOSE(10)

                    WRITE(*,'(A, I, A, A)') "States Population in histogram format for occupation ", occupation, " written in", outputFile

                CASE DEFAULT

                    WRITE(*,*) "No valid option for obligatory occupation parameter in population output method"

            END SELECT

        ELSE
            WRITE(*,*) "Plotbool is 0, therefore no population output in histogram format is generated"
        END IF
    end subroutine write_population_histo_tofile


    subroutine write_population_tofile(self, string, occupation)

        CLASS(observables), intent(in)                          ::self !< self object
        CHARACTER(len=*), intent(in)                            :: string!< String to specifiy the output file exactly
        INTEGER, intent(in)                                     :: occupation

        INTEGER                                                 :: state_number

        !Outputfile

        CHARACTER(len=1000)                          ::outputFile

        !DO LOOP
        INTEGER                                     ::i,j

        state_number = self%pop_number

        IF(self%pop_bool.eq.1) THEN
            SELECT CASE(occupation)

                CASE(0)

                    outputFile =trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//".p0"
                    OPEN(UNIT=9, FILE= outputFile, ACTION="write", STATUS="replace")

                    DO i = 1,self%parameter_grid_length, 1
                        WRITE (9,'(E, 10000E)')  self%parameter_grid(i), (self%state_population(j , i, 0), j = 1, state_number)
                    END DO

                    CLOSE(9)

                    WRITE(*,'(A, I1, A, A)') "States Population for occupation ", occupation, " written in", outputFile

                CASE(1)

                    outputFile =trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//".p1"

                    !Occupied
                    OPEN(UNIT=10, FILE= outputFile, ACTION="write", STATUS="replace")

                    DO i = 1,self%parameter_grid_length, 1
                        WRITE (10,'(E, 10000E)')  self%parameter_grid(i), (self%state_population(j , i, 1), j = 1, state_number)
                    END DO

                    CLOSE(10)

                    WRITE(*,'(A, I, A, A)') "States Population for occupation ", occupation, " written in", outputFile

                CASE DEFAULT

                    WRITE(*,*) "No valid option for obligatory occupation parameter in population output method"

            END SELECT

        ELSE
            WRITE(*,*) "Plotbool is 0, therefore no population output is generated"
        END IF
    end subroutine write_population_tofile

    subroutine write_summary_tofile(self, sys, cou, number_wave, string)

        CLASS(observables), intent(in):: self
        CLASS(system), intent(in) ::sys
        CLASS(coupling), intent(in)::cou

        INTEGER, intent(in)                                    :: number_wave
        CHARACTER(len=*), intent(in)                            :: string

        !Outputfile
        CHARACTER(len=1000)                          ::outputFile
        !DO LOOP
        INTEGER                                     ::i,j

        IF(self%summary_bool.eq.1) THEN

            outputFile =trim(adjustl(self%sOutput_summary))//"_"//trim(adjustl(string))//".sum"

            WRITE(*,'(A,A)') 'Systemdata in summary format has been written to file: ', outputFile

            OPEN(UNIT=9, FILE=outputfile  , ACTION="write", STATUS="replace")

            !Loop over the dvr position grid points
            DO i = 1,sys%N
                WRITE(9, '(10000E)') sys%grid(i), sys%hsEN(i,0)*au2eV,  sys%hsEN(i,1)*au2eV, &
                    sys%potential(sys%grid(i), sys%id)*au2ev, sys%potential(sys%grid(i), sys%id + 1)*au2eV,&
                    cou%switch_table(i), 0.0d0, (sys%hsOV(i,j, 1) + sys%hsEN(j,1)*au2ev , j=1,number_wave)
            END DO

            CLOSE(9)

        ELSE

            WRITE(*,'(A)')  "Plotbool is 0, thus no summary file of the observable object is written"

        END IF
    end subroutine write_summary_tofile

    subroutine write_performance_tofile(self, string)

        CLASS(observables)::self

        CHARACTER(len=*), intent(in)                            :: string

        !Outputfile
        CHARACTER(len=1000)                          ::outputFile
        !DO LOOP
        INTEGER                                     ::i,j

        REAL(8)                                     :: complete_time

        IF(self%performance_bool.eq.1) THEN

            complete_time = SUM(self%performance)

            outputFile =trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//".per"

            WRITE(*,'(A,A)') "Performance File is written to: ", outputFile

            OPEN(UNIT=9, FILE=outputfile  , ACTION="write", STATUS="replace")
            WRITE(9,'(A,E)') '?', complete_time
            WRITE(9,'(A,E)') '?'
            !Loop over the dvr position grid points
            DO i = 1, self%parameter_grid_length
                WRITE(9, '(10000E)') self%parameter_grid(i), self%performance(i)
            END DO

            CLOSE(9)

        ELSE

            WRITE(*,'(A)') "Plotbool is 0, no performance file is written"

        END IF

    end subroutine write_performance_tofile


end module class_observables
