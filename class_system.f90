module class_system

    USE UnitConv
    USE class_inputdata
    USE class_grid
    USE class_density

    implicit none
    private

    type, public :: system

        INTEGER                                            :: id
        INTEGER                                            :: N
        INTEGER                                            :: solve_method
        INTEGER                                            :: abstol

        INTEGER                                            :: dimension

        CHARACTER(len=1000)                                   :: sOutput !< Output path and prefix for outputf iles
        INTEGER                                     :: plot_bool !< Pseudo Boolean for outputfiles fetched from the inputdata object in the intialization. 0 - no output 1 - output.
        INTEGER                                     :: surf_bool !< Pseudo Boolean for outputfiles fetched from the inputdata object in the intialization. 0 - no output 1 - output.

        !GridFlag
        LOGICAL                                            :: first_call_flag = .false.


        !Double Well Parameters
        REAL(8),private                                    :: l !< Distance beteween local minima
        REAL(8),private                                   :: delta  !< Energy shift, here for the right minima.
        REAL(8),private                                    :: Vb     !<Value V(0)=Vb, barrier height
        REAL(8),private                                    :: shiftEnergy !< constant shifting of the occupied system against the unoccupied
        REAL(8),private                                    :: x_shift !< vertical shift, position shift of the occupied system against the unoccupoied

        !2nd Set of Double Well Parameters
        REAL(8),private                                    :: l_2 !< Distance beteween local minima
        REAL(8),private                                   :: delta_2  !< Energy shift, here for the right minima.
        REAL(8),private                                    :: Vb_2     !<Value V(0)=Vb, barrier height
        REAL(8),private                                    :: shiftEnergy_2 !< constant shifting of the occupied system against the unoccupied
        REAL(8),private                                    :: x_shift_2 !< vertical shift, position shift of the occupied system against the unoccupoied

        !Vibronic
        REAL(8),private                                    :: Omega !< parameter omega for the harmonic potential
        REAL(8),private                                    :: lambda !< parameter lambda for the harmonic potential

        REAL(8),private                                    :: parabola_shift !< shift parameter for the harmonic potential
        !Parameters
        !---Electrical Field
        REAL(8),private                                      :: d !< distance the electrial field is applied, distance betwenn the two leads
        REAL(8),private                                    :: Gate_Voltage !< Gate voltage applied to the molecule

        !---Leads
        REAL(8),private                                     :: bias_voltage
        REAL(8), private                                    :: Lj
        REAL(8), private                                    :: angle

        !Copy of the primary paramter grid
        REAL(8), ALLOCATABLE, DIMENSION(:)        :: parameter_grid !< this is the paramter the calculating loop varies
        INTEGER                                   :: parameter_grid_length


        !Eigenvalues and Eigenvectors
        !---Eigenvalues (Eigenvalues, 0:1 (unoccupied, occupied))
        REAL(8), public, ALLOCATABLE, DIMENSION(:,:)                :: hsEN

        !---Eigenvectors ( Matrix of Eigenvectors (n-th Column = n -th eigenvector), 0:1 (unoccupied, occupied) 2:Additional Information)
        REAL(8), public, ALLOCATABLE, DIMENSION(:,:,:)              :: hsOV

        !---Own Copy of corresponding grid points
        REAL(8), public, ALLOCATABLE, DIMENSION(:)                 ::grid

        !X in Energyspace (for Calculation of mean Value <x>)
        REAL(8), ALLOCATABLE, DIMENSION(:, :, :)                :: X_energybase

        !Transition Frequencies
        REAL(8), ALLOCATABLE, DIMENSION(:, :, :)                :: tran_freq


    contains
        procedure, public :: init_system
        procedure, public :: potential
        procedure, private :: set_system
        procedure, private :: compute_system
        procedure, private :: set_x_in_energybase
        procedure, public :: write_system_to_file
        procedure, public :: change_gate
        procedure, public :: change_gate_and_bias_field
        procedure, public :: potential_surface
        procedure, private :: put_console_info
        procedure, private :: compute_transition_frequencies
        procedure, public :: write_x_energy
        final :: destructor
    end type system

contains

    subroutine init_system(self, input, grid_points, density, gate_in)

        IMPLICIT NONE

        CLASS(system), intent(inout)        :: self
        CLASS(densityvector), intent(in)    :: density
        CLASS(inputdata), intent(in)        :: input
        CLASS(grid), intent(in)             :: grid_points

        INTEGER, intent(in), optional       :: gate_in

        IF(present(gate_in)) THEN

            SELECT CASE(gate_in)

                CASE(1)

                    self%Gate_Voltage = input%start_gate_voltage

                CASE(2)

                    self%Gate_Voltage = input%end_gate_voltage

                CASE DEFAULT

                    WRITE(*,*) "No valid value for the gate_in argument. Default value is taken from the start_gate_voltage variable"

            END SELECT

        !Default is start_gate_voltage
        ELSE

            self%Gate_Voltage = input%start_gate_voltage

        END IF

        self%N = grid_points%N
        self%dimension = density%maxDiagValue
        self%id = input%id

        !1st Set
        self%l = input%l
        self%delta = input%delta
        self%Vb = input%Vb

        !2st Set
        self%delta_2 = input%delta_2
        self%Vb_2 = input%Vb_2
        self%l_2 = input%l_2

        self%shiftEnergy = input%shiftEnergy
        self%x_shift = input%x_shift

        !Harmonic Set
        self%Omega = input%Omega
        self%parabola_shift = input%parabola_shift

        !Electric Field of the Leads
        self%d = input%d
        self%Lj = input%Lj
        self%angle = input%angle

        !Computational
        self%abstol = input%abstol
        self%solve_method=input%solve_method

        ALLOCATE(self%hsEN(1:self%N, 0:1))
        ALLOCATE(self%HsOV(1:self%N,1:self%N,0:1))
        ALLOCATE(self%grid(1:self%N))
        ALLOCATE(self%X_energybase(1:self%dimension, 1:self%dimension, 0:1))
        ALLOCATE(self%tran_freq(1:self%dimension, 1:self%dimension, 0:1))

        !Copy the grid from grid_points object
        self%grid = grid_points%dvrX

        !Get output file path with prefix
        self%sOutput = input%system_output_filepath

        self%plot_bool = input%plot_bool
        self%surf_bool = input%surf_bool

        !First system set
        CALL self%change_gate(self%Gate_Voltage)

        !Bias Grid for potential surface Method
        ALLOCATE(self%parameter_grid(input%parameter_grid_length)) !< this is the paramter the calculating loop varies
        self%parameter_grid = input%parameter_grid
        self%parameter_grid_length = input%parameter_grid_length

        self%first_call_flag = .true.

        CALL self%compute_transition_frequencies


    end subroutine init_system

        !---unoccupied
    subroutine set_system(self)
        IMPLICIT NONE
        CLASS(system):: self

        !Eigensystem Calculation
        CALL self%compute_system(0)
        !---occupied
        CALL self%compute_system(1)

        CALL self%put_console_info

    END SUBROUTINE set_system

    SUBROUTINE set_x_in_energybase(self)

        CLASS(system)::self

        REAL(8), ALLOCATABLE, DIMENSION(:,:,:) ::hsOV_copy

        !Summation
        REAL(8)                                         :: summation
        INTEGER                                         :: occupation
        !Loop
        INTEGER             ::i,j
        INTEGER             ::k

        ALLOCATE(hsOV_copy(1:self%N, 1:self%N, 0:1))

        !Create a working copy of systems hsOV for faster processing
        hsOV_copy = self%hsOV

        DO occupation = 0,1
            DO i = 1,self%dimension
                DO j = 1, self%dimension
                    summation = 0.0d0
                    DO k = 1, self%N
                        summation = summation + self%grid(k)*self%hsOV(k,i,occupation)*hsOV_copy(k,j,occupation)
                    END DO
                    self%X_energybase(i,j, occupation) = summation
                END DO
            END DO
        END DO

    END SUBROUTINE

    subroutine compute_system(self, occupation)
        IMPLICIT NONE

        CLASS(system):: self
         !---Variables
        INTEGER, intent(in)                             :: occupation
        INTEGER                                         :: N                    !Dimension

        !Hamilton Matrix
        REAL(8), ALLOCATABLE, DIMENSION(:,:)            :: Hmatrix
        REAL(8), ALLOCATABLE, DIMENSION(:)              :: eigenValueH         !Calculated eigenvalues of X which equals the Gauss Pivots
        REAL(8), ALLOCATABLE, DIMENSION(:,:)            :: eigenVectorH

        !Lapack Variables
        REAL(8), DIMENSION(:,:), ALLOCATABLE              :: Z                 !Work Array for Lapack Routine
        INTEGER                                         :: ldz
        INTEGER                                         :: m

        REAL(8), DIMENSION(:), ALLOCATABLE              :: work                 !Work Array for Lapack Routine
        INTEGER                                         :: lwork                !Dimensions for work and iwork

        REAL(8), DIMENSION(:), ALLOCATABLE              :: iwork
        INTEGER                                         :: liwork

        INTEGER, DIMENSION(:), ALLOCATABLE              :: isuppz                 !Work Array for Lapack Routine
        INTEGER :: nisuppz

        INTEGER                                         :: info


        !Eigenvalue Boundarys
        INTEGER :: vl, vu
        INTEGER :: il, iu


        !Loop Variables
        INTEGER                                         ::i,j

        !Initialization
        !---Dimension of System
        N = self%N
        !---Calculation
        ALLOCATE(Hmatrix(N,N), eigenValueH(N))
        Hmatrix = 0

        !---Diagonal Elements Hermite
        DO i = 1, N
            Hmatrix(i,i) = Minverse*((4*N -1 - 2* self%grid(i) * self%grid(i)) / 6 * 0.5d0 )&
                + self%potential(self%grid(i), self%id + occupation)
        END DO


        !Upper Triangle Block
        DO j = 2,N
            DO i = 1 , (j-1)
                Hmatrix(i,j) = Minverse*(((-1)**(i-j))*(2/(self%grid(i)-self%grid(j))**2 - 0.5d0)*0.5d0)
            END DO
        END DO


        SELECT CASE(self%solve_method)

            CASE(1)

                eigenValueH = 0
                lwork = 5*N

                !---MKL Workarrays
                ALLOCATE(work(lwork))

                CALL DSYEV('V', 'U', N,  Hmatrix, N, eigenValueH, work, lwork, info)

                !Copy and normcheck
                DO j = 1, N
                    self%hsEN(j, occupation) = eigenValueH(j)
                    DO i = 1, N
                        self%hsOV(i,j, occupation) = Hmatrix(i,j)
                    END DO
                END DO


            CASE(2)
                m = N
                lwork = 40*N
                liwork = 40*N
                nisuppz = 3*N
                ldz = N

                ALLOCATE(eigenVectorH(ldz,ldz))

                eigenValueH = 0; eigenVectorH = 0

                !---Allocation
                ALLOCATE(work(lwork), iwork(liwork))
                !---MKL Workarrays
                ALLOCATE(isuppz(nisuppz))

                !EigenValue Calc
                CALL DSYEVR('V', 'A', 'U', N, Hmatrix, N, vl, vu, il, iu, self%abstol, m, eigenValueH, eigenVectorH, ldz, isuppz, work, lwork, iwork, liwork, info)

                DO j = 1, N
                    self%hsEN(j, occupation) = eigenValueH(j)
                    DO i = 1, N
                        self%hsOV(i,j, occupation) = eigenVectorH(i,j)
                    END DO
                END DO


            CASE DEFAULT

                WRITE(*,*) "****************************************************"
                WRITE(*,*) "No correct solve method selected in computation file"
                WRITE(*,*) "****************************************************"


        END SELECT

    END SUBROUTINE compute_system

    subroutine normcheck(self)
        IMPLICIT NONE

        CLASS(system)::self
    !
    !         !NormChecking
    !        REAL(8), ALLOCATABLE, DIMENSION(:)              :: normCheck
    !        REAL(8)                                         :: normSum
    !
    !       !Norm Checking
    !        normSum = 0
    !        !---Software Unrolling if N is even
    !        IF(mod(N,2).eq.0)  THEN
    !            DO i = 1, N, 2
    !                normSum = normSum + normCheck(i)
    !                normSum = normSum + normCheck(i+1)
    !            END DO
    !        ELSE
    !            DO i = 1, N
    !                normSum = normSum + normCheck(i)
    !            END DO
    !        END IF
    !
    !                        normCheck(i) = normCheck(i) +  Abs(hsOV(i,j, occupation))* Abs(hsOV(i,j, occupation))
    !                        normCheck(i) = normCheck(i) +  Abs(hsOV(i,j, occupation))* Abs(hsOV(i,j, occupation))
    !        !Check if normSum is between 98% and 102% of N
    !        IF(((normSum).GE.(normSum*0.98d0)).AND.((normSum).LE.(normSum*1.02d0))) THEN
    !            !WRITE(*,'(A, E , T55, A)') 'NormCheck', normSum, 'POSITIV'
    !
    !        ELSE
    !
    !            WRITE(*,'(A)') 'Attention please. Caclulation of Eigensytem is wrong.'
    !            WRITE(*,'(A)') 'Program Exit'
    !            WRITE(*,'(A, E , T55, A)') 'NormCheck', normSum, 'negativ'
    !            CALL Exit(0)
    !        END IF


    end subroutine normcheck

    subroutine compute_transition_frequencies(self)

        CLASS(system) :: self

        !Loop Variables
        INTEGER n, m
        INTEGER v, w

        !energyDiff Matrix
        !00 electronic subspac
        DO m = 1, self%dimension
            DO n = 1, self%dimension
                self%tran_freq(n ,m, 0) = self%hsEN(n, 0) -  self%hsEN(m, 0)
            ENDDO
        ENDDO

        !11 electronic subspace
        DO w = 1, self%dimension
            DO v = 1, self%dimension
                self%tran_freq(v ,w, 1) = self%hsEN(v, 1) -  self%hsEN(w, 1)
            ENDDO
        ENDDO

    end subroutine


    subroutine potential_surface(self, string, occupation)

        CLASS(system), intent(inout)                   :: self
        CHARACTER(len=*), intent(in)                :: string
        INTEGER, intent(in)                         :: occupation
        CHARACTER(len=1000)                         :: outputFile

        INTEGER                                     :: i, j !< Loop Variable

        outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"potenial.surf"

        IF(self%surf_bool.eq.1) THEN
            OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

            ASSOCIATE(bias_end => self%parameter_grid_length, &
                bias => self%parameter_grid)

                DO  i = 1, bias_end

                    self%bias_voltage = bias(i)*Vm2au
                    DO j = 1, self%N
                        WRITE(16, '(E,E,E)') self%grid(j), bias(i), self%potential(self%grid(j), self%id + occupation)
                    END DO
                    WRITE(16, '(A)') ' '

                END DO

            END ASSOCIATE
            CLOSE(16)

            WRITE(*,*) "Potential Surface is written in: ", outputFile

        ELSE

            WRITE(*,*) "Surface bool is 0. No data file is created"

        END IF

    end subroutine potential_surface

    REAL(8) FUNCTION potential(self, x, id)
        IMPLICIT NONE

        CLASS(system) :: self
        !---Position
        REAL(8), intent(in)                           ::x
        INTEGER, intent(in)                           ::id

        !Better Readbility
        ASSOCIATE(  &

            !1st Set:Double Well
            Vb => self%Vb,            &
            l => self%l,              &
            delta => self%delta,      &
            x_shift => self%x_shift,   &
            shiftEnergy => self%shiftEnergy, &

            !2nd Set:Double Well
            Vb_2 => self%Vb_2,            &
            l_2 => self%l_2,              &
            delta_2 => self%delta_2,      &

            !Electric Field Gate
            U => self%gate_voltage, &
            d => self%d, &
            !Electric Field Biash
            V => self%bias_voltage, &
            Lj => self%Lj, &
            phi => self%angle, &
            !Vibronic
            Omega => self%Omega, &
            lambda => self%lambda, &
            parabola_shift => self%parabola_shift)

            SELECT CASE(id)

                !Double Well, regular where the occupied potential is shifted in energy and position
                CASE(0)
                    potential = ((0.5d0*(x + l)*delta)) / l + ( (Vb - 0.5d0 * delta)*(x + l)*(x + l)*(x - l)*(x - l)/ (l*l*l*l)) - ((U/d))*x
                CASE(1)
                    potential = (((0.5d0*((x-x_shift) + l)*delta)) / l + &
                        ( (Vb - 0.5d0 * delta)*((x-x_shift) + l)*((x-x_shift) + l)*((x-x_shift) - l)*((x-x_shift) - l)/ (l*l*l*l))) &
                        - ((U/d))*(x-x_shift) + shiftEnergy
                !Harmonic
                CASE(2)
                    !Ivans Model DELTA = 0.5, lambdda 0.3, epsilon0 = 0.5
                    potential = 0.5d0*(protonMass)*Omega*Omega*x**2
                CASE(3)
                    !potential = 0.5d0*(protonMass)*Omega*Omega*x**2 + lambda*DSQRT(2*protonMass*Omega)*x + parabola_shift
                    potential = 0.5d0*(protonMass)*Omega*Omega*(x-x_shift)**2 + parabola_shift

                !Double Well Diode, means that here the bias voltage V is included. Moreover
                CASE(10)
                    potential = ((0.5d0*(x + l)*delta)) / l + ( (Vb - 0.5d0 * delta)*(x + l)*(x + l)*(x - l)*(x - l)/ (l*l*l*l)) &
                        - ((U/d))*x*Sin(phi) - (V/Lj)*x*Cos(phi)
                CASE(11)
                    potential = (((0.5d0*((x-x_shift) + l)*delta)) / l + ( (Vb - 0.5d0 * delta)*((x-x_shift) + l)*((x-x_shift) + l)*((x-x_shift) - l)*((x-x_shift) - l)/ (l*l*l*l))) &
                        - ((U/d))*(x-x_shift)*Sin(phi) -  (V/Lj)*(x-x_shift)*Cos(phi) + shiftEnergy

                !Double Well where two different Sets for the unoccupied and the occupied are used.
                CASE(20)
                    potential = ((0.5d0*(x + l)*delta)) / l + ( (Vb - 0.5d0 * delta)*(x + l)*(x + l)*(x - l)*(x - l)/ (l*l*l*l)) &
                        - ((U/d))*x*Sin(phi) - (V/Lj)*x*Cos(phi)
                CASE(21)
                    potential = (((0.5d0*((x-x_shift) + l_2)*delta_2)) / l_2 + ( (Vb_2 - 0.5d0 * delta_2)*((x-x_shift) + l_2)*((x-x_shift) + l_2)*((x-x_shift) - l_2)*((x-x_shift) - l_2)/ (l_2*l_2*l_2*l_2))) &
                        - ((U/d))*(x-x_shift)*Sin(phi) -  (V/Lj)*(x-x_shift)*Cos(phi) + shiftEnergy

                !Unoccupied: Double Well| Occupied: Harmonic Oscillator
                CASE(30)
                    potential = ((0.5d0*(x + l)*delta)) / l + ( (Vb - 0.5d0 * delta)*(x + l)*(x + l)*(x - l)*(x - l)/ (l*l*l*l)) &
                        - ((U/d))*x*Sin(phi) - (V/Lj)*x*Sin(phi)
                CASE(31)
                    potential = 0.5d0*(protonMass)*Omega*Omega*x**2 + lambda*DSQRT(2*protonMass*Omega)*x + parabola_shift

                CASE(100)
                    potential = -((U/d))*x
               
                !This potentials are shifted simultanously by x_shift to study the effect of non-zero pointed electrical field 
                CASE(200)
                    potential = ((0.5d0*(x + l)*delta)) / l + ( (Vb - 0.5d0 * delta)*(x + l)*(x + l)*(x - l)*(x - l)/ (l*l*l*l)) - ((U/d))*(x-x_shift)
                CASE(201)
                    potential = ((0.5d0*(x + l)*delta)) / l + ( (Vb - 0.5d0 * delta)*(x + l)*(x + l)*(x - l)*(x - l)/ (l*l*l*l)) - &
                    ((U/d))*(x-x_shift) - ((U/d))*(x-x_shift) + shiftEnergy

            END SELECT

        END ASSOCIATE

    END FUNCTION potential



    subroutine destructor(self)
        IMPLICIT NONE
        type(system), intent(in) :: self




    end subroutine


    SUBROUTINE write_x_energy(self, occupation, string, number)

        IMPLICIT NONE

        CLASS(system):: self

        INTEGER, intent(in) :: number
        INTEGER, intent(in) :: occupation

        CHARACTER(len=*) :: string

        !Additional Strings to Filename
        CHARACTER(len=1000)                          ::outputFile

        !Loop Variables
        INTEGER                                     ::m, v

        IF((occupation.eq.0).or.(occupation.eq.1)) THEN

            IF(self%plot_bool.eq.1) THEN

                !OutputFile
                outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//".xen"

                OPEN(UNIT=9, FILE= outputFile, ACTION="write", STATUS="replace")
                DO m = 1, number
                    WRITE (9, '(1000E)') (Abs(self%X_energybase(m,v,occupation))*Abs(self%X_energybase(m,v,occupation)), v=1, number)
                END DO
                CLOSE(9)

            END IF

        END IF

    END SUBROUTINE write_x_energy

    subroutine write_system_to_file(self, number_wave, string)
        IMPLICIT NONE

        CLASS(system), intent(in) ::self

        INTEGER, intent(in)                                    :: number_wave
        CHARACTER(len=*), intent(in)                            :: string

        !Outputfile
        CHARACTER(len=1000)                          ::outputFile_occ
        CHARACTER(len=1000)                          ::outputFile_unocc
        !DO LOOP
        INTEGER                                     ::i,j

        IF(self%plot_bool.eq.1) THEN

            outputFile_occ =trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_occ.sys"
            outputFile_unocc =trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_unocc.sys"

            !WRITE unoccupied File
            OPEN(UNIT=9, FILE=outputfile_unocc  , ACTION="write", STATUS="replace")

                DO i = 1,self%N
                 WRITE(9, '(10000E)') self%grid(i), self%hsEN(i,0)*au2eV,  self%hsEN(i,1)*au2eV, &
                       self%potential(self%grid(i), self%id)*au2ev, self%potential(self%grid(i), self%id + 1)*au2eV, (self%hsOV(i,j, 0)*self%hsOV(i,j, 0)*0.5d0 + self%hsEN(j,0)*au2ev , j=1,number_wave)
        !          WRITE(9, '(10000E)') self%grid(i), self%hsEN(i,0)*au2eV,  self%hsEN(i,1)*au2eV, &
        !               self%potential(self%grid(i), self%id)*au2ev, self%potential(self%grid(i), self%id + 1)*au2eV, (self%hsOV(i,j, 0), j=1,number_wave)

                END DO

            CLOSE(9)
            WRITE(*,'(A,A)') 'Systemdata for UNOCCUPIED has been written to file: ', outputFile_unocc

            !WRITE Occupied File
!           OPEN(UNIT=19, FILE=outputfile_occ  , ACTION="write", STATUS="replace")

!               DO i = 1,self%N
!                   WRITE(19, '(10000E)') self%grid(i), self%hsEN(i,0)*au2eV,  self%hsEN(i,1)*au2eV, &
!                       self%potential(self%grid(i), self%id)*au2ev, self%potential(self%grid(i), self%id + 1)*au2eV, (self%hsOV(i,j, 1) + self%hsEN(j,1)*au2ev , j=1,number_wave)
!               END DO

!           CLOSE(19)
!           WRITE(*,'(A,A)') 'Systemdata for OCCUPIED has been written to file: ', outputFile_occ


        ELSE

            WRITE(*,'(A)') 'Plotbool is 0, thus no sytem summary file is generated'

        END IF

    end subroutine write_system_to_file

    subroutine change_gate(self, gate_voltage)

        CLASS(system), intent(inout)::self
        REAL(8), intent(in) :: gate_voltage

        self%gate_voltage = gate_voltage*Vm2au

        CALL self%set_system
        CALL self%set_x_in_energybase

        CALL self%compute_transition_frequencies
    
    end subroutine change_gate

    subroutine change_gate_and_bias_field(self,gate_voltage, bias_voltage)

        CLASS(system), intent(inout)::self
        REAL(8), intent(in) :: bias_voltage
        REAL(8), intent(in) :: gate_voltage

        self%bias_voltage = bias_voltage*Vm2au
        self%gate_voltage = gate_voltage*Vm2au
                   
        CALL self%set_system
        CALL self%set_x_in_energybase

        CALL self%compute_transition_frequencies

    end subroutine change_gate_and_bias_field

    subroutine put_console_info(self)

        CLASS(system), intent(in) :: self

        WRITE(*,*)
        WRITE(*,'(A)') 'SYSTEM'
        WRITE(*,'(A)')  '------------------------------------------------------------'
        WRITE(*,'(A, T60)') '---Parameters in SI'
        WRITE(*,'(A, F6.3)') 'l =', self%l*au2ang
        WRITE(*,'(A, F6.3)') 'delta = ', self%delta*au2ev
        WRITE(*,'(A, F6.3)') 'Vb = ', self%Vb*au2ev
        WRITE(*,'(A, F6.3)') 'shiftEnergy= ', self%shiftEnergy*au2ev
        WRITE(*,'(A, F6.3)') 'x_shift = ', self%x_shift
        WRITE(*,'(A, F6.3)') 'Current Gate Voltage = ', self%gate_voltage*au2Vm
        WRITE(*,'(A)')       '------------------------------------------------------------'

    end subroutine put_console_info

end module class_system
