module class_density

    IMPLICIT NONE
    PRIVATE

    TYPE, PUBLIC :: densityvector

        COMPLEX(8), ALLOCATABLE, DIMENSION(:)               :: rho

        !Tulp Array for the indices
        INTEGER, ALLOCATABLE, DIMENSION(:,:,:)       :: tulpVector

        INTEGER                                             :: meDim00           !< Number of the Diagonalelements 00
        INTEGER                                             :: meDim11           !< Number of the Diagonalelements 11
        INTEGER                                             :: meBND             !< Width of the band centered by the main diagonal

        INTEGER                                             :: offBand00
        INTEGER                                             :: offBand11

        INTEGER                                             :: NumDim00 !< Number of all band elements in the 00 subspace
        INTEGER                                             :: NumDim11 !< Number of all band elements in the 11 subsapce

        INTEGER                                             :: meDim  !< sum of NumDim00 + NumDim11

        !Diagonal Elements in rho
        INTEGER                                             :: maxValue

        !---00
        INTEGER                                             :: lower_diag_bound00
        INTEGER                                             :: upper_diag_bound00

        !---11
        INTEGER                                             :: lower_diag_bound11
        INTEGER                                             :: upper_diag_bound11

        !Non-Diagonal
        !---00
        INTEGER                                             :: lower_non_diag_bound00
        INTEGER                                             :: upper_non_diag_bound00

        !---11
        INTEGER                                             :: lower_non_diag_bound11
        INTEGER                                             :: upper_non_diag_bound11

        !Matrix Constructor
        INTEGER                                             :: seperator
        INTEGER                                             :: LTensor_dimension

        !Max value in diag
        INTEGER                                             :: maxDiagValue


        !Computation
        REAL(8)                                             :: LTensor_storage

        !Initial State Preparation 
        INTEGER                                             :: index_of_coherences 
        INTEGER                                             :: state1 
        INTEGER                                             :: state2 
 
    CONTAINS

        PROCEDURE, private           :: set => set_elements
        PROCEDURE, private           :: print => print_set
        PROCEDURE, public            :: init_densityvector
        PROCEDURE, public            :: set_rho_pure
        PROCEDURE           :: check_trace => check_trace
        PROCEDURE, public   :: set_rho_superposition_pure
        PROCEDURE, public   :: set_rho_superposition_symmetric
    END TYPE

contains

    SUBROUTINE init_densityvector(self, input, selector)
        USE class_inputdata
        CHARACTER(len=1) :: selector
        CLASS(densityvector), intent(out) :: self

        CLASS(inputdata), intent(in) :: input

        self%meDim00 = input%meDim0
        self%meDim11 = input%meDim1

        SELECT CASE(selector)

            CASE("L")
                self%meBND = input%meBND

            CASE("S")
                self%meBND = input%meBND_small
            CASE DEFAULT
                WRITE(*,*) "No proper selection for 'init_densityvector' method"

        END SELECT

        CALL self%set

        CALL self%print

    END SUBROUTINE init_densityvector

    SUBROUTINE set_elements(this)

        CLASS(densityvector), intent(out) :: this

        INTEGER                                         :: minValue
        INTEGER                                         :: maxValue

        !Kommentar
        !Initialzation
        minValue = min(this%medim11-1, this%medim00-1)
        maxValue = max(this%medim00, this%medim11)

        IF( this%meBND.le.minValue) THEN

            !Calculated
            this%offband00 = -0.5d0*this%meBND*(1 + this%meBND - 2 *this%medim00)
            this%offband11 = -0.5d0*this%meBND*(1 + this%meBND - 2 *this%medim11)

            this%NumDim00 = 2*this%offband00 + this%medim00
            this%NumDim11 = 2*this%offband11 + this%medim11

            this%maxValue = max(this%NumDim00, this%NumDim11)

            this%maxDiagValue = maxValue

            !Diagonal Elements
            !---00
            this%lower_diag_bound00 = 1
            this%upper_diag_bound00 = this%meDim00
            !---11
            this%lower_diag_bound11 = 1
            this%upper_diag_bound11 = this%meDim11

            !Non-Diagonal Elements
            !---00
            this%lower_non_diag_bound00 = this%meDim00 + 1
            this%upper_non_diag_bound00 = this%NumDim00

            !---11
            this%lower_non_diag_bound11 = this%meDim11 + 1
            this%upper_non_diag_bound11 = this%NumDim11

            !Total number
            this%meDim = this%NumDim00 + this%NumDim11

            !Matrix Constructors
            this%seperator =  this%NumDim00 + 1

            !Set the Elements
            CALL tulp_set(this)

            !Calculate Tensor Dimension
            this%LTensor_dimension = this%medim*this%medim

            !Calculate storage
            this%LTensor_storage = (REAL(this%medim*this%meDim,4)*16) / Real(1.0d06,4)


            ALLOCATE(this%rho(this%medim))
            this%rho= (0.0d0,0.0d0)

        ELSE

            WRITE(*,*)
            WRITE(*,*) '--------------------------------------------------'
            WRITE(*,*) 'No valid values for medim0, medim1 and this%meBND given'
            WRITE(*,*) 'They have to obey this%meBND < min(medim1-1, medim0-1).'
            WRITE(*,*) 'Program exit'
            WRITE(*,*) '--------------------------------------------------'
            WRITE(*,*)
            CALL EXIT(0)

        END IF

    END SUBROUTINE set_elements

    SUBROUTINE tulp_set(this)

        CLASS(densityvector), intent(inout) :: this

        INTEGER         ::diag
        INTEGER         ::tulpPointer

        !Loop
        INTEGER         :: i


        !Initialization
        !---Allocation (is optimal if both Dimension of 00 and 11 are equal)
        ALLOCATE(this%tulpVector(1:this%maxValue, 1:2 , 0:1))

        tulpPointer = 0

        !00 - Subspace
        ! (+|*)
        ! (*|*)

        DO diag = 0, this%meBND

            DO i = (diag + 1) ,  this%meDim00

                tulpPointer = tulpPointer + 1

                !lower band
                this%tulpVector(tulpPointer,1,0) = i
                this%tulpVector(tulpPointer,2,0) = i-diag

                !upper band
                this%tulpVector(tulpPointer + this%offband00 ,1,0) = i-diag
                this%tulpVector(tulpPointer + this%offband00 ,2,0) = i

            END DO

        END DO

        !        !11 - Subspace
        !        ! (*|*)
        !        ! (*|+)

        tulpPointer = 0

        DO diag = 0, this%meBND
            Do i = (diag + 1) ,  this%meDim11

                tulpPointer = tulpPointer + 1

                !lower band
                this%tulpVector(tulpPointer,1,1) = i
                this%tulpVector(tulpPointer,2,1) = i-diag

                !upper band
                this%tulpVector(tulpPointer + this%offband11 , 1,1) = i-diag
                this%tulpVector(tulpPointer + this%offband11 , 2,1) = i

            END DO
        END DO

    END SUBROUTINE tulp_set

    SUBROUTINE print_set(this)

        CLASS(densityvector), intent(in) :: this

        WRITE(*,*)
        WRITE(*,'(A,T45,I)') 'Bandwidth meBND: ', this%meBND
        WRITE(*,*)
        WRITE(*,'(A,T45,I)') '---00-Subspace'
        WRITE(*,'(A,T45,I)') 'Dimension', this%meDim00
        WRITE(*,'(A,T45,I)') 'Number of Bandelements', 2*this%offband00
        WRITE(*,'(A,T45,I)') 'Number of total matrix elements', this%NumDim00
        WRITE(*,*)
        WRITE(*,'(A,T45,I)') '---11-Subspace'
        WRITE(*,'(A,T45,I)') 'Dimension', this%meDim11
        WRITE(*,'(A,T45,I)') 'Number of Bandelements', 2*this%offband11
        WRITE(*,'(A,T45,I)') 'Number of total matrix elements', this%NumDim11
        WRITE(*,*)
        WRITE(*, '(A,T45,I)') 'Total Number: ', this%meDim

        WRITE(*,'(A)')
        WRITE(*,'(A)') "--- rho-vector"
        WRITE(*, '(A,T45, E10.3)') 'Size in MB: ', (REAL(this%medim,4)*16) / Real(1.0d06,4)
        WRITE(*,'(A)')

        WRITE(*, '(A,T45, E10.3)') 'Matrix Elements: ', REAL(this%medim*this%meDim,4)
        WRITE(*, '(A,T45, E10.3)') 'Bytes: ', REAL(this%medim*this%meDim,4)*16
        WRITE(*, '(A,T45, F10.3)') 'MBytes: ', (REAL(this%medim*this%meDim,4)*16) / Real(1.0d06,4)
        WRITE(*,*)

    END SUBROUTINE print_set

    SUBROUTINE check_trace(this)

        CLASS(densityvector), intent(in) :: this !< this object

        REAL(8) :: trace !< calculated trace

        REAL(8) :: accuracy = 0.01000d0 !< Accuracy of the density matrix trace
        INTEGER :: i

        ASSOCIATE(  all_m => this%NumDim00,              &
            all_v => this%NumDim11,              &
            diag_m => this%medim00,          &
            diag_v => this%medim11 )

            trace = 0.d0

            !Trace computation
            !00 - Subspace
                  DO i = 1, diag_m
                      trace = trace+ REAL(this%rho(i),8)
                  END DO

!           trace = SUM(this%rho(1:diag_m))

            !11 - Subspace
                      DO i = 1, diag_v
                          trace = trace + REAL(this%rho(i + all_m),8)
                     END DO

!           trace = trace + SUM(this%rho((all_m + 1):(all_m + diag_v)))
            !Check
            IF(Abs(1.0d0-trace).ge.accuracy) THEN
                WRITE(*,*) '**************************************************************************'
                WRITE(*,'(A, E)') 'Stop. Trace of checked density matrix is not equal 1: ', trace
                WRITE(*,*) '**************************************************************************'

                CALL EXIT(0)
            ELSE

            WRITE(*,'(A,E)') "Trace is: ",trace

            END IF

        END ASSOCIATE

    END SUBROUTINE check_trace

    SUBROUTINE set_rho_pure(this, state, electronic)

        CLASS(densityvector), intent(inout) :: this
        INTEGER, intent(in)               :: state
        INTEGER, intent(in)               :: electronic

        !Loop Variables
        INTEGER                           :: i, j

        ASSOCIATE(  all_m => this%NumDim00,              &
            all_v => this%NumDim11,              &
            diag_m => this%medim00,          &
            diag_v => this%medim11,          &

            !Indize
            m1 => this%tulpVector(:,1,0),        &
            m2 => this%tulpVector(:,2,0),        &
            v1 => this%tulpVector(:,1,1),        &
            v2 => this%tulpVector(:,2,1))

        WRITE(*,*) '++++++++++++++++++++++++++++++++++'
        WRITE(*,*) 'Pure State is generated'
        WRITE(*,*) '++++++++++++++++++++++++++++++++++'

            WRITE(*,'(A)'), state, electronic
            WRITE(*,*) SIZE(this%rho,1)
            !00 electronic Subspace
            !------------------------------------

            IF(electronic.eq.0) THEN

                this%rho(state)= (1.0d0,0.0d0)

                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"
                WRITE(*,'(A)') "Pure State is generated"
                WRITE(*,'(A)') "Electronic subspace: ", electronic
                WRITE(*,'(A,I,I)') "Indices: ", m1(state), m2(state)
                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"

                ELSE

                !11 electronic Subspace
                !------------------------------------
                this%rho(state + all_m)= (1.0d0,0.0d0)

                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"
                WRITE(*,'(A)') "Pure State is generated"
                WRITE(*,'(A)') "Electronic subspace: ", electronic
                WRITE(*,'(A,I,I)') "Indices: ", v1(state + all_m), v2(state + all_m)
                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"


            END IF



        END ASSOCIATE

    END SUBROUTINE set_rho_pure

    SUBROUTINE set_rho_superposition_arbitrary(this, state1, state2, amplitude_1, amplitude_2, phase, electronic)

        CLASS(densityvector), intent(inout) :: this
        INTEGER, intent(in)               :: state1
        INTEGER, intent(in)               :: state2
        INTEGER, intent(in)               :: electronic

        REAL(8), intent(in)               :: amplitude_1
        REAL(8), intent(in)               :: amplitude_2
        REAL(8), intent(in)               :: phase
        
        !Loop Variables
        INTEGER                           :: i, j

        ASSOCIATE(  all_m => this%NumDim00,              &
            all_v => this%NumDim11,              &
            diag_m => this%medim00,          &
            diag_v => this%medim11,          &

            !Indize
            m1 => this%tulpVector(:,1,0),        &
            m2 => this%tulpVector(:,2,0),        &
            v1 => this%tulpVector(:,1,1),        &
            v2 => this%tulpVector(:,2,1))

        WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++'
        WRITE(*,*) 'Pure State (Superposition) is generated'
        WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++'

            WRITE(*,'(A, I, I)'), "State 1: ",  state1, electronic
            WRITE(*,'(A, I, I)'), "State 2: ",  state2, electronic

            !00 electronic Subspace
            !------------------------------------

            IF(electronic.eq.0) THEN

                !Populations
                this%rho(state1)= Abs(amplitude_1)*Abs(amplitude_1)
                this%rho(state2)= Abs(amplitude_2)*Abs(amplitude_2)
                
                !memorize the state number in the object for getting the
                !elements for output 
                this%state1 = state1 
                this%state2 = state2 

                !Coherences
                DO i = 1, all_m 
                       
                        !Check Indizes 
                        IF((m1(i).eq.state1).and.(m2(i).eq.state2)) THEN
                            this%rho(i) = Abs(amplitude_1)*Abs(amplitude_2)*EXP((0.0d0,1d0)*phase)
                            this%index_of_coherences = i
                        END IF
                        
                        IF((m1(i).eq.state2).and.(m2(i).eq.state1)) THEN
                            this%rho(i) = Abs(amplitude_1)*Abs(amplitude_2)*EXP((0.0d0,-1d0)*phase)
                        END IF
        
                END DO

                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"
                WRITE(*,'(A)') "Pure State is generated"
                WRITE(*,'(A)') "Electronic subspace: ", electronic
                WRITE(*,'(A,I,I)') "Indices 1: ", m1(state1), m2(state1)
                WRITE(*,'(A,I,I)') "Indices 2: ", m1(state2), m2(state2)
                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"

                ELSE

                !11 electronic Subspace
                !------------------------------------
                this%rho(state1 + all_m)= Abs(amplitude_1)*Abs(amplitude_1)
                this%rho(state2 + all_m)= Abs(amplitude_2)*Abs(amplitude_2)
                
                this%state1 = state1 + all_m
                this%state2 = state2 + all_m 

                !Coherences
                DO i = 1, all_v
                        !Check Indizes 
                        IF((v1(i).eq.state1).and.(v2(i).eq.state2)) THEN
                            this%rho(i + all_m) = Abs(amplitude_1)*Abs(amplitude_2)*EXP((0.0d0,1d0)*phase)
                            this%index_of_coherences = i + all_m
                        END IF
                        
                        IF((v1(i).eq.state2).and.(v2(i).eq.state1)) THEN
                            this%rho(i + all_m) = Abs(amplitude_1)*Abs(amplitude_2)*EXP((0.0d0,-1d0)*phase)
                        END IF
                END DO

                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"
                WRITE(*,'(A)') "Pure Superposition State is generated"
                WRITE(*,'(A)') "Electronic subspace: ", electronic
                WRITE(*,'(A,I,I)') "Indices 1: ", v1(state1), v2(state1)
                WRITE(*,'(A,I,I)') "Indices 2: ", v1(state2), v2(state2)
                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"


            END IF


        END ASSOCIATE

    END SUBROUTINE set_rho_superposition_arbitrary

    SUBROUTINE set_rho_superposition_symmetric(this, symmetry, state1, state2, electronic)

        CLASS(densityvector), intent(inout) :: this
        INTEGER, intent(in)               :: state1
        INTEGER, intent(in)               :: state2
        INTEGER, intent(in)               :: electronic
        INTEGER, intent(in)               :: symmetry 

        COMPLEX(8)                           :: coherence

        
        !Loop Variables
        INTEGER                           :: i, j

        ! Set Coherences according to the symmtry paramter
        SELECT CASE (symmetry)
            CASE(0)
                coherence = (0.5d0, 0.0d0)
            CASE(1)
                coherence = (-0.5d0, 0.0d0)
            CASE DEFAULT 
                WRITE(*,'(A)') 'No valid symmetry paramter given in setrhosuperpositionsymmtric'
         END SELECT
        
        ASSOCIATE(  all_m => this%NumDim00,              &
            all_v => this%NumDim11,              &
            diag_m => this%medim00,          &
            diag_v => this%medim11,          &

            !Indize
            m1 => this%tulpVector(:,1,0),        &
            m2 => this%tulpVector(:,2,0),        &
            v1 => this%tulpVector(:,1,1),        &
            v2 => this%tulpVector(:,2,1))

        WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++'
        WRITE(*,*) 'Pure State (Superposition) is generated'
        WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++'

            WRITE(*,'(A, I, I)'), "State 1: ",  state1, electronic
            WRITE(*,'(A, I, I)'), "State 2: ",  state2, electronic

            !00 electronic Subspace
            !------------------------------------

            IF(electronic.eq.0) THEN

                !Populations
                this%rho(state1)= (0.5d0,0.0d0)
                this%rho(state2)= (0.5d0,0.0d0)
                
                this%state1 = state1 
                this%state2 = state2 

                !Coherences
                DO i = 1, all_m 
                       
                        !Check Indizes 
                        IF((m1(i).eq.state1).and.(m2(i).eq.state2)) THEN
                            this%rho(i) = coherence
                            this%index_of_coherences = i
                        END IF
                        
                        IF((m1(i).eq.state2).and.(m2(i).eq.state1)) THEN
                            this%rho(i) = coherence 
                        END IF
        
                END DO

                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"
                WRITE(*,'(A)') "Pure State is generated"
                WRITE(*,'(A)') "Electronic subspace: ", electronic
                WRITE(*,'(A,I,I)') "Indices 1: ", m1(state1), m2(state1)
                WRITE(*,'(A,I,I)') "Indices 2: ", m1(state2), m2(state2)
                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"

                ELSE

                !11 electronic Subspace
                !------------------------------------
                this%rho(state1 + all_m)= (0.5d0,0.0d0)
                this%rho(state2 + all_m)= (0.5d0,0.0d0)
                
                this%state1 = state1 + all_m
                this%state2 = state2 + all_m 

                !Coherences
                DO i = 1, all_v
                       
                        !Check Indizes 
                        IF((v1(i).eq.state1).and.(v2(i).eq.state2)) THEN
                            this%rho(i + all_m) = coherence 
                            this%index_of_coherences = i + all_m
                        END IF
                        
                        IF((v1(i).eq.state2).and.(v2(i).eq.state1)) THEN
                            this%rho(i + all_m) = coherence 
                        END IF
        
                END DO

                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"
                WRITE(*,'(A)') "Pure Superposition State is generated"
                WRITE(*,'(A)') "Electronic subspace: ", electronic
                WRITE(*,'(A,I,I)') "Indices 1: ", v1(state1), v2(state1)
                WRITE(*,'(A,I,I)') "Indices 2: ", v1(state2), v2(state2)
                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"


            END IF


        END ASSOCIATE

    END SUBROUTINE set_rho_superposition_symmetric


    SUBROUTINE set_rho_superposition_pure(this, state1, state2, electronic)

        CLASS(densityvector), intent(inout) :: this
        INTEGER, intent(in)               :: state1
        INTEGER, intent(in)               :: state2
        INTEGER, intent(in)               :: electronic

        !Loop Variables
        INTEGER                           :: i, j

        ASSOCIATE(  all_m => this%NumDim00,              &
            all_v => this%NumDim11,              &
            diag_m => this%medim00,          &
            diag_v => this%medim11,          &

            !Indize
            m1 => this%tulpVector(:,1,0),        &
            m2 => this%tulpVector(:,2,0),        &
            v1 => this%tulpVector(:,1,1),        &
            v2 => this%tulpVector(:,2,1))

        WRITE(*,*) '++++++++++++++++++++++++++++++++++'
        WRITE(*,*) 'Pure State (Superposition) is generated'
        WRITE(*,*) '++++++++++++++++++++++++++++++++++'

            WRITE(*,'(A, I, I)'), "State 1: ",  state1, electronic
            WRITE(*,'(A, I, I)'), "State 2: ",  state2, electronic

            WRITE(*,*) SIZE(this%rho,1)
            !00 electronic Subspace
            !------------------------------------

            IF(electronic.eq.0) THEN

                !Populations
                this%rho(state1)= (0.5d0,0.0d0)
                this%rho(state2)= (0.5d0,0.0d0)
                
                !Coherences
                DO i = 1, all_m 
                       
                        !Check Indizes 
                        IF((m1(i).eq.state1).and.(m2(i).eq.state2)) THEN
                            this%rho(i) = (0.5d0,0.0d0)
                        END IF
                        
                        IF((m1(i).eq.state2).and.(m2(i).eq.state1)) THEN
                            this%rho(i) = (0.5d0,0.0d0)
                        END IF
        
                END DO

                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"
                WRITE(*,'(A)') "Pure State is generated"
                WRITE(*,'(A)') "Electronic subspace: ", electronic
                WRITE(*,'(A,I,I)') "Indices 1: ", m1(state1), m2(state1)
                WRITE(*,'(A,I,I)') "Indices 2: ", m1(state2), m2(state2)
                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"

                

                ELSE

                !11 electronic Subspace
                !------------------------------------
                this%rho(state1 + all_m)= (0.5d0,0.0d0)
                this%rho(state2 + all_m)= (0.5d0,0.0d0)

                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"
                WRITE(*,'(A)') "Pure Superposition State is generated"
                WRITE(*,'(A)') "Electronic subspace: ", electronic
                WRITE(*,'(A,I,I)') "Indices: ", v1(state1 + all_m), v2(state1 + all_m)
                WRITE(*,'(A,I,I)') "Indices: ", v1(state2 + all_m), v2(state2 + all_m)
                WRITE(*,'(A)') "++++++++++++++++++++++++++++++++++"


            END IF



        END ASSOCIATE

    END SUBROUTINE set_rho_superposition_pure

    subroutine destructor(self)
        type(densityvector), intent(in) :: self
    end subroutine

end module class_density
