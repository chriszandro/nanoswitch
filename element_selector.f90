!* Module element_selecter
!> @brief: Constructing the band of the density matrix.
!> @detail: It is object-oriented written.
!> @param rho the denstiy matrix
MODULE element_selector

    IMPLICIT NONE
    PRIVATE

    TYPE, PUBLIC :: matrix_elements

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

        !Tulp Array for the indices
        INTEGER, ALLOCATABLE, DIMENSION(:,:,:)       :: tulpVector

    CONTAINS

        PROCEDURE           :: set => set_elements
        PROCEDURE           :: print => print_set
        PROCEDURE           :: tulpel_out => out_tulpelmatrix
        PROCEDURE           :: set_copy => set_copy
        PROCEDURE           :: set_matrix => set_matrix
        PROCEDURE           :: check_trace => check_trace
        PROCEDURE           :: set_ltensor => set_ltensor

    END TYPE


CONTAINS


    SUBROUTINE set_elements(this, medim00, medim11, meBND)

        CLASS(matrix_elements), intent(out) :: this

        !Parameters
        INTEGER, intent(in)                             :: medim00
        INTEGER, intent(in)                             :: medim11
        INTEGER, intent(in)                             :: meBND

        !Density Vector
        !Stationary
        ! COMPLEX(8), ALLOCATABLE, DIMENSION(:), intent(inout)            ::rho_stat
        !Time evolution
        ! COMPLEX(8), ALLOCATABLE, DIMENSION(:), intent(inout)            ::rho_t
        !Initial state before switching process
        ! COMPLEX(8), ALLOCATABLE, DIMENSION(:), intent(inout)            ::rho_initial
        !Initial state after switching process
        ! COMPLEX(8), ALLOCATABLE, DIMENSION(:), intent(inout)            ::rho_0

        INTEGER                                         :: minValue
        INTEGER                                         :: maxValue

!Kommentar
        !Initialzation
        minValue = min(medim11-1, medim00-1)
        maxValue = max(medim00, medim11)

        IF( meBND.le.minValue) THEN

            this%medim00 = medim00
            this%medim11 = medim11
            this%meBND = meBND

            !Calculated
            this%offband00 = -0.5d0*meBND*(1 + meBND - 2 *this%medim00)
            this%offband11 = -0.5d0*meBND*(1 + meBND - 2 *this%medim11)

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

            !Set Density Matrix
            !--- stationary
            !      ALLOCATE(rho_stat(1:this%meDim))
            !--- time evaluating
            !       ALLOCATE(rho_t(1:this%meDim))
            !--- Initial
            !       ALLOCATE(rho_initial(1:this%meDim))
            !--- Initial
            !       ALLOCATE(rho_0(1:this%meDim))

            !Console output of the structur data
            CALL this%print

            !For printing the data tuple structure into a file then uncomment:
            !CALL this%tulpel_out

        ELSE

            WRITE(*,*)
            WRITE(*,*) '--------------------------------------------------'
            WRITE(*,*) 'No valid values for medim0, medim1 and meBND given'
            WRITE(*,*) 'They have to obey meBND < min(medim1-1, medim0-1).'
            WRITE(*,*) 'Program exit'
            WRITE(*,*) '--------------------------------------------------'
            WRITE(*,*)
            CALL EXIT(0)

        END IF

    END SUBROUTINE set_elements

    SUBROUTINE set_ltensor(this, Ltensor)

        CLASS(matrix_elements), intent(in) :: this
        !Stationary
        COMPLEX(8), ALLOCATABLE, DIMENSION(:,:), intent(inout)            ::Ltensor

        ALLOCATE(Ltensor(1:this%medim,1:this%medim))


    END SUBROUTINE set_ltensor


    !*SUBROUTINE set_matrix
    !> @brief Subroutine allocates density vectors according to the defined structure of %this
    SUBROUTINE set_matrix(this, rho)

        CLASS(matrix_elements), intent(in) :: this
        !Stationary
        COMPLEX(8), ALLOCATABLE, DIMENSION(:), intent(inout)            ::rho

        ALLOCATE(rho(1:this%meDim))

    END SUBROUTINE set_matrix


    SUBROUTINE print_set(this)

        CLASS(matrix_elements), intent(in) :: this

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



    !Subroutine set_copy
    !@brief Copying the corresponding matrix elements from the objects matrix to a smaller one
    !@note The subroutine first creates map between the indeces of the large and the small matrix.
    !@todo zcopy nutzen

    SUBROUTINE set_copy(this, rho_origin, rho_target, density_small)

        COMPLEX(8), intent(in)               ::rho_origin(1:)
        COMPLEX(8), intent(inout)            ::rho_target(1:)

        CLASS(matrix_elements), intent(in) :: this
        CLASS(matrix_elements), intent(in) :: density_small

        INTEGER, ALLOCATABLE, DIMENSION(:)  :: map

        INTEGER                             :: i,j

        ALLOCATE(map(1:density_small%medim))

        ASSOCIATE(  low00 => density_small%lower_non_diag_bound00,  &
            up00 => density_small%upper_non_diag_bound00,   &

            low11 => density_small%lower_non_diag_bound11,  &
            up11 => density_small%upper_non_diag_bound11,   &

            allm_target => density_small%NumDim00,              &
            allv_target => density_small%NumDim11,              &

            allm_origin => this%NumDim00, &
            allv_origin => this%NumDim11,  &

            seperator => density_small%seperator,         &

            meDim00 => density_small%meDim00,             &
            meDim11 => density_small%meDim11,             &
            meDim => density_small%meDim,                 &

            !Indizes
            !00 Subspace
            uM1 => this%tulpVector(:,1,0),        &
            uM2 => this%tulpVector(:,2,0),        &
            m1 => density_small%tulpVector(:,1,0),       &
            m2 => density_small%tulpVector(:,2,0),       &

            !11 Subspace
            uV1 => this%tulpVector(:,1,1),        &
            uV2 => this%tulpVector(:,2,1),         &
            v1 => density_small%tulpVector(:,1,1),       &
            v2 => density_small%tulpVector(:,2,1) )

            !Mapping Algorithm
            IF(this%medim.eq.density_small%medim) THEN

                rho_target = rho_origin

            ELSE

                !Mapping 00
                !If the indeces of both tulpvectors coincides the mapping stores:
                ! i-th component of small density matrix -> j-th compoment of large density matrix
                smalldensity: DO i = 1, allm_target

                    largedensity: DO j = 1, allm_origin

                        IF((m1(i).eq.uM1(j)).and.(m2(i).eq.uM2(j))) THEN

                            map(i) = j

                            CYCLE smalldensity

                        END IF

                    END DO largedensity

                END DO smalldensity

                !Mapping 11
                smalldensity: DO i = 1, allv_target

                    largedensity: DO j = 1, allv_origin

                        IF((v1(i).eq.uV1(j)).and.(v2(i).eq.uV2(j))) THEN

                            map(i+ allm_target) = j + allm_origin

                            CYCLE smalldensity

                        END IF

                    END DO largedensity

                END DO smalldensity

            END IF

            !Copy Loops

            DO i = 1, allm_origin

                rho_target(i) = rho_origin(map(i))

            END DO

            DO i = 1, allv_target

                rho_target(i + allm_target) = rho_origin(map(i))

            END DO

        END ASSOCIATE

    END SUBROUTINE set_copy





    SUBROUTINE tulp_set(this)

        CLASS(matrix_elements), intent(inout) :: this

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

        !Output of TulpVector; Optional
        !CALL out_tulpelmatrix(this)

    END SUBROUTINE tulp_set


    SUBROUTINE out_tulpelmatrix(this)

        IMPLICIT NONE

        CLASS(matrix_elements), intent(in) :: this

        !tulpvector( medim x 2, 0,1 )
        INTEGER, DIMENSION(:,:), ALLOCATABLE                 ::tempTulpVector

        !Loop Variables
        INTEGER                                             :: i

        !Initialization
        !---Allocation
        ALLOCATE(tempTulpVector(1:this%meDim, 1:2))


        !Copying
        !---00
        tempTulpVector(1:this%NumDim00,1) = this%tulpVector(:, 1, 0)
        tempTulpVector(1:this%NumDim00,2) = this%tulpVector(:, 2, 0)

        !---11
        tempTulpVector((this%NumDim00+1):this%meDim ,1) = this%tulpVector(: , 1, 1)
        tempTulpVector((this%NumDim00+1):this%meDim ,2) = this%tulpVector(: , 2, 1)

        OPEN(13, file = './output/mastereq/tulpelVector.out',  ACTION="write", STATUS="replace")

        DO i = 1, this%meDim

            IF(i.EQ.1) THEN
                WRITE(13,'(A)') 'M -------------'
            END IF

            WRITE(13, '(I4,I4)'), tempTulpVector(i,1),tempTulpVector(i,2)

            IF(i.EQ.this%NumDim00) THEN
                WRITE(13,'(A)') 'V -------------'
            END IF

        END DO
        CLOSE(13)

        WRITE(*,*) "Tulpelvector has been written"

    END SUBROUTINE out_tulpelmatrix


!*SUBROUTINE check_trace
!> @brief Subroutine to check if the trace of the density vector is one
!> @note The density vector in the input have to be allocated by the same object otherwise it leads to an error .

    SUBROUTINE check_trace(this, rho)

        CLASS(matrix_elements), intent(in) :: this !< this object
        COMPLEX(8), intent(in)               ::rho(1:) !< incoming density vector to check
        REAL(8) :: trace !< calculated trace

        INTEGER :: i

        trace = 0.d0

        ASSOCIATE(  all_m => this%NumDim00,              &
            all_v => this%NumDim11,              &
            diag_m => this%medim00,          &
            diag_v => this%medim11 )

            !Trace computation
            !00 - Subspace
            DO i = 1, diag_m
                trace = trace+ REAL(rho(i),8)
            END DO
            !11 - Subspace
            DO i = 1, diag_v
                trace = trace + REAL(rho(i + all_m),8)
            END DO

            !Check
            IF(Abs(1-trace).ge.0.1d0) THEN
                WRITE(*,*) '****************************************************'
                WRITE(*,'(A, E)') 'Stop. Trace of checked density matrix is not equal 1', trace
                WRITE(*,*) '****************************************************'
                CALL EXIT(0)
            ELSE
                WRITE(*,'(A, E)') 'Normalization of density matrix is preserved', trace
            END IF


        END ASSOCIATE

    END SUBROUTINE check_trace

END MODULE element_selector
