module class_ltensor

    use class_system
    use class_density
    use class_coupling
    use class_junction
    use class_bosonic_bath
    use class_inputdata

    implicit none
    private

    type, public :: ltensor
        private
        COMPLEX(8), ALLOCATABLE, DIMENSION(:,:),public         :: LTensor
        LOGICAL                                                 :: init_flag = .false.
        INTEGER                                                 :: fermion_bath
        INTEGER                                                 :: boson_bath
        INTEGER, public                                         :: dimension

        !Theoretical Validity Check:
        REAL(8), ALLOCATABLE, DIMENSION(:,:,:), private             :: redfield_check
        REAL(8), ALLOCATABLE, DIMENSION(:,:,:), private             :: hbath_checker
        REAL(8), private                                        :: eta

        INTEGER                                                   :: check_dimension
        ! Absolute Values of Frequencies, Tensorelements

        !Output of Tensor
        CHARACTER(len=1000)                                     :: sOutput

    contains
        procedure, private :: compute_stationary_rho
        procedure, public :: get_stationary_rho
        procedure, public  :: get_and_alloc_tensor
        procedure, public :: set_tensor
        procedure, public :: init_tensor
        procedure, public :: check_bath_ltensor
        procedure, public:: analyse_bath_ltensor
        procedure, public  :: write_checked_bath_ltensor
        final :: destructor
    end type ltensor

contains

    subroutine init_tensor(self, input, density)

        CLASS(ltensor):: self
        CLASS(densityvector)::density
        CLASS(inputdata):: input

        IF(self%init_flag.eq..false.) THEN

            self%check_dimension = density%maxDiagValue
            self%dimension = density%medim
            ALLOCATE(self%LTensor(1:self%dimension, 1:self%dimension))
            ALLOCATE(self%hbath_checker(1:self%check_dimension, 1:self%check_dimension, 0:1))
            ALLOCATE(self%redfield_check(1:self%check_dimension, 1:self%check_dimension, 0:1))

            self%sOutput = input%system_output_filepath

            self%init_flag = .true.
            self%eta = input%eta

        ELSE
            WRITE(*,*) 'Ltensor already initialized'
        END IF

    end subroutine init_tensor

    subroutine get_stationary_rho(self, density, sys, junction, coupl, hbath, input)

        class(ltensor), intent(inout) :: self
        class(system), intent(in)::sys
        class(leads), intent(in)::junction
        class(densityvector), intent(inout)::density
        class(coupling), intent(in)::coupl
        class(bath), intent(in):: hbath
        class(inputdata), intent(in) :: input

        self%fermion_bath = input%fermion_bath
        self%boson_bath = input%boson_bath

        CALL self%set_tensor(density, sys, junction, coupl, hbath)
        CALL self%compute_stationary_rho(density)

    end subroutine get_stationary_rho

    subroutine get_and_alloc_tensor(self, matrix)

        class(ltensor), intent(out) :: self
        COMPLEX(8), intent(out)            :: matrix(1:,1:)

        matrix = self%ltensor

    end subroutine get_and_alloc_tensor

    subroutine set_tensor(self, density, sys, junction, coupl, hbath)

        class(ltensor), intent(inout) :: self
        class(system), intent(in)::sys
        class(leads), intent(in)::junction
        class(densityvector), intent(in)::density
        class(coupling), intent(in)::coupl
        class(bath), intent(in):: hbath

        !Loop
        INTEGER                                      ::i,j
        !Loop - Summation
        INTEGER                                     ::m, v

        COMPLEX(8)  :: sumLeft, sumRight    !! Left and Right Fermionic Lead summation variables
        COMPLEX(8)  :: sumBosonic !! Bosonic Summation Variables


        IF(self%init_flag.eq..true.) THEN

            !---Initialization
            self%LTensor = dcmplx(0.0d0,0.0d0)

            !WRITE(*,*)
            !WRITE(*,'(A)') 'MASTER EQUATION '
            !WRITE(*,'(A)')  '------------------------------------------------------------'
            !Tulpvector indizes of column = mapped ColumnIndizes
            !Indize of Row = actual Indize

            ASSOCIATE(  low00 => density%lower_non_diag_bound00,  &
                up00 => density%upper_non_diag_bound00,   &
                low11 => density%lower_non_diag_bound11,  &
                up11 => density%upper_non_diag_bound11,   &
                all_m => density%NumDim00,              &
                all_v => density%NumDim11,              &
                seperator => density%seperator,         &
                meDim00 => density%meDim00,             &
                meDim11 => density%meDim11,             &
                meDim => density%meDim,                 &
                !Indize
                m1 => density%tulpVector(:,1,0),        &
                m2 => density%tulpVector(:,2,0),        &
                uM1 => density%tulpVector(:,1,0),       &
                uM2 => density%tulpVector(:,2,0),       &
                v1 => density%tulpVector(:,1,1),        &
                v2 => density%tulpVector(:,2,1),         &
                uV1 => density%tulpVector(:,1,1),       &
                uV2 => density%tulpVector(:,2,1), &
                hsEN => sys%hsEN, &
                jnVDM => coupl%jnVDM,&
                jnSEM => junction%jnSEM, &
                x_energy => sys%x_energybase, &
                bath => hbath%bathGAM)

            !m-subspace

            !L0^00
            ! (+|*)
            ! (*|*)

            !Loop along sub-diagonal

            !Elements i = 1, medim0-1 are of the form (1,1), (2,2)...(medim0,medim0)
            !* i = (medim0-1, NumDim0)


            !WRITE(*, '(A)') '---Set LTensor'

            !L_0^00
            ! (+|*)
            ! (*|*)
            !WRITE(*, '(A)') '---L_0^00'

            DO i= 1, all_m
                self%LTensor(i,i) = dcmplx(0.0d0 , hsEN(m2(i),0) - hsEN(m1(i),0))
            END DO


            !L_1^00
            ! (+|*)
            ! (*|*)
            !WRITE(*, '(A)') '---L_1^00'

            DO j=1, all_m

                DO i=1, all_m

                    !IF(uM1(j).eq.m1(i)) LTensor(i,j) = LTensor(i,j) + dcmplx(0.0d0 , 11.0d0)
                    !IF(uM2(j).eq.m2(i)) LTensor(i,j) = LTensor(i,j) + dcmplx(0.0d0 , -11.0d0)

                    IF(uM1(j).eq.m1(i)) THEN

                        sumLeft = 0.0d0
                        DO v = 1, meDim11
                            sumLeft = sumLeft + jnSEM(uM2(j), v, 1, 1)*jnVDM(uM2(j),v, 1)*jnVDM(m2(i),v, 1)
                        END DO

                        sumRight = 0.0d0
                        DO v = 1, meDim11
                            sumRight = sumRight + jnSEM(uM2(j), v, 2, 1)*jnVDM(uM2(j),v, 2)*jnVDM(m2(i),v, 2)
                        END DO

                        self%LTensor(i,j) = self%LTensor(i,j) + dcmplx(0.0d0, 1.0d0)*(sumLeft + sumRight)

                    END IF

                    IF(uM2(j).eq.m2(i)) THEN

                        sumLeft = 0.0d0
                        DO v = 1, meDim11
                            sumLeft = sumLeft + dconjg(jnSEM(uM1(j), v, 1, 1))*jnVDM(uM1(j),v, 1)*jnVDM(m1(i), v, 1)
                        END DO

                        sumRight = 0.0d0
                        DO v = 1, meDim11
                            sumRight = sumRight +  dconjg(jnSEM(uM1(j), v, 2, 1))*jnVDM(uM1(j),v, 2)*jnVDM(m1(i),v, 2)
                        END DO

                        self%LTensor(i,j) = self%LTensor(i,j) + dcmplx(0.0d0, -1.0d0)*(sumLeft + sumRight)

                    END IF

                END DO

            END DO

            !L_1^01
            ! (*|+)
            ! (*|*)
            !WRITE(*, '(A)') '---L_1^01'

            DO j= 1, all_v !Columns
                DO i= 1, all_m !Rows

                    ! LTensor(i, j + all_m) = dcmplx(0.0d0, 1.0d0)

                    !Left
                    self%LTensor(i, j + all_m) = self%LTensor(i, j + all_m) + dcmplx(0.0d0,1.0d0)*dconjg(jnSEM(m1(i),uV1(j), 1, 2))*jnVdM(m1(i),uV1(j),1)*jnVdM(m2(i),uV2(j),1) &
                        -  dcmplx(0.0d0,1.0d0)*jnSEM(m2(i),uV2(j),1,2)*jnVdM(m1(i),uV1(j),1)*jnVdM(m2(i),uV2(j),1)

                    !Right
                    self%LTensor(i, j + all_m) = self%LTensor(i, j + all_m) + dcmplx(0.0d0,1.0d0) * dconjg(jnSEM(m1(i),uV1(j),2,2))*jnVdM(m1(i),uV1(j),2)*jnVdM(m2(i),uV2(j),2) &
                        - dcmplx(0.0d0,1.0d0) * jnSEM(m2(i),uV2(j),2,2) * jnVdM(m1(i),uV1(j),2)*jnVdM(m2(i),uV2(j),2)

                END DO
            END DO

            !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


            !v-subspace
            !L0^11

            !WRITE(*, '(A)') '---L_0^11'

            DO i = 1, all_v
                self%LTensor(i + all_m, i + all_m) = dcmplx(0.0d0 , hsEN(v2(i),1) - hsEN(v1(i),1))
            END DO

!            !L_1^11
!            ! (*|*)
!            ! (*|+)
!
            !WRITE(*, '(A)') '---L_1^11'



            DO j= 1, all_v
                DO i= 1, all_v
                    !                    !IF(uV1.eq.v1) LTensor(i,j) = LTensor(i,j) + dcmplx(0.0d0 , 23.0d0)
                    !                    !IF(uV2.eq.v2) LTensor(i,j) = LTensor(i,j) + dcmplx(0.0d0 , -23.0d0)

                    IF(uV1(j).eq.v1(i)) THEN

                        sumLeft = 0.0d0
                        DO m = 1, meDim00
                            sumLeft = sumLeft + jnSEM(m, uV2(j), 1, 2)*jnVDM(m,uV2(j),1)*jnVDM(m, v2(i),1)
                        END DO

                        sumRight = 0.0d0
                        DO m = 1, meDim00
                            sumRight = sumRight + jnSEM(m, uV2(j), 2, 2)*jnVDM(m, uV2(j), 2)*jnVDM(m, v2(i), 2)
                        END DO

                        self%LTensor(i + all_m, j + all_m) = self%LTensor(i + all_m, j + all_m) + dcmplx(0.0d0, 1.0d0)*(sumLeft + sumRight)

                    END IF

                    IF(uV2(j).eq.v2(i)) THEN

                        sumLeft = 0.0d0
                        DO m = 1, meDim00
                            sumLeft = sumLeft + dconjg(jnSEM(m, uV1(j), 1, 2))*jnVDM(m, uV1(j),1)*jnVDM(m, v1(i), 1)
                        END DO

                        sumRight = 0.0d0
                        DO m = 1, meDim00
                            sumRight = sumRight + dconjg(jnSEM(m, uV1(j), 2, 2))*jnVDM(m, uV1(j), 2)*jnVDM(m, v1(i), 2)
                        END DO

                        self%LTensor(i + all_m,j + all_m) = self%LTensor(i + all_m,j + all_m) + dcmplx(0.0d0, -1.0d0)*(sumLeft + sumRight)

                    END IF

                END DO
            END DO

                        !L_1^10
                        ! (*|*)
                        ! (+|*)


            !WRITE(*, '(A)') '---L_1^10'

            DO j=1, all_m
                DO i= 1, all_v

                    !LTensor(i + all_m, j) = dcmplx(0.0d0, 55.0d0) zum Testen

                    !Left
                    self%LTensor(i + all_m,j) = self%LTensor(i + all_m, j) + dcmplx(0.0d0,1.0d0) * dconjg(jnSEM(uM1(j),v1(i),1,1))*jnVdM(uM1(j),v1(i),1)*jnVdM(uM2(j),v2(i),1) &

                        - dcmplx(0.0d0,1.0d0) * jnSEM(uM2(j),v2(i),1,1) * jnVdM(uM1(j),v1(i),1)*jnVdM(uM2(j),v2(i),1)

                    !Right
                    self%LTensor(i + all_m, j) = self%LTensor(i + all_m, j) + dcmplx(0.0d0,1.0d0) * dconjg(jnSEM(uM1(j),v1(i),2,1))*jnVdM(uM1(j),v1(i),2)*jnVdM(uM2(j),v2(i),2) &
                        - dcmplx(0.0d0,1.0d0) * jnSEM(uM2(j),v2(i),2,1) * jnVdM(uM1(j),v1(i),2)*jnVdM(uM2(j),v2(i),2)

                END DO
            END DO
                 
                 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                 !FERMIONIC PART END
                 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                 !**************************************************************
                 !**************************************************************
                 !BOSONIC PART START
                 !**************************************************************
                 !**************************************************************
                 !UNOCCUPIED SYSTEM
                 !L_1^00
                 ! (+|*)
                 ! (*|*)
                 !WRITE(*, '(A)') '---L_1^00'

                IF(self%boson_bath.eq.1) THEN

                    DO j=1, all_m
                        DO i=1, all_m

                            !Gamma- part
                            IF(uM1(j).eq.m1(i)) THEN

                                sumBosonic = 0.0d0
                                DO v = 1, meDim11
                                    sumBosonic = sumBosonic + bath(v,uM2(j),0,1) * x_energy(uM2(j),v,0)*x_energy(v,m2(i),0)*(-1.d0)
                                END DO
                                self%LTensor(i,j) = self%LTensor(i,j) + sumBosonic

                            END IF

                            self%LTensor(i,j) = self%LTensor(i,j) + bath(m2(i),uM2(j),0,1) * x_energy(uM2(j),m2(i),0)*x_energy(m1(i),uM1(j),0)

                            !Gamm+ part
                            IF(uM2(j).eq.m2(i)) THEN

                                sumBosonic = 0.0d0
                                DO v = 1, meDim11
                                    sumBosonic = sumBosonic + bath(uM1(j),v,0,0) * x_energy(m1(i),v,0)*x_energy(v,uM1(j),0)*(-1.d0)
                                END DO

                                self%LTensor(i,j) = self%LTensor(i,j)+ sumBosonic

                            END IF

                            self%LTensor(i,j) = self%LTensor(i,j) + bath(uM1(j), m1(i),0,0) * x_energy(uM2(j),m2(i),0)*x_energy(m1(i),uM1(j),0)

!                            IF((m1(i).eq.m2(i)).and.(uM1(j).eq.uM2(j))) THEN
!                                WRITE(*,*) "CHECK"
!                            END IF
!                            WRITE(*,'(A,I,A,I)') "m1: ", m1(i), "m2: ", m2(i)
!                            WRITE(*,'(A,I,A,I)') "M1: ", uM1(j), "M2: ", uM2(j)
!                            WRITE(*,'(A,E)') "MatrixElement: ", self%Ltensor(i,j)

                        END DO

                    END DO

                    !OCCUPIED SYSTEM
                    !            !L_1^11
                    !(*|*)
                    !(*|+)
                    DO j= 1, all_v
                        DO i= 1, all_v

                            !Gamma -
                            IF(uV1(j).eq.v1(i)) THEN

                                sumBosonic = 0.0d0
                                DO m = 1, meDim00
                                    sumBosonic = sumBosonic + bath(m,uV2(j),1,1) * x_energy(uV2(j),m,1)*x_energy(m,v2(i),1)*(-1.d0)
                                END DO

                                self%LTensor(i + all_m, j + all_m) = self%LTensor(i + all_m, j + all_m) + sumBosonic

                            END IF

                            self%LTensor(i + all_m, j + all_m) = self%LTensor(i + all_m, j + all_m) +  bath(v2(i),uV2(j),1,1) * x_energy(uV2(j),v2(i),1)*x_energy(v1(i),uV1(j),1)

                            !Gamma +
                            IF(uV2(j).eq.v2(i)) THEN

                                sumBosonic = 0.0d0
                                DO m = 1, meDim00
                                    sumBosonic = sumBosonic + bath(uV1(j),m,1,0)*x_energy(v1(i),m,1)* x_energy(m,uV1(j),1)*(-1.d0)
                                END DO

                                self%LTensor(i + all_m,j + all_m) = self%LTensor(i + all_m,j + all_m) + sumBosonic

                            END IF

                            self%LTensor(i + all_m, j + all_m) = self%LTensor(i + all_m, j + all_m) +  bath(uV1(j), v1(i),1,0) * x_energy(uV2(j),v2(i),1)*x_energy(v1(i),uV1(j),1)

!                            IF((v1(i).eq.v2(i)).and.(uV1(j).eq.uV2(j))) THEN
!                                WRITE(*,*) "CHECK"
!                            END IF
!
!                            WRITE(*,'(A,I,A,I)') "v1:", v1(i), "v2:", v2(i)
!                            WRITE(*,'(A,I,A,I)') "V1:", uV1(j), "V2:", uV2(j)
!                            WRITE(*,'(A,E)') "MatrixElement: ", self%Ltensor(i + all_m,j + all_m)

                        END DO
                    END DO

                END IF

            END ASSOCIATE


        ELSE

            WRITE(*,*) 'Ltensor is not initialized yet'

        END IF



    end subroutine set_tensor

    subroutine check_bath_ltensor(self, density)

        class(ltensor), intent(inout) :: self
        class(densityvector), intent(in)::density

        !Loop
        !Loop
        INTEGER            :: i,j

             !Get Elements
        ASSOCIATE(  low00 => density%lower_non_diag_bound00,  &
            up00 => density%upper_non_diag_bound00,   &
            low11 => density%lower_non_diag_bound11,  &
            up11 => density%upper_non_diag_bound11,   &
            all_m => density%NumDim00,              &
            all_v => density%NumDim11,              &
            seperator => density%seperator,         &
            meDim00 => density%meDim00,             &
            meDim11 => density%meDim11,             &
            meDim => density%meDim,                 &
            !Indize
            m1 => density%tulpVector(:,1,0),        &
            m2 => density%tulpVector(:,2,0),        &
            uM1 => density%tulpVector(:,1,0),       &
            uM2 => density%tulpVector(:,2,0),       &
            v1 => density%tulpVector(:,1,1),        &
            v2 => density%tulpVector(:,2,1),         &
            uV1 => density%tulpVector(:,1,1),       &
            uV2 => density%tulpVector(:,2,1) )

            !00 Subspace
            DO j=1, all_m

                IF(uM1(j).eq.uM2(j)) THEN

                    DO i=1, all_m

                        IF(m1(i).eq.m2(i)) THEN
                            self%redfield_check(m1(i),uM1(j),0) = self%ltensor(i,j)
                        END IF

                    END DO

                END IF
            END DO

            !11 Subspace
            DO j=1, all_v

                IF(uV1(j).eq.uV2(j)) THEN

                    DO i=1, all_v

                        IF(v1(i).eq.v2(i)) THEN

                            self%redfield_check(v1(i),uV1(j),1) = self%ltensor(i + all_m,j + all_m)

                        END IF

                    END DO

                END IF

            END DO

        END ASSOCIATE

    end subroutine check_bath_ltensor


    subroutine analyse_bath_ltensor(self,sys, test_pass)

        class(ltensor), intent(inout) :: self
        class(system), intent(in) :: sys

        INTEGER, intent(inout)    :: test_pass

        INTEGER                   :: i,j
        INTEGER                   :: occupation
        REAL(8)                   :: factor = 10.0d0
 
        test_pass = 0
        self%hbath_checker = 0.0d0

        WRITE(*,*) "analyse_bath_ltensor"

        !Analysis
        DO occupation = 0,1
            DO i = 1, self%check_dimension
                DO j = 1, self%check_dimension

                    IF(self%redfield_check(i,j,occupation).le.(factor*Abs(sys%tran_freq(i,j,occupation)))) THEN
                        self%hbath_checker(i,j,occupation) = 1.0d0
                    ELSE
                        test_pass = 1
                        WRITE(*,*) 'No Pass'
                    END IF

                END DO
            END DO
        END DO

        IF(test_pass.eq.1) THEN
                WRITE(*,*) "******************************"
                WRITE(*,*) "******************************"
                WRITE(*,*) "Hbath Test not passed"
                WRITE(*,*) "******************************"
                WRITE(*,*) "******************************"
        END IF

    end subroutine
!
!Die Spalten entsprechen den Großbuchstaben bzw. den Zielzuständen
!
    subroutine write_checked_bath_ltensor(self, sys)

        class(ltensor), intent(inout) :: self
        class(system), intent(in) :: sys

        CHARACTER(len=1000)                :: outputFile11Red
        CHARACTER(len=1000)                :: outputFile00Red

        CHARACTER(len=1000)                :: outputFile11trans_freq
        CHARACTER(len=1000)                :: outputFile00trans_freq

        CHARACTER(len=1000)                :: outputFile11Check
        CHARACTER(len=1000)                :: outputFile00Check

        INTEGER                            :: i,j
        CHARACTER(len=20)                  :: eta_string
        CHARACTER(len=30)                   :: string_suffix

        WRITE(eta_string, '(F10.3)' ) self%eta
        string_suffix = "eta"//trim(adjustl(eta_string))

        outputFile11Red = trim(adjustl(self%sOutput))//trim(adjustl(string_suffix))//"_Red00.hbath"
        outputFile00Red = trim(adjustl(self%sOutput))//trim(adjustl(string_suffix))//"_Red11.hbath"

        outputFile11Check = trim(adjustl(self%sOutput))//trim(adjustl(string_suffix))//"_Check00.hbath"
        outputFile00Check = trim(adjustl(self%sOutput))//trim(adjustl(string_suffix))//"_Check11.hbath"

        outputFile00trans_freq = trim(adjustl(self%sOutput))//trim(adjustl(string_suffix))//"_trans_freq00.hbath"
        outputFile11trans_freq = trim(adjustl(self%sOutput))//trim(adjustl(string_suffix))//"_trans_freq11.hbath"

        !!REDFIELD ELEMENTS
        !!*********************************************
        !!00 subspace
        OPEN(16, file = outputFile00Red, ACTION="write", STATUS="replace")
                DO i = 1, self%check_dimension
                    WRITE(16,'(1000E)') (self%redfield_check(i,j,0), j = 1, self%check_dimension)
                END DO
        CLOSE(16)

        !!11 subspace
        OPEN(17, file = outputFile11Red, ACTION="write", STATUS="replace")
        DO i = 1, self%check_dimension
                    WRITE(17,'(1000E)') (self%redfield_check(i,j,1), j = 1, self%check_dimension)
        END DO
        CLOSE(17)
        !!*********************************************


        !!CHECKER ELEMENTS
        !!*********************************************
        !!00 subspace
        OPEN(16, file = outputFile00Check, ACTION="write", STATUS="replace")
        DO i = 1, self%check_dimension
                   WRITE(16,'(1000E)') (self%hbath_checker(i,j,0), j=1,self%check_dimension)
        END DO
        CLOSE(16)

        !!11 subspace
        OPEN(17, file = outputFile11Check, ACTION="write", STATUS="replace")
        DO i = 1, self%check_dimension
                   WRITE(17,'(1000E)') (self%hbath_checker(i,j,1), j=1,self%check_dimension)
        END DO
        CLOSE(17)

        !!TRANSITION FREQUENCIES
        !!*********************************************
        !!00 subspace
        OPEN(16, file = outputFile00trans_freq, ACTION="write", STATUS="replace")
        DO i = 1, self%check_dimension
                   WRITE(16,'(1000E)') (Abs(sys%tran_freq(i,j,0)), j=1,self%check_dimension)
        END DO
        CLOSE(16)

        !!11 subspace
        OPEN(17, file = outputFile11trans_freq, ACTION="write", STATUS="replace")
        DO i = 1, self%check_dimension
                   WRITE(17,'(1000E)') (Abs(sys%tran_freq(i,j,1)), j=1,self%check_dimension)
        END DO
        CLOSE(17)

        WRITE(*,*) "Results of test are written to: "

        WRITE(*,*) trim(adjustl(outputFile11Red))
        WRITE(*,*) trim(adjustl(outputFile00Red))

        WRITE(*,*) trim(adjustl(outputFile11Check))
        WRITE(*,*) trim(adjustl(outputFile00Check))

        WRITE(*,*) trim(adjustl(outputFile11Check))
        WRITE(*,*) trim(adjustl(outputFile00Check))

    end subroutine write_checked_bath_ltensor

    
    subroutine compute_stationary_rho(self, density)

        class(ltensor), intent(inout) :: self
        class(densityvector), intent(inout)::density

        !Temporary Ltensor
        COMPLEX(8), ALLOCATABLE, DIMENSION(:,:)         :: Ltensor_temp

        !Variables for MKL Subroutines
        !---LU Factorization
        INTEGER                             :: infoFactorize                !Information if Factorization done successfull
        INTEGER                             :: infoLinSolution              !Information if Solving done successfull
        INTEGER, DIMENSION(1:density%medim)         :: ipiv                         !Array what rows interchanged during LU Factorization
        !Error Code
        INTEGER                            ::error_code

        !Loop
        INTEGER                                         ::i,j

        !Initializiation
        density%rho = dcmplx(0.0d0, 0.0d0);
        !VERY IMPORTANT
        density%rho(1) = dcmplx(1.0d0, 0.0d0)

        !Make a working copy of the ltensor
        ALLOCATE(Ltensor_temp(self%dimension, self%dimension))
        Ltensor_temp = self%Ltensor

        ASSOCIATE(  diag_m => density%medim00,          &
            diag_v => density%medim11,          &
            all_m => density%NumDim00,          &
            medim => density%medim)

            !WRITE(*,'(A)') '---Solving master Equation'

            !Including Tr(rho) = 1
            !---Diagonalelemente in 00 - Erste Zeile
            DO j = 1, diag_m
                Ltensor_temp(1,j) = Ltensor_temp(1,j) + 1
            END DO
            !---Diagonalelemente in 11 - Erste Zeile
            DO j = 1, diag_v
                Ltensor_temp(1,j + all_m) = Ltensor_temp(1,j+ all_m) + 1
            END DO


            !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            !Solving the linear equation system
            !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            !LU Factorization
            CALL zgetrf( self%dimension,  self%dimension, Ltensor_temp, self%dimension, ipiv, infoFactorize )

            IF(infoFactorize.ne.0) THEN

                WRITE(*,*) 'Attention: Infor variable of LU Factorization not equals 0'
                WRITE(*,*) 'Something is wrong with the Matrix'
                WRITE(*,*) 'INFO Variable of Factorization: ', infoFactorize

                error_code = 1

            ELSE

                !Matrix Solving
                CALL zgetrs('N', self%dimension, 1, LTensor_temp, self%dimension, ipiv, density%rho, meDim, infoLinSolution )

                IF(infoLinSolution.ne.0) THEN

                    WRITE(*,*) 'Attention: Info variable of Linear System Solver not equals 0'
                    WRITE(*,*) 'Something is wrong with the Matrix'
                    WRITE(*,*) 'INFO Variable of Linear System Solving: ', infoLinSolution

                    error_code = 1

                ELSE
                    !NORMCHECK
                    CALL density%check_trace
                END IF

            END IF

        END ASSOCIATE


    end subroutine compute_stationary_rho

    subroutine destructor(self)
        type(ltensor), intent(in) :: self
    end subroutine

end module class_ltensor
