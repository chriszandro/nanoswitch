module class_junction

    USE class_system
    USE UnitConv

    implicit none

    TYPE, private :: single_lead

        REAL(8) :: mu              !!< Chemical Potential

        REAl(8) :: beta1           !!< Beta 1
        REAL(8) :: beta2           !!< Beta 2

    END TYPE single_lead

    TYPE, public :: leads

        TYPE(single_lead):: left_lead
        TYPE(single_lead):: right_lead

        REAL(8) :: temperature     !!< Temparature
        REAL(8) :: fermi_level     !!< Fermi-level
        REAL(8) :: voltage

        INTEGER :: dimension       !<dimension given by the density matrix

        CHARACTER(len=1000) :: sOutput

        !SELF ENERGY
        !( left:right lead, function 1 or 2, m-space, v-space, left:right, function 1 or 2)
        COMPLEX(8), ALLOCATABLE, DIMENSION(:,:,:,:)    :: jnSEM


        !Ancilla Matrix for tabulating Selfenergy containig energy differences
        REAL(8), ALLOCATABLE, DIMENSION(:,:)                 :: enDiff
    contains
        procedure, public:: init_leads
        procedure, private :: set_leads
        procedure, private:: set_endiff
        procedure, private ::set_chem_potentials
        procedure, public :: change_voltage
        procedure, public :: change_temperature
        procedure, private :: put_console_info
        procedure, public :: write_leads
        procedure, public :: update

    END TYPE

CONTAINS

    SUBROUTINE init_leads(self, input, density, sys, voltage_in)

        USE class_density
        USE class_inputdata

        CLASS(leads)        ::self
        CLASS(inputdata)    ::input
        CLASS(densityvector)::density
        CLASS(system)       ::sys

        INTEGER(8), intent(in), optional :: voltage_in
        ! Get data from corresponding density-object
        self%dimension = density%maxDiagValue
        ! Get data from inputfile
        self%left_lead%beta1 = input%beta1L
        self%left_lead%beta2 = input%beta2L
        self%right_lead%beta1 = input%beta1R
        self%right_lead%beta2 = input%beta2R
        self%temperature = input%T
        self%fermi_level = input%fermi_level
        self%sOutput = input%system_output_filepath

        ALLOCATE(self%jnSEM(1:self%dimension,1:self%dimension, 1:2, 1:2))
        !---Allocation
        ALLOCATE(self%enDiff(self%dimension,self%dimension))

         IF(present(voltage_in)) THEN

            SELECT CASE(voltage_in)

                CASE(1)

                    self%voltage = input%start_bias_voltage

                CASE(2)

                    self%voltage = input%end_bias_voltage

                CASE DEFAULT

                    WRITE(*,*) "No valid value for the gate_in argument. Default value is taken from the start_gate_voltage variable"

                END SELECT
            END IF

        CALL self%update(sys)

        CALL self%put_console_info

    END SUBROUTINE init_leads

    subroutine set_endiff(self, sys)

        CLASS(leads)::self
        CLASS(system)::sys

        INTEGER v, m

        !energyDiff Matrix
        DO v = 1, self%dimension          !v-subspace
            DO m = 1, self%dimension      !m-subspace
                self%enDiff(m ,v) = sys%hsEN(v, 1) -  sys%hsEN(m, 0)
            ENDDO
        ENDDO

    end subroutine set_endiff

    subroutine change_temperature(self, temperature)

        REAL(8), intent(in) :: temperature 
        CLASS(leads), intent(inout)::self

        self%temperature = temperature*Kb
        CALL self%set_leads

    end subroutine

    subroutine change_voltage(self, voltage)

        REAL(8), intent(in) :: voltage
        CLASS(leads), intent(inout)::self

        self%voltage = voltage
        CALL self%set_chem_potentials
        CALL self%set_leads

    end subroutine

    subroutine set_chem_potentials(self)
        CLASS(leads), intent(inout) ::self

        self%left_lead%mu = self%voltage*ev2au*0.5d0 + self%fermi_level
        self%right_lead%mu = -self%voltage*ev2au*0.5d0 + self%fermi_level


    end subroutine set_chem_potentials


    subroutine update(self, sys)
        CLASS(leads), intent(inout)::self
        CLASS(system), intent(in)::sys

        CALL self%set_endiff(sys)
        CALL self%change_voltage(self%voltage)

    end subroutine update


    SUBROUTINE set_leads(self)

        CLASS(leads), intent(inout):: self

        !Loop-Subspaces
        INTEGER                 :: m,v

        ASSOCIATE(  dimension => self%dimension, &
            muL => self%left_lead%mu, &
            muR => self%right_lead%mu, &
            beta1R => self%right_lead%beta1, &
            beta2R => self%right_lead%beta2, &
            beta1L => self%left_lead%beta1, &
            beta2L => self%left_lead%beta2, &
            enDiff => self%enDiff, &
            T => self%temperature)

            DO v = 1,dimension          !m-subspace
                DO m = 1,dimension      !v-subspace

                    IF ((enDiff(m,v).GE.(muL-2.0d0*beta1L)).AND.(enDiff(m,v).LE.(muL + 2.0d0*beta1L)))THEN
                        self%jnSEM(m,v,1,1) = dcmplx(0.0d0, hlFD(enDiff(m,v), T, muL)*0.5d0*DSQRT(4.0d0*beta1L*beta1L - (enDiff(m,v)-muL)**2) * (beta2L/beta1L)**2)
                    END IF

                ENDDO
            ENDDO


            !epsilon2
            DO v = 1,dimension              !m-subspace
                DO m = 1,dimension          !v-subspace

                    IF ((enDiff(m,v).GE.(muL-2.0d0*beta1L)).AND.(enDiff(m,v).LE.(muL + 2.0d0*beta1L)))THEN

                        self%jnSEM(m,v,1,2) = dcmplx(0.0d0, (1.0d0 - hlFD(enDiff(m,v), T, muL))*0.5d0*DSQRT(4.0d0*beta1L*beta1L - (enDiff(m,v)-muL)**2) * (beta2L/beta1L)**2)

                    END IF

                ENDDO
            END DO

            !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            !Right
            !----------------------
            !epsilon1
            do v = 1,dimension              !m-subspace
                do m = 1,dimension          !v-subspace

                    IF ((enDiff(m,v).GE.(muR-2.0d0*beta1R)).AND.(enDiff(m,v).LE.(muR + 2.0d0*beta1R)))THEN

                        self%jnSEM(m,v,2,1) = dcmplx(0.0d0, hlFD(enDiff(m,v), T, muR)*0.5d0*DSQRT(4.0d0*beta1R*beta1R - (enDiff(m,v)-muR)**2) * (beta2R/beta1R)**2)

                    END IF

                ENDDO
            ENDDO

            !epsilon2
            DO v = 1,dimension              !m-subspace
                DO m = 1, dimension         !v-subspace

                    IF ((enDiff(m,v).GE.(muR-2.0d0*beta1R)).AND.(enDiff(m,v).LE.(muR + 2.0d0*beta1R)))THEN

                        self%jnSEM(m,v,2,2) = dcmplx(0.0d0, (1.0d0 - hlFD(enDiff(m,v), T, muR))*0.5d0*DSQRT(4.0d0*beta1R*beta1R - (enDiff(m,v)-muR)**2) * (beta2R/beta1R)**2)

                    END IF

                ENDDO
            ENDDO

        END ASSOCIATE
    END SUBROUTINE set_leads

    subroutine put_console_info(self)

        CLASS(leads), intent(in) :: self

        WRITE(*,*)
        WRITE(*,'(A)') 'LEADS | SELF-ENERGY'
        WRITE(*,'(A)')  '------------------------------------------------------------'
        WRITE(*,'(A)') '---Left'
        WRITE(*,'(A, T20, F10.3)') 'beta1: ', self%left_lead%beta1*au2ev
        WRITE(*,'(A,T20, F10.3)') 'beta2: ', self%left_lead%beta2*au2ev
        WRITE(*,'(A)') '---Right'
        WRITE(*,'(A, T20, F10.3)') 'beta1: ',self%right_lead%beta1*au2ev
        WRITE(*,'(A, T20, F10.3)') 'beta2: ',self%right_lead%beta2*au2ev
        WRITE(*,*)
        WRITE(*,'(A,T20, F10.3)') 'Temperature: ',self%temperature*(1/Kb)
        WRITE(*,'(A,T20, F10.3)') 'Fermi-Level: ', self%fermi_level*au2ev
        WRITE(*,*)
        WRITE(*,'(A,T20, F10.3)') 'Voltage: ', self%voltage
        WRITE(*,'(A)')  '------------------------------------------------------------'
        WRITE(*,*)

    end subroutine put_console_info

    subroutine write_leads(self)
        CLASS(leads)::self

        !Loop Variables
        INTEGER                             :: i,j

        CHARACTER(len=1000)                 :: outputfile

        !Left Lead----------------------------------------------------------------------------

        OPEN(21, file = trim(adjustl(self%sOutput))//'jnSEM_IM_1.out', ACTION="write", STATUS="replace")
        DO i=1,self%dimension
            WRITE(21,'(1000E)') (AIMAG(self%jnSEM(i,j,1,1)), j=1,self%dimension)
        END DO
        CLOSE(21)

        OPEN(22, file = trim(adjustl(self%sOutput))//'jnSEM_REAL_1.out', ACTION="write", STATUS="replace")
        DO i=1,self%dimension
            WRITE(22,'(1000E)') (REAL(self%jnSEM(i,j,1,1)), j=1,self%dimension)
        END DO
        CLOSE(22)

        OPEN(23, file = trim(adjustl(self%sOutput))//'jnSEM_IM_2.out', ACTION="write", STATUS="replace")
        DO i=1,self%dimension
            WRITE(23,'(1000E)') (AIMAG(self%jnSEM(i,j,1,2)), j=1,self%dimension)
        END DO
        CLOSE(23)

        OPEN(23, file =trim(adjustl(self%sOutput))//'jnSEM_REAL_2.out', ACTION="write", STATUS="replace")
        DO i=1,self%dimension
            WRITE(23,'(1000E)') (REAL(self%jnSEM(i,j,1,2)), j=1,self%dimension)
        END DO
        CLOSE(23)

        !Right Lead----------------------------------------------------------------------------

        OPEN(21, file =trim(adjustl(self%sOutput))//'jnSEM_IM_1.out', ACTION="write", STATUS="replace")
        DO i=1,self%dimension
            WRITE(21,'(1000E)') (AIMAG(self%jnSEM(i,j,2,1)), j=1,self%dimension)
        END DO
        CLOSE(21)

        OPEN(22, file = trim(adjustl(self%sOutput))//'jnSEM_REAL_1.out', ACTION="write", STATUS="replace")
        DO i=1,self%dimension
            WRITE(22,'(1000E)') (REAL(self%jnSEM(i,j,2,1)), j=1,self%dimension)
        END DO
        CLOSE(22)

        OPEN(23, file = trim(adjustl(self%sOutput))//'jnSEM_IM_2.out', ACTION="write", STATUS="replace")
        DO i=1,self%dimension
            WRITE(23,'(1000E)') (AIMAG(self%jnSEM(i,j,2,2)), j=1,self%dimension)
        END DO
        CLOSE(23)

        OPEN(23, file = trim(adjustl(self%sOutput))//'jnSEM_REAL_2.out', ACTION="write", STATUS="replace")
        DO i=1,self%dimension
            WRITE(23,'(1000E)') (REAL(self%jnSEM(i,j,2,2)), j=1,self%dimension)
        END DO
        CLOSE(23)

        WRITE(*,'(A, A)') 'Lead matrices are written to', trim(adjustl(self%sOutput))
    end subroutine write_leads

    ! ----------------------------------------- Fermi-Dirac destribution --------------------------------------------------------
    !****************************************************************************************************
    !Subroutine
    !
    !****************************************************************************************************

    REAL(8) FUNCTION hlFD(energy, T, mu)

        REAL(8), intent(in)                 :: energy       !Energy
        REAL(8), intent(in)                 :: T            !Temperatur
        REAL(8), intent(in)                 :: mu           !Chemical Potential

        !For T>0 Fermi Function
        IF (T.GT.0.0D0) THEN
            hlFD = ( (1.0D0)/ (EXP((energy-mu)/T) +1.0D0) )
        ELSE
            !Stepfunction if T==0
            IF(energy.LT.mu) THEN
                hlFD = 1.0D0
            ELSE IF(energy.GT.mu) THEN
                hlFD = 0.0D0
            ELSE IF(energy.EQ.mu) THEN
                hlFD = 0.5d0
            END IF
        END IF

    END FUNCTION hlFD

    ! ---------------------------------------------------------- Density of states ------------------------------------------------------------------
    !****************************************************************************************************
    !Subroutine
    !
    !****************************************************************************************************
    REAL(8) FUNCTION hlDS(en, mu, beta1)

        REAL(8), intent(in)             :: en       !Energy
        REAL(8), intent(in)             :: mu       !chemical potential
        REAL(8), intent(in)             :: beta1       !band-width

        IF ((en.gt.(mu-2.0d0*beta1)).and.(en.lt.(mu + 2.0d0*beta1))) THEN
            hlDS = (4.0d0*beta1*beta1-(en-mu)**2 )**(-0.5d0)/Pi
        ELSE
            hlDS = 0.0d0
        ENDIF

    END FUNCTION



end module class_junction
