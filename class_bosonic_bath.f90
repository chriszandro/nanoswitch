module class_bosonic_bath

    USE class_system
    USE UnitConv


    implicit none
    TYPE, public :: bath

        INTEGER :: dimension       !<dimension given by the density matrix
        REAL(8) :: temperature     !!< Temparature

        CHARACTER(len=1000) :: sOutput

        !bathGAM
        !( Matrix : electronic subspaces, 0:1, +- )
        REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:)    :: bathGAM

        REAL(8) :: eta !! Coupling paramter system-harmonic bath in the spectral density
        REAL(8) :: wcut !! Cut-off frequency in the spectral density

        !Ancilla Matrix for tabulating Selfenergy containig energy differences
        !( Matrix(n,m) of energy differences, last indice labels the electronic subspace)
        REAL(8), ALLOCATABLE, DIMENSION(:,:,:)                 :: enDiff

    contains
        procedure, public :: init_bath
        procedure, private :: set_endiff
        procedure, private :: set_bath
        procedure, public :: update_bath
        procedure, private :: put_console_info
        procedure, public :: write_bath
        procedure, public :: write_bath_functions
        procedure, public :: update_and_change_temperature

    END TYPE

CONTAINS

    SUBROUTINE init_bath(self, input, density, sys, voltage_in)

        USE class_density
        USE class_inputdata

        CLASS(bath)        ::self
        CLASS(inputdata)    ::input
        CLASS(densityvector)::density
        CLASS(system)       ::sys

        INTEGER(8), intent(in), optional :: voltage_in
        ! Get data from corresponding density-object
        self%dimension = density%maxDiagValue
        self%sOutput = input%system_output_filepath

        WRITE(*,*) input%system_output_filepath

        !---Allocation
        ALLOCATE(self%enDiff(self%dimension,self%dimension, 0:1))
        ALLOCATE(self%bathGAM(1:self%dimension,1:self%dimension,0:1, 0:1))


        ! Get data from inputfile
        self%eta = input%eta
        self%temperature = input%tharmonic
        self%wcut = input%wcut

        CALL self%put_console_info

        !---Initialization of the bath
        CALL self%set_endiff(sys)
        CALL self%set_bath

    END SUBROUTINE init_bath

    subroutine update_and_change_temperature(self, sys, temperature)

        CLASS(bath)        ::self
        CLASS(system)       ::sys
        REAL(8) :: temperature     !!< Temparature....yap....clearly
        
        !---Recalculation of the transition frequencies
        CALL self%set_endiff(sys)
        self%temperature = temperature*Kb
        CALL self%set_bath

    end subroutine update_and_change_temperature
    
    subroutine update_bath(self, sys)

        CLASS(bath)        ::self
        CLASS(system)       ::sys
        
        !---Recalculation of the transition frequencies
        CALL self%set_endiff(sys)
        CALL self%set_bath

    end subroutine

    subroutine set_endiff(self, sys)

        CLASS(bath)::self
        CLASS(system)::sys

        !Loop Variables
        INTEGER n, m
        INTEGER v, w

        !energyDiff Matrix
        !00 electronic subspace
        DO m = 1, self%dimension
            DO n = 1, self%dimension
                self%enDiff(n ,m, 0) = sys%hsEN(n, 0) -  sys%hsEN(m, 0)
            ENDDO
        ENDDO

        !11 electronic subspace
        DO w = 1, self%dimension
            DO v = 1, self%dimension
                self%enDiff(v ,w, 1) = sys%hsEN(v, 1) -  sys%hsEN(w, 1)
            ENDDO
        ENDDO


    end subroutine set_endiff

    SUBROUTINE set_bath(self)

        CLASS(bath), intent(inout):: self

        !Loop-Subspaces
        !00 electronic subspace
        INTEGER                 :: n,m
        !11 electronic subspace
        INTEGER                 :: w,v

        ASSOCIATE(  dimension => self%dimension, &
            enDiff => self%enDiff, &
            T => self%temperature)

            !Gamma Plus
            !---00 electronic subspace
            DO m = 1,dimension
                DO n = 1,dimension
                    self%bathGAM(n,m,0,0) = pi*(1 + hlBE(self%endiff(n,m,0), self%temperature))*bosonic_spectral_density( self%enDiff(n,m,0) , self%wcut, self%eta)
                ENDDO
            ENDDO
            !---11 electronic subspace
            DO w = 1,dimension
                DO v = 1,dimension
                    self%bathGAM(v,w,1,0) = pi*(1 + hlBE(self%endiff(v,w,1), self%temperature))*bosonic_spectral_density( self%enDiff(v,w,1) , self%wcut, self%eta)
                ENDDO
            ENDDO

            !Gamma Minus
            !---00 electronic subspace
            DO m = 1,dimension
                DO n = 1,dimension
                    self%bathGAM(n,m,0,1) = pi* hlBE(self%endiff(n,m,0), self%temperature)*bosonic_spectral_density( self%enDiff(n,m,0) , self%wcut, self%eta)
                ENDDO
            ENDDO

            !---11 electronic subspace
            DO w = 1,dimension
                DO v = 1,dimension
                    self%bathGAM(v,w,1,1) = pi * hlBE(self%endiff(v,w,1), self%temperature)*bosonic_spectral_density( self%enDiff(v,w,1) , self%wcut, self%eta)
                ENDDO
            ENDDO

        END ASSOCIATE
    END SUBROUTINE set_bath

    subroutine put_console_info(self)

        CLASS(bath), intent(in) :: self

        WRITE(*,*)
        WRITE(*,'(A)') 'HARMONIC BATH'
        WRITE(*,'(A)')  '------------------------------------------------------------'
        WRITE(*,'(A,T20, F10.3)') 'Coupling Eta: ',self%eta
        WRITE(*,'(A,T20, F10.3)') 'Cut-Off Frequency: ', self%wcut
        WRITE(*,'(A,T20, F10.3)') 'Temperature: ',self%temperature*(1/Kb)
        WRITE(*,'(A)')  '------------------------------------------------------------'
        WRITE(*,*)

    end subroutine put_console_info

    SUBROUTINE write_bath_functions(self)

        CLASS(bath), intent(in) :: self
        !Loop Variables
        INTEGER                             :: i

        !Omega Grid
        REAL(8)                             :: omega
        REAL(8)                             :: omega_start= 0
        REAL(8)                             :: omega_end = 1
        REAL(8)                             :: omega_step
        INTEGER                             :: grid=1000

        CHARACTER(len=1000)                 :: outputfile

        omega_start = omega_start*ev2au
        omega_end = omega_end*ev2au
        omega_step = ABS(omega_end -omega_start) / grid

        outputfile = trim(adjustl(self%sOutput))//"harmonic_bath.out"

        !Write File
        OPEN(50, file = outputfile, ACTION="write", STATUS="replace")

        DO i = 0, grid

            omega = omega_start + i*omega_step

            WRITE(50, '(E,E,E)') omega*au2eV, hlBE(omega, self%temperature), bosonic_spectral_density(omega, self%wcut, self%eta)

        END DO

        CLOSE(50)

        WRITE(*,'(A,A)') 'Harmonic Bath functions written in:', outputfile

    END SUBROUTINE write_bath_functions


    subroutine write_bath(self)
        CLASS(bath)::self

        !Loop Variables
        INTEGER                             :: i,j

        CHARACTER(len=1000)                 :: outputfile


        !Matrix Output

        !00 subspace----------------------------------------------------------------------------
        OPEN(21, file = trim(adjustl(self%sOutput))//'bathGAM+00.out', ACTION="write", STATUS="replace")
        DO i=1,self%dimension
            WRITE(21,'(10000E)') ( self%bathGAM(i,j,0,0), j=1,self%dimension)
        END DO
        CLOSE(21)

        OPEN(22, file = trim(adjustl(self%sOutput))//'bathGAM-00.out', ACTION="write", STATUS="replace")
        DO i=1,self%dimension
            WRITE(22,'(10000E)') ( self%bathGAM(i,j,0,1), j=1,self%dimension)
        END DO
        CLOSE(22)


        !11 subspace----------------------------------------------------------------------------
        OPEN(21, file = trim(adjustl(self%sOutput))//'bathGAM+11.out', ACTION="write", STATUS="replace")
        DO i=1,self%dimension
            WRITE(21,'(10000E)') ( self%bathGAM(i,j,1,0), j=1,self%dimension)
        END DO
        CLOSE(21)

        OPEN(22, file = trim(adjustl(self%sOutput))//'bathGAM-11.out', ACTION="write", STATUS="replace")
        DO i=1,self%dimension
            WRITE(22,'(10000E)') ( self%bathGAM(i,j,1,1), j=1,self%dimension)
        END DO
        CLOSE(22)

        WRITE(*,'(A,A)') "Harmonic Heath Bath is written in: ", self%sOutput
    end subroutine write_bath

    ! ----------------------------------------- Boson destribution --------------------------------------------------------
    !****************************************************************************************************
    !Subroutine
    !Bose-Einstein distribution in canonical ensemble
    !****************************************************************************************************

    REAL(8) FUNCTION hlBE(energy, T)

        REAL(8), intent(in)                 :: energy       !Energy
        REAL(8), intent(in)                 :: T            !Temperatur

        !For T>0 Bose 
        IF (T.GT.0.0D0) THEN
            IF(energy.GT.0.0D0) THEN
                hlBE = ((1.0D0)/(EXP((energy)/T) -1.0D0))
            ELSE
                hlBE = 1.0D30
            END IF
        ELSE

        END IF

    END FUNCTION hlBE

    ! ---------------------------------------------------------- Spectral Density ------------------------------------------------------------------
    !****************************************************************************************************
    !Subroutine
    !
    !****************************************************************************************************
    REAL(8) FUNCTION bosonic_spectral_density(en, wcut, eta) result(J)

        REAL(8), intent(in)             :: en       !Energy/Frequency
        REAL(8), intent(in)             :: wcut     !Cut-off Frequency
        REAL(8), intent(in)             :: eta     !Coupling parameter

        IF(en.gt.0.0d0) THEN
            J = eta*en*exp(-en/wcut)
        ELSE
            J = 0.d0
        END IF

    END FUNCTION

end module class_bosonic_bath
