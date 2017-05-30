module class_coupling

    use class_system
    IMPLICIT NONE

    type coupling

        !Switch Function with Tanh[]
        REAL(8)                                    :: A                         !upper step
        REAL(8)                                    :: B                         !low step
        REAL(8)                                    :: F                         !Smoothing Factor
        REAL(8)                                    :: switchShift               !Horizontal Shift


        !COUPLING MATRIX
        !<m| V_d | v >
        !Is Symmetric.
        !( left:right lead , m-space, v-space  )
        REAL(8), ALLOCATABLE, DIMENSION(:, :, :), public            :: jnVdM            !(m,v, ln)           Coupling tabulated as matrix

        CHARACTER(len=1000)                                 :: sOutput

        !Grid Copy
        INTEGER                                             :: N !< Number of points of associated grid
        INTEGER                                             :: dimension!< Dimension of the densityvector
        REAL(8), ALLOCATABLE, DIMENSION(:)                  :: grid!<Copy of the associated grid received of a grid object in the Instantation

        !Switch Function Table
        REAL(8), ALLOCATABLE, DIMENSION(:)                  :: switch_table!< Table of the switch function on the i-th grid point of self%grid

        !Plot Bool
        INTEGER                                             :: coupling_bool

    contains

        procedure, public :: init_coupling
        procedure, private :: get_switch_table
        procedure, private :: set_coupling
        procedure, public :: switch_func
        procedure, private :: put_console_info
        procedure, public :: write_couplings
        procedure, public :: update

    END TYPE


CONTAINS

    subroutine init_coupling(self, sys, grid_points, density, input)

        USE class_system
        USE class_grid
        USE class_density
        USE class_inputdata

        CLASS(coupling):: self
        CLASS(system)::sys
        CLASS(grid)::grid_points
        CLASS(densityvector)::density
        CLASS(inputdata)::input

        !Get parameters from input
        self%A = input%A
        self%B = input%B
        self%F = input%F
        self%switchShift = input%switchShift
        self%sOutput = input%system_output_filepath


        self%coupling_bool = input%coupling_bool

        !Copy the grid
        self%N = grid_points%N
        ALLOCATE(self%grid(self%N))
        self%grid = grid_points%dvrX

        !Get the information about the density matrix
        self%dimension = density%maxDiagValue

        !Allocate the switch table and calculate it
        ALLOCATE(self%switch_table(self%N))
        CALL self%get_switch_table

        !First Calc of the couplings
        ALLOCATE(self%jnVdM(1:self%dimension, 1:self%dimension, 1:2))
        CALL self%set_coupling(sys)

        CALL self%put_console_info

    end subroutine init_coupling

    subroutine update(self, sys)

        CLASS(coupling)::self
        CLASS(system)::sys

        CALL self%set_coupling(sys)

    end subroutine update

    subroutine set_coupling(self, sys)

        CLASS(coupling)::self
        CLASS(system)::sys

        REAL(8)                                  :: summation

        REAL(8), DIMENSION(:,:,:), ALLOCATABLE   ::hsOV_copy

        !Loop Subspaces
        INTEGER                                ::m,v      !subspace v, m
        INTEGER                                ::i        !summation

        ALLOCATE(hsOV_copy(sys%N,sys%N,0:1))

        !Copy
        hsOV_copy = sys%hsOV

        !Left Side
        DO v = 1, self%dimension
            DO m = 1, self%dimension
                summation = 0.0d0
                !Loop Unrolling
                DO i = 1, self%N
                    summation = summation + sys%hsOV(i,m,0)*hsOV_copy(i,v,1)*self%switch_table(i)
                END DO
                self%jnVdM(m,v,1) = summation
            ENDDO
        ENDDO

        !Just Copy Left to Right Side
        self%jnVDM(:,:,2) = self%jnVDM(:,:,1)

    end subroutine

    subroutine put_console_info(self)
        CLASS(coupling):: self

        WRITE(*,*)
        WRITE(*,'(A)') 'COUPLING'
        WRITE(*,'(A)')  '------------------------------------------------------------'
        WRITE(*,'(A, E)') 'A : ', self%A
        WRITE(*,'(A, E)') 'B : ', self%B
        WRITE(*,'(A, E)') 'F : ', self%F
        WRITE(*,'(A, E)') 'X_Shift : ', self%switchShift
        WRITE(*,'(A)')  '------------------------------------------------------------'
        WRITE(*,*)

    end subroutine

    subroutine get_switch_table(self)

        CLASS(coupling)::self

        !Loop Variable
        INTEGER i

        !Filling the Coupling_Array
        DO i = 1, self%N
            self%switch_table(i) = self%switch_func(self%grid(i))
        END DO

    end subroutine get_switch_table

    !****************************************************************************************************
    !Subroutine
    !Writes Coupling Factors |<m|V|v>|^2 into a file for better analysis.
    !
    !File extension: .cou
    !File organisation: As the orginal matrix   jnvdmwith dimension number x number
    !Rows:      Unoccupation Number m
    !Columns:   Occupation Number v
    !
    !Attention: If the coupling to the leads is asymatric the subroutine has to be modified
    !****************************************************************************************************

    SUBROUTINE write_couplings(self, occupation, string, number)

        IMPLICIT NONE

        CLASS(coupling):: self

        INTEGER, intent(in) :: number
        CHARACTER(len=*) :: string
        INTEGER, intent(in) :: occupation

        !Additional Strings to Filename
        CHARACTER(len=1000)                          ::outputFile
        CHARACTER(len=1000)                          ::outputFile_log

        !Loop Variables
        INTEGER                                     ::m, v


        IF((occupation.eq.0).or.(occupation.eq.1)) THEN

            IF(self%coupling_bool.eq.1) THEN

                !OutputFile
                outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//".cou"

                OPEN(UNIT=9, FILE= outputFile, ACTION="write", STATUS="replace")
                DO m = 1, number
                    WRITE (9, '(1000E)') (Abs(self%jnVDM(m,v,1))*Abs(self%jnVdM(m,v,1)), v=1, number)
                END DO
                CLOSE(9)

                !OutputFile_log
                outputFile_log = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_log.cou"
                OPEN(UNIT=10, FILE= outputFile_log, ACTION="write", STATUS="replace")
                DO m = 1, number
                    WRITE (10, '(1000E)') ( log10(Abs(self%jnVDM(m,v,1))*Abs(self%jnVdM(m,v,1))), v=1, number)
                END DO
                CLOSE(10)

            END IF

        END IF

        WRITE(*,*) "Frank - Condon Matrix written to: ", outputFile

    END SUBROUTINE write_couplings


    !****************************************************************************************************
    !Function
    !
    !****************************************************************************************************

    REAL(8) FUNCTION switch_func(p, x)

        IMPLICIT NONE

        !Parameters
        REAL(8)                 ::x    !Position
        CLASS(coupling) :: p

        switch_func = ((p%A+p%B)*0.5d0)+((p%A-p%B)*0.5d0)*tanh(p%F*x + p%switchShift)

    END FUNCTION switch_func




end module class_coupling
