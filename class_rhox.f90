module class_rhox
    USE UnitConv
    USE class_inputdata
    USE class_density
    USE class_grid
    USE class_system
   ! USE vtk

    implicit none
    private

    type, public :: rhox

        REAL(8), ALLOCATABLE, DIMENSION(:,:) :: rhox !< Array to store the density matrix elements for different parameter values
        
        !Probability Integral 
        REAL(8), ALLOCATABLE, DIMENSION(:,:) :: prob_integral !< Tunneling and Survival Probability
        REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: prob_integral_matrix !< Tunneling and Survival Probability Matrix for Heatmaps
        
        REAL(8), ALLOCATABLE, DIMENSION(:) :: x_grid !< Grid x position copied from the grid object
        INTEGER                            :: x_grid_length
        REAL(8), ALLOCATABLE, DIMENSION(:, :) :: dvr_subset_weights

        REAL(8), ALLOCATABLE, DIMENSION(:) :: parameter_grid !< Paramtergrid, either time or bias/gate voltage
        REAL(8)                            :: parameter_grid_length

        REAL(8), ALLOCATABLE, DIMENSION(:) :: sec_parameter_grid !< Paramtergrid, either time or bias/gate voltage
        REAL(8)                             ::sec_parameter_grid_length
        
        INTEGER                             :: rhox_bool
        INTEGER                             :: plot_bool

        COMPLEX(8), ALLOCATABLE, DIMENSION(:, :) :: transform_matrix !< transformation matrix

        !---Integral
        REAL(8), ALLOCATABLE, DIMENSION(:,:,:)      :: wave_integral
!       REAL(8), ALLOCATABLE, DIMENSION(:,:)      :: wave_matrix
        INTEGER                                     :: wave_number=15
        INTEGER                                     :: wave_matrix_elements
        
        !Error Management
        INTEGER                             :: error_negative=0
        REAL(8)                             :: smallest_value=0.0d0
        REAL(8)                             :: tolerance

        !Error entitys
        INTEGER, ALLOCATABLE, DIMENSION(:)  :: negative_elements

        CHARACTER(len=1000)                 :: sOutput
        CHARACTER(len=1)                    :: mode_specifier

        LOGICAL                             :: init_flag = .false.
    contains
        procedure, public::init_rhox
        procedure, public::write_prob_integral
        procedure, public::get_and_store_rhox
        procedure, private::set_transform_matrix
        procedure, public::write_rhox_tofile_xyz
        procedure, public::update_rhox
        procedure, private::check_negative_elements
        procedure, public::write_rhox_negative_elements
        procedure, public:: write_rhox_tofile_mformat
        procedure, public::tunnel_probability
        procedure, public::put_tunnel_probability_matrix
        procedure, public::write_tunnel_probability_matrix
        procedure, public::get_wave_integral
        procedure, public::write_wave_integral
        
    end type rhox

contains

    subroutine init_rhox(self, input, density, grid_points, sys)

        class(rhox), intent(inout) :: self
        class(inputdata), intent(in)::input
        class(densityvector), intent(in)::density
        class(grid), intent(in)::grid_points
        class(system), intent(in):: sys

        INTEGER i

        IF(self%init_flag.eq..false.) THEN

            ALLOCATE(self%rhox(grid_points%subgrid%index_range_x, input%parameter_grid_length))
            !Probability Integral 
            ALLOCATE(self%prob_integral(input%parameter_grid_length, 0:1))
            ALLOCATE(self%prob_integral_matrix(input%parameter_grid_length, input%sec_parameter_grid_length, 0:1))
            
            !Wave Matrix
            ALLOCATE(self%wave_integral(self%wave_number, 0:1, input%parameter_grid_length))
            !ALLOCATE(self%wave_matrix(self%wave_number, self%parameter_grid_length))
             
            !Error Handling
            ALLOCATE(self%negative_elements(input%parameter_grid_length))
            ALLOCATE(self%transform_matrix(density%medim, grid_points%subgrid%index_range_x))

            !Get a position grid copy
            ALLOCATE(self%x_grid(grid_points%subgrid%index_range_x))
            ALLOCATE(self%parameter_grid(input%parameter_grid_length))
            ALLOCATE(self%dvr_subset_weights(grid_points%subgrid%index_range_x, 1:3))
            self%x_grid_length = grid_points%subgrid%index_range_x
            self%x_grid=grid_points%dvrX

            self%dvr_subset_weights=grid_points%subgrid%dvrX_subset_weights

            !Get output paths + präfix
            self%sOutput = input%result_output_filepath
            self%mode_specifier = input%mode_specifier
            self%rhox_bool = input%rhox_bool
            self%plot_bool = input%plot_bool

            self%parameter_grid_length = input%parameter_grid_length
            self%parameter_grid = input%parameter_grid
            self%sec_parameter_grid_length = input%sec_parameter_grid_length
            self%sec_parameter_grid = input%sec_parameter_grid
            
            !Error Handling
            self%negative_elements = 0
            self%tolerance = input%rhox_tolerance

            CALL self%set_transform_matrix(density, grid_points, sys)
            self%init_flag=.true.

        ELSE

            WRITE(*,*) "Rhox Object is already initialized"

        END IF
    end subroutine init_rhox

    subroutine set_transform_matrix(self, density, grid_points, sys)

        CLASS(rhox)::self
        CLASS(densityvector)::density
        CLASS(system), intent(in):: sys
        CLASS(grid), intent(in)::grid_points

        INTEGER x !< Integer Variable to loop over the indiczes position x in the subgrid
        INTEGER i !< Integer Variable to loop over the indices of the densityvector

        ASSOCIATE(  all_m => density%NumDim00,              &
            all_v => density%NumDim11,              &
            diag_m => density%medim00,          &
            diag_v => density%medim11,          &

            !Indize
            m1 => density%tulpVector(:,1,0),        &
            m2 => density%tulpVector(:,2,0),        &
            v1 => density%tulpVector(:,1,1),        &
            v2 => density%tulpVector(:,2,1), &
            x_left => grid_points%subgrid%x_left, &
            x_right => grid_points%subgrid%x_right, &
            index_map => grid_points%subgrid%index_map)

            DO x = x_left, x_right

                !--- Matrixelements in 00
                DO i = 1, all_m
                    self%transform_matrix(i, index_map(x)) =  sys%hsOV(x,m1(i),0)*sys%hsOV(x,m2(i),0)
                END DO
                !--- Matrixelements in 11
                DO i = 1, all_v
                    self%transform_matrix(i + all_m, index_map(x)) = sys%hsOV(x,v1(i),1)*sys%hsOV(x,v2(i),1)
                END DO

            END DO

        END ASSOCIATE

    end subroutine set_transform_matrix


    subroutine update_rhox(self, density, sys, grid_points)

        CLASS(rhox)::self
        CLASS(densityvector)::density
        CLASS(system), intent(in):: sys
        CLASS(grid), intent(in)::grid_points

        CALL self%set_transform_matrix(density, grid_points, sys)

    end subroutine update_rhox

    subroutine tunnel_probability(self, indexx, grid_points)

        CLASS(rhox), intent(inout) :: self
        INTEGER, intent(in) :: indexx
        CLASS(grid)::grid_points

        INTEGER x!<Loop over the indices of the position x
        REAL(8) summationL
        REAL(8) summationR


        ASSOCIATE( x_left => grid_points%subgrid%x_left, &
            x_right => grid_points%subgrid%x_right, &
            index_map => grid_points%subgrid%index_map, &
            x_zero =>  grid_points%subgrid%x_index_zero, &
            index_range => grid_points%subgrid%index_range_x)

            !Right
            summationR = 0.0d0
            DO x = x_zero, x_right
                summationR = summationR + ( grid_points%dvrX(x+1) -  grid_points%dvrX(x))*(self%rhox(index_map(x+1), indexx)+self%rhox(index_map(x), indexx))
            END DO
            summationR = 0.5d0*summationR
            self%prob_integral(indexx, 1)  = summationR

            !Left
!           summationL = 0.0d0
!           DO x = x_left, x_zero 
!               summationL = summationL + ( grid_points%dvrX(x+1) -  grid_points%dvrX(x))*(self%rhox(index_map(x+1), indexx)+self%rhox(index_map(x), indexx))
!           END DO
!           summationL = 0.5d0*summationL

            summationL = 1.d0 - summationR
            self%prob_integral(indexx, 0) = summationL

        END ASSOCIATE

    end subroutine tunnel_probability


    subroutine put_tunnel_probability_matrix(self, indexx)

        CLASS(rhox), intent(inout) :: self
        INTEGER, intent(in) :: indexx

        self%prob_integral_matrix(:,indexx, 0) = self%prob_integral(:, 0) 
        self%prob_integral_matrix(:,indexx, 1) = self%prob_integral(:, 1) 

    end subroutine put_tunnel_probability_matrix

    subroutine write_tunnel_probability_matrix(self, string)

        CLASS(rhox), intent(inout) :: self
        CHARACTER(len=*), intent(in)                 :: string
        
        CHARACTER(len=1000)                         :: outputFileR
        CHARACTER(len=1000)                         :: outputFileL

        INTEGER                                     :: i, j !< Loop Variable

        outputFileL = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_left_heatmap.prob"
        outputFileR = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_right_heatmap.prob"


        ASSOCIATE(sec_end => self%sec_parameter_grid_length, &
            sec_param => self%sec_parameter_grid,&

            prim_end => self%parameter_grid_length, &
            prim_param => self%parameter_grid)

        OPEN(16, file = outputFileL, ACTION="write", STATUS="replace")
        
            DO  i = 1, prim_end
                WRITE(16, '(10000E)') (self%prob_integral_matrix(i,j,0), j=1,sec_end)
            END DO

        CLOSE(16)
        
        OPEN(17, file = outputFileR, ACTION="write", STATUS="replace")

            DO  i = 1, prim_end
                WRITE(17, '(10000E)') (self%prob_integral_matrix(i,j,1), j=1,sec_end)
            END DO
        
        CLOSE(17)
        END ASSOCIATE

        WRITE(*,*) "Left Heatmap for Tunnel Probability", outputFileL
        WRITE(*,*) "Right Heatmap for Tunnel Probability", outputFileR

    end subroutine write_tunnel_probability_matrix

    subroutine get_wave_integral(self, grid_points, sys, indexx)

        CLASS(rhox), intent(inout) :: self
        CLASS(grid)::grid_points
        CLASS(system), intent(in):: sys
        INTEGER, intent(in) :: indexx
        
        INTEGER x!<Loop over the indices of the position x
        INTEGER j!<Loop over the indices of the position x

        REAL(8) summationR
        REAL(8) summationL

        ASSOCIATE( x_left => grid_points%subgrid%x_left, &
            x_right => grid_points%subgrid%x_right, &
            index_map => grid_points%subgrid%index_map, &
            x_zero =>  grid_points%subgrid%x_index_zero, &
            index_range => grid_points%subgrid%index_range_x)

           DO j = 1, self%wave_number 
                
               !Right
               summationR = 0.0d0
               DO x = x_zero, x_right
                
                  summationR = summationR + ( grid_points%dvrX(x+1) - grid_points%dvrX(x))*&
                  (Abs(sys%hsOV(x+1,j,0)*sys%hsOV(x+1,j,0)) *self%dvr_subset_weights(index_map(x+1),2)  + & 
                   Abs(sys%hsOV(x, j,0)*sys%hsOV(x,j,0)) *self%dvr_subset_weights(index_map(x),2) ) 

               END DO
               summationR = 0.5d0*summationR
               self%wave_integral(j, 0, indexx)  = summationR
               self%wave_integral(j, 1, indexx)  = 1.0d0 - summationR
                
               !Left
            !  summationL = 0.0d0
            !  DO x = x_left, x_right
            !     summationL = summationL + ( grid_points%dvrX(x+1) - grid_points%dvrX(x))*&
            !     (Abs(sys%hsOV(x+1,j,0))*sys%hsOV(x+1,j,0) *self%dvr_subset_weights(index_map(x+1),2)  + & 
            !      Abs(sys%hsOV(x, j,0)*sys%hsOV(x,j,0) *self%dvr_subset_weights(index_map(x),2) )) 

            !  END DO
            !  summationL = 0.5d0*summationL
            !  self%wave_integral(j, 1)  = summationL

           END DO
          
        END ASSOCIATE
    end subroutine get_wave_integral

    subroutine write_wave_integral(self, string)

        CLASS(rhox), intent(inout) :: self
        CHARACTER(len=*), intent(in)                 :: string

        INTEGER                                     ::i,j

        CHARACTER(len=1000)                         :: outputFile
        
        IF(self%plot_bool.eq.0) THEN

            WRITE(*,*) "No wave integral file from rhox class is generated since rhox bool is 0"

        ELSE

        outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_0_rhox.wave"

        OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

            DO i = 1, self%parameter_grid_length
                WRITE(16, '(1000E)') self%parameter_grid(i), (self%wave_integral(j,0,i) , j=1, self%wave_number) 
            END DO

        CLOSE(16)

        WRITE(*,'(A,A)') "Wave matrix outputfile has been written to: ", outputFile

        END IF


    end subroutine write_wave_integral

    subroutine get_and_store_rhox(self, indexx, density, grid_points, sys)

        CLASS(rhox), intent(inout) :: self
        INTEGER, intent(in) :: indexx
        CLASS(densityvector), intent(in)::density
        CLASS(grid)::grid_points
        CLASS(system)::sys

        INTEGER x!<Loop over the indices of the position x
        COMPLEX(8) summation

        COMPLEX(8), external :: zdotu

        INTEGER i,j

        ASSOCIATE(  all_m => density%NumDim00,              &
            all_v => density%NumDim11,              &
            diag_m => density%medim00,          &
            diag_v => density%medim11,          &

            !Indize
            m1 => density%tulpVector(:,1,0),        &
            m2 => density%tulpVector(:,2,0),        &
            v1 => density%tulpVector(:,1,1),        &
            v2 => density%tulpVector(:,2,1), &
            x_left => grid_points%subgrid%x_left, &
            x_right => grid_points%subgrid%x_right, &
            index_map => grid_points%subgrid%index_map, &
            index_range => grid_points%subgrid%index_range_x)

            DO j = 1, index_range
                summation = zdotu(density%medim, density%rho,1, self%transform_matrix(:,j), 1)

                self%rhox(j, indexx) = REAL(summation, 8)*self%dvr_subset_weights(j,2)
            END DO

            !Calculate the tunnel and the survival probability
            CALL self%tunnel_probability(indexx, grid_points)
            !Check Negative Elements in calculated rhox
            self%negative_elements(indexx) = self%check_negative_elements(indexx)
        END ASSOCIATE

    end subroutine

    INTEGER FUNCTION check_negative_elements(self, indexx) RESULT(number_of_elements)

        CLASS(rhox), intent(inout)::self
        INTEGER, intent(in) :: indexx

        INTEGER i

        number_of_elements = 0

        DO i = 1, self%x_grid_length
            IF(self%rhox(i, indexx).le.0.d0) THEN

                IF(Abs(self%rhox(i, indexx)).le.self%tolerance) THEN
                    self%rhox(i, indexx) = 0.0d0

                ELSE
                    self%smallest_value = max(abs(self%rhox(i,indexx)), self%smallest_value)
                    number_of_elements = number_of_elements + 1
                END IF

            END IF
        END DO

        !Counter
        IF(number_of_elements.ge.1) THEN
            self%error_negative = 1
            WRITE(*, *) '**************************************'
            WRITE(*,'(A)') 'Negative Values in rhox: '
            WRITE(*,'(I,I)') number_of_elements, self%x_grid_length
            WRITE(*, *) '**************************************'
        ELSE

            ! WRITE(*,*) 'Check is negative'

        END IF
    end function check_negative_elements

    subroutine write_prob_integral(self, string)

        CLASS(rhox), intent(in)::self
        CHARACTER(len=*), intent(in)                 :: string


        CHARACTER(len=1000)                         :: outputFile

        INTEGER                                     ::i,j

        IF(self%plot_bool.eq.0) THEN

            WRITE(*,*) "No probability integral file from rhox class is generated since rhox bool is 0"

        ELSE

        outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_rhox.prob"

        OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

            DO i = 1, self%parameter_grid_length
                WRITE(16, '(E,E, E)') self%parameter_grid(i), self%prob_integral(i, 0),  self%prob_integral(i, 1)
            END DO

        CLOSE(16)

        END IF

    end subroutine write_prob_integral

    subroutine write_rhox_negative_elements(self, string)

        CLASS(rhox), intent(in)::self
        CHARACTER(len=*), intent(in)                 :: string


        CHARACTER(len=1000)                         :: outputFile

        INTEGER                                     ::i,j

        IF(self%rhox_bool.eq.0) THEN

            WRITE(*,*) "No eper file from rhox class is generated since rhox bool is 0"

        ELSE

        outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_rhox.eper"

        OPEN(16, file = outputFile, ACTION="write", STATUS="replace")

        WRITE(16, '(A, E)') "# Smalles Value in computation: ", self%smallest_value

        DO i = 1, self%parameter_grid_length
            WRITE(16, '(E,E)') self%parameter_grid(i), REAL(self%negative_elements(i),8)
        END DO

        CLOSE(16)

        END IF


    end subroutine write_rhox_negative_elements

    subroutine write_rhox_tofile_xyz(self, string, grid_points)

        CLASS(rhox), intent(in)::self
        CLASS(grid), intent(in)::grid_points
        CHARACTER(len=*), intent(in)                 :: string

        CHARACTER(len=1000)                         :: outputFile
        INTEGER                                     ::i,j

        INTEGER                                     :: x, dt
        REAL(8)                                     :: parameter_converter

        !---Time Grid
        SELECT CASE(self%mode_specifier)

            CASE("V", "G", "B", "W")
                parameter_converter = 1
            CASE("T", "K")
                parameter_converter = au2femto
        END SELECT

        IF(self%rhox_bool.eq.0) THEN

            WRITE(*,*) "No eper file from rhox class is generated since rhox bool is 0"

        ELSE

        outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_xyz.rhox"

        ASSOCIATE( x_left =>  grid_points%subgrid%x_left,  x_right =>  grid_points%subgrid%x_right, index_map => grid_points%subgrid%index_map)

                OPEN(16, file = outputFile, ACTION="write", STATUS="replace")
                !Time Loop
                DO dt = 1, self%parameter_grid_length
                    !Important Blank for gnuplot
                    WRITE(16,'(A)') ' '
                    !Loop over the position on the grid
                    DO x = x_left, x_right
                        WRITE(16, '(E, A, E, A, E)') grid_points%dvrX(x)*au2ang,",", self%parameter_grid(dt)*parameter_converter,",", self%rhox(index_map(x), dt)*ang2au
                    END DO
                END DO
                CLOSE(16)
        END ASSOCIATE

        WRITE(*,*) "rho(x) written in csv Format for gnuplot e.g.: ", outputFile

        END IF
    end subroutine write_rhox_tofile_xyz

    subroutine write_rhox_tofile_mformat(self, string, grid_points)

        CLASS(rhox), intent(in)::self
        CLASS(grid), intent(in)::grid_points
        CHARACTER(len=*), intent(in)                 :: string


        CHARACTER(len=1000)                         :: outputFile
        CHARACTER(len=1000)                         :: outputFile_xgrid
        CHARACTER(len=1000)                         :: outputFile_paragrid_log
        CHARACTER(len=1000)                         :: outputFile_paragrid

        REAL(8)                                     :: parameter_converter

        INTEGER                                     ::i,j

        INTEGER                                     :: x, dt

        !Write Axes Vector


        IF(self%rhox_bool.eq.0) THEN

            WRITE(*,*) "No eper file from rhox class is generated since rhox bool is 0"

        ELSE

        !---Position Grid
        outputFile_xgrid = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//".rhomx"
        ASSOCIATE( x_left =>  grid_points%subgrid%x_left,  x_right =>  grid_points%subgrid%x_right, index_map => grid_points%subgrid%index_map)

            OPEN(14, file = outputFile_xgrid, ACTION="write", STATUS="replace")
            DO x = x_left, x_right
                WRITE(14, '(E)') grid_points%dvrX(x)*au2ang
            END DO
            CLOSE(14)

        END ASSOCIATE

        !---Time Grid
        SELECT CASE(self%mode_specifier)

            CASE("V", "G")
                parameter_converter = 1.0d0
            CASE("T", "K")
                parameter_converter = au2femto
        END SELECT

        outputFile_paragrid= trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//".rhot"
        ASSOCIATE( x_left =>  grid_points%subgrid%x_left,  x_right =>  grid_points%subgrid%x_right, index_map => grid_points%subgrid%index_map)
            OPEN(15, file = outputFile_paragrid, ACTION="write", STATUS="replace")
            DO i = 1, self%parameter_grid_length
                WRITE(15, '(E)') self%parameter_grid(i)*parameter_converter
            END DO
            CLOSE(15)
        END ASSOCIATE

        outputFile_paragrid_log= trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_log.rhot"
        ASSOCIATE( x_left =>  grid_points%subgrid%x_left,  x_right =>  grid_points%subgrid%x_right, index_map => grid_points%subgrid%index_map)
            OPEN(15, file = outputFile_paragrid_log, ACTION="write", STATUS="replace")
            DO i = 1, self%parameter_grid_length
                WRITE(15, '(E)') log(self%parameter_grid(i)*parameter_converter)
            END DO
            CLOSE(15)
        END ASSOCIATE

        !Main File
        ASSOCIATE( x_left =>  grid_points%subgrid%x_left,  x_right =>  grid_points%subgrid%x_right, index_map => grid_points%subgrid%index_map)
            outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//".rhob"
            OPEN(16, file = outputFile, ACTION="write", STATUS="replace")
            DO dt = 1, self%parameter_grid_length
                WRITE(16, '(10000E)') (self%rhox(index_map(x), dt)*ang2au, x = x_left, x_right)
            END DO
            CLOSE(16)
        END ASSOCIATE

        WRITE(*,*) "rho(x) written in m Format: ", outputFile
        WRITE(*,*) "with respective grids"

        END IF

    end subroutine write_rhox_tofile_mformat


    subroutine write_rhox_tofile_vtk(self, string, grid_points)

        CLASS(rhox), intent(in)::self
        CLASS(grid), intent(in)::grid_points
        CHARACTER(len=*), intent(in)                 :: string

        CHARACTER(len=1000)                         :: outputFile
        INTEGER                                     ::i,j

        INTEGER                                     :: x, dt
        REAL(8)                                     :: parameter_converter

     !  TYPE(VTK_file_handle)::fd

     !  CALL VTK_open_file(PREFIX="rhox", FD=fd)
     !  CALL VTK_close_file(FD=fd)

!       CALL VTK_collect_file(FD=fd)


        !---Time Grid
        SELECT CASE(self%mode_specifier)

            CASE("V", "G", "B", "W")
                parameter_converter = 1
            CASE("T", "K")
                parameter_converter = au2femto
        END SELECT

        IF(self%rhox_bool.eq.0) THEN

            WRITE(*,*) "No eper file from rhox class is generated since rhox bool is 0"

        ELSE



        outputFile = trim(adjustl(self%sOutput))//"_"//trim(adjustl(string))//"_xyz.rhox"

        ASSOCIATE( x_left =>  grid_points%subgrid%x_left,  x_right =>  grid_points%subgrid%x_right, index_map => grid_points%subgrid%index_map)

                OPEN(16, file = outputFile, ACTION="write", STATUS="replace")
                !Time Loop
                DO dt = 1, self%parameter_grid_length
                    !Important Blank for gnuplot
                    WRITE(16,'(A)') ' '
                    !Loop over the position on the grid
                    DO x = x_left, x_right
                        WRITE(16, '(E, A, E, A, E)') grid_points%dvrX(x)*au2ang,",", self%parameter_grid(dt)*parameter_converter,",", self%rhox(index_map(x), dt)*ang2au
                    END DO
                END DO
                CLOSE(16)
        END ASSOCIATE

        WRITE(*,*) "rho(x) written in vtk: ", outputFile

        END IF
    end subroutine write_rhox_tofile_vtk
end module class_rhox
