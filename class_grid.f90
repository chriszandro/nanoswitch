module class_grid

    USE class_inputdata

    implicit none

    private
    real(8), parameter    :: Pi = 3.1415926535897932385d0
    real(8), parameter    :: gaussianIntegral = DSQRT(Pi)

    TYPE, public :: subgrid_structure
        !Integer for left and right side of the Intervall
        INTEGER                                         :: x_left !< index of left x
        INTEGER                                         :: x_right !< index of right index

        !Integreation & Tunnel/Survival Probability
        INTEGER                                         ::x_integrate_left
        INTEGER                                         ::x_integrate_right
        REAL(8)                                         ::x_range_integrate=4.0d0
        INTEGER                                         ::x_index_zero !< index of right index
        
        INTEGER                                         :: index_range_x !< range of indices from x_left to x_right
        REAL(8)                                         :: xrange !< x range

        REAL(8), ALLOCATABLE, DIMENSION(:,:)            :: dvrX_subset_weights !> Weightsh
        REAL(8), ALLOCATABLE, DIMENSION(:)              :: position
        REAL(8), ALLOCATABLE, DIMENSION(:)              :: variable
        INTEGER, ALLOCATABLE, DIMENSION(:)              :: index_map
        INTEGER, ALLOCATABLE, DIMENSION(:)              :: index_map_integrate

    END TYPE

    type, public :: grid

        integer, public :: N

        !Grid
        !---Hermitian Grid Points Calculation
        REAL(8), ALLOCATABLE, DIMENSION(:), public                  :: dvrX !< Contains the grid points
        REAL(8), ALLOCATABLE, DIMENSION(:,:), public                  ::eigenvectorX !< Eigenvektors of X for the calculations of the weights
        !---Weights
        REAL(8), ALLOCATABLE, DIMENSION(:,:), public                  ::dvrX_weights !> Weightsh

        type(subgrid_structure) subgrid

    contains
        procedure, private :: set_grid
        procedure, private :: get_weights
        procedure, private :: set_subgrid
        procedure, public :: init_grid
        final :: destructor
    end type grid

contains

    subroutine init_grid(self, input)
        class(grid), intent(inout) :: self
        class(inputdata), intent(in) :: input

        !Fetch data of inputfile
        INTEGER :: N
        N = input%N

        self%N = N
        self%subgrid%xrange = input%xrange

        !Allocate Arrays
        !---Global Grid
        ALLOCATE(self%dvrX(N),self%eigenvectorX(N,N), self%dvrX_weights(N,1:3))
        !---Sub Grid
        ALLOCATE(self%subgrid%index_map(N))
        ALLOCATE(self%subgrid%index_map_integrate(N))


        CALL set_grid(self)
        !Get Weights
        CALL get_weights(self)
        CALL set_subgrid(self)

    end subroutine init_grid

    subroutine set_grid(self)
        class(grid), intent(inout) :: self

        !X Matrix
        REAL(8), ALLOCATABLE, DIMENSION(:)                ::off_diagX             !first off-diagonal of X

        !Lapack MKL Arrays and Variables
        REAL(8), DIMENSION(:), ALLOCATABLE                ::work                !Work Array for Lapack Routine
        INTEGER, DIMENSION(:), ALLOCATABLE                ::iwork
        INTEGER                                           ::info

        !Loop Variables
        INTEGER                                            :: i

        !Initialization
        !---Allocation
        ALLOCATE (work(2*self%N*self%N + 4*self%N + 1),iwork(2*5*self%N*self%N+3))         !Dimension according to MKL reference
        ALLOCATE(off_diagX(self%N - 1))

        !---Array
        off_diagX = 0; self%dvrX = 0; self%eigenVectorX = 0

        !Setting Matrix X in MKL Shape
        !Optimization: Set Matrix M(i,j+1) =0.5d0*I afterwards MKL Routine with DSQRT on every Element)
        DO i = 1, self%N-1
            off_diagX(i) = DSQRT(0.5d0*i)
        END DO

        !Diagonalization of X
        CALL DSTEVD("V", self%N, self%dvrX, off_diagX, self%eigenVectorX, self%N, work, SIZE(work), iwork, SIZE(iwork), info)


        !Correctness Checks
        WRITE(*,*)
        WRITE(*,'(5X,A)'), 'Hermitian Grid Points (one-time initialization)'
        WRITE(*,'(5X,A)') '--------------------------------------------------------'

        IF(info.eq.0) THEN
            WRITE(*,'(5X,A,I6,A,T60,A)') 'Calculation of ', self%N ,' Hermitian grid points', 'ok'
            WRITE(*,'(5X,A,F10.3,A,F10.3)') 'From', self%dvrX(1), ' to ',self%dvrX(self%N)
        ELSE
            WRITE(*,'(5X,A,I6,A,T60,A)') 'Calculation of ',self%N ,' Hermitian grid points', 'wrong'
            WRITE(*,*)  'Something wrong in Module systemdvr1 Subroutine Init_HermitePivots'
            WRITE(*, *) 'Info variable of diagonalization: ',self%N
            WRITE(*, *) 'Program Exit'
            CALL EXIT(0)
        END IF

    end subroutine set_grid

    subroutine get_weights(self)
        class(grid), intent(inout) :: self

        INTEGER i

        DO i = 1,self%N
            self%dvrX_weights(i,1) = (self%eigenvectorX(1,i)**2 * gaussianIntegral) * exp((self%dvrX(i)**2))
            self%dvrX_weights(i,2) = 1/(self%dvrX_weights(i,1))
            self%dvrX_weights(i,3) = DSQRT(self%dvrX_weights(i,2))
        END DO

    END SUBROUTINE get_weights

    subroutine set_subgrid(self)
        class(grid), intent(inout) :: self

        !Loop Variables
        INTEGER                               :: i
        INTEGER                               :: indexpointer

        self%subgrid%index_map = 0

        !Search for the range for the position of rho
        DO i = 1,self%N
            !rhox_grid
            IF(self%dvrX(i).le.((-1)*self%subgrid%xrange)) self%subgrid%x_left = i
            IF(self%dvrX(i).le.self%subgrid%xrange) self%subgrid%x_right = i
            !integration_grid 
            IF(self%dvrX(i).le.((-1)*self%subgrid%x_range_integrate)) self%subgrid%x_integrate_left = i
            IF(self%dvrX(i).le.self%subgrid%x_range_integrate) self%subgrid%x_integrate_right = i
            !zero point
            IF(self%dvrX(i).le.0.0d0) self%subgrid%x_index_zero = i
        
        END DO
        
        
        !This is a very subtile + 1 addition ever.
        self%subgrid%index_range_x = ABS(self%subgrid%x_left - self%subgrid%x_right)  + 1

        ALLOCATE(self%subgrid%position(self%subgrid%index_range_x))
        ALLOCATE(self%subgrid%dvrX_subset_weights(self%subgrid%index_range_x, 1:3))


        !Create index map for rhox. It maps a position in the dvrX grid to the subgrid
        indexpointer = 0
        DO i = self%subgrid%x_left, self%subgrid%x_right
            indexpointer = indexpointer + 1
            self%subgrid%index_map(i) = indexpointer
        END DO
        !Create index map for integration. It maps a position in the dvrX grid to the subgrid
        indexpointer = 0
        DO i = self%subgrid%x_integrate_left, self%subgrid%x_integrate_right
            indexpointer = indexpointer + 1
            self%subgrid%index_map_integrate(i) = indexpointer
        END DO
       

        !Write Index Entries for the corresponding positions
        DO i =  self%subgrid%x_left, self%subgrid%x_right
            self%subgrid%position(self%subgrid%index_map(i)) = self%dvrX(i)

            self%subgrid%dvrX_subset_weights(self%subgrid%index_map(i), 1)= self%dvrX_weights(i,1)
            self%subgrid%dvrX_subset_weights(self%subgrid%index_map(i), 2)= self%dvrX_weights(i,2)
            self%subgrid%dvrX_subset_weights(self%subgrid%index_map(i), 3)= self%dvrX_weights(i,3)

        END DO
        !Output
        WRITE(*,*)
        WRITE(*,'(5X, A)') 'Search Routine'
        WRITE(*,'(5X, A)') '--------------------------------------------------------'
        WRITE(*,'(A, T20, F10.3)') "Range : ", self%subgrid%xrange
        WRITE(*,'(A, T20, I)') "Indize Range: ", self%subgrid%index_range_x
        WRITE(*,'(A, T20, F10.3, A, T50, I)') "x_left: ", self%dvrX(self%subgrid%x_left), " at index: ", self%subgrid%x_left
        WRITE(*,'(A, T20, F10.3, A, T50, I)') "x_right: ", self%dvrX(self%subgrid%x_right), " at index: ", self%subgrid%x_right
        WRITE(*,'(5X, A)') '--------------------------------------------------------'
        WRITE(*,*)

    end subroutine set_subgrid

    subroutine destructor(self)
        type(grid), intent(in) :: self
    end subroutine

end module class_grid
