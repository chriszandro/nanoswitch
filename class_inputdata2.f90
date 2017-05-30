!****************************************************************************************************
!Module datain
!
!Module reads in all data needed for computation
!****************************************************************************************************

MODULE datain

    use UnitConv
  ! use junction
  ! use systemdvr1
  ! use couplings

    IMPLICIT NONE

    CHARACTER(len=100) :: sInputparameters             !Name of the inputfile file for physical parameters
    CHARACTER(len=300) :: sInputparameters_copy       !A Copy: Name of the inputfile file for physical parameters
    CHARACTER(len=400) :: sInputpath                  !Inputpath
    CHARACTER(len=400) :: sInputfile


    CHARACTER(len=400) :: sComputationparameters       !Name of inputfile handling computational parameters
    CHARACTER(len=400) :: sComputationparameters_copy !Copy:



    INTEGER                        ::inputNumber   !Number of inputfile data


    !Input File
    !----------------------------------------------------------------------------------------
    !Variable Set for file management
    CHARACTER(len=60)                           :: sInputfiles               !Path of the directory containing the input-files
    CHARACTER(len=60)                           :: sOutputpath               !output Path

    !Potential Parameters
    REAL(8)                                    :: l              !Value V(-l)=0 and V(l)=0 for Vb = 0
    REAL(8)                                    :: delta   !Shift of the right minima
    REAL(8)                                    :: Vb                     !Value V(0)=Vb
    REAL(8)                                    :: shiftEnergy
    REAL(8)                                    :: x_shift

    !Switch Function
    REAL(8)                                    :: A                         !upper step
    REAL(8)                                    :: B                         !low step
    REAL(8)                                    :: F                         !Smoothing Factor
    REAL(8)                                    :: switchshift

    !Vibronic
    REAL(8)                                    :: Omega
    REAL(8)                                    :: lambda
    REAL(8)                                    :: parabola_shift

    !Voltage Grid
    INTEGER                                    :: gridV
    REAL(8)                                    :: voltage1
    REAL(8)                                    :: voltage2
    !Field Grid
    INTEGER                                    :: gridV_gate
    REAL(8)                                    :: voltage_gate_1
    REAL(8)                                    :: voltage_gate_2
    REAL(8)                                    :: d
    !Junctions
    !---Left
    REAL(8)                                    ::beta1L
    REAL(8)                                    ::beta2L
    !---Right
    REAL(8)                                    ::beta1R
    REAL(8)                                    ::beta2R
    !---Both
    REAL(8)                                    ::T
    REAL(8)                                    ::fermi_level
    !---Time Evolution
    INTEGER                                    ::time_grid
    REAL(8)                                    ::time_end
    REAL(8)                                    ::time_start
    !--- Various Parameter Grid
    INTEGER                                    :: grid_vparameter
    REAL(8)                                    :: parameter1
    REAL(8)                                    :: parameter2

    !Computation File
    !----------------------------------------------------------------------------------------
    !Variable set for program flow control
    !---Program Mode
    INTEGER                                     :: mode
    INTEGER                                     :: plot_bool
    !---Dimension
    INTEGER                :: medim1 !< dimension of the subspace 11
    INTEGER                :: medim0 !< dimension of the subspace 00
    INTEGER                :: meBND !< band width for the density matrix
    INTEGER                :: meBnd_small !< band width for the large density matrix for the time evolution calculations

    CHARACTER(len=10)               :: temp1 !<temporal Integer 1 for input from console
    CHARACTER(len=10)                :: temp2 !<temporal Integer 2 for input from console
    CHARACTER(len=10)                :: temp3 !<temporal Integer 3 for input from console
    CHARACTER(len=10)                :: temp4 !<temporal Integer 4 for input from console

    !---Paths
    CHARACTER(len=400)       ::system_output_filepath
    CHARACTER(len=400)       ::summary_output_filepath
    CHARACTER(len=400)       ::result_output_filepath
    !---Logicals
    INTEGER                 ::evo_method
    INTEGER                 ::solve_method
    REAL(8)                 ::abstol
    INTEGER                 ::ideg

    INTEGER                 ::timebool
    REAL(8)                 ::xrange



CONTAINS

    !****************************************************************************************************
    !Subroutine readConsoleInput
    !
    !Subroutine reads input from console
    !****************************************************************************************************

    SUBROUTINE readConsoleInput

        IMPLICIT NONE

        inputNumber = iargc()             !iarg() reading the number of inputfiles

        SELECT CASE (inputNumber)
            CASE (2)
                CALL getarg(1, sInputparameters)
                CALL getarg(2, sComputationparameters)
            CASE(3)
                CALL getarg(1, sInputparameters)
                CALL getarg(2, sComputationparameters)
                CALL getarg(3, sInputpath)
            CASE(5)
                CALL getarg(1, sInputparameters)
                CALL getarg(2, sComputationparameters)
                CALL getarg(3, sInputpath)
                CALL getarg(4, temp1)
                CALL getarg(5, temp2)
             CASE(7)
                CALL getarg(1, sInputparameters)
                CALL getarg(2, sComputationparameters)
                CALL getarg(3, sInputpath)
                CALL getarg(4, temp1)
                CALL getarg(5, temp2)
                CALL getarg(6, temp3)
                CALL getarg(7, temp4)
            CASE DEFAULT
                PRINT *, 'No valid number of files. Program exit'
                CALL exit(0)
        END SELECT

    END SUBROUTINE readConsoleInput

    !****************************************************************************************************
    !Subroutine readConsoleInput
    !
    !Subroutine reads inputfile
    !****************************************************************************************************

    SUBROUTINE readInputPhysics

        IMPLICIT NONE

        sInputfile = trim(adjustl(sInputparameters))

        IF(inputNumber.eq.3) THEN
            sInputfile = trim(adjustl(sInputpath))//trim(adjustl(sInputparameters))
        END IF

        OPEN(unit=100,file=sInputfile, status='old', action='read')  !INPUT-FILME

        READ(100, *)                                            !Seperatora
        READ(100,'(1D20.10)') l                                    !
        READ(100,'(1D20.10)') delta                                !
        READ(100,'(1D20.10)') Vb                                   !
        READ(100,'(1D20.10)') shiftEnergy
        READ(100,'(1D20.10)') x_shift
        READ(100,*)
        READ(100,'(1D20.10)') Omega                                 !
        READ(100,'(1D20.10)') lambda                                 !
        READ(100,'(1D20.10)') parabola_shift
        READ(100,*)
        READ(100,'(1D20.10)') A                                    !
        READ(100,'(1D20.10)') B                                 !
        READ(100,'(1D20.10)') F                                 !
        READ(100,'(1D20.10)') switchShift
        READ(100,*)
        READ(100,'(1I20)') gridV
        READ(100,'(1D20.10)') voltage1                                    !
        READ(100,'(1D20.10)') voltage2                                    !
        READ(100,*)
        READ(100,'(1I20)') gridV_gate
        READ(100,'(1D20.10)') Voltage_Gate_1
        READ(100,'(1D20.10)') Voltage_Gate_2
        READ(100, '(1D20.10)') d
        READ(100,*)
        READ(100,*)
        READ(100,'(1D20.10)') beta1L                                 !
        READ(100,'(1D20.10)') beta2L
        READ(100,*)
        READ(100,*)
        READ(100,'(1D20.10)') beta1R                                 !
        READ(100,'(1D20.10)') beta2R
        READ(100,*)
        READ(100,*)
        READ(100,'(1D20.10)') T
        READ(100,'(1D20.10)') fermi_level
        READ(100,*)
        READ(100,'(1I20)') time_grid
        READ(100,'(1D20.10)') time_start
        READ(100,'(1D20.10)') time_end
        READ(100,*)
        READ(100,'(1I20)') grid_vparameter
        READ(100,'(1D20.10)') parameter1
        READ(100,'(1D20.10)') parameter2
        CLOSE(100)

        !Unit Conversion
        !--Double Well Potential
        l = l*ang2au
        delta = delta*ev2au
        Vb = Vb*ev2au
        shiftEnergy = shiftEnergy*ev2au
        !---Parabolas
        Omega = Omega*ev2au
        lambda = lambda*ev2au
        parabola_shift = parabola_shift*ev2au

        !--Switchfunction
        switchShift = switchShift

        !---Gate Voltage
        !-- conversion in nano meters
        d = d * 1d-9

        !Self-energy
        !---Left
        beta1L = beta1L * ev2au
        beta2L = beta2L * ev2au
        !---Right
        beta1R = beta1R * ev2au
        beta2R = beta2R * ev2au
        !---Left + Right
        T = T * kb
        fermi_level = fermi_level * ev2au


    END SUBROUTINE readInputPhysics


    SUBROUTINE readInputComputation

        IMPLICIT NONE

        IF(inputNumber.eq.3) THEN
            sComputationparameters_copy = trim(adjustl(sInputpath))//trim(adjustl(sComputationparameters))
            sComputationparameters = sComputationparameters_copy
        END IF

        OPEN(unit=101,file=sComputationparameters, status='old', action='read')  !INPUT-FILME

        !System
        READ(101,*)
        READ(101,'(4I20)') mode
        READ(101,'(4I20)') plot_bool
        READ(101,*)                                                   !program-modes
   !     READ(101,'(10I20)') N
   !     READ(101,'(4I20)') id
        READ(101,*)
        READ(101,'(10I20)') medim0
        READ(101,'(10I20)') medim1
        READ(101,'(10I20)') meBnd
        READ(101,'(10I20)') meBnd_small
        READ(101, *)
        READ(101, *)
        READ(101,'(A)') system_output_filepath
        READ(101, *)
        READ(101,'(A)') summary_output_filepath
        READ(101, *)
        READ(101,'(A)') result_output_filepath
        READ(101, *)
        READ(101,'(10I20)') evo_method
        READ(101,'(10I20)') ideg !!< For Pade Approximation.
        READ(101, *)
        READ(101,'(4I20)') solve_method !!< Method switch for Diagonalization of the Schroedinger Eqauation in systemdvr.f90 module
        READ(101,'(1D20.10)') abstol
        READ(101, *)
        READ(101,'(4I20)') timebool !!<Boolean for timestamp.
        READ(101,'(1D20.10)') xrange !!<Range parameter for the position
        CLOSE(101)

        IF(inputNumber.ge.5) THEN

            READ( temp1, '(I10)' ) meBnd
            READ( temp2, '(I10)') meBnd_small

        END IF


    END SUBROUTINE readInputComputation

END MODULE datain
