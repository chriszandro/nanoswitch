module UnitConv

    IMPLICIT NONE

    !Time
    real(8), parameter :: au2second = 2.4188843d-17
    real(8), parameter :: second2au = 1.0d0/au2second
    real(8), parameter :: au2femto = au2second*1.0d15

    ! Energy
    real(8), parameter :: au2ev = 27.2113834d0         ! Hartree energy (au)
    real(8), parameter :: ev2au = 1.0d0/au2ev          ! and eV
    ! Length
    real(8), parameter :: au2ang = 0.5291772083d0      ! Bohr radius (au)
    real(8), parameter :: ang2au = 1.0d0/au2ang        ! and angstrom
    ! Mass
    real(8), parameter :: au2kg = 9.10938188d-31       ! Electron rest mass (au)
    real(8), parameter :: kg2au = 1.0d0/au2kg          ! and kg (si)
    real(8), parameter :: amu2kg = 1.66053873d-27      ! Unified atomic mass unit
    real(8), parameter :: kg2amu = 1.0d0/amu2kg        ! and kg (si)
    real(8), parameter :: au2amu = au2kg*kg2amu        ! Electron rest mass (au)
    real(8), parameter :: amu2au = 1.0d0/au2amu        ! and unified atomic mass unit
    ! Current
    real(8), parameter :: au2mua = 0.00662362108d6     ! Atomic unit of current (au)
                                                       ! and microAmpers
    ! Mass Proton
    REAL(8), parameter :: protonMass = 1836.1524d0     ! m_p = 1836.1524*me
    REAL(8), parameter :: Minverse =  1/protonMass                      !0.0005446171d0

    !Field
    real(8), parameter :: au2Vm = 5.14220652d11
    real(8),  parameter :: Vm2au=1.0d0/au2Vm

    ! Fundamental constants
    real(8), parameter    :: Pi = 3.1415926535897932385d0
    real(8), parameter    :: gaussianIntegral = DSQRT(Pi)
    real(8), parameter    :: Kb = 8.61734315d-5*ev2au
    complex(8), parameter :: Im = dcmplx(0.0d0,1.0d0)

    real(8), parameter    :: au2kelvin = 1.0/Kb 
    
    ! Angle Conversion
    REAL, PARAMETER :: Degree180 = 180.0
    REAL, PARAMETER :: R_to_D    = Degree180/PI
    REAL, PARAMETER :: D_to_R    = PI/Degree180

end module UnitConv
