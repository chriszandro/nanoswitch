################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../UnitConv.f90 \
../VTR_mod.f90 \
../class_bosonic_bath.f90 \
../class_coupling.f90 \
../class_density.f90 \
../class_grid.f90 \
../class_inputdata.f90 \
../class_inputdata2.f90 \
../class_junction.f90 \
../class_ltensor.f90 \
../class_observables.f90 \
../class_rhox.f90 \
../class_system.f90 \
../element_selector.f90 \
../module_base_transformation.f90 \
../module_evolution.f90 \
../module_run_stationary.f90 \
../module_snapshots.f90 \
../program.f90 \
../timestamper.f90 

OBJS += \
./UnitConv.o \
./VTR_mod.o \
./class_bosonic_bath.o \
./class_coupling.o \
./class_density.o \
./class_grid.o \
./class_inputdata.o \
./class_inputdata2.o \
./class_junction.o \
./class_ltensor.o \
./class_observables.o \
./class_rhox.o \
./class_system.o \
./element_selector.o \
./module_base_transformation.o \
./module_evolution.o \
./module_run_stationary.o \
./module_snapshots.o \
./program.o \
./timestamper.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Intel(R) Intel(R) 64 Fortran Compiler'
	ifort -g -O3 -c -mkl -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

UnitConv.o: ../UnitConv.f90

VTR_mod.o: ../VTR_mod.f90

class_bosonic_bath.o: ../class_bosonic_bath.f90 UnitConv.o class_density.o class_inputdata.o class_system.o

class_coupling.o: ../class_coupling.f90 class_density.o class_grid.o class_inputdata.o class_system.o

class_density.o: ../class_density.f90 class_inputdata.o

class_grid.o: ../class_grid.f90 class_inputdata.o

class_inputdata.o: ../class_inputdata.f90 UnitConv.o timestamper.o

class_inputdata2.o: ../class_inputdata2.f90 UnitConv.o

class_junction.o: ../class_junction.f90 UnitConv.o class_density.o class_inputdata.o class_system.o

class_ltensor.o: ../class_ltensor.f90 class_bosonic_bath.o class_coupling.o class_density.o class_inputdata.o class_junction.o class_system.o

class_observables.o: ../class_observables.f90 class_coupling.o class_density.o class_grid.o class_inputdata.o class_junction.o class_system.o

class_rhox.o: ../class_rhox.f90 UnitConv.o class_density.o class_grid.o class_inputdata.o class_system.o

class_system.o: ../class_system.f90 UnitConv.o class_density.o class_grid.o class_inputdata.o

element_selector.o: ../element_selector.f90

module_base_transformation.o: ../module_base_transformation.f90 class_density.o class_system.o

module_evolution.o: ../module_evolution.f90 class_bosonic_bath.o class_coupling.o class_density.o class_grid.o class_inputdata.o class_junction.o class_ltensor.o class_observables.o class_rhox.o class_system.o module_base_transformation.o

module_run_stationary.o: ../module_run_stationary.f90 class_bosonic_bath.o class_coupling.o class_density.o class_grid.o class_inputdata.o class_junction.o class_ltensor.o class_observables.o class_rhox.o class_system.o

module_snapshots.o: ../module_snapshots.f90 class_bosonic_bath.o class_coupling.o class_density.o class_grid.o class_inputdata.o class_junction.o class_ltensor.o class_observables.o class_rhox.o class_system.o

program.o: ../program.f90 class_inputdata.o module_evolution.o module_run_stationary.o module_snapshots.o

timestamper.o: ../timestamper.f90


