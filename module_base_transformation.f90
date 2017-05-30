module module_base_transformation

    USE class_system
    USE class_density

    implicit none


contains

    subroutine transform_to_newbase_eq(density_old, density_new, old_sys, new_sys)

        CLASS(densityvector), intent(inout)::density_new
        CLASS(densityvector), intent(inout)::density_old
        CLASS(system), intent(in):: old_sys
        CLASS(system), intent(in):: new_sys

        COMPLEX(8)                         :: summation
        INTEGER                            :: i, j

        ASSOCIATE(  all_m => density_new%NumDim00,              &
            all_v => density_new%NumDim11,              &
            diag_m => density_new%medim00,          &
            diag_v => density_new%medim11,          &

            !Indize
            m1 => density_new%tulpVector(:,1,0),        &
            m2 => density_new%tulpVector(:,2,0),        &
            m1new => density_new%tulpVector(:,1,0),        &
            m2new => density_new%tulpVector(:,2,0),        &
            v1 => density_new%tulpVector(:,1,1),        &
            v2 => density_new%tulpVector(:,2,1),        &
            v1new => density_new%tulpVector(:,1,1),        &
            v2new => density_new%tulpVector(:,2,1))

            !> Check for norm preservation
            CALL density_old%check_trace()

            WRITE(*,*) "Performing base transformation"
            !--- Matrixelements in 00
            DO i = 1, all_m
                summation = (0.0d0,0.0d0)

                !00
                DO j = 1, all_m
                    summation = summation + dot_product(new_sys%hsOV(:, m1(i), 0), old_sys%hsOV(:, m1new(j), 0) )*density_old%rho(j)*dot_product(new_sys%hsOV(:, m2(i), 0), old_sys%hsOV(:, m2new(j), 0))
                END DO

                density_new%rho(i) = summation

            END DO

            !--- Matrixelements in 11
            DO i = 1, all_v
                summation = (0.0d0,0.0d0)

                !11
                DO j = 1, all_v
                    summation = summation + dot_product(new_sys%hsOV(:, v1(i), 1), old_sys%hsOV(:, v1new(j), 1) )*density_old%rho(j+all_m)*dot_product(new_sys%hsOV(:, v2(i), 1), old_sys%hsOV(:, v2new(j), 1))
                END DO

                density_new%rho(i + all_m) = summation
            END DO

        END ASSOCIATE


        !> Check for norm preservation
        CALL density_new%check_trace()

    end subroutine transform_to_newbase_eq



end module module_base_transformation
