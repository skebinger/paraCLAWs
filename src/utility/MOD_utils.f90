! ====================== FORTRAN ======================= !
!                         ____ _        ___        __    !
!  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  !
! | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| !
! | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ !
! | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ !
! |_|                                                    !
! ====================================================== !

module MOD_utils
    use omp_lib

    double precision :: t1,t2

contains
    subroutine tic()
      implicit none
#ifdef _OPENMP
      !$OMP PARALLEL
      !$OMP SINGLE
      t1=omp_get_wtime()
      !$OMP END SINGLE
      !$OMP END PARALLEL
#else
      call cpu_time(t1)
#endif

    end subroutine tic
    
    subroutine toc()
      implicit none
#ifdef _OPENMP
      !$OMP PARALLEL
      !$OMP SINGLE
      t2=omp_get_wtime()
      !$OMP END SINGLE
      !$OMP END PARALLEL
#else
      call cpu_time(t2)
#endif
      print*,"Time Taken -->", real(t2-t1)
    end subroutine toc
end module