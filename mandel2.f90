program mandelbrot
  
  implicit none

  integer, parameter :: nx = 70
  integer, parameter :: ny = 30
  integer, parameter :: max_iters = 1000000
  
  integer, parameter :: qp = kind(0.q0)
  real(qp), parameter :: D = 0.035_qp
  !real(qp), parameter :: xmin = -0.743643887037157704752191506114774_qp
  !real(qp), parameter :: xmax = -0.743643887037159704752191506114774_qp
  !real(qp), parameter :: ymin =  0.131825904205310970493132056385139_qp
  !real(qp), parameter :: ymax =  0.131825904205312970493132056385139_qp
  real(qp), parameter :: xmin = -1.9_qp
  real(qp), parameter :: xmax = 0.9_qp
  real(qp), parameter :: ymin = -0.9_qp
  real(qp), parameter :: ymax = 0.9_qp

  integer :: iters(0:ny, 0:nx)
  integer :: ix, iy

  !$omp parallel default( none ) shared( iters )
  call do_mandelbrot(xmin, xmax, ymin, ymax, iters, max_iters)
  !$omp end parallel


  do iy=0,ny
     do ix=0,nx
        write(*,'(i0)',advance="no") (9*iters(iy,ix))/(max_iters)
     end do
     write(*,*)
  end do


contains

  subroutine do_mandelbrot(xmin, xmax, ymin, ymax, iters, max_iters)

    !$ use omp_lib, only : omp_get_thread_num
    !$ use omp_lib, only : omp_get_num_threads

    implicit none

    real(qp), intent(in) :: xmin, xmax, ymin, ymax
    integer, intent(in) :: max_iters
    integer, intent(out) :: iters(0:, 0:)

    
    complex(qp) :: c, z
    real(qp) :: stepx, stepy
    integer :: nx, ny, ix, iy, k
    logical escaped

    ny = ubound(iters, dim=1)
    nx = ubound(iters, dim=2)
    
    stepy = ( ymax - ymin ) / ny
    stepx = ( xmax - xmin ) / nx

!    !$ if( omp_get_thread_num() == 0 ) then
!    !$    write( *, * ) 'On ', omp_get_num_threads(), ' threads'
!    !$ end if
    
    !$omp do collapse(2) schedule(dynamic)
    do ix=0,nx
       do iy=0,ny
          c = cmplx(xmin+stepx*ix, ymin+stepy*iy, qp)
          z = 0
          escaped = .false.
          do k=1,max_iters
             z = z**2 + c
             escaped = abs(z) > 2
             if (escaped) exit
          end do
          iters(iy, ix) = k      
       end do
    end do
    !$omp end do

  end subroutine do_mandelbrot
  
end program mandelbrot
