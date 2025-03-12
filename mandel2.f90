program mandelbrot
implicit none
integer, parameter :: qp = kind(0.q0)
real(qp), parameter :: D = 0.035_qp
real(qp), parameter :: xmin = -0.743643887037157704752191506114774_qp
real(qp), parameter :: xmax = -0.743643887037159704752191506114774_qp
real(qp), parameter :: ymin =  0.131825904205310970493132056385139_qp
real(qp), parameter :: ymax =  0.131825904205312970493132056385139_qp
!real(qp), parameter :: xmin = -1.9_qp
!real(qp), parameter :: xmax = 0.9_qp
!real(qp), parameter :: ymin = -0.9_qp
!real(qp), parameter :: ymax = 0.9_qp

integer :: ix, iy, k, max_iters
integer :: iters(0:30, 0:70)
complex(qp) :: c, z
logical escaped
real(qp) :: stepx, stepy

stepx = (xmax - xmin) / 70.0_qp
stepy = (ymax - ymin) / 30.0_qp

max_iters = 25000
iters(:,:) = 0

!$omp parallel do reduction(+:iters)
do iy=0,30
   do ix=0,70
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
!$omp end parallel do

do iy=0,30
   do ix=0,70
      write(*,'(i0)',advance="no") (9*iters(iy,ix))/(max_iters)
   end do
   write(*,*)
end do

end program mandelbrot
