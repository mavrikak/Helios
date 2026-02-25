module regridpack_c
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use regridpack_module, only: regrid   ! <- use the public generic
  implicit none
contains

! C signature:
!   void rgrd2_eval_c(int mx, int my,
!                     const double* x, const double* y, const double* f_colmajor,
!                     int m, int n,
!                     const double* xq, const double* yq,
!                     int method, /* 0=linear, 1=cubic */
!                     double* fq_out, int* ier);
subroutine rgrd2_eval_c(mx, my, x, y, f, m, n, xq, yq, method, fq, ier) bind(C, name="rgrd2_eval_c")
  integer(c_int),              value :: mx, my, m, n, method
  real(c_double), intent(in)           :: x(mx), y(my)
  real(c_double), intent(in),  target  :: f(mx*my)     ! TARGET so we can remap
  real(c_double), intent(in)           :: xq(m), yq(n)
  real(c_double), intent(out), target  :: fq(m*n)      ! TARGET so we can remap
  integer(c_int), intent(out)          :: ier

  integer :: intpol(2)
  real(c_double), pointer :: p2d(:,:), q2d(:,:)

  ! Remap flat arrays to 2D (column-major) without copies
  p2d(1:mx,1:my) => f
  q2d(1:m,1:n)   => fq

  ! 0 -> linear (1), 1 -> cubic (3) for both dims
  if (method == 0) then
     intpol = [1,1]
  else
     intpol = [3,3]
  end if

  ! Call the public generic; it will dispatch to the 2D wrapper
  call regrid(x, y, p2d, xq, yq, q2d, intpol, ier)
end subroutine rgrd2_eval_c

end module regridpack_c
