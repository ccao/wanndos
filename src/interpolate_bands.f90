include 'lapack.f90'

SUBROUTINE interpolate_bands(kv)
  !
  use lapack95,  only : heev
  use constants, only : dp, twopi, cmplx_0, cmplx_i
  use wanndata,  only : rvec, ham, weight, nrpt, norb
  use dosdata,   only : eig, work
  !
  implicit none
  !
  real(dp) kv(1:3)
  real(dp) rdotk
  complex(dp) fact
  !
  integer ir, info
  !
  work(:,:)=cmplx_0
  do ir=1, nrpt
    rdotk=SUM(kv(:)*rvec(:, ir))
    fact=exp(-cmplx_i*twopi*rdotk)/weight(ir)
    work(:,:)=work(:,:)+ham(:,:,ir)*fact
  enddo ! ir
  !
  call heev(work, eig, 'V', 'U', info)
  !
  ! work(ii, io): ii\th element eigenvector of eig(io)
  !
END SUBROUTINE
