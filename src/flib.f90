FUNCTION gaussian(x)
  !
  use constants,  only : dp
  !
  implicit none
  !
  real(dp) gaussian, x
  !
  gaussian=exp(-x*x)/1.77245385090552
  !
END FUNCTION
