program main

  
  implicit none
 
  integer, parameter                    :: n = 11
  integer, parameter                    :: m = 100
  real(kind=8)                          :: x1(n), v1(n)
  real(kind=8)                          :: xq1(m), vq1(m)

  real(kind=8)                          :: dx, a, b 
  integer                               :: i
  a = 1.0
  b = 10.0
  dx = (b-a) / real(m-1, kind=8)

  x1 = (/1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 7.0, 8.0, 9.0, 9.5, 10.0/)
  v1 = (/0.0, 0.0, 0.0, 0.5, 0.4, 1.2, 1.2, 0.1, 0.0, 0.3, 0.6 /)

  xq1(1) = a
  do i=2, m
   xq1(i) = xq1(i-1) + dx
   vq1(i) = 0.0
  enddo

  call makima(x1, v1, n, xq1, vq1, m)

  !! open file 
  open( unit=100, file='data', status='unknown' )
  do i=1, m
    write(100, *) xq1(i), vq1(i)
    !!write(* , *) xq1(i), vq1(i)
  enddo
  close(100)

end program main


subroutine diff(x, n, res)
!!
!! Calculates the difference between adjascent element 
!!

  implicit none 


  integer, intent(in)           :: n               !! number of elements in x
  real(kind=8), intent(in)      :: x(n)            !! array of values
  real(kind=8), intent(out)     :: res(n-1)        !! output array  

  integer                       :: i
  
  do i=1, n-1
    res(i) = x(i+1)-x(i)
  enddo

end subroutine


subroutine makimaSlopes(delta_in , m, s)
!!
!! Derivative values for modified Akima cubic Hermite interpolation
!!
!! Akima's derivative estimate at grid node x(i) requires the four finite
!! differences corresponding to the five grid nodes x(i-2:i+2).
!!
!! For boundary grid nodes x(1:2) and x(n-1:n), append finite differences
!! which would correspond to x(-1:0) and x(n+1:n+2) by using the following
!! uncentered difference formula correspondin to quadratic extrapolation
!! using the quadratic polynomial defined by data at x(1:3)
!! (section 2.3 in Akima's paper):
!!

  implicit none
  
  integer, intent(in)           :: m
  real(kind=8), intent(out)      :: s(m+1)
  real(kind=8), intent(in)      :: delta_in(m)

  integer                       :: n, i
  real(kind=8)                  :: delta_0
  real(kind=8)                  :: delta_m1
  real(kind=8)                  :: delta_n
  real(kind=8)                  :: delta_n1
  real(kind=8)                  :: delta(m+4)
  real(kind=8)                  :: weights(m+3)
  real(kind=8)                  :: weights1(m+1)
  real(kind=8)                  :: weights2(m+1)
  real(kind=8)                  :: weights12(m+1)
  real(kind=8)                  :: delta1(m+1)
  real(kind=8)                  :: delta2(m+1)
  real(kind=8)                  :: eps

  eps = 1.0e-20
  n = m+1
  delta_0  = 2.0*delta_in(1)   - delta_in(2)
  delta_m1 = 2.0*delta_0       - delta_in(1)
  delta_n  = 2.0*delta_in(n-1) - delta_in(n-2)
  delta_n1 = 2.0*delta_n      - delta_in(n-1)

  delta(1) = delta_m1
  delta(2) = delta_0
  delta(3:n+1) = delta_in
  delta(n+2) = delta_n
  delta(n+3) = delta_n1

  !! Akima's derivative estimate formula (equation (1) in the paper):
  !!
  !!       H. Akima, "A New Method of Interpolation and Smooth Curve Fitting
  !!       Based on Local Procedures", JACM, v. 17-4, p.589-602, 1970.
  !!
  !! s(i) = (|d(i+1)-d(i)| * d(i-1) + |d(i-1)-d(i-2)| * d(i))
  !!      / (|d(i+1)-d(i)|          + |d(i-1)-d(i-2)|)
  !!
  !! To eliminate overshoot and undershoot when the data is constant for more
  !! than two consecutive nodes, in MATLAB's 'makima' we modify Akima's
  !! formula by adding an additional averaging term in the weights.
  !! s(i) = ( (|d(i+1)-d(i)|   + |d(i+1)+d(i)|/2  ) * d(i-1) +
  !!          (|d(i-1)-d(i-2)| + |d(i-1)+d(i-2)|/2) * d(i)  )
  !!      / ( (|d(i+1)-d(i)|   + |d(i+1)+d(i)|/2  ) +
  !!          (|d(i-1)-d(i-2)| + |d(i-1)+d(i-2)|/2)
  call diff( delta, m+4, weights)
  weights = abs(weights) + abs((delta(1:m+3)+delta(2:m+4))*0.5) 

  
  weights1 = weights(1:n);   !! |d(i-1)-d(i-2)|
  weights2 = weights(3:n+2); !! |d(i+1)-d(i)|
  delta1 = delta(2:n+1);     !! d(i-1)
  delta2 = delta(3:n+2);     !! d(i)

  weights12 = weights1 + weights2;

  do i=1, n
    !!If the data is constant for more than four consecutive nodes, then the
    !!denominator is zero and the formula produces an unwanted NaN result.
    !!Replace this NaN with 0
    if(abs(weights12(i)) < eps) then
       s(i) = 0.00
    else
       s(i) = (weights2(i)/weights12(i)) * delta1(i) &
            + (weights1(i)/weights12(i)) * delta2(i);
    endif
  enddo

end subroutine


subroutine makima(x, v, n, xq, vq, m)
!!
!!
!!
  implicit none

  integer, intent(in)           :: n
  integer, intent(in)           :: m
  real(kind=8), intent(in)      :: x(n)
  real(kind=8), intent(in)      :: v(n)
  real(kind=8), intent(in)      :: xq(m)
  real(kind=8), intent(out)     :: vq(m)

  integer                       :: i,j
  real(kind=8)                  :: h(n-1)
  real(kind=8)                  :: delta(n-1)
  real(kind=8)                  :: slopes(n)
  real(kind=8)                  :: d(m)
  real(kind=8)                  :: s(m)
  real(kind=8)                  :: t(m)
  

  !! Calculate modified Akima slopes
  call diff(x, n, h)
  call diff(v, n, delta)
  print*, 'h=', h

  do i=1, n-1
    delta(i) = delta(i) / h(i)
  enddo 

  call makimaSlopes(delta, n-1, slopes)

  print*, 'slopes = ', slopes
  !! Evaluate piece wise cubic Hermite interpolants

  do i=1, n-1
    do j=1, m
      if( xq(j) .eq. x(i) ) then
        vq(j) = v(i)
      else if( xq(j) .eq. x(i+1) ) then
        vq(j) = v(i+1)
      else if( x(i) < xq(j) .and. xq(j) < x(i+1) ) then
        call hermite_cubic_value( x(i), v(i), slopes(i), x(i+1), v(i+1), slopes(i+1), &
                                 1, xq(j), vq(j), d(j), s(j), t(j) )
      endif
   enddo
  enddo

end subroutine

subroutine hermite_cubic_value ( x1, f1, d1, x2, f2, d2, n, x, f, d, s, t )

!*****************************************************************************80
!
!! HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.
!
!  Discussion:
!
!    The input arguments can be vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2011
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, F1, D1, the left endpoint, function value
!    and derivative.
!
!    Input, real ( kind = 8 ) X2, F2, D2, the right endpoint, function value
!    and derivative.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the points at which the Hermite cubic
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), D(N), S(N), T(N), the value and first
!    three derivatives of the Hermite cubic at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) df
  real ( kind = 8 ) dx(n)
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) h
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  h =    x2 - x1
  df = ( f2 - f1 ) / h

  c2 = - ( 2.0 * d1 - 3.0 * df + d2 ) / h
  c3 =   (       d1 - 2.0 * df + d2 ) / h / h

  dx(1:n) = x(1:n) - x1
  f(1:n) = f1 + dx(1:n) &
             * ( d1 + dx(1:n) * (           c2 + dx(1:n) *           c3 ) )
  d(1:n) =       d1 + dx(1:n) * ( 2.0D+00 * c2 + dx(1:n) * 3.0D+00 * c3 )
  s(1:n) =                        2.0D+00 * c2 + dx(1:n) * 6.0D+00 * c3
  t(1:n) =                                                 6.0D+00 * c3

  return
end
