subroutine calck0 ( arg, result_, jint )

!*****************************************************************************80
!
!! CALCK0 computes various K0 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order zero, K0(X) and EXP(X)*K0(X), for real
!    arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind=4 ) ARG, the argument.  0 < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real ( kind=4 ) result_, the value of the function,
!    which depends on the input value of JINT:
!    1, result_ = K0(x);
!    2, result_ = exp(x) * K0(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, K0(x);
!    2, exp(x) * K0(x);
!
  implicit none

  integer i
  integer jint
  real arg
  real ( kind=8 ) f(4)
  real ( kind=8 ) g(3)
  real ( kind=8 ) p(6)
  real ( kind=8) pp(10)
  real ( kind=8 ) q(2)
  real ( kind=8 ) qq(10)
  real result_
  real ( kind=8 ) sumf
  real ( kind=8 ) sumg
  real ( kind=8 ) sump
  real ( kind=8) sumq
  real ( kind=8) temp
  real ( kind=8) x
  real ( kind=8) xinf
  real ( kind=8 ) xmax
  real ( kind=8) xsmall
  real ( kind=8 ) xx
!
!  Machine-dependent constants
!
  data xsmall /1.11d-16/
  data xinf /1.79d+308/
  data xmax /705.342d0/
!
!  Coefficients for XSMALL <= ARG <= 1.0
!
  data   p/ 5.8599221412826100000d-04, 1.3166052564989571850d-01, &
            1.1999463724910714109d+01, 4.6850901201934832188d+02, &
            5.9169059852270512312d+03, 2.4708152720399552679d+03/
  data   q/-2.4994418972832303646d+02, 2.1312714303849120380d+04/
  data   f/-1.6414452837299064100d+00,-2.9601657892958843866d+02, &
           -1.7733784684952985886d+04,-4.0320340761145482298d+05/
  data   g/-2.5064972445877992730d+02, 2.9865713163054025489d+04, &
           -1.6128136304458193998d+06/
!
!  Coefficients for  1.0 < ARG
!
  data  pp/ 1.1394980557384778174d+02, 3.6832589957340267940d+03, &
            3.1075408980684392399d+04, 1.0577068948034021957d+05, &
            1.7398867902565686251d+05, 1.5097646353289914539d+05, &
            7.1557062783764037541d+04, 1.8321525870183537725d+04, &
            2.3444738764199315021d+03, 1.1600249425076035558d+02/
  data  qq/ 2.0013443064949242491d+02, 4.4329628889746408858d+03, &
            3.1474655750295278825d+04, 9.7418829762268075784d+04, &
            1.5144644673520157801d+05, 1.2689839587977598727d+05, &
            5.8824616785857027752d+04, 1.4847228371802360957d+04, &
            1.8821890840982713696d+03, 9.2556599177304839811d+01/

  x = arg
!
!  0.0 < ARG <= 1.0.
!
  if ( 0.0D+00 < x ) then

    if ( x <= 1.0D+00 ) then

      temp = log ( x )

      if ( x < xsmall ) then
!
!  Return for small ARG.
!
        result_ = p(6) / q(2) - temp

      else

        xx = x * x

        sump = (((( &
                 p(1) &
          * xx + p(2) ) &
          * xx + p(3) ) &
          * xx + p(4) ) &
          * xx + p(5) ) &
          * xx + p(6)

        sumq = ( xx + q(1) ) * xx + q(2)
        sumf = ( ( &
                 f(1) &
          * xx + f(2) ) &
          * xx + f(3) ) &
          * xx + f(4)

        sumg = ( ( xx + g(1) ) * xx + g(2) ) * xx + g(3)

        result_ = sump / sumq - xx * sumf * temp / sumg - temp

        if ( jint == 2 ) then
          result_ = result_ * dexp(x)
        end if

      end if

    else if ( jint == 1 .and. xmax < x ) then
!
!  Error return for XMAX < ARG.
!
      result_ = 0.0D+00

    else
!
!  1.0 < ARG.
!
      xx = 1.0D+00 / x
      sump = pp(1)
      do i = 2, 10
        sump = sump * xx + pp(i)
      end do

      sumq = xx
      do i = 1, 9
        sumq = ( sumq + qq(i) ) * xx
      end do
      sumq = sumq + qq(10)
      result_ = sump/sumq/dsqrt(x)

      if ( jint == 1 ) then
        result_ = result_*dexp (-x)
      end if

    end if

  else
!
!  Error return for ARG <= 0.0.
!
    result_= xinf

  end if

  return
end
subroutine calck1 ( arg, result_, jint )

!*****************************************************************************80
!
!! CALCK1 computes various K1 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order one, K1(X) and EXP(X)*K1(X), for real arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind=4 ) ARG, the argument.  XLEAST < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real ( kind=4 ) result__, the value of the function,
!    which depends on the input value of JINT:
!    1, result__ = K1(x);
!    2, result__ = exp(x) * K1(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, K1(x);
!    2, exp(x) * K1(x);
!
  implicit none

  real ( kind=4 ) arg
  real ( kind=4 ) f(5)
  real ( kind=4 ) g(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jint
  real ( kind=4 ) p(5)
  real ( kind=4 ) pp(11)
  real ( kind=4 ) q(3)
  real ( kind=4 ) qq(9)
  real result_
  real ( kind=4 ) sumf
  real ( kind=4 ) sumg
  real ( kind=4 ) sump
  real ( kind=4 ) sumq
  real ( kind=4 ) x
  real ( kind=8 ) xinf
  real ( kind=4 ) xmax
  real ( kind=8 ) xleast
  real ( kind=4 ) xsmall
  real ( kind=4 ) xx
!
!  Machine-dependent constants
!
  data xleast /2.23d-308/
  data xsmall /1.11d-16/
  data xinf /1.79d+308/
  data xmax /705.343d+0/
!
!  Coefficients for  XLEAST <=  ARG  <= 1.0
!
  data   p/ 4.8127070456878442310d-1, 9.9991373567429309922d+1, &
            7.1885382604084798576d+3, 1.7733324035147015630d+5, &
            7.1938920065420586101d+5/
  data   q/-2.8143915754538725829d+2, 3.7264298672067697862d+4, &
           -2.2149374878243304548d+6/
  data   f/-2.2795590826955002390d-1,-5.3103913335180275253d+1, &
           -4.5051623763436087023d+3,-1.4758069205414222471d+5, &
           -1.3531161492785421328d+6/
  data   g/-3.0507151578787595807d+2, 4.3117653211351080007d+4, &
           -2.7062322985570842656d+6/
!
!  Coefficients for  1.0 < ARG
!
  data  pp/ 6.4257745859173138767d-2, 7.5584584631176030810d+0, &
            1.3182609918569941308d+2, 8.1094256146537402173d+2, &
            2.3123742209168871550d+3, 3.4540675585544584407d+3, &
            2.8590657697910288226d+3, 1.3319486433183221990d+3, &
            3.4122953486801312910d+2, 4.4137176114230414036d+1, &
            2.2196792496874548962d+0/
  data  qq/ 3.6001069306861518855d+1, 3.3031020088765390854d+2, &
            1.2082692316002348638d+3, 2.1181000487171943810d+3, &
            1.9448440788918006154d+3, 9.6929165726802648634d+2, &
            2.5951223655579051357d+2, 3.4552228452758912848d+1, &
            1.7710478032601086579d+0/

  x = arg
!
!  Error return for ARG < XLEAST.
!
  if ( x < xleast ) then

    result_ = xinf
!
!  XLEAST <= ARG <= 1.0.
!
  else if ( x <= 1.0D+00 ) then

    if ( x < xsmall ) then
!
!  Return for small ARG.
!
      result_ = 1.0D+00 / x

    else

      xx = x * x

      sump = (((( &
               p(1) &
        * xx + p(2) ) &
        * xx + p(3) ) &
        * xx + p(4) ) &
        * xx + p(5) ) &
        * xx + q(3)

      sumq = (( &
          xx + q(1) ) &
        * xx + q(2) ) &
        * xx + q(3)

      sumf = ((( &
               f(1) &
        * xx + f(2) ) &
        * xx + f(3) ) &
        * xx + f(4) ) &
        * xx + f(5)

      sumg = (( &
          xx + g(1) ) &
        * xx + g(2) ) &
        * xx + g(3)

      result_ = ( xx * log ( x ) * sumf / sumg + sump / sumq ) / x

      if ( jint == 2 ) then
        result_ = result_ * exp ( x )
      end if

    end if

  else if ( jint == 1 .and. xmax < x ) then
!
!  Error return for XMAX < ARG.
!
    result_ = 0.0D+00

  else
!
!  1.0 < ARG.
!
    xx = 1.0D+00 / x

    sump = pp(1)
    do i = 2, 11
      sump = sump * xx + pp(i)
    end do

    sumq = xx
    do i = 1, 8
      sumq = ( sumq + qq(i) ) * xx
    end do
    sumq = sumq + qq(9)

    result_ = sump / sumq / sqrt ( x )

    if ( jint == 1 ) then
      result_ = result_ * exp ( -x )
    end if

  end if
  return
end
   
real function besek0 ( x )

!*****************************************************************************80
!
!! BESEK0 evaluates the exponentially scaled Bessel K0(X) function.
!!  Parameters:
!
!    Input, real  X, the argument of the function.
!    0 < X.
!
!    Output,  BESK0, the value of the function.
!
  implicit none
  integer jint
  real x,result_
  jint = 1
  call calck0(x,result_,jint )
  besek0 = result_
  return
end
real function besek1(x)
!! BESEK1 evaluates the exponentially scaled Bessel K1(X) function.
!    This routine computes approximate values for the
!    modified Bessel function of the second kind of order one
!    multiplied by the exponential function, for arguments
!    XLEAST <= ARG <= XMAX.
!  Parameters:
!    Input, real X, the argument of the function.
!    Output, real BESEK1, the value of the function.
  implicit none
  integer jint
  real result_
  real  x
  jint = 1
  call calck1(x,result_, jint)
  besek1 = result_
  return
end
