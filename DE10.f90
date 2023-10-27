
subroutine DE10(DS, A, B, N, Y) !решение уравнения сверху методом РК 4 порядка
interface
  subroutine DS(X, Y, DY)
  real, intent(in):: X, Y
  real:: DY
  end subroutine DS
end interface
real, intent(in):: A, B
integer, intent(in):: N
real, intent(inout):: Y
real, parameter:: C(3)=(/0.5, 0.5, 1.0/),   &
  P(3)=(/0.166666666666667, 0.33333333333333, 0.33333333333333/)
real:: CH(3)
real:: Y0, YH, K
real:: X, XH, H
integer:: i, j, l, M
!begin
  if (A==B) return
  H=(B-A)/float(N)
  do j=1, 3
    CH(j)=H*C(j)
  end do
  Y0=Y
  do l=0, N-1
    X=A+H*float(l)
    call DS(X, Y, K)
    do j=1, 3
      XH=X+CH(j)
      K =H*K 
      YH =Y0 +C(j)*K 
      Y =Y +P(j)*K 
      call DS(XH, YH, K)
    end do
    Y =Y +P(1)*H*K 
    Y0=Y
  end do
  return
end subroutine DE10
