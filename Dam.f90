!���������� ����������� � ������������ ��������
!��� y � �����, ��� � ������, ��� z ����
! AB ������� ��� ���� �������� �����, ������� ���������� ������
! ������: 1 - ������ ������, 2 - ���� �����, ��������� ������
! 3 - ���� � ���� ������, ����� ������ � ����������
! 4- ���� � ���� �����,���� ���������
! 5 - ���� ������, ���� ���������, ����� ����. �� ������� � ����
! �������� ������ ���� ������ - ������ ���� � ������ - ������� ����
! ����������� ����� �� �������
! ��������� ������� �������� �� ������� ������ KL.
! ��������������� ����� �� ������� � ��������
! ����� ����� ������� N1  �  N1+1

program dambawater5
implicit none

integer kol,jtab,Nk,Ny,it,N,N1,N2,N3,N4,N5,N5_1,N5_2,N6,Model !��� ���������� �� ��� �����
integer  i,j,k,kk,nks,Nnodes
real AB,ABdim,MN,deltaH,width,widthDim,beta,StartPoint,Apoint,Ampl
real left,right,leftH,rightH,leftHeight, stime,etime,rightHeight,xC,zC,xE,zE,xK,zK,xL,zL
real eps,tmp,tmp1,tmp2,yy,hx,a_ij,a_zar,f_1,f_2,k1_
real h,h0,h1,h2,h3,h4,kymax,klimit,ymax,h_kln,h_yln,delta,delta1
real DC1(1:3),DC2(1:3),xzar,yzar,zzar,ro1,ro2,ro3,ro4,kappa
real,allocatable::RelX(:),RelZ(:),ReldZ(:),AM(:),AN(:),RoKA(:),xMN(:),zMN(:)
real, allocatable,dimension(:):: lambda,y,rhand,dlambda
real,dimension(:),allocatable :: x,z,ksi,ksi1 ! ���� �����
real,dimension(:),allocatable :: nu,nu1,U_,B,B1,k_y
real, allocatable,dimension(:,:):: AA,BB,U,UT
real, allocatable,dimension(:,:):: A,A1,nu_transf,U_transf,nu2D,nuxy
!*************************************************************
real,parameter:: PI=3.14159265358979
real, external:: besek0,besek1,integr1D
namelist /param/ widthDim,StartPoint,beta,jtab ! ������ ��� ���������� ����� �������
!*************************************************************
!*****             PROGRAM                       ***********
!*************************************************************
open(4,file='Inputdata1.txt') ! ��������� ������� - ������ ���. �����, ����, ����� �����
read(4,NML=param)
close(4)
write(*,NML=param)
print *, 'Enter AB (meters)'
read *,ABDim ! ��������� �����
print *, 'Enter position of A point from relief start point (meters)'
read *,Apoint ! �� ����� ���������� �� x �� ����� ������� ������� ������ �������� � ������
AB=2. ! ������������ �����, ������� ��������� ����� ������ ����� �������� ��� �����
kol=45
MN=AB/kol
Ampl=tan(beta*pi/180)/pi*widthDim ! ��������� ��������� ������ ��� ���������� �������
print *, 'amplitude(m)=',  Ampl
width=widthDim/ABDim*2 ! ��������������� ������ �������
Apoint=Apoint/ABDim*2.
! ������ ������� ����� � ������
open(10, file='Surface1.dat',form='unformatted')
read(10) jtab !total points
allocate(RelX(1:jtab),RelZ(1:jtab))
do i=1,jtab
    read(10) RelX(i),RelZ(i)  ! tabulated x and z coordinates
enddo
close(10)

open(10, file='DifSurface1.dat',form='unformatted')
read(10) j !total points
allocate(ReldZ(1:jtab))
do i=2,jtab-1
    read(10) RelX(i),ReldZ(i)  ! tabulated x and z coordinates
enddo
ReldZ(1)=ReldZ(2)
ReldZ(jtab)=ReldZ(jtab-1)
close(10)

! ������� � ������������ ������ � ������� ������ ��������� �� StartPoint � ����� A:
! x1=x-StartPoint , x2=x1-Apoint=x-Spoint-Apoint
do i=1,jtab
    RelX(i)=(RelX(i)-Startpoint)/ABDim*2-Apoint
	RelZ(i)=RelZ(i)/ABDim*2!  dimesionless x and z coordinates
enddo
!********************************************************
!                           R B F
!********************************************************

call cpu_time(stime)

eps=1/MN!for 100 parazit occilations for fi1, smaller values does not work for fi1

allocate(lambda(1:jtab),dlambda(1:jtab))
allocate(AA(1:jtab,1:jtab),U(1:jtab,1:jtab))
allocate(y(1:jtab),rhand(1:jtab))
allocate(UT(1:jtab,1:jtab),BB(1:jtab,1:jtab))

do i=1,jtab
do j=1,jtab
    AA(i,j)=0.
    do k=1,jtab
    AA(i,j)=AA(i,j)+fi1(abs(RelX(i)-RelX(k)))*fi1(abs(RelX(j)-RelX(k)))
    enddo
enddo
enddo
 !Reshenie SLAU metodom kvadratnih kornei (Ax=UTUx=b)
U(1,1)=sqrt(AA(1,1))
do j=1,jtab
   U(1,j)=AA(1,j)/U(1,1)
enddo

do i=2,jtab
   tmp1=0.
   do k=1,i-1
      tmp2=0.
      do kk=1,k-1
         tmp2=tmp2+U(kk,k)*U(kk,i)
      enddo
      U(k,i)=(AA(k,i)-tmp2)/U(k,k)
      tmp1=tmp1+U(k,i)**2
   enddo
   U(i,i)=sqrt(AA(i,i)-tmp1)
enddo
! (((((((((((((((  ��� �������� ������� )))))))))))))))))
do j=1,jtab
   do i=1,jtab
      if (j<i) then
        U(i,j)=0.
      endif
   enddo
enddo
do k=1,jtab
   do j=1,jtab
      UT(j,k)=U(k,j)
   enddo
enddo

do i=1,jtab
    do j=1,jtab
       tmp2=0.
       do k=1,jtab
           tmp2=tmp2+UT(i,k)*U(k,j)
       enddo
       BB(i,j)=tmp2
    enddo
enddo
tmp=0.
do i=1,jtab
   do j=1,jtab
       tmp1=abs(AA(i,j)-BB(i,j))
       if (tmp1>tmp) tmp=tmp1
   enddo
enddo
print *, 'error in matrice inversion', tmp
! ����� �������� ))))))))))))))))))))))
! ������ ����� � ��������� Ay=rhand
do i=1,jtab
    rhand(i)=0.
    do k=1,jtab
        rhand(i)=rhand(i)+RelZ(k)*fi1(abs(RelX(i)-RelX(k)))
    enddo
enddo

do i=1,jtab
   tmp1=0.
   do k=1,i-1
      tmp1=tmp1+U(k,i)*y(k)
   enddo
   y(i)=(rhand(i)-tmp1)/U(i,i)
enddo
do i=jtab,1,-1
   tmp1=0.d+0
   do k=i+1,jtab
      tmp1=tmp1+lambda(k)*U(i,k)
   enddo
   lambda(i)=(y(i)-tmp1)/U(i,i)
enddo
! ������� ������������ ��� ������������� �����������
! ������ ����� � ��������� Ay=rhand
do i=1,jtab
    rhand(i)=0.
    do k=1,jtab
        rhand(i)=rhand(i)+ReldZ(k)*fi1(abs(RelX(i)-RelX(k)))
    enddo
enddo

do i=1,jtab
   tmp1=0.
   do k=1,i-1
      tmp1=tmp1+U(k,i)*y(k)
   enddo
   y(i)=(rhand(i)-tmp1)/U(i,i)
enddo
do i=jtab,1,-1
   tmp1=0.d+0
   do k=i+1,jtab
      tmp1=tmp1+dlambda(k)*U(i,k)
   enddo
   dlambda(i)=(y(i)-tmp1)/U(i,i)
enddo

! free memory space
deallocate(y,rhand); deallocate(AA,U);deallocate(BB,UT)
! check the relief form
hx=(RelX(jtab)-RelX(1))/500
open(11,file='relieftest.txt')
do k=0,500
    tmp=RelX(1)+hx*k
    tmp1=relief(tmp)
    write(11,10) tmp,tmp1
enddo
close(11)

open(12,file='normalComp.txt')
do k=0,500
    tmp=RelX(1)+hx*k
    tmp1=F1(tmp)
    tmp2=F2(tmp)
    write(12,10) tmp,tmp1,tmp2
enddo
close(12)
!pause
!stop ! ��� ������������ RBF

ro1=10. ! ���� �������
ro2=500. ! water
!ro3=ro1 ! � ������ ���� ������� ����. ��������� ��������� � ����� �������

print *, ' Model: 1 - relief, 2 - water on the left, 3 - two water surfaces,  '
print *, '4- two water surfaces, damb on the base , 5- water on the rigth and base'
print *, '6- one water surfaces, damb on the base , 7- two water spaces, leaking and base'
read *, Model

if (Model>=4 ) then ! �������� ���������
ro3=500. ! ������������� ���������
!ro3=ro1 ! ��� ������������
endif
if (Model>5) then ! �������� ��������
ro4=30. !
!ro4=ro1 ! ��� ������������
endif
print *, 'Height of the A =',-relief(0.)*ABDim*0.5
print *, 'ro_damb, ro_water, ro_base,Ro_leaking', ro1,ro2,ro3,ro4
! ******************************************************************
! **                 � � � � �                                    **
! ******************************************************************
print *, 'Enter nks - number of segments between MN; 1<nks<=3'
read *, nks
h=MN/nks ! ��� �����
N=nks*kol ! ����� ������� �� ��
N1=N! �����  ����� �� �
N2=N+N/2 ! �����  ������ �� �
N3=N ! ����� �� ���� �����
if (Model.eq.5 .or. Model .eq.1 .or. Model.eq.6) N3=0
N4=N ! ����� �� ���� ������
if (Model.eq.2 .or. Model.eq.1 ) N4=0
N5=widthDim/ABDim*N ! ����� �� ��������� ������� ����� ����� ��������������� ������
if (ro3.eq.ro1) N5=0 ! ��������� ��������� ��� ������ ��������������
xC=-100.
if (Model>1 .and. Model.ne.5 .and. Model.ne.6 ) then ! ���� ������ ����
    ! ��������� ���� ������ ���� ���� ���������
    print *, ' Heght of the water at left side(m) <' , -relief(0.)*(ABdim*0.5)
    read *, leftHeight! ������ ���� �����
    if (leftHeight>-relief(0.)*(ABdim*0.5)) then
        print *, 'leftHeight>A position',-relief(0.)*(ABdim*0.5)
        !pause
       ! stop
    endif
    print *, 'leftHeight = ', leftHeight
    leftH=leftHeight/(ABdim*0.5)
    ! ������� �� ����� ���������� �� ������ ����� ��� ������
    tmp=-Apoint
    tmp1=0.
    call DE10(DS1,tmp1,-leftH,30000,tmp)
    xC=tmp
    zC=relief(xC)
    print *, 'x water left', xC, 'exact height for that x', -zC*(ABdim*0.5)
! ��������� ���������� �� � �� ����
    left=integr1D(curve_element,xC,0.);
endif
xE=100.
if (Model>=3 ) then ! ���� ������� ���� ������
    if (Apoint>width/2) then !  ���� �������� ������ �������
        print *, ' Heght of the water at the right side(m)>',leftHeight,'and<', -relief(0.)*ABDim*0.5
    else
        print *, ' Heght of the water at the right side(m)>',leftHeight,'and<', Ampl
    endif
    read *, rightHeight! ������ ���� �����
    if (rightHeight>Ampl) then
        print *, 'rightHeight > total height',Ampl
       ! pause
       ! stop
    endif
    print *, ' rightHeight = ', rightHeight
    rightH=rightHeight/(ABdim*0.5)
    ! ������� �� ����� ���������� �� ������ ����� ��� ������
    tmp=0.
    tmp1=relief(0.)
    call DE10(DS1,tmp1,-rightH,30000,tmp)
    print *, 'x water right', tmp, 'exact height for that x', -relief(tmp)*(ABdim*0.5)
    if (tmp<0.) then
        tmp=width-Apoint*2-tmp
        print *, 'x', tmp, 'height', -relief(tmp)*(ABdim*0.5)
    endif
    xE=tmp
    zE=relief(xE)
    right=integr1D(curve_element,0.,xE);
endif
if (Model>=6) then ! ��������
    print *, 'Height of water leaking (m):'
    read *, zL
    zL=zL/(ABDim*0.5)
    print *, 'Enter left distance to the water leaking (m)'
    read *, xK
    xK=xK/(ABDim*0.5)
    xK=xK-Apoint
    zK=0.
! ������ ��������� ������ �������� KL
    tmp=width-Apoint ! �������� ����� �����
    tmp1=0. ! z
    call DE10(DS1,tmp1,-zL,5000,tmp)
    if (tmp<0.) then ! ������������ ����� � ������ ������� �����
        tmp=width-Apoint*2-tmp
    endif
    print *, 'xK,zK', xK*(ABdim*0.5),zK*(ABdim*0.5),'xL(m),zL(m)',tmp*(ABdim*0.5), relief(tmp)*(ABdim*0.5)
    xL=tmp
    zL=relief(tmp)
    N6=dist_a(xL,zL,xK,zK)/AB*N ! ��������������� ����� KL
    ! ���������� ���������� ����� �� ���� ������ ���������
    tmp=dist_a(-Apoint,0.,xK,zK)
    N5_1=tmp/h/2 !  ����� ����� �� ��������-���������, ����� � ��� ���� ����
    h2=tmp/N5_1 ! ��� ����� �� �������� �� ���������
     ! ����� �� � �� ������ ����� ���������
    tmp=dist_a(xK,zK,width-Apoint,0.)
    N5_2=tmp/h/2 ! ����� ����� �� ��������� � ��� ���� ����
    h3=tmp/N5_2

endif

if (Model.eq.1)  Nnodes=N1+N2  ! ������ ������ � ������
if (Model.eq.2)  Nnodes=N1+N2+N3  ! ������ � ������ ����
if (Model.eq.3)  Nnodes=N1+N2+N3+N4 ! ������, ������� � ������ �����
if (Model.eq.4)  Nnodes=N1+N2+N3+N4+N5 !! ������, ������� � ������ �����, + ���������
if (Model.eq.5)  Nnodes=N1+N2+N4+N5 ! ������, ������� ����, ���������
if (Model.eq.6)  Nnodes=N1+N2+N4+N5_1+N5_2+N6 ! ������, ������� ����, ���������, ��������
if (Model.eq.7)  Nnodes=N1+N2+N3+N4+N5_1+N5_2+N6 ! ������, ������� � ������ ����, ���������, ��������

if (Model<1 .or. Model>7 ) then
    print *, 'there is no such a model'
   ! pause
   ! stop
endif
!********  ������������ �������� ���������� DirectCurrent (DC)  ********
 DC1(1)=0.!��������� ��������� ��������� A � ������ ���������
 DC1(2)=0. !��������
 DC1(3)=relief(DC1(1)) !������ (�������)

!******* ������������� ���������      ******************************
allocate(AM(1:kol),AN(1:kol),xMN(1:kol),zMN(1:kol))
!**������ x,z ��������� �������� ����������
! ���������� ������������� ����������
tmp=0. ! ������ ��������� � �
do i=1,kol
      call DE10(DS, MN*(i-1),MN*i,500*nks,tmp) !
      xMN(i)=tmp
      zMN(i)=relief(tmp)
      AM(i)=MN*i
      AN(i)=AM(i)+MN
enddo
!��������� ������� ��������� ��������� � ����� ������������� �����
   DC2(1)=xMN(kol)  ! ��������� ��������� ��������� B
   DC2(2)=0.
   DC2(3)=zMN(kol)
allocate(x(1:Nnodes),z(1:Nnodes))
tmp=0.  ! ���������� ����� ����� �� �, N1 ��������
call DE10(DS, 0.,-h*0.5,1000,tmp)   !������ x
x(N1)=tmp !�������� ����� A - ������ ���������
z(N1)=relief(tmp) ! ������
do j = 1,N1-1
    call DE10(DS, -h*(j-0.5),-h*(j+0.5),500,tmp) !������ x
    x(N1-j)=tmp !��������
    z(N1-j)=relief(tmp) ! ������
enddo
tmp=0.  ! xA
call DE10(DS,0.,h*0.5,500,tmp)   !������ x
x(N1+1)=tmp !
z(N1+1)=relief(tmp) ! ������
do j = 1,N2-1
    call DE10(DS, h*(j-0.5),h*(j+0.5),500,tmp) !������ x
    x(N1+1+j)=tmp !��������
    z(N1+1+j)=relief(tmp) ! ������
enddo

if (Model>=2 .and. Model<=4 .or. Model.eq.7) then !  ����� ����� �� ����������� ���� �.�.
    do i=1,N3 !
        x(N1+N2+i)= xC-h*(i-0.5) ! ������������ �� ����� �, �� ������� �� �
        z(N1+N2+i)= zC     ! ������ ����� ���� �
    enddo
endif
if (Model>=3 ) then ! ����� ����� �� ����������� ���� �������� �����
    do i=1,N4 !
        x(N1+N2+N3+i)= xE+h*(i-0.5) ! ������������ �� ����� E, �� ������� �� E
        z(N1+N2+N3+i)= zE     ! ������ ����� ���� E
    enddo
endif
! ��������� ���� �� ������� ���������-����� ��� ��������
if (Model >=4 .and. Model<=5 ) then ! ��������� � ������ ��������������
do i=1,N5
    x(N1+N2+N3+N4+i)=-Apoint+(i-0.5)*width/N5
    z(N1+N2+N3+N4+i)=0.
enddo
endif

! ��������� ���� �� ������� �������� - �����
if (Model >5) then !
    tmp=dist_a(xL,zL,xK,zK)
    h4=tmp/N6
    tmp1=atan((zK-zL)/(xL-xK))
    do i=1,N6
        x(N1+N2+N3+N4+i)=xL-(i-0.5)*h4*cos(tmp1)
        z(N1+N2+N3+N4+i)=zL+(i-0.5)*h4*sin(tmp1)
    enddo
! ����� �� ����� ����� ��������� �� K
    do i=1,N5_1
        x(N1+N2+N3+N4+N6+i)=-Apoint+(i-0.5)*h2
        z(N1+N2+N3+N4+N6+i)=0.
    enddo
    do i=1,N5_2
        x(N1+N2+N3+N4+N6+N5_1+i)=xK+(i-0.5)*h3
        z(N1+N2+N3+N4+N6+N5_1+i)=0.
    enddo
endif

! ��������� ��������� ������ - ������ ����� �
!����� ��������� ��� ���������� � �
open(11, file='Electrodes.csv')
write(11, '(*(G0.7,:,","))') 'ksi,x,z'
write(11,'(*(G0.7,:,","))') 0.,0.,relief(0.)*(ABDim*0.5)
do i=1,kol
    write(11,'(*(G0.7,:,","))') tmp, xMN(i)*(ABDim*0.5),zMN(i)*(ABDim*0.5) !���������� ������������� ����������
    tmp=tmp+h
enddo
close(11)
open(10, file='XZsurface.csv')
write(10, '(*(G0.7,:,","))') 'ksi,x,z'
tmp=0.
do i=1,N1+N2
    write(10, '(*(G0.7,:,","))') tmp,x(i)*(ABDim*0.5),z(i)*(ABDim*0.5) !���������� ������ �������
    tmp=tmp+h
enddo
close(10)
if (Model>=2 .and. Model<=4 .or. Model.eq.7) then
open(11, file='XZwater1.csv')
write(11, '(*(G0.7,:,","))') 'x,z'
do i=N1+N2+1,N1+N2+N3  !
    write(11, '(*(G0.7,:,","))') x(i)*(ABDim*0.5),z(i)*(ABDim*0.5)! ���������� �� ���� �����
enddo
close(11)
endif
if (Model>=2) then
open(11, file='XZwater2.csv')
write(11,'(*(G0.7,:,","))') 'x,z'
do i=N1+N2+N3+1,N1+N2+N3+N4
    write(11,'(*(G0.7,:,","))') x(i)*(ABDim*0.5),z(i)*(ABDim*0.5)! ����������� ���� ������
enddo
close(11)
endif
if (Model.eq.4 .or. Model .eq.5) then
open(11, file='XZbase.txt')
write(11, *) '  x              z  '
do i=N1+N2+N3+N4+1,N1+N2+N3+N4+N5
    write(11,*) x(i)*(ABDim*0.5),z(i)*(ABDim*0.5)
enddo
close(11)
endif
if (Model>5) then
open(11, file='XZleaking.txt')
write(11, *) '  x              z  '
do i=N1+N2+N3+N4+1,N1+N2+N3+N4+N6
    write(11,*) x(i)*(ABDim*0.5),z(i)*(ABDim*0.5)
enddo
close(11)
open(10, file='XZbase1.txt')
write(10, *) '  x              z  '
do i=N1+N2+N3+N4+N6+1,N1+N2+N3+N4+N6+N5_1
    write(10,*) x(i)*(ABDim*0.5),z(i)*(ABDim*0.5)
enddo
close(11)
open(11, file='XZbase2.txt')
write(11, *) '  x              z  '
do i=N1+N2+N3+N4+N6+N5_1+1,N1+N2+N3+N4+N6+N5_1+N5_2
    write(11,*) x(i)*(ABDim*0.5),z(i)*(ABDim*0.5)
enddo
close(11)
endif

open(11, file='Parameters.txt')
write (11,*) 'ABDim=', ABDim
write (11,*) ' number of measuring electrodes=',kol
write (11,*) ' MN=',MN
write (11,*) 'max inclination angle in grad,beta=', beta
write (11,*) 'damba width=', widthDim
write (11,*) 'Apoint', Apoint*ABdim/2.
write (11,*) 'leftHeight, rightHeight',LeftHeight,RightHeight
write(11,*)  'N=',N,'N1=',N1, 'N2=',N2,'Model ',Model
write(11,*) 'total number of points Nnodes=', Nnodes,'N3,N4,N5,N6', N3,N4,N5,N6
if (Model>1) write(11,*) 'coordinates (x,z) of lower water points ', xC*ABdim*0.5,zC*ABdim*0.5
if (Model>2 ) write(11,*) 'coordinates (x,z) of upper water point', xE*ABdim*0.5,zE*ABdim*0.5
write(11,*) 'ro1,ro2,ro3,ro4 - damb, water,base,leaking:', ro1,ro2,ro3,ro4
write(11,*) 'positions of A ', 0.,relief(0.),' and B ', xMN(kol),zMN(kol)
if (Model>5) then
    write(11,*) 'positions of K ', xK*ABDim*0.5,zK*ABDim*0.5,' and L ',xL*ABDim*0.5,zL*ABDim*0.5
    write(11,*) 'N5_1, N5_2', N5_1,N5_2
endif

close(11)
allocate(B(1:Nnodes))
allocate(A(1:Nnodes,1:Nnodes))
allocate(B1(1:Nnodes),A1(1:Nnodes,1:Nnodes))! ��� �� � ��������� �����
allocate(nu(1:Nnodes),nu1(1:Nnodes))

xzar=DC1(1); yzar=DC1(2); zzar=DC1(3)
!kymax=20. ! ������� ������� ���� �������������� � ��������� ����������
kymax=8.
klimit=1./MN/4. !
print *,'klimit',klimit
Nk=64
Ny=64
ymax=4.
h_kln=log(1.+kymax)/Nk ! ��� �� k_�
h_yln=(log(1.+ymax))/Ny
allocate(k_y(0:Nk))
allocate(nu2D(1:Nnodes,0:Nk))

! ******************************************************
!
!              SOLVING OF INTEGRAL EQUATION
!  �� 1  �� N1+jw1  ����� ����� �� � lambda = , ������� F1,F2
!  �� N1+1      �� N1+N2+jw2    ������ �� �     ����������� ����� - ������      lambda =1  ������� F1,F2
!  �� N1+N2+jw2+1 �� N1+N2+jw2+N3  ����������� ���� lambda =1, ������� (0,-1)
!  �� N1+N2+jw2+N3+1 �� N1+N2+jw2+N3+N4 ����������� ����
!  o� N1+N2+jw2+N3+N4+1 �� N1+N2+jw2+N3+N4+N5 ��������� �����
!  ���� ��� ������� �����, ��   N3=0.
!   ���� ��� �������� �����, ��   N4=0.
!******************************************************
! Nk=1  ! test ��� ����� �������
do k=0,Nk
!���������� ������������� ���� � ������ �����
k_y(k)=exp(h_kln*k)-1.

do i=1,Nnodes
    if (i<=N1+N2  ) then  ! ����� �� ������� ECABD ��� - ����� - ���
        f_1=F1(x(i))!������� �� � �� �������
        f_2=F2(x(i))!������� �� z �� �������
        kappa=1. ! ������ - �����
        if (Model>=3 .and. Model<=4 ) then ! ������ �������� ���� 2,3,4
             if (x(i)<xC .or. x(i)>xE) then ! ��������� ����� �������
                kappa=(ro2-ro1)/(ro1+ro2)  ! ���� - �����
             endif
        endif
        if (Model.eq.2 .and. x(i)<xC ) kappa=(ro2-ro1)/(ro1+ro2)  ! ���� - �����
        if (Model.eq.5 .and. x(i)>xE) kappa=(ro2-ro1)/(ro1+ro2)
        if ((x(i)<-Apoint .or. x(i)>width-Apoint) .and. Model>3) kappa=(ro2-ro3)/(ro2+ro3) ! ������� ��������� - ����
        if ((x(i)<-Apoint) .and. (Model.eq.5 .or. Model.eq.6)) kappa=1.  ! ��� 5� � 6� ������ ����� ��� ����
        if ( x(i)>xL .and. x(i)<width-Apoint .and. (Model>=6)) kappa=0. !�������� - ���� ������

        h0=h

    else  ! i>N1+N2 ����������� ����, ��������-�����, ��������-���������, ���������-�����, ���������-����
        kappa=1. ! ��� ����������� ����
        h0=h
        f_1=0.
        f_2=-1.
        if (i>N1+N2+N3+N4 .and.(Model.eq.4 .or. Model .eq.5))  then
            kappa=(ro1-ro3)/(ro1+ro3) !  ��������� - �����
            h0=width/N5 ! ��� �� ���������
        endif
        if (i>N1+N2+N3+N4+1 .and. i<=N1+N2+N3+N4+N6 .and.(Model>=6))  then !
            kappa=(ro1-ro2)/(ro1+ro2) !  ����� (ro1) - ��������(ro2)
            h0=h4 ! ��� �� ����� ��������� �� ��������
            f_1=nx(xL,zL,xK,zK)
            f_2=nz(xL,zL,xK,zK)
        endif
        if (i>N1+N2+N3+N4+N6 .and. i<=N1+N2+N3+N4+N6+N5_1 .and.(Model>=6))  then !
            kappa=(ro1-ro3)/(ro1+ro3) !  ��������� - �����
            h0=h2 ! ��� �� ����� ��������� �� ��������
        endif
        if (i>N1+N2+N3+N4+N6+N5_1 .and. i<=N1+N2+N3+N4+N6+N5_1+N5_2 .and. Model>=6 )  then !
            kappa=(ro2-ro3)/(ro2+ro3) !  ��������� (ro3) - ��������(ro2)
            h0=h3 ! ��� �� ����� ��������� �� ��������
        endif

    endif

    if (x(i) .eq. xzar .and. z(i) .eq. zzar) then
        B(i)=0.
    else
!�������� ���������� �� ������ ������ �� ������
        a_zar=dist_a(x(i),z(i),xzar,zzar)
        if (k.eq.0)then
            B(i)=((x(i)-xzar)*f_1+(z(i)-zzar)*f_2)/a_zar/a_zar
        else
            k1_=besek1(a_zar*k_y(k))
            B(i)=((x(i)-xzar)*f_1+(z(i)-zzar)*f_2)*k1_*(k_y(k)/2./a_zar) ! ������ ����� ��� ���������
        endif
        B(i)=B(i)*2./pi*kappa !
   endif
do j=1,Nnodes
A(i,j)=0.
a_ij=dist_a(x(i),z(i),x(j),z(j))
if (a_ij==0.) cycle ! ��������� �������� �����  ��� ���������� ���� ����� ��� ����������� �����
if (k.eq.0)then
    A(i,j)=((x(i)-x(j))*f_1+(z(i)-z(j))*f_2)/a_ij/a_ij
else
    k1_=besek1(a_ij*k_y(k))
    A(i,j)=((x(i)-x(j))*f_1+(z(i)-z(j))*f_2)*k1_*(k_y(k)/2./a_ij)
endif
A(i,j)=A(i,j)*kappa*h0/pi

enddo !��������� ���� �� j
enddo !��������� ���� �� i
! ������ ������������ ������ ����� � ����� ������ ����� �������� �����
!B(N+1)=(B(N+2)+B(N))*2/3-(B(N+3)+B(N-1))/6
!B(N1+1)=(B(N1+2)+B(N1))/2
print *,' frequency', k_y(k)
! ������ ������ ����� � ������� ��� �� � ��������� �����
B1=matmul(A,B)+B
!B1=matmul(A,B1)+B1
A1=matmul(A,A)
!A1=matmul(A,A1)
! ���������� ������ �����
B1(N1)=(B1(N1-2)*2.+B1(N1-1)*2.+B1(N1)+B1(N1+1))/6.
B1(N1+1)=(B1(N1+1)+B1(N1+2)*2.+B1(N1+3)*2.+B1(N1))/6.
! ���������� ��������
B1(N1)=(B1(N1-2)*2.+B1(N1-1)*2.+B1(N1)+B1(N1+1))/6.
B1(N1+1)=(B1(N1+1)+B1(N1+2)*2.+B1(N1+3)*2.+B1(N1))/6.
!������������� ��������
do i=1,Nnodes
   nu(i)=B1(i)
enddo
   it=0
   delta=100000.;delta1=200000.
   do while((delta>0.00000001 .and. delta<delta1))
   delta1=delta
   nu1=matmul(A1,nu)+B1
   delta=dot_product(nu-nu1,nu-nu1)
   delta=sqrt(delta); it=it+1
   print *, 'iteration number=', it, ' delta=', delta
   nu=nu1
   if(delta>delta1) print *, 'Iterations diverge'
   enddo !��������� ���� �� while

    do i=1,Nnodes
       nu2D(i,k)=nu(i)
    enddo

enddo !��������� ���� �� k

! ������� ��������� ���������� ��� ������� �������
nu2D(:,0)=nu2D(:,1)*2-nu2D(:,2)
! ������ ������� ������ ��� ������ ���� klimit
do k=0,Nk
    nu2D(:,k)=nu2D(:,k)*exp(-(k_y(k)/klimit)**2);
    !nu2D(:,k)=nu2D(:,k)*FrequencyFilter(k*h_k,klimit,klimit*0.1);
enddo

deallocate(A,B) !������������ ������
deallocate(A1,B1) !������������ ������
deallocate(lambda,RelX,RelZ)
deallocate(dlambda,ReldZ)
allocate(nuxy(1:Nnodes,0:Ny))
! �������� �������������� ����� ��� ���������������� �������
do i=1,Nnodes
    do j=0,Ny
        tmp1=exp(h_yln*j)-1. !��������
        tmp=0.
        !do k=1,Nk
        !    tmp=tmp+nu2D(i,k)*(k_y(k)+1)*cos(k_y(k)*tmp1) !������� ���������������
        !enddo
        ! ������� ��������
        tmp=nu2D(i,0)*(k_y(0)+1.)*cos(k_y(0)*tmp1)
        tmp=(tmp+nu2D(i,Nk)*(k_y(Nk)+1.)*cos(k_y(Nk)*tmp1))/2
        do k=1,Nk-1
            tmp=tmp+nu2D(i,k)*(k_y(k)+1.)*cos(k_y(k)*tmp1)
        enddo
        nuxy(i,j)=tmp*h_kln/pi
    enddo
enddo

!*******************************************************
!*****************���������� ����������*****************
!*******************************************************
! ������ ���������� � ���������������� �����������
allocate(U_(1:kol))
do i=1,kol   ! ��������� ����������� � ������ � �������� xMN,zMN(i)
   tmp=0.
do j=0,Ny-1
   tmp1=h_yln*j
   yy=exp(tmp1)-1.
    do k=1,Nnodes
        if (xMN(i).ne. x(k) .and. zMN(i).ne. z(k)) then
        a_ij=(xMN(i)-x(k))*(xMN(i)-x(k))
        a_ij=a_ij+(zMN(i)-z(k))*(zMN(i)-z(k))
        a_ij=sqrt(a_ij+yy*yy)
        ! ���� ����� ��� ������ ��������
        h0=h
        if (i>N1+N2+N3+N4 .and.(Model.eq.4 .or. Model .eq.5))  then
            h0=width/N5 ! ��� �� ���������
        endif
        if (i>N1+N2+N3+N4+1 .and. i<=N1+N2+N3+N4+N6 .and.(Model>=6))  then !
             h0=h4 ! ��� �� ������� �������� - �����
        endif
        if (i>N1+N2+N3+N4+N6 .and. i<=N1+N2+N3+N4+N6+N5_1 .and.(Model>=6))  then !
             h0=h2 ! ��� �� ����� ��������� �� ��������
        endif
        if (i>N1+N2+N3+N4+N6+N5_1 .and. i<=N1+N2+N3+N4+N6+N5_1+N5_2 .and. Model>=6 )  then !
            h0=h3 ! ��� �� ����� ��������� �� ��������
        endif
        tmp=tmp+nuxy(k,j)*h0/a_ij*(yy+1.)
        endif
    enddo
enddo
U_(i)=tmp*h_yln/(pi*2.)
enddo

open (11,file="nuxy1.csv")
    do j=0,Ny
        write(11,'(*(G0.7,:,","))'), (nuxy(i,j),i=1,N1+N2)
    enddo
close(11)
open (10,file="nuxk.csv")
    do j=0,Nk
        write(10,'(*(G0.7,:,","))'), (nu2D(i,j),i=1,N1+N2)
    enddo
close(10)
open (11,file="x1.txt")
    do j=1,N1+N2
        write(11,10), x(j)
    enddo
close(11)
open (11,file="y1.txt")
    do j=0,Ny
        write(11,10), exp(h_yln*j)-1.
    enddo
close(11)
if (Model>1 .and. Model.ne.5 .and. Model .ne.6 ) then
open (11,file="nuWaterLeft.csv") ! ����������� ���� ����� N3 �����
    do j=0,Ny
        write(11,'(*(G0.7,:,","))'), (nuxy(i,j),i=N1+N2+1,N1+N2+N3)
    enddo
close(11)

open (11,file="x2.txt")
    do j=N1+N2+1,N1+N2+N3
        write(11,10), x(j)
    enddo
close(11)
endif
if (Model>2) then
open (10,file="nuWaterRight.csv") ! ����������� ���� ������ N4 �����
    do j=0,Ny
        write(10,'(*(G0.7,:,","))'), (nuxy(i,j),i=N1+N2+N3+1,N1+N2+N3+N4)
    enddo
close(10)
open (11,file="x3.txt")
    do j=N1+N2+N3+1,N1+N2+N3+N4
        write(11,10), x(j)
    enddo
close(11)
endif
if (Model.eq.4 .or. Model .eq.5) then
open (10,file="nuBase.txt") ! ��������� N5 �����
    do j=0,Ny
    write(10,10), (nuxy(i,j),i=N1+N2+N3+N4+1,N1+N2+N3+N4+N5)
    enddo
close(10)
endif

if (Model.eq.6 .or. Model .eq.7) then
open (10,file="nuBase1.txt") ! ��������� N5_1 �����
    do j=0,Ny
    write(10,10), (nuxy(i,j),i=N1+N2+N3+N4+N6+1,N1+N2+N3+N4+N6+N5_1)
    enddo
close(10)
open (11,file="nuBase2.txt") ! ��������� N5_2 �����
    do j=0,Ny
    write(11,10), (nuxy(i,j),i=N1+N2+N3+N4+N6+N5_1+1,N1+N2+N3+N4+N6+N5_1+N5_2)
    enddo
close(11)
open (10,file="nuLeaking.txt") ! �������� ����� N6 �����
    do j=0,Ny
    write(10,10), (nuxy(i,j),i=N1+N2+N3+N4+1,N1+N2+N3+N4+N6)
    enddo
close(10)
endif

allocate(RoKA(1:kol-1))
open(10, file='Roka1.csv')
write(10,'(*(G0.7,:,","))') 'x,M+0.5MN,roK'
!���������� RoK(i) �� �������
do i=1,kol-1
tmp=2.*pi*AM(i)*AN(i)/MN !����������� ��������� �� ��������� ��������� � ����� ������ ��� ������������������ ���������, �� ���� ������������ �� ����� ������ �
RoKA(i) = 1.+tmp*(U_(i)-U_(i+1)) !�������� ��������� ������������� ���������� �� ro1 � ������� ��
write (10,'(*(G0.7,:,","))'), (xMN(i)+xMN(i+1))/2.*ABDim*0.5, (AM(i)+AN(i))/2.*ABDim*0.5, RoKA(i)
enddo

call cpu_time(etime)
      print *, ' time for grid calculation ', etime-stime;


deallocate(RoKA,U_,AM,AN,k_y,xMN,zMN,x,z)
deallocate(nu,nu1)
deallocate(nu2D,nuxy)
10 format(2048(e15.7))

    CONTAINS
!*******************************************************
!************       ������      �� ������ RBF     ******
!*******************************************************
real function fi1(r)
real r
if (r>=0.) then
   fi1=exp(-eps*r*eps*r)   !Gaussian
else
   fi1=0.
endif
return
end function fi1

real function relief(xx)
integer i
real xx,s1
s1=0.
do i=1,jtab
!s1=s1+lambda(i)*fi2(abs(xx-RelX(i)))
s1=s1+lambda(i)*fi1(abs(xx-RelX(i)))
enddo
relief=s1
return
end function relief
!
real function drelief(xx)
integer i
real xx,s1
s1=0.
do i=1,jtab
!s1=s1+lambda(i)*fi2(abs(xx-RelX(i)))
s1=s1+dlambda(i)*fi1(abs(xx-RelX(i)))
enddo
drelief=s1
return
end function drelief

real function curve_element(x_)
real x_,t
t=drelief(x_)
curve_element=sqrt(1+t*t)
return
end function curve_element
!
SUBROUTINE DS(X_, Y_, DY)!��.����� ��� y'=1/dsqrt(1.+f(y)**2)
      real, intent(in):: X_, Y_
      real DY, t_
      t_=drelief(Y_)
      DY=1/sqrt(1+t_*t_)
end SUBROUTINE DS
SUBROUTINE DS1(X_, Y_, DY)
      real, intent(in):: X_, Y_
      real DY, t_
      t_=drelief(Y_)
      DY=1/t_
end SUBROUTINE DS1
!������� ���������� � ������� �������
!������� �� �
      real function F1(x)
      real x,t
      t=drelief(x)
      F1=t/sqrt(1+t*t)
      return
      end function F1

!������� �� z
      real function F2(x)
      real x,t_
      t_=drelief(x)
      F2=-1./sqrt(1.+t_*t_) ! ��� z ���������� ����
      return
      end function F2

!������� � ������� �������� ���������� � ������� �����
!������� �� �
      real function nx(xL,zL,xK,zK)
      real xL,zL,xK,zK,t
      t=(zL-zK)*(zL-zK)+(xL-xK)*(xL-xK) ! ������
      nx=(zL-zK)/sqrt(t)
      return
      end function nx

!������� �� z
      real function nz(xL,zL,xK,zK)
      real xL,zL,xK,zK,t
      t=(zL-zK)*(zL-zK)+(xL-xK)*(xL-xK) ! ������
      nz=(xK-xL)/sqrt(t)
      return
      end function nz

      real function dist_a(xi,zi,xj,zj) !
      real xi,zi,xj,zj,t1,t3
      t1=xi-xj;t3=zi-zj
      dist_a=sqrt(t1*t1+t3*t3)
      return
      end function dist_a

      real*8 function func_R_U(xi,zi,xj,zj)
      real xi,zi,xj,zj,t,t1,t2,t3
      t1=xi-xj; t3=zi-zj
      t=t1*t1+t3*t3
      t=sqrt(t)
      func_R_U=1.d+0/t
      return
      end function func_R_U


end program dambawater5

!**************	E X T E R N A L       F U N C T I O N S  ******
      real function integr1D(f,x1,x2)
      integer m,i
      real x1,x2,h,s1,s2,x
      real,external:: f
      m=10000
      s1=0.; h=(x2-x1)/m; x=x1+h/2
      do i =1, m*2-1,2
      s1 =  f(x)+ s1
      x = x+h
      end do
      S2 = 0;x =h
      do i =2, m*2-2,2
      S2 = f(x)+S2
      x = x + h
      end do
       integr1D = ( f(x1) + f(x2) + S1 * 4.+ S2 * 2.) * h/ 6.
      end function integr1D

