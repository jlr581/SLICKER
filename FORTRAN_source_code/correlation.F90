module correlation_mod

use data_mod

implicit none

real (kind=r8_kind), allocatable :: timebase_for_integration(:,:,:)
real (kind=r8_kind), allocatable ::  delta_x(:,:)
real (kind=r8_kind), allocatable ::  recip_delta_x(:,:)
real (kind=r8_kind), allocatable ::  interp_factor(:,:,:)
real (kind=r8_kind), allocatable ::  orig_data_for_integ(:,:,:)
#ifdef ifort
real (kind=r8_kind), allocatable ::  sum_x(:)
#endif
real (kind=r8_kind), allocatable ::  denom1(:)
integer, allocatable :: index_for_integration(:,:,:),timebase_length(:)
real (kind=r8_kind) :: one_third=1d0/3d0

contains

function r_critical(dof)

real :: r_critical
integer, intent(in) :: dof
integer, parameter :: n_crit=38

real :: r_crit(n_crit,2)

integer :: i
real :: tmp

r_crit(:,1)=(/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0, &
            15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,30.0,&
            35.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,150.0,300.0,500.0,1000.0/)
r_crit(:,2)=(/0.9969,0.9500,0.8783,0.8114,0.7545,0.7067,0.6664,0.6319,&
            0.6021,0.5760,0.5529,0.5324,0.5140,0.4973,0.4821,0.4683,0.4555,&
            0.4438,0.4329,0.4227,0.4132,0.4044,0.3961,0.3882,0.3809,0.3494,&
            0.3246,0.3044,0.2732,0.2500,0.2319,0.2172,0.2050,0.1946,0.1593,& 
            0.1129,0.0875,0.0619/)

if (dof.lt.1) then
  tmp=r_crit(1,2)
elseif (dof.gt.r_crit(n_crit,1)) then
  tmp=r_crit(n_crit,2)
else
  do i=1,n_crit-1
    if (dof.lt.r_crit(i,1)) exit
  enddo
  tmp=r_crit(i-1,2)+(r_crit(i,2)-r_crit(i-1,2))*real(dof-r_crit(i-1,1))/real(r_crit(i,1)-r_crit(i-1,1))
endif
r_critical=tmp

end function

recursive function quick_select_r8(list,k,piv,n) result(kth_element)

integer, intent(in) :: k,n
real (kind=r8_kind), intent(in) :: list(n),piv
real (kind=r8_kind) :: kth_element
real (kind=r8_kind) :: left(n),middle(n),right(n)
integer :: nl,nm,nr,i

nl=0
nm=0
nr=0
do i=1,n
  if (list(i).lt.piv) then
    nl=nl+1
    left(nl)=list(i)
  elseif (list(i).gt.piv) then
    nr=nr+1
    right(nr)=list(i)
  else
    nm=nm+1
    middle(nm)=list(i)
  endif
enddo

if (k.le.nl) then
  kth_element=quick_select_r8(left(1:nl),k,left(nl/2+1),nl)
elseif (k.gt.nl+nm) then
  kth_element=quick_select_r8(right(1:nr),k-nl-nm,right(nr/2+1),nr)
else
  kth_element=middle(1)
endif

end function quick_select_r8

subroutine merge_double(l,r,u,orig,n)
integer, intent(in) :: l,r,u,n
real (kind=r8_kind), intent(inout) :: orig(n)
integer i,j,k
real (kind=r8_kind) :: copy(n)
i=l
j=r
k=l
do 
  if (.not.((i.lt.r).and.(j.lt.u))) exit
  if (orig(i).lt.orig(j)) then
    copy(k)=orig(i) 
    i=i+1
  else 
    copy(k)=orig(j) 
    j=j+1
  endif
  k=k+1
enddo 
do
  if (.not.(i.lt.r)) exit  
  copy(k)=orig(i)
  i=i+1
  k=k+1
enddo 
do
  if (.not.(j.lt.u)) exit 
  copy(k)=orig(j) 
  j=j+1
  k=k+1
enddo
orig(l:k-1)=copy(l:k-1) 
end subroutine

subroutine merge_sort_double(a,n)

integer, intent(in) :: n
real (kind=r8_kind), intent(inout) :: a(n)
integer :: k,u,i
real (kind=r8_kind) :: orig(n)
real (kind=r8_kind) :: copy(n)

orig(1:n)=a(1:n)

k=1
do
  if (.not.(k.lt.n)) exit
  i=1
  do
  if (.not.(i+k.lt.n+1)) exit
    u=i+k*2
    if (u.gt.n) u=n+1
    call merge_double(i,i+k,u,orig,n)
    i=i+k*2
  enddo
  k=k*2
enddo
a(1:n)=orig(1:n)
end subroutine

function correlation_full(y1,y2,x1,x2,n,m,hc_optional,dof)

real (kind=r8_kind) :: correlation_full
integer, intent(in) :: n,m
real, intent(in) :: y1(n),y2(m),x1(n),x2(m)
real, optional, intent(inout) :: hc_optional
integer, optional, intent(out) :: dof
integer :: i,j,k,max_samples
integer :: start_1,start_2,end_1,end_2,n_xbase
integer :: idx1,idx2
real (kind=r8_kind) :: hc
real (kind=r8_kind) :: mean1,mean2,xmin,xmax
real (kind=r8_kind) :: num,den1,den2
real (kind=r8_kind) :: dx(max(n,m)-1)
real (kind=r8_kind) :: median1,median2,iqr1,iqr2,temp1,temp2
real (kind=r8_kind) :: xbase(3*(n+m))
logical :: valid1,valid2
real (kind=r8_kind) :: a1,a2,b1,b2
real (kind=r8_kind) :: delta_x
real (kind=r8_kind) :: data_for_integration(3*(n+m),6)

hc=0.4d0
if(present(hc_optional)) hc=hc_optional

! find valid domain for both series, discard any data more than one point
! outside range of other data series
xmin=max(x1(1),x2(1))
xmax=min(x1(n),x2(m))
if (xmin.eq.x1(1)) then
  start_1=1
  do start_2=2,m
    if (x2(start_2).gt.xmin) exit
  enddo
  start_2=start_2-1
else
  start_2=1
  do start_1=2,n
    if (x1(start_1).gt.xmin) exit
  enddo
  start_1=start_1-1
endif
if (xmax.eq.x1(n)) then
  end_1=n
  do end_2=m-1,1,-1
    if (x2(end_2).lt.xmax) exit
  enddo
  end_2=end_2+1
else
  end_2=m
  do end_1=n-1,1,-1
    if (x1(end_1).lt.xmax) exit
  enddo
  end_1=end_1+1
endif

if ((end_1-start_1.le.2).or.(end_2-start_2.le.2)) then
  correlation_full=-huge(correlation_full)
  return
endif

! find median and interquartile range of data spacing
dx(1:end_1-start_1)=x1(start_1+1:end_1)-x1(start_1:end_1-1)
if (mod(end_1-start_1,2).eq.1) then ! odd number of points
  median1=quick_select_r8(dx(1:end_1-start_1),(end_1-start_1)/2+1,dx((end_1-start_1)/2),end_1-start_1)
  i=floor((end_1-start_1)/4d0)
  if (mod(end_1-start_1,4).eq.1) then
    temp1=(quick_select_r8(dx(1:end_1-start_1),i,dx((end_1-start_1)/2),end_1-start_1)+ &
       3d0*quick_select_r8(dx(1:end_1-start_1),i+1,dx((end_1-start_1)/2),end_1-start_1))/4d0
    temp2=(3d0*quick_select_r8(dx(1:end_1-start_1),3*i+1,dx((end_1-start_1)/2),end_1-start_1)+ &
        quick_select_r8(dx(1:end_1-start_1),3*i+2,dx((end_1-start_1)/2),end_1-start_1))/4d0
  else
    temp1=(3d0*quick_select_r8(dx(1:end_1-start_1),i+1,dx((end_1-start_1)/2),end_1-start_1)+ & 
       quick_select_r8(dx(1:end_1-start_1),i+2,dx((end_1-start_1)/2),end_1-start_1))/4d0
    temp2=(quick_select_r8(dx(1:end_1-start_1),3*i+2,dx((end_1-start_1)/2),end_1-start_1)+ &
       3d0*quick_select_r8(dx(1:end_1-start_1),3*i+3,dx((end_1-start_1)/2),end_1-start_1))/4d0
  endif
else
  median1=(quick_select_r8(dx(1:end_1-start_1),(end_1-start_1)/2,dx((end_1-start_1)/2),end_1-start_1)+ &
      quick_select_r8(dx(1:end_1-start_1),(end_1-start_1)/2+1,dx((end_1-start_1)/2),end_1-start_1))/2d0
  temp1=(quick_select_r8(dx(1:end_1-start_1),(end_1-start_1)/4,dx((end_1-start_1)/2),end_1-start_1)+ &
      quick_select_r8(dx(1:end_1-start_1),(end_1-start_1)/4+1,dx((end_1-start_1)/2),end_1-start_1))/2d0
  temp2=(quick_select_r8(dx(1:end_1-start_1),3*(end_1-start_1)/4,dx((end_1-start_1)/2),end_1-start_1)+ &
      quick_select_r8(dx(1:end_1-start_1),3*(end_1-start_1)/4,dx((end_1-start_1)/2),end_1-start_1))/2d0
endif
iqr1=temp2-temp1

dx(1:end_2-start_2)=x2(start_2+1:end_2)-x2(start_2:end_2-1)
if (mod(end_2-start_2,2).eq.1) then ! odd number of points
  median2=quick_select_r8(dx(1:end_2-start_2),(end_2-start_2)/2+1,dx((end_2-start_2)/2),end_2-start_2)
  i=floor((end_2-start_2)/4d0)
  if (mod(end_2-start_2,4).eq.1) then
    temp1=(quick_select_r8(dx(1:end_2-start_2),i,dx((end_2-start_2)/2),end_2-start_2)+ &
       3d0*quick_select_r8(dx(1:end_2-start_2),i+1,dx((end_2-start_2)/2),end_2-start_2))/4d0
    temp2=(3d0*quick_select_r8(dx(1:end_2-start_2),3*i+1,dx((end_2-start_2)/2),end_2-start_2)+ &
        quick_select_r8(dx(1:end_2-start_2),3*i+2,dx((end_2-start_2)/2),end_2-start_2))/4d0
  else
    temp1=(3d0*quick_select_r8(dx(1:end_2-start_2),i+1,dx((end_2-start_2)/2),end_2-start_2)+ & 
       quick_select_r8(dx(1:end_2-start_2),i+2,dx((end_2-start_2)/2),end_2-start_2))/4d0
    temp2=(quick_select_r8(dx(1:end_2-start_2),3*i+2,dx((end_2-start_2)/2),end_2-start_2)+ &
       3d0*quick_select_r8(dx(1:end_2-start_2),3*i+3,dx((end_2-start_2)/2),end_2-start_2))/4d0
  endif
else
  median2=(quick_select_r8(dx(1:end_2-start_2),(end_2-start_2)/2,dx((end_2-start_2)/2),end_2-start_2)+ &
      quick_select_r8(dx(1:end_2-start_2),(end_2-start_2)/2+1,dx((end_2-start_2)/2),end_2-start_2))/2d0
  temp1=(quick_select_r8(dx(1:end_2-start_2),(end_2-start_2)/4,dx((end_2-start_2)/2),end_2-start_2)+ &
      quick_select_r8(dx(1:end_2-start_2),(end_2-start_2)/4+1,dx((end_2-start_2)/2),end_2-start_2))/2d0
  temp2=(quick_select_r8(dx(1:end_2-start_2),3*(end_2-start_2)/4,dx((end_2-start_2)/2),end_2-start_2)+ &
      quick_select_r8(dx(1:end_2-start_2),3*(end_2-start_2)/4,dx((end_2-start_2)/2),end_2-start_2))/2d0
endif
iqr2=temp2-temp1

! set interval around points
hc=hc*max(median1,median2,iqr1,iqr2)

if(present(hc_optional)) hc_optional=hc

! construct points where kernel starts, midpoint, stops
j=1
do i=start_1,end_1
  if (x1(i)-hc.gt.x1(start_1)) then
    xbase(j)=x1(i)-hc
    j=j+1
  endif
  xbase(j)=x1(i)
  j=j+1
  if (x1(i)+hc.lt.x1(end_1)) then
    xbase(j)=x1(i)+hc
    j=j+1
  endif
enddo
do i=start_2,end_2
  if (x2(i)-hc.gt.x2(start_2)) then
    xbase(j)=x2(i)-hc
    j=j+1
  endif
  xbase(j)=x2(i)
  j=j+1
  if (x2(i)+hc.lt.x2(end_2)) then
    xbase(j)=x2(i)+hc
    j=j+1
  endif
enddo
n_xbase=j-1
call merge_sort_double(xbase(1:n_xbase),n_xbase)
! remove duplicates
i=2
do
  if (xbase(i).eq.xbase(i-1)) then
    xbase(i:n_xbase-1)=xbase(i+1:n_xbase)
    n_xbase=n_xbase-1
  else
    i=i+1
  endif
  if (i.gt.n_xbase) exit
enddo

! generate points of valid data for integration
idx1=2
idx2=2
j=1
do i=2,n_xbase
  if (xbase(i).le.max(x1(start_1),x2(start_2))) cycle
  if (xbase(i).gt.min(x1(end_1),x2(end_2))) cycle
  ! find both series have a data point within valid interval
  valid1=.false.
  valid2=.false.
  idx1=max(1,idx1-1)
  idx2=max(1,idx2-1)
  do k=idx1,n
    if (abs(x1(k)-(xbase(i-1)+xbase(i))/2d0).le.hc) then
      valid1=.true.
      idx1=max(1,k-1)
      exit
    endif
  enddo
  do k=idx2,m
    if (abs(x2(k)-(xbase(i-1)+xbase(i))/2d0).le.hc) then
      valid2=.true.
      idx2=max(1,k-1)
      exit
    endif
  enddo
  if (valid1.and.valid2) then
    data_for_integration(j,1:2)=xbase(i-1:i)
    ! linear interpolate
    do k=max(1,idx1-1),n-1
      if (x1(k).ge.xbase(i-1)) exit
    enddo
    k=max(1,k-1)
    data_for_integration(j,3)=y1(k)+(y1(k+1)-y1(k))*(xbase(i-1)-x1(k))/(x1(k+1)-x1(k))
    do k=max(1,idx1-1),n-1
      if (x1(k).ge.xbase(i)) exit
    enddo
    k=max(1,k-1)
    data_for_integration(j,4)= &
      (y1(k)+(y1(k+1)-y1(k))*(xbase(i)-x1(k))/(x1(k+1)-x1(k))- &
      data_for_integration(j,3))/(data_for_integration(j,2)- &
      data_for_integration(j,1))
    do k=max(1,idx2-1),m-1
      if (x2(k).ge.xbase(i-1)) exit
    enddo
    k=max(1,k-1)
    data_for_integration(j,5)=y2(k)+(y2(k+1)-y2(k))*(xbase(i-1)-x2(k))/(x2(k+1)-x2(k))
    do k=max(1,idx2-1),m-1
      if (x2(k).ge.xbase(i)) exit
    enddo
    k=max(1,k-1)
    data_for_integration(j,6)= &
      (y2(k)+(y2(k+1)-y2(k))*(xbase(i)-x2(k))/(x2(k+1)-x2(k))- &
      data_for_integration(j,5))/(data_for_integration(j,2)- &
      data_for_integration(j,1))
    j=j+1
  endif
enddo

j=j-1
if (present(dof)) dof=j

! calculate means by integration
mean1=0d0
mean2=0d0
num=0d0
do i=1,j
  delta_x=data_for_integration(i,2)-data_for_integration(i,1)
  a1=data_for_integration(i,3)
  b1=data_for_integration(i,4)
  mean1=mean1+delta_x*(a1+delta_x*b1*0.5d0)
  a2=data_for_integration(i,5)
  b2=data_for_integration(i,6)
  mean2=mean2+delta_x*(a2+delta_x*b2*0.5d0)
  num=num+delta_x
enddo
mean1=mean1/num
mean2=mean2/num

! calculate correlation by integration
data_for_integration(1:j,3)=data_for_integration(1:j,3)-mean1
data_for_integration(1:j,5)=data_for_integration(1:j,5)-mean2
num=0d0
den1=0d0
den2=0d0
do i=1,j
  delta_x=data_for_integration(i,2)-data_for_integration(i,1)
  a1=data_for_integration(i,3)
  b1=data_for_integration(i,4)
  a2=data_for_integration(i,5)
  b2=data_for_integration(i,6)
  num=num+delta_x*(a1*a2+delta_x*(b1*a2/2d0+b2*a1/2d0+delta_x*b1*b2*one_third))
  den1=den1+delta_x*(a1**2+delta_x*(b1*a1+delta_x*b1**2*one_third))
  den2=den2+delta_x*(a2**2+delta_x*(b2*a2+delta_x*b2**2*one_third))
enddo

correlation_full=num/sqrt(den1*den2)

end function

function correlation(y2,p,n,m)

real (kind=r8_kind) :: correlation
integer, intent(in) :: p,n,m
real, intent(in) :: y2(m)
integer :: i,k
real (kind=r8_kind) :: mean2
real (kind=r8_kind) :: num,den2
real (kind=r8_kind) :: a1,a2,b1,b2
real (kind=r8_kind) :: data_for_integration(n,3:6)
#ifdef ifort
real :: y5(n+m),dy5(n+m),y6(n+m),dy6(n+m)

do i=1,timebase_length(p)
  k=index_for_integration(i,5,p)
  y5(i)=y2(k)
  dy5(i)=y2(k+1)-y5(i)
  k=index_for_integration(i,6,p)
  y6(i)=y2(k)
  dy6(i)=y2(k+1)-y6(i)
enddo

do i=1,timebase_length(p)
  data_for_integration(i,5)=y5(i)+dy5(i)*interp_factor(i,5,p)
  data_for_integration(i,6)=(y6(i)+dy6(i)*interp_factor(i,6,p)-data_for_integration(i,5))*recip_delta_x(i,p)
enddo

#else

do i=1,timebase_length(p)
  k=index_for_integration(i,5,p)
  data_for_integration(i,5)=y2(k)+(y2(k+1)-y2(k))*interp_factor(i,5,p)
  k=index_for_integration(i,6,p)
  data_for_integration(i,6)=(y2(k)+(y2(k+1)-y2(k))*interp_factor(i,6,p)-data_for_integration(i,5))*recip_delta_x(i,p)
enddo

#endif


! calculate means by integration
mean2=0d0
#ifndef ifort
num=0d0
#endif
do i=1,timebase_length(p)
  a2=data_for_integration(i,5)
  b2=data_for_integration(i,6)
  mean2=mean2+delta_x(i,p)*(a2+delta_x(i,p)*b2*0.5d0)
#ifndef ifort
  num=num+delta_x(i,p)
#endif
enddo
#ifdef ifort
mean2=mean2/sum_x(p)
#else
mean2=mean2/num
#endif


! calculate correlation by integration
data_for_integration(1:timebase_length(p),5)=data_for_integration(1:timebase_length(p),5)-mean2
num=0d0
den2=0d0
do i=1,timebase_length(p)
  a1=orig_data_for_integ(i,3,p)
  b1=orig_data_for_integ(i,4,p)
  a2=data_for_integration(i,5)
  b2=data_for_integration(i,6)
  num=num+delta_x(i,p)*(a1*a2+delta_x(i,p)*(b1*a2*0.5d0+b2*(a1*0.5d0+delta_x(i,p)*b1*one_third)))
  den2=den2+delta_x(i,p)*(a2**2+delta_x(i,p)*b2*(a2+delta_x(i,p)*b2*one_third))
enddo

correlation=num/sqrt(denom1(p)*den2)

end function

subroutine correlation_setup(x1,x2,y1,p,n,m,hc_optional)

integer, intent(in) :: p,n,m
real, intent(in) :: x1(n),x2(m),y1(n)
real, optional, intent(in) :: hc_optional
integer :: i,j,k,max_samples
integer :: start_1,start_2,end_1,end_2,n_xbase
integer :: idx1,idx2
real (kind=r8_kind) :: hc
real (kind=r8_kind) :: xmin,xmax
real (kind=r8_kind) :: dx(max(n,m)-1)
real (kind=r8_kind) :: median1,median2,iqr1,iqr2,temp1,temp2
real (kind=r8_kind) :: xbase(3*(n+m))
real (kind=r8_kind) :: a1,b1,mean1,den
logical :: valid1,valid2

hc=0.4d0
if(present(hc_optional)) hc=hc_optional

! find valid domain for both series, discard any data more than one point
! outside range of other data series
xmin=max(x1(1),x2(1))
xmax=min(x1(n),x2(m))
if (xmin.eq.x1(1)) then
  start_1=1
  do start_2=2,m
    if (x2(start_2).gt.xmin) exit
  enddo
  start_2=start_2-1
else
  start_2=1
  do start_1=2,n
    if (x1(start_1).gt.xmin) exit
  enddo
  start_1=start_1-1
endif
if (xmax.eq.x1(n)) then
  end_1=n
  do end_2=m-1,1,-1
    if (x2(end_2).lt.xmax) exit
  enddo
  end_2=end_2+1
else
  end_2=m
  do end_1=n-1,1,-1
    if (x1(end_1).lt.xmax) exit
  enddo
  end_1=end_1+1
endif

! construct points where kernel starts, midpoint, stops
j=1
do i=start_1,end_1
  if (x1(i)-hc.gt.x1(start_1)) then
    xbase(j)=x1(i)-hc
    j=j+1
  endif
  xbase(j)=x1(i)
  j=j+1
  if (x1(i)+hc.lt.x1(end_1)) then
    xbase(j)=x1(i)+hc
    j=j+1
  endif
enddo
do i=start_2,end_2
  if (x2(i)-hc.gt.x2(start_2)) then
    xbase(j)=x2(i)-hc
    j=j+1
  endif
  xbase(j)=x2(i)
  j=j+1
  if (x2(i)+hc.lt.x2(end_2)) then
    xbase(j)=x2(i)+hc
    j=j+1
  endif
enddo
n_xbase=j-1
call merge_sort_double(xbase(1:n_xbase),n_xbase)
! remove duplicates 
i=2
do
  if (xbase(i).eq.xbase(i-1)) then
    xbase(i:n_xbase-1)=xbase(i+1:n_xbase)
    n_xbase=n_xbase-1
  else
    i=i+1
  endif
  if (i.gt.n_xbase) exit
enddo


! generate points of valid data for integration
idx1=2
idx2=2
j=1
do i=2,n_xbase
  if (xbase(i).le.max(x1(start_1),x2(start_2))) cycle
  if (xbase(i).gt.min(x1(end_1),x2(end_2))) cycle
  ! find both series have a data point within valid interval
  valid1=.false.
  valid2=.false.
  idx1=max(1,idx1-1)
  idx2=max(1,idx2-1)
  do k=idx1,n
    if (abs(x1(k)-(xbase(i-1)+xbase(i))/2d0).le.hc) then
      valid1=.true.
      idx1=max(1,k-1)
      exit
    endif
  enddo
  do k=idx2,m
    if (abs(x2(k)-(xbase(i-1)+xbase(i))/2d0).le.hc) then
      valid2=.true.
      idx2=max(1,k-1)
      exit
    endif
  enddo
  if (valid1.and.valid2) then
    timebase_for_integration(j,1:2,p)=xbase(i-1:i)
    do k=max(1,idx1-1),n-1
      if (x1(k).ge.xbase(i-1)) exit
    enddo
    k=max(1,k-1)
    index_for_integration(j,3,p)=k
    do k=max(1,idx1-1),n-1
      if (x1(k).ge.xbase(i)) exit
    enddo
    k=max(1,k-1)
    index_for_integration(j,4,p)=k
    do k=max(1,idx2-1),m-1
      if (x2(k).ge.xbase(i-1)) exit
    enddo
    k=max(1,k-1)
    index_for_integration(j,5,p)=k
    do k=max(1,idx2-1),m-1
      if (x2(k).ge.xbase(i)) exit
    enddo
    k=max(1,k-1)
    index_for_integration(j,6,p)=k
    j=j+1
  endif
enddo

timebase_length(p)=j-1

do i=1,timebase_length(p)
  delta_x(i,p)=timebase_for_integration(i,2,p)-timebase_for_integration(i,1,p)
enddo
do i=1,timebase_length(p)
  recip_delta_x(i,p)=1d0/delta_x(i,p)
enddo

do i=1,timebase_length(p)
  k=index_for_integration(i,3,p)
  interp_factor(i,3,p)=(timebase_for_integration(i,1,p)-x1(k))/(x1(k+1)-x1(k))
  k=index_for_integration(i,4,p)
  interp_factor(i,4,p)=(timebase_for_integration(i,2,p)-x1(k))/(x1(k+1)-x1(k))
  k=index_for_integration(i,5,p)
  interp_factor(i,5,p)=(timebase_for_integration(i,1,p)-x2(k))/(x2(k+1)-x2(k))
  k=index_for_integration(i,6,p)
  interp_factor(i,6,p)=(timebase_for_integration(i,2,p)-x2(k))/(x2(k+1)-x2(k))
enddo

do i=1,timebase_length(p)
  k=index_for_integration(i,3,p)
  orig_data_for_integ(i,3,p)=y1(k)+(y1(k+1)-y1(k))*interp_factor(i,3,p)
  k=index_for_integration(i,4,p)
  orig_data_for_integ(i,4,p)=(y1(k)+(y1(k+1)-y1(k))*interp_factor(i,4,p)-orig_data_for_integ(i,3,p))*recip_delta_x(i,p)
enddo

! calculate means by integration
mean1=0d0
#ifdef ifort
sum_x(p)=0d0
#else
den=0d0
#endif
do i=1,timebase_length(p)
  a1=orig_data_for_integ(i,3,p)
  b1=orig_data_for_integ(i,4,p)
  mean1=mean1+delta_x(i,p)*(a1+delta_x(i,p)*b1*0.5d0)
#ifdef ifort
  sum_x(p)=sum_x(p)+delta_x(i,p)
enddo
mean1=mean1/sum_x(p)
#else
  den=den+delta_x(i,p)   
enddo
mean1=mean1/den
#endif

! calculate correlation by integration
orig_data_for_integ(1:timebase_length(p),3,p)=orig_data_for_integ(1:timebase_length(p),3,p)-mean1

denom1(p)=0d0
do i=1,timebase_length(p)
  a1=orig_data_for_integ(i,3,p)
  b1=orig_data_for_integ(i,4,p)
  denom1(p)=denom1(p)+delta_x(i,p)*(a1**2+delta_x(i,p)*(b1*a1+delta_x(i,p)*b1**2*one_third))
enddo

end subroutine

end module
