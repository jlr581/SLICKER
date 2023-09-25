module med_stats_mod

implicit none

contains

subroutine Hoare_partition(list,init_left,init_right,k,n,k_prime) 

integer, intent(in) :: n,k,init_left,init_right
integer, intent(out) :: k_prime
real, intent(inout) :: list(n)

real :: swap
integer :: left,right

swap=list(init_left)
list(init_left)=list(k)
list(k)=swap

left=init_left+1
right=init_right

do
  do
    if (left.gt.right) exit
    if (list(left).gt.list(init_left)) exit
    left=left+1
  enddo
  if (left.gt.right) exit
  do
    if (list(init_left).ge.list(right)) exit
    right=right-1
  enddo
  if (left.ge.right) exit
  swap=list(left)
  list(left)=list(right)
  list(right)=swap
  left=left+1
  right=right-1
enddo
swap=list(init_left)
list(init_left)=list(left-1)
list(left-1)=swap
k_prime=left-1
                                                                                
end subroutine Hoare_partition

function quick_select(list,k,n) result(kth_element)

integer, intent(in) :: k,n
real, intent(in) :: list(n)
real :: kth_element
real :: copy(n)

integer :: left,right,k_prime,k_new

copy(:)=list(:)

left=1
right=n
k_prime=k

do
  if (left.gt.right) exit
  call Hoare_partition(copy(1:n),left,right,k_prime,n,k_new)
  if (k_new.lt.k_prime) then
    left=k_new+1
  elseif (k_new.gt.k_prime) then
    right=k_new-1
  else
    kth_element=copy(k_prime)
    exit
  endif
enddo

end function quick_select

function med_stats(a,n)
real :: med_stats(2)

integer, intent(in) :: n
real, intent(in) :: a(n)
integer :: pl,pu,i
real :: b(n),med,mad,al,au,sm

if (mod(n,2).eq.1) then ! odd
  med=quick_select(a,floor(real(n)/2)+1,n)
else ! even
  med=(quick_select(a,floor(real(n)/2),n)+quick_select(a,floor(real(n)/2)+1,n))/2.0e0
endif
pl=floor(real((n+2)/4))
pu=floor(real((3*n+2)/4))
al=real(mod(n+2,4))/4.0
au=real(mod(3*n+2,4))/4.0
sm=(1.0-au)*quick_select(a,pu,n)+au*quick_select(a,pu+1,n)+ &
  ((1.0-al)*quick_select(a,pl,n)+al*quick_select(a,pl+1,n))

med_stats(1)=(2e0*med+sm)/4e0

do i=1,n
  b(i)=abs(med_stats(1)-a(i))
enddo

if (mod(n,2).eq.1) then ! odd
  mad=quick_select(b,floor(real(n)/2)+1,n)
else ! even
  mad=(quick_select(b,floor(real(n)/2),n)+quick_select(b,floor(real(n)/2)+1,n))/2.0e0
endif

med_stats(2)=mad

end function

function MAD(a,n)
real :: MAD(2)

integer, intent(in) :: n
real, intent(in) :: a(n)
integer :: i
real :: b(n),med

if (mod(n,2).eq.1) then ! odd
  med=quick_select(a,floor(real(n)/2)+1,n)
else ! even
  med=(quick_select(a,floor(real(n)/2),n)+quick_select(a,floor(real(n)/2)+1,n))/2.0e0
endif

MAD(1)=med

do i=1,n
  b(i)=abs(a(i)-med)
enddo

if (mod(n,2).eq.1) then ! odd
  med=quick_select(b,floor(real(n)/2)+1,n)
else ! even
  med=(quick_select(b,floor(real(n)/2),n)+quick_select(b,floor(real(n)/2)+1,n))/2.0e0
endif

MAD(2)=med

end function

function M_estimator(a,n)
real :: M_estimator

integer, intent(in) :: n
real, intent(in) :: a(n)
integer :: i,u,l
real :: mad_stats(2),sm

mad_stats=MAD(a,n)

mad_stats(2)=mad_stats(2)/0.6745e0

sm=0e0
u=0
l=0
do i=1,n
  if (abs(a(i)-mad_stats(1)).gt.1.28e0*mad_stats(2)) then
    if (a(i).lt.mad_stats(1)) then
      l=l+1
    else
      u=u+1
    endif
  else
    sm=sm+a(i)
  endif
enddo

M_estimator=(1.28e0*mad_stats(2)*(u-l)+sm)/(n-l-u)

end function

function Qn(a,n)
real :: Qn

integer, intent(in) :: n
real, intent(in) :: a(n)
integer :: i,j,l,pl
real :: b(n*(n-1)/2),lq,al,d

l=1
do i=2,n
  do j=1,i-1
    b(l)=abs(a(i)-a(j))
    l=l+1
  enddo
enddo
l=l-1

pl=floor(real((l+2)/4))
al=real(mod(l+2,4))/4.0
lq=(1.0-al)*quick_select(b(1:l),pl,l)+al*quick_select(b(1:l),pl+1,l)

d=2.2219e0
if (mod(n,2).eq.1) then ! odd
  d=d*real(n)/(n+1.4e0)
else ! even
  d=d*real(n)/(n+3.8e0)
endif

Qn=d*lq

end function

function cube_root(a,n)

integer, intent(in) :: n
real, intent(in) :: a(n)
real :: cube_root(n)
integer :: i

do i=1,n
  if (a(i).ne.0) then
    cube_root(i)=sign(exp(log(abs(a(i)))/3e0),a(i))
  else
    cube_root(i)=0e0
  endif
enddo
  
end function

function alt_med_stats(a,n)
real :: alt_med_stats(2)

integer, intent(in) :: n
real, intent(in) :: a(n)
integer :: i
real :: mean,std_dev

mean=a(1)
do i=2,n
  mean=mean+a(i)
enddo
mean=mean/n
std_dev=0e0
do i=1,n
  std_dev=std_dev+(a(i)-mean)**2
enddo
alt_med_stats(1)=mean
alt_med_stats(2)=sqrt(std_dev/n)

end function

end module
