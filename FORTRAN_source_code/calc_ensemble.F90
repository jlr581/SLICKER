module ensemble_mod

use med_stats_mod
use grad_solve_mod
use correlation_mod
#ifdef openmp
use omp_lib, only: omp_get_num_threads
#endif

implicit none

integer :: n_subset

contains

function calc_ensemble(n_members,recon_length,n_series,series_length,recon_time_base, &
  orig_med,var,proxy,cal_index,hc,verbose,time_out)

integer, intent(in) :: n_members,recon_length
integer, intent(in) :: n_series,series_length(:),cal_index(2)
real, intent(in) :: recon_time_base(:),orig_med(:),proxy(:,:,:),hc(:)
real (kind=r8_kind), intent(in) :: var(:)
logical, intent(in) :: verbose,time_out

real :: calc_ensemble(recon_length,3),cal_time
real :: jackknife(n_members),jack_temp(n_members-1),jack_mean

integer :: member,i,j,k,l
real :: med_qr(2),hc_temp
real, allocatable :: ensemble(:,:)
real :: recon(1:recon_length)
real (kind=r8_kind) :: min_cost,station(n_members)
real (kind=r8_kind) :: eps,full_var(n_series),swap
real, allocatable :: initial_values(:,:)
real :: random(recon_length),dis,min_dis,slope,qn_proxy1
integer :: seed(33),extra

#ifdef openmp
!$OMP parallel default(shared)
extra=omp_get_num_threads()
!$OMP end parallel 
#else
extra=0
#endif

eps=epsilon(eps)

allocate (timebase_for_integration(3*(max_length+recon_length),1:2,2:n_series))
allocate (delta_x(0:3*(max_length+recon_length)+1,2:n_series))
allocate (recip_delta_x(0:3*(max_length+recon_length)+1,2:n_series))
allocate (interp_factor(0:3*(max_length+recon_length)+1,3:6,2:n_series))
allocate (orig_data_for_integ(0:3*(max_length+recon_length)+1,3:4,2:n_series))
allocate (index_for_integration(3*(max_length+recon_length),3:6,2:n_series))
#ifdef ifort
allocate (sum_x(2:n_series))
#endif
allocate (denom1(2:n_series))
allocate (timebase_length(2:n_series))
allocate (ensemble(n_members+extra,1:recon_length))
allocate (initial_values(n_members+extra,recon_length))

seed=1
number_valid=0

call random_seed(put=seed)

do i=1,n_members+extra
  ensemble(i,1)=huge(ensemble(i,1))
enddo

initial_values=0d0
qn_proxy1=Qn(proxy(1:series_length(1),1,2),series_length(1))
do j=start_series,n_series
  slope=var(j)*Qn_proxy1/Qn(proxy(1:series_length(j),j,2),series_length(j))
  do i=1,recon_length
    min_dis=huge(min_dis)
    do k=1,series_length(j)
      dis=abs(proxy(k,j,1)-recon_time_base(i))
      if (dis.lt.min_dis) then
        min_dis=dis
        l=k
      endif
    enddo
    if ((proxy(1,j,1).le.recon_time_base(i)).and.(proxy(series_length(j),j,1).ge.recon_time_base(i))) &
       initial_values(1,i)=initial_values(1,i)+slope*proxy(l,j,2)
  enddo
enddo
med_qr(1)=M_estimator(initial_values(1,:),recon_length)
med_qr(2)=Qn(initial_values(1,:),recon_length)

initial_values(1,:)=(initial_values(1,:)-med_qr(1))*orig_med(2)/med_qr(2)+orig_med(1)

do i=n_members+extra,1,-1
  call random_number(random(:))
  initial_values(i,:)=initial_values(1,:)+(random(:)-0.5)*orig_med(2)
enddo

! pre-calculate SLICK correlation linear interpolation coefficients and indices
!$OMP parallel do schedule(dynamic) default(shared), private(i)  
do i=start_series,n_series
  call correlation_setup(proxy(1:series_length(i),i,1), &
    recon_time_base,proxy(1:series_length(i),i,2), &
    i,series_length(i),recon_length,hc_optional=hc(i))
enddo
!$OMP end parallel do

if (verbose)  print *,'     solution_number         tolerance        time         valid'
!$OMP parallel do schedule(dynamic) default(shared), private(member,recon,min_cost)  
  do member=1,n_members+extra
    recon(:)=initial_values(member,:)
    call grad_solve(recon,recon_length,n_series,var,min_cost, &
      orig_med,cal_index,n_members,verbose,time_out)
    ensemble(member,:)=recon
    if (verbose.and.(mod(member,100).eq.0)) &
      write(11,'(a,i7,a,g15.6)')'!ensemble member',member,' solution residual',min_cost
  enddo
!$OMP end parallel do

do member=1,n_members+extra
  if (ensemble(member,1).ne.ensemble(member,1)) ensemble(member,1)=huge(ensemble(member,1))
enddo

#ifdef openmp
do member=1,n_members
  if (ensemble(member,1).eq.huge(ensemble(member,1))) then
    do i=n_members+1,n_members+extra
      if (ensemble(i,1).ne.huge(ensemble(i,1))) then
        ensemble(member,:)=ensemble(i,:)
        ensemble(i,:)=huge(ensemble(i,1))
        exit
      endif
    enddo
  endif
enddo
#endif

cal_time=recon_time_base(cal_index(2))-recon_time_base(cal_index(1))
station(:)=0d0
!$OMP parallel do schedule(dynamic) default(shared), private(i,l,hc_temp,full_var,j,k,min_cost)  
do i=1,n_members
  do l=start_series,n_series
    hc_temp=hc(l)
    full_var(l)=correlation_full( &
      ensemble(i,:), &
      proxy(1:series_length(l),l,2), &
      recon_time_base, &
      proxy(1:series_length(l),l,1), &
      recon_length,series_length(l),hc_optional=hc_temp)
  enddo
  j=1
  do k=recon_length,j+1,-1
    if (recon_time_base(k)-recon_time_base(j).le.cal_time/2) exit
  enddo
  do l=start_series,n_series
    hc_temp=hc(l)
    min_cost=correlation_full( &
      ensemble(i,j:k), &
      proxy(1:series_length(l),l,2), &
      recon_time_base(j:k), &
      proxy(1:series_length(l),l,1), &
      k-j+1,series_length(l),hc_optional=hc_temp)
    if (min_cost.gt.-10e0) station(i)=station(i)+(full_var(l)-min_cost)**2
  enddo
  do
    do k=recon_length,j+1,-1
      if (recon_time_base(k)-recon_time_base(j).le.cal_time) exit
    enddo
    do l=start_series,n_series
      hc_temp=hc(l)
      min_cost=correlation_full( &
        ensemble(i,j:k), &
        proxy(1:series_length(l),l,2), &
        recon_time_base(j:k), &
        proxy(1:series_length(l),l,1), &
        k-j+1,series_length(l),hc_optional=hc_temp)
      if (min_cost.gt.-10e0) station(i)=station(i)+(full_var(l)-min_cost)**2
    enddo
    do k=recon_length,j+1,-1
      if (recon_time_base(k)-recon_time_base(j).le.cal_time/2) exit
    enddo
    j=k
    if (j.ge.recon_length) exit
    if (recon_time_base(recon_length)-recon_time_base(j).le.cal_time/4) exit
  enddo
enddo
!$OMP end parallel do

min_cost=quick_select_r8(station,n_subset,station(n_subset),n_members)

j=1
k=n_members
do
  if (station(j).gt.min_cost) then
    swap=station(j)
    station(j)=station(k)
    station(k)=swap
    recon=ensemble(j,:)
    ensemble(j,:)=ensemble(k,:)
    ensemble(k,:)=recon
    k=k-1
  else
    j=j+1
  endif
  if (j.eq.k) exit
enddo

!$OMP parallel do schedule(dynamic) default(shared), private(i,med_qr)  
do i=1,recon_length
  med_qr(1)=M_estimator(ensemble(1:n_subset,i),n_subset)
  med_qr(2)=Qn(ensemble(1:n_subset,i),n_subset)

  calc_ensemble(i,1)=med_qr(1)
  calc_ensemble(i,2)=med_qr(2)
enddo
!$OMP end parallel do

med_qr(1)=M_estimator(calc_ensemble(cal_index(1):cal_index(2),1),cal_index(2)-cal_index(1)+1)
med_qr(2)=orig_med(2)/Qn(calc_ensemble(cal_index(1):cal_index(2),1),cal_index(2)-cal_index(1)+1)

do i=1,recon_length
  ensemble(1:n_members,i)=(ensemble(1:n_members,i)-med_qr(1))*med_qr(2)+orig_med(1)
enddo

!$OMP parallel do schedule(dynamic) default(shared), private(i,med_qr,j,l,k,jack_temp,jackknife,jack_mean)  
do i=1,recon_length
  med_qr(1)=M_estimator(ensemble(1:n_subset,i),n_subset)
  med_qr(2)=Qn(ensemble(1:n_subset,i),n_subset)
  do j=1,n_subset
    l=1
    do k=1,n_subset
      if (k.eq.j) cycle
      jack_temp(l)=ensemble(k,i)
      l=l+1
    enddo
    jackknife(j)=M_estimator(jack_temp(1:n_subset-1),n_subset-1)
  enddo
  jack_mean=0e0
  do j=1,n_subset
    jack_mean=jack_mean+jackknife(j)
  enddo
  jack_mean=jack_mean/(n_subset)
  calc_ensemble(i,3)=0e0
  do j=1,n_subset
    calc_ensemble(i,3)=calc_ensemble(i,3)+(jackknife(j)-jack_mean)**2
  enddo
  calc_ensemble(i,3)=sqrt(calc_ensemble(i,3)/(n_subset))* &
    1.95569e0*sqrt(real(n_subset))

  calc_ensemble(i,1)=med_qr(1)
  calc_ensemble(i,2)=med_qr(2)
enddo
!$OMP end parallel do

deallocate(timebase_for_integration,index_for_integration,timebase_length)
deallocate(delta_x,interp_factor,recip_delta_x,ensemble,initial_values)
deallocate(orig_data_for_integ,denom1)
#ifdef ifort
deallocate(sum_x)
#endif

end function

end module
