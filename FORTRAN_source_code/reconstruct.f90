program reconstuction

use ensemble_mod
use data_mod
use grad_solve_mod, only : max_length, max_time, tol, start_series
use setup_mod

implicit none

integer :: n_members
integer :: i,j,k
integer :: recon_length,dof
integer, allocatable :: series_length_test(:)
real, allocatable :: proxy_test(:,:,:)
real, allocatable :: recon_time_base(:)
real, allocatable :: hc_test(:),proxy_time(:)
real, allocatable :: ensemble_summary(:,:)
real (kind=r8_kind) :: corr,neg_corr,corr_zero
real :: orig_med(2),tmp1,tmp2,med_qr(2),r_crit
real :: tol_orig,opt_thresh,mn_time,mx_time,full_time
character(len=255) :: ft,f_temp,str
character(len=255), allocatable :: fp(:)
character :: yn
logical :: test_inversion

print *,'Enter filename for calibration target data series'
read(*,'(a)')ft
print *,'Enter number of proxy_orig data series'
read *,n_series
n_series=n_series+1
allocate (fp(n_series),series_length(3*n_series-1))
allocate (series_length_test(4))
allocate (proxy_time(n_series))

fp(1)=ft

do i=2,n_series
  print *,'Enter filename for proxy_orig data series ',i-1,' <optional max time for sign checking of this proxy>'
  read(*,'(a)')fp(i)
  j=index(fp(i),'<')
  k=index(fp(i),'>')
  if ((j.gt.0).and.(k.gt.0)) then
    f_temp=fp(i)
    read(f_temp(j+1:k-1),'(f8.3)')proxy_time(i)
    fp(i)=f_temp(1:j-1)
  else
    proxy_time(i)=-999
  endif
enddo

do i=1,n_series
  open(11,file=trim(fp(i)),form="formatted")
  series_length(i)=0
  do 
    read (11,*,end=10,err=10)tmp1,tmp2
    series_length(i)=series_length(i)+1
  enddo
  10 close(11)
enddo

max_length=maxval(series_length(1:n_series))
allocate (proxy(max_length,3*n_series-1,2))
allocate (proxy_test(max_length,4,2))
allocate (hc(3*n_series-1))
allocate (var(3*n_series-1))
allocate (hc_test(4))

mn_time=huge(mn_time)
mx_time=-mn_time
do i=1,n_series
  open(11,file=trim(fp(i)),form="formatted")
  do j=1,series_length(i)
    read (11,*,end=20,err=20)proxy(j,i,:)
  enddo
20 close(11)
  if (i.gt.1) then
    mn_time=min(mn_time,proxy(1,i,1))
    mx_time=max(mx_time,proxy(series_length(i),i,1))
  endif
enddo

print *,'Enter filename for time base of reconstruction'
read(*,'(a)')fp(1)

open(11,file=trim(fp(1)),form="formatted")
recon_length=0
do 
  read (11,*,end=30,err=30)tmp1
  if ((tmp1.ge.mn_time).and.(tmp1.le.mx_time)) recon_length=recon_length+1
enddo
30 close(11)

allocate (recon_time_base(1:recon_length))

open(11,file=trim(fp(1)),form="formatted")
j=1
do 
  read (11,*,end=40,err=40)tmp1
  if ((tmp1.ge.mn_time).and.(tmp1.le.mx_time)) then
   recon_time_base(j)=tmp1
   j=j+1
 endif
enddo
40 close(11)

print *,'Enter number of ensemble members (default 4096)'
read(*,FMT='(a)')str
if (len_trim(str).eq.0) then
  n_members=4096
else
  read(str,*)n_members
endif

print *,'Enter value for SLACK width parameter (h) (default 0.4,1.6)'
read(*,FMT='(a)')str
if (len_trim(str).eq.0) then
  hc(1)=0.4e0
  hc(n_series+1)=1.6e0
else
  read(str,*)hc(1),hc(n_series+1)
endif
hc(2:n_series)=hc(1)
hc(n_series+2:2*n_series-1)=hc(n_series+1)
hc(2*n_series:)=hc(1)

print *,'Enter maximum CPU time per ensemble member'
read *,full_time
do i=2,n_series
  if (proxy_time(i).lt.0) proxy_time(i)=full_time
enddo
print *,'Enter solution tolerance'
read *,tol

print *,'Enter filename for output'
read (*,'(a)')fp(1)

print *,'Attempt non-linear reconstruction (y/n)?'
read *,yn
non_linear=.false.
if ((yn.eq.'y').or.(yn.eq."Y")) non_linear=.true.

print *,'Enter percentile for stationarity optimisatipn (0.5 default)'
read(*,FMT='(a)')str
if (len_trim(str).eq.0) then
  tmp1=0.5e0
else
  read(str,*)tmp1
endif
tmp1=min(max(0.0,tmp1),1.0)
n_subset=n_members*tmp1

print *,'check for proxy inversion (y/n)?'
read *,yn
test_inversion=.false.
if ((yn.eq.'y').or.(yn.eq."Y")) test_inversion=.true.

open(11,file=trim(fp(1)),form="formatted")
write(11,'(a)')'!target data'
do j=1,series_length(1)
  write(11,'(a,2g15.6)')'!',proxy(j,1,:)
enddo
write(11,'(a)')'!'

orig_med(1)=M_estimator(proxy(1:series_length(1),1,2),series_length(1))
orig_med(2)=Qn(proxy(1:series_length(1),1,2),series_length(1))

proxy(1:series_length(1),1,2)=(proxy(1:series_length(1),1,2)-orig_med(1))/orig_med(2)


start_series=2
tol_orig=tol

do i=2,n_series
  proxy(:,i+n_series-1,:)=proxy(:,i,:)
  proxy(:,i+2*n_series-2,:)=proxy(:,i,:)
  series_length(i+n_series-1)=series_length(i)
  series_length(i+2*n_series-2)=series_length(i)
enddo

if (test_inversion) then
  proxy_test(:,1,:)=proxy(:,1,:)
  series_length_test(1)=series_length(1)
  hc_test(1)=hc(1)

  print *,'chekcing proxy sign (* denotes selected sign)'
  print *,'proxy   +ve proxy correlation   -ve proxy correlation    critical_value'
  write(11,'(a)')'!checking proxy sign (* denotes selected sign)'
  write(11,'(a)')'!proxy   +ve proxy correlation   -ve proxy correlation    critial_value'
  do i=2,n_series
  
    tol=tol_orig**2
    hc_test(2)=hc(i)
    hc_test(3)=hc(n_series+1)
    hc_test(4)=hc(i)
    series_length_test(2)=series_length(i)
    series_length_test(3)=series_length(i)
    series_length_test(4)=series_length(i)
    proxy_test(:,2,1)=proxy(:,i,1)
    proxy_test(:,3,1)=proxy(:,i,1)
    proxy_test(:,4,1)=proxy(:,i,1)
  
    proxy_test(:,2,2)=-proxy(:,i,2)
    proxy_test(:,3,2)=-proxy(:,i,2)
    proxy_test(:,4,2)=-proxy(:,i,2)

    max_time=proxy_time(i)
  
    call setup(2,proxy_test,series_length_test,hc_test, &
      recon_length,recon_time_base,.false.)
  
    ensemble_summary=calc_ensemble(n_members,recon_length,3, &
      series_length,recon_time_base, &
      orig_med,var,proxy,cal_index,hc,.false.,.true.) 
  
    neg_corr=correlation_full( &
      proxy(1:series_length(1),1,2), &
      ensemble_summary(1:recon_length,1), &
      proxy(1:series_length(1),1,1), &
      recon_time_base, &
      series_length(1),recon_length,dof=dof)
  
    proxy_test(:,2,2)=proxy(:,i,2)
    proxy_test(:,3,2)=proxy(:,i,2)
    proxy_test(:,4,2)=proxy(:,i,2)
  
    call setup(2,proxy_test,series_length_test,hc_test, &
      recon_length,recon_time_base,.false.)
  
    ensemble_summary=calc_ensemble(n_members,recon_length,3, &
      series_length,recon_time_base, &
      orig_med,var,proxy,cal_index,hc,.false.,.true.) 

    corr=correlation_full( &
      proxy(1:series_length(1),1,2), &
      ensemble_summary(1:recon_length,1), &
      proxy(1:series_length(1),1,1), &
      recon_time_base, &
      series_length(1),recon_length,dof=dof)
  
    r_crit=+r_critical(dof)/2e0
  
    if (neg_corr.gt.corr+r_crit) then
      proxy(:,i,2)=-proxy(:,i,2)
      write(*,'(i5,2f15.6,a,f15.6)')i-1,corr,neg_corr,'*',r_crit
      write(11,'(a,i5,2f15.6,a,f15.6)')'!',i-1,corr,neg_corr,'*',r_crit
    else
      write(*,'(i5,f15.6,a,2f15.6)')i-1,corr,'*',neg_corr,r_crit
      write(11,'(a,i5,f15.6,a,2f15.6)')'!',i-1,corr,'*',neg_corr,r_crit
    endif
  enddo
endif

write(11,'(a)')'!proxy_orig data'
do i=2,n_series
  write(11,'(a,i5)')'!series ',i-1
  do j=1,series_length(i)
    write(11,'(a,2g15.6)')'!',proxy(j,i,:)
  enddo
  write(11,'(a)')'!'
enddo

call setup(n_series,proxy,series_length,hc, &
  recon_length,recon_time_base,.true.)

if (non_linear) then
  n_series=3*n_series-2
else
  n_series=2*n_series-1
endif

tol=tol_orig**2*(n_series-start_series+1)

max_time=full_time

ensemble_summary=calc_ensemble(n_members,recon_length,n_series,series_length,recon_time_base, &
  orig_med,var,proxy,cal_index,hc,.true.,.false.) 

write(11,'(a)')"! time base     M-estimator   +/- for_95%_CI      Qn"
do i=1,recon_length
  write(11,*)recon_time_base(i),ensemble_summary(i,1), &
    ensemble_summary(i,3),ensemble_summary(i,2)
enddo
close(11)

deallocate(proxy,series_length,recon_time_base,fp)
deallocate(hc,ensemble_summary,var,proxy_time)

end program
