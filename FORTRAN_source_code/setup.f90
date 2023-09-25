module setup_mod

use correlation_mod
use data_mod
use med_stats_mod
use grad_solve_mod, only : max_length,start_series

implicit none

logical :: non_linear
integer :: cal_index(2),n_series
integer, allocatable :: series_length(:)
real (kind=r8_kind), allocatable :: var(:)
real, allocatable :: hc(:)
real, allocatable :: proxy(:,:,:)

contains

subroutine setup(n_series_orig,proxy_orig,series_length_orig, &
    hc_orig,recon_length,recon_time_base,verbose)

integer, intent(in) :: n_series_orig,series_length_orig(:),recon_length
real, intent(inout) :: proxy_orig(:,:,:)
real, intent(inout) :: hc_orig(:)
real, intent(in) :: recon_time_base(:)
logical, intent(in) :: verbose

integer :: i,j,start_index,length,k
integer :: si,ei
integer :: total_windows
real, allocatable :: test2(:),prox(:)
real (kind=r8_kind) :: alt_var
real :: alp
real :: cal_time,start_time,hc_temp,hc_test,pmin,pmax

allocate (test2(max_length))
allocate (prox(max_length))

! find cal epoch for rescaling
do i=1,recon_length
  if (proxy_orig(1,1,1).le.recon_time_base(i)) exit
enddo
cal_index(1)=i

do i=recon_length,1,-1
  if (proxy_orig(series_length_orig(1),1,1).ge.recon_time_base(i)) exit
enddo
cal_index(2)=i

if (verbose) then
  write(11,'(a)')'!correlations'
  print *,'       proxy    correlation    hc'   
endif

!$OMP parallel do schedule(dynamic) default(shared), private(i)
do i=2,2*n_series_orig-1
  var(i)=correlation_full( &
    proxy_orig(1:series_length_orig(1),1,2), &
    proxy_orig(1:series_length_orig(i),i,2), &
    proxy_orig(1:series_length_orig(1),1,1), &
    proxy_orig(1:series_length_orig(i),i,1), &
    series_length_orig(1),series_length_orig(i),hc_optional=hc_orig(i))
enddo
!$OMP end parallel do 

if (verbose) then
  do i=2,2*n_series_orig-1
    write(11,'(a,i5,a,g15.6)')'!target & proxy_orig ',i-1,' : ',var(i)
    print *,i,var(i),hc_orig(i)
  enddo
endif

! test for non linearities
if (non_linear) then
!$OMP parallel do schedule(dynamic) default(shared), private(i,test2,hc_temp,prox,k,alp,hc_test,alt_var,pmin,pmax)
  do i=2*n_series_orig,3*n_series_orig-2
    hc_temp=hc_orig(i)
    var(i)=correlation_full( &
      proxy_orig(1:series_length_orig(1),1,2), &
      proxy_orig(1:series_length_orig(i),i,2), &
      proxy_orig(1:series_length_orig(1),1,1), &
      proxy_orig(1:series_length_orig(i),i,1), &
      series_length_orig(1),series_length_orig(i),hc_optional=hc_orig(i))
    prox(1:series_length_orig(i))=proxy_orig(1:series_length_orig(i),i,2)
    if (verbose) write(11,'(a,i5,a,g15.6)')'!target & proxy_orig ',i-1,' : ',var(i)
    pmin=minval(proxy_orig(1:series_length_orig(i),i,2))
    pmax=maxval(proxy_orig(1:series_length_orig(i),i,2))
    do k=-400,400
      alp=(pmin+pmax)/2d0+(pmax-pmin)*k/100e0   
      test2(1:series_length_orig(i))=(proxy_orig(1:series_length_orig(i),i,2)-alp)**2
      hc_test=hc_temp
      alt_var=correlation_full( &
        proxy_orig(1:series_length_orig(1),1,2), &
        test2(1:series_length_orig(i)), &
        proxy_orig(1:series_length_orig(1),1,1), &
        proxy_orig(1:series_length_orig(i),i,1), &
        series_length_orig(1),series_length_orig(i),hc_optional=hc_test)
      if (abs(alt_var).gt.abs(var(i))) then
        prox(1:series_length_orig(i))=test2(1:series_length_orig(i))
        var(i)=alt_var
        hc_orig(i)=hc_test
        if (verbose) write(11,'(a,i5,a,g15.6,a,g15.6)')'!target & (proxy_orig ',i-1,' -',alp,')^2 : ',alt_var
      endif
    enddo
    do k=-400,400
      alp=(pmin+pmax)/2d0+(pmax-pmin)*k/100e0   
      test2(1:series_length_orig(i))= &
        abs(proxy_orig(1:series_length_orig(i),i,2)-alp)*(proxy_orig(1:series_length_orig(i),i,2)-alp)
      hc_test=hc_temp
      alt_var=correlation_full( &
        proxy_orig(1:series_length_orig(1),1,2), &
        test2(1:series_length_orig(i)), &
        proxy_orig(1:series_length_orig(1),1,1), &
        proxy_orig(1:series_length_orig(i),i,1), &
        series_length_orig(1),series_length_orig(i),hc_optional=hc_test)
      if (abs(alt_var).gt.abs(var(i))) then
        prox(1:series_length_orig(i))=test2(1:series_length_orig(i))
        hc_orig(i)=hc_test
        var(i)=alt_var
        if (verbose) write(11,'(a,i5,a,g15.6,a,i5,a,g15.6,a,g15.6)') &
          '!target & |proxy_orig ',i-1,' -',alp,'|*(proxy_orig ',i-1,' -',alp,') : ',alt_var
      endif
    enddo
    do k=-400,400
      alp=(pmin+pmax)/2d0+(pmax-pmin)*k/100e0   
      test2(1:series_length_orig(i))=(proxy_orig(1:series_length_orig(i),i,2)-alp)**3
      hc_test=hc_temp
      alt_var=correlation_full( &
        proxy_orig(1:series_length_orig(1),1,2), &
        test2(1:series_length_orig(i)), &
        proxy_orig(1:series_length_orig(1),1,1), &
        proxy_orig(1:series_length_orig(i),i,1), &
        series_length_orig(1),series_length_orig(i),hc_optional=hc_test)
      if (abs(alt_var).gt.abs(var(i))) then
        prox(1:series_length_orig(i))=test2(1:series_length_orig(i))
        var(i)=alt_var
        hc_orig(i)=hc_test
        if (verbose) write(11,'(a,i5,a,g15.6,a,g15.6)')'!target & (proxy_orig ',i-1,' -',alp,')^3 : ',alt_var
      endif
    enddo
    do k=-400,400
      alp=(pmin+pmax)/2d0+(pmax-pmin)*k/100e0   
      test2(1:series_length_orig(i))=sqrt(abs(proxy_orig(1:series_length_orig(i),i,2)-alp))
      hc_test=hc_temp
      alt_var=correlation_full( &
        proxy_orig(1:series_length_orig(1),1,2), &
        test2(1:series_length_orig(i)), &
        proxy_orig(1:series_length_orig(1),1,1), &
        proxy_orig(1:series_length_orig(i),i,1), &
        series_length_orig(1),series_length_orig(i),hc_optional=hc_test)
      if (abs(alt_var).gt.abs(var(i))) then
        prox(1:series_length_orig(i))=test2(1:series_length_orig(i))
        hc_orig(i)=hc_test
        var(i)=alt_var
        if (verbose) write(11,'(a,i5,a,g15.6,a,g15.6)') &
          '!target & sqrt(|proxy_orig ',i-1,' -',alp,'|) : ',alt_var
      endif
    enddo
    do k=-400,400
      alp=(pmin+pmax)/2d0+(pmax-pmin)*k/100e0   
      test2(1:series_length_orig(i))= &
        sign(sqrt(abs(proxy_orig(1:series_length_orig(i),i,2)-alp)),proxy_orig(1:series_length_orig(i),i,2)-alp)
      hc_test=hc_temp
      alt_var=correlation_full( &
        proxy_orig(1:series_length_orig(1),1,2), &
        test2(1:series_length_orig(i)), &
        proxy_orig(1:series_length_orig(1),1,1), &
        proxy_orig(1:series_length_orig(i),i,1), &
        series_length_orig(1),series_length_orig(i),hc_optional=hc_test)
      if (abs(alt_var).gt.abs(var(i))) then
        prox(1:series_length_orig(i))=test2(1:series_length_orig(i))
        hc_orig(i)=hc_test
        var(i)=alt_var
        if (verbose) write(11,'(a,i5,a,g15.6,a,i5,a,g15.6,a,g15.6)') &
          '!target & sign(sqrt(|proxy_orig ',i-1,' -',alp,'|),(proxy_orig ',i-1,' -',alp,') : ',alt_var
      endif
    enddo
    do k=-400,400
      alp=(pmin+pmax)/2d0+(pmax-pmin)*k/100e0   
      test2(1:series_length_orig(i))= &
        cube_root(proxy_orig(1:series_length_orig(i),i,2)-alp,series_length_orig(i))
      hc_test=hc_temp
      alt_var=correlation_full( &
        proxy_orig(1:series_length_orig(1),1,2), &
        test2(1:series_length_orig(i)), &
        proxy_orig(1:series_length_orig(1),1,1), &
        proxy_orig(1:series_length_orig(i),i,1), &
        series_length_orig(1),series_length_orig(i),hc_optional=hc_test)
      if (abs(alt_var).gt.abs(var(i))) then
        prox(1:series_length_orig(i))=test2(1:series_length_orig(i))
        hc_orig(i)=hc_test
        var(i)=alt_var
        if (verbose) write(11,'(a,i5,a,g15.6,a,g15.6)') &
          '!target & cube_root(proxy_orig ',i-1,' -',alp,') : ',alt_var
      endif
    enddo
    proxy_orig(1:series_length_orig(i),i,2)=prox(1:series_length_orig(i))
  enddo
!$OMP end parallel do 
endif

deallocate(test2,prox)

end subroutine

end module
