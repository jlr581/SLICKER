module grad_solve_mod

use correlation_mod
use data_mod
use med_stats_mod

implicit none

logical :: influence_check=.false.
real :: max_time,tol
integer :: number_valid,start_series,max_length
real, parameter :: cn=1e-5

contains

subroutine grad_solve(recon,recon_length,n_series,var,min_cost, &
  orig_med,cal_index,n_members,verbose,time_out) 

integer, intent(in) :: recon_length,n_series
real, intent(inout) :: recon(1:recon_length)
real, intent(in) :: orig_med(:)
real (kind=r8_kind), intent(in) :: var(:)
real (kind=r8_kind), intent(out) :: min_cost
integer, intent(in) :: cal_index(:),n_members
logical, intent(in) :: verbose,time_out

integer :: iter,i,j,loop
real :: grad(1:recon_length),best_recon(1:recon_length)
real :: delta_recon(1:recon_length)
real (kind=r8_kind) :: cost,cost_delta,cost_relax,last_relax,reg_cost
real :: last_recon(1:recon_length)
real :: relax_recon(1:recon_length)
real :: relax,med_qr(2)
integer*8 :: time1,time2,rate,targ_time
logical :: valid

do
  min_cost=huge(min_cost)
  med_qr(:)=alt_med_stats(recon(cal_index(1):cal_index(2)),cal_index(2)-cal_index(1)+1)

  best_recon=(recon-med_qr(1))/med_qr(2)
   
  min_cost=0d0
  do i=start_series,n_series
    min_cost=min_cost+(var(i)-correlation(recon,i,timebase_length(i),recon_length))**2
  enddo

  cost=min_cost

  valid=.false.

  call system_clock(time1,rate)
  targ_time=time1+(max_time)*rate
  iter=0
  do 
    iter=iter+1
  
    if (mod(iter,600).eq.0) then
      if (number_valid.ge.n_members) exit
      if (iter.eq.600) then
        last_recon=best_recon
        last_relax=min_cost
      else
        last_relax=min_cost
        grad=best_recon-last_recon
        last_recon=best_recon
        relax_recon=best_recon+grad
        cost_relax=0d0
        do i=start_series,n_series
          cost_relax=cost_relax+(var(i)-correlation(relax_recon,i,timebase_length(i),recon_length))**2
        enddo
        if (cost_relax.lt.min_cost) then
          min_cost=cost_relax
          med_qr(:)=alt_med_stats(relax_recon(cal_index(1):cal_index(2)),cal_index(2)-cal_index(1)+1)
          best_recon=(relax_recon-med_qr(1))/med_qr(2)
        endif
      endif
    endif

    call random_number(grad)
    do j=1,recon_length
      if (grad(j).lt.0.5e0) then
        grad(j)=-cn
      else
        grad(j)=cn
      endif
    enddo

    delta_recon=best_recon+grad
    cost_delta=0d0
    do i=start_series,n_series
      cost_delta=cost_delta+(var(i)-correlation(delta_recon,i,timebase_length(i),recon_length))**2
    enddo

    grad=(cost_delta-cost)/grad

    relax=1e3*epsilon(relax)/log10(iter+1e0)
    do loop=1,10
      recon=best_recon-min_cost/grad*relax
      cost=0d0
      do i=start_series,n_series
        cost=cost+(var(i)-correlation(recon,i,timebase_length(i),recon_length))**2
      enddo
      if (cost.gt.min_cost) exit
      min_cost=cost
      med_qr(:)=alt_med_stats(recon(cal_index(1):cal_index(2)),cal_index(2)-cal_index(1)+1)
      best_recon=(recon-med_qr(1))/med_qr(2)
      if (min_cost.lt.tol) then
        valid=.true.
        exit
      endif
      relax=relax*3.0e0
    enddo
    if (cost_delta.lt.min_cost) then
      min_cost=cost_delta
      med_qr(:)=alt_med_stats(delta_recon(cal_index(1):cal_index(2)),cal_index(2)-cal_index(1)+1)
      best_recon=(delta_recon-med_qr(1))/med_qr(2)
    endif
    if (min_cost.lt.tol) then
      valid=.true.
      exit
    endif
    call system_clock(time2)
    if (time2.gt.targ_time) exit
  enddo
  if (time_out) valid=.true.
  if (best_recon(1).ne.best_recon(1)) valid=.FALSE.
  if (verbose.or.((number_valid.eq.0).and.(.not.valid))) &
    print *,number_valid,sqrt(min_cost/(n_series-start_series+1)),real(time2-time1)/real(rate),valid
  if (valid) then
    number_valid=number_valid+1
    exit
  endif
  if (number_valid.ge.n_members) exit
  call random_number(recon)
enddo

if (valid) then
  med_qr(1)=M_estimator(best_recon(cal_index(1):cal_index(2)),cal_index(2)-cal_index(1)+1)
  med_qr(2)=Qn(best_recon(cal_index(1):cal_index(2)),cal_index(2)-cal_index(1)+1)
  if (abs(med_qr(2)).gt.1e-7) then
    recon=(best_recon-med_qr(1))*orig_med(2)/med_qr(2)+orig_med(1)
  else
    recon=best_recon-med_qr(1)+orig_med(1)
  endif

else
  recon=huge(recon)
endif

end subroutine

end module
