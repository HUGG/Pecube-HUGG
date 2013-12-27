!%%%%%%%%%%%%%%%%%%%%%%%%
subroutine calc_raman (temperature,ntime,raman)

     implicit none
     
     integer ntime
     real temperature(ntime),raman
     
     ! raman is the peak temperature
     raman=maxval(temperature)

end subroutine calc_raman
!%%%%%%%%%%%%%%%%%%%