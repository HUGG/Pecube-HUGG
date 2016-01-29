      subroutine calculate_ages (nrec,nsurf,nz,istep,age1,age2,age3,age4,age5,age6,age7,age8, &
                                 ftdist,tprevious,ftlflag,ageflag,iproc,nd)
!__________________________________________________________________________
!This file has been changed to use ketcham program. 
!original version is found in calculate_ages_old.f90 in this same folder
!by VKP
!__________________________________________________________________________

      implicit none

      integer nsurf,nrec,nz,istep,i,ij,irec,k,ftlflag,ageflag(9)
      real*4 age1(nsurf),age2(nsurf),age4(nsurf),age5(nsurf),age3(nsurf) 
      real*4 age6(nsurf),age7(nsurf),age8(nsurf)
      real*4,dimension(:),allocatable::depth,ztime,ztemp,ztime_3,ztemp_3    
      real*4 zjunk,timenow,ftld(17),ftldmean,ftldsd
      double precision fdist(200),alo,oldest_age,final_age,fmean
      double precision ftdist(17,nsurf),tprevious
      integer aftmodel
      integer iproc,nd

      aftmodel=ftlflag

      allocate (ztime(nrec),ztemp(nrec),depth(nrec),ztime_3(nrec),ztemp_3(nrec))
  
      age1=0.
      age2=0.
      age3=0.
      age4=0.
      age5=0.
      age6=0.
      age7=0.
      age8=0.

        do i=1,nsurf
        if (iproc.eq.0.and.nd.eq.0) call screen_counter (i,nsurf,0)
        ij=(i-1)*nz+1
            do irec=1,nrec
            read (100+istep,rec=irec) ztime(irec),(zjunk,k=1,i-1),ztemp(irec), &
                                    (zjunk,k=i+1,nsurf), &
                                    (zjunk,k=1,i-1),depth(irec)
            enddo
        timenow=ztime(nrec)
        ztime=ztime-timenow

        if (aftmodel.eq.1) then
        do irec=1,nrec
            ztime_3(irec)=ztime(nrec-irec+1) !MOD VKP
            if (ztemp(nrec-irec+1).gt.500) ztemp(nrec-irec+1)=500
            ztemp_3(irec)=ztemp(nrec-irec+1) !MOD VKP
        enddo
        alo = 16.0
        final_age = 0
        fmean = 0
        oldest_age  = 0
          if (ageflag(3).eq.1) then
          call ketch_main(nrec,ztime_3,ztemp_3,alo,final_age,oldest_age,fmean,fdist)
          age3(i)=real(final_age,4)
          endif
          if (ageflag(9).eq.1) then
          ftdist(:,i)=fdist(1:17)/100.d0
          else
          ftdist(:,i)=0.
          endif

        else
        if (ageflag(3).eq.1) call Mad_Trax (ztime,ztemp,nrec,1,2, &
                                            age3(i),ftld,ftldmean,ftldsd)
          if (ageflag(9).eq.1) then
          ftdist(:,i)=ftld(1:17)/100.d0
          else
          ftdist(:,i)=0.
          endif
        endif
    
        if (ageflag(4).eq.1) call Mad_Zirc (ztime,ztemp,nrec,0,2, &
                                            age4(i),ftld,ftldmean,ftldsd)

        if (ageflag(1).eq.1) call Mad_He (ztime,ztemp,nrec,age1(i),1)
        if (ageflag(2).eq.1) call Mad_He (ztime,ztemp,nrec,age2(i),2)
        if (ageflag(5).eq.1) call Mad_He (ztime,ztemp,nrec,age5(i),3)
        if (ageflag(6).eq.1) call Mad_He (ztime,ztemp,nrec,age6(i),4)
        if (ageflag(7).eq.1) call Mad_He (ztime,ztemp,nrec,age7(i),5)
        if (ageflag(8).eq.1) call Mad_He (ztime,ztemp,nrec,age8(i),6)
        age1(i)=age1(i)+timenow
        age2(i)=age2(i)+timenow
        age3(i)=age3(i)+timenow
        age4(i)=age4(i)+timenow
        age5(i)=age5(i)+timenow
        age6(i)=age6(i)+timenow
        age7(i)=age7(i)+timenow
        age8(i)=age8(i)+timenow
        enddo

      if (iproc.eq.0.and.nd.eq.0) print*,''

      deallocate (ztime,ztime_3,ztemp,ztemp_3,depth)

      return
      end
