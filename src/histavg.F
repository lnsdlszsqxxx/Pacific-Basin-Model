c=======================================================================
c CVS $Id: histavg.F,v 1.6 1998/08/25 20:46:13 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/histavg.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE histavg ()

      implicit none

      !---- global ----
      include 'control.h'      ! global time step & control info

      !----- local -----
      integer      i,j         ! character string index
      character*64 mss_fn      ! mass-store file name
      character*32 loc_fn      ! local      file name (current)
      character*32 old_fn      ! local      file name (previous)
      integer      ier         ! mswrite error code
      integer      nsamples    ! # of samples in t-avg data

      DATA         old_fn      / "null" /
      DATA         nsamples    / -1     /

      SAVE     

c-----------------------------------------------------------------------
c PURPOSE:
c   if hist_tavg = TRUE, then create monthly average history files
c
c NOTES:
c o history file directory, password, retention period, and mswrite 
c   options are the same as for restart files
c o if the model is not run uninterupted for the entire month, no
c   monthly avg data will be created.  nsamples = -1 is used to flag a
c   model restart during the middle of a month.
c-----------------------------------------------------------------------


      if ( hist_tavg.eq.0 ) RETURN

      if ( day.eq.1 .and. sec.eq.0 ) then

        if ( nsamples.gt.0 ) then
          call time_avg("avg")
  
          !oc_fn='123456789+1.nc '
          loc_fn='hyyyy-mm-av.nc '
          if (month.eq.1) then
            write(loc_fn( 2: 5),'(i4.4)') year-1
            write(loc_fn( 7: 8),'(i2.2)') 12
          else
            write(loc_fn( 2: 5),'(i4.4)') year
            write(loc_fn( 7: 8),'(i2.2)') month-1
          endif
          write(*,9002) '(histavg) date, sec, nstep: ',cdate,sec,nstep
          call nc_histavg(loc_fn,ier)

          if (ier.ne.0) then
            write(*,9000) '(histavg) ERROR creating netcdf history file'
          else
            i      = index(mss_dir,' ')-1
            j      = index(loc_fn  ,' ')-1
            mss_fn = mss_dir(1:i)//loc_fn(1:j)//' '
            call MSwrite(ier, loc_fn,mss_fn, mss_pass,mss_rtpd,mss_opts)
            i=index(mss_fn,' ')-1
            j=index(loc_fn,' ')-1
            write(*,9004) ier,loc_fn(1:j),mss_fn(1:i)
            if ( old_fn(1:4) .ne. "null") call systemcall('rm '//old_fn)
            if ( mss_rmlf    .ne. 0     ) old_fn = loc_fn
          end if
        end if
        call time_avg("zero")
        nsamples = 0
      end if

      if ( nsamples.lt.0 .and. day.lt.5 ) then
        call time_avg("zero")
        nsamples = 0
      end if

      if ( sec.eq.0 .and. nsamples.ge.0 ) then
        call time_avg("sum")
        nsamples = nsamples + 1
      end if

      RETURN
 9000 format(2a)
 9002 format(a,i8.8,1x,i5,1x,i7)
 9004 format('(histavg) return code ',i4, ' on  mswrite ',a,' ',a)
      END

c=======================================================================

      SUBROUTINE time_avg (cFlag)

      implicit      none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields
      include 'histavg.h'      ! time-avg data

      !----- arguments -----
      character*(*) cFlag           ! selects zero, sum, or average

      !----- local -----
      real     T2 (nxp1,nyp1,nz2p1) ! T: data on 2*nz+1 sigma grid
      real     T3 (nxp1,nyp1,nzta ) ! T: data on fixed depth grid
      real     U2 (nxp1,nyp1,nz2p1) ! u: data on 2*nz+1 sigma grid
      real     U3 (nxp1,nyp1,nzta ) ! u: data on fixed depth grid
      real     V2 (nxp1,nyp1,nz2p1) ! v: data on 2*nz+1 sigma grid
      real     V3 (nxp1,nyp1,nzta ) ! v: data on fixed depth grid
      integer  ns                   ! number of time samples
      real     fact                 ! 1/(number of time samples)
      integer  i,j,k                ! generic loop indicies
      logical  first_call           ! initialization flag
      data     first_call /.true./

      SAVE

c-----------------------------------------------------------------------
c PURPOSE:
c   form time average data
c-----------------------------------------------------------------------
      
      IF ( first_call ) THEN
         first_call = .false.
         !--- define z-coordinate array ---
         if (nzta.eq.11 .and. z0.eq.0.0 .and. z1.eq.400.0) then
           z( 1) =   0.0
           z( 2) =  33.3
           z( 3) =  66.7
           z( 4) = 100.0
           z( 5) = 133.3
           z( 6) = 166.7
           z( 7) = 200.0
           z( 8) = 233.3
           z( 9) = 266.7
           z(10) = 300.0
           z(11) = 400.0
         else
           call UnifGrid(z0,z1,nzta,z)
         end if
         !--- 3D fields ---
         do k=1,nzta
         do j=1,nyp1
         do i=1,nxp1
            aT(i,j,k) = SPVAL
            aU(i,j,k) = SPVAL
            aV(i,j,k) = SPVAL
         end do
         end do
         end do
         !--- surface fields ---
         do j=1,nyp1
         do i=1,nxp1
            aTAUx(i,j) = SPVAL
            aTAUy(i,j) = SPVAL
            aQtot(i,j) = SPVAL
         end do
         end do
      END IF

      IF (cFlag .eq. 'zero') THEN
         ns=0
         !--- 3D fields ---
         do k=1,nzta
         do j=1,nyp1
         do i=it1(j),it2(j)
            aT(i,j,k) = 0.0
            aU(i,j,k) = 0.0
            aV(i,j,k) = 0.0
         end do
         end do
         end do
         !--- surface fields ---
         do j=1,nyp1
         do i=it1(j),it2(j)
            aTAUx(i,j) = 0.0
            aTAUy(i,j) = 0.0
            aQtot(i,j) = 0.0
         end do
         end do
      ELSE IF (cFlag .eq. 'sum' ) THEN
         ns=ns+1
         !--- re-grid data ---
         call calcZZZ(H,ZZZ,it1,it2,nxp1,nyp1,nz,SPVAL)
         call Expand (T ,T2,'T',Tbot,sigma,it1,it2,nxp1,nyp1,nz,SPVAL)
         call Expand (U ,U2,'U', 0.0,sigma,it1,it2,nxp1,nyp1,nz,SPVAL)
         call Expand (V ,V2,'V', 0.0,sigma,it1,it2,nxp1,nyp1,nz,SPVAL)
         do k=1,nzta
          call sliceXY(T2,T3(1,1,k),z(k),nxp1,nyp1,nz2p1,ZZZ,tbot,spval)
          call sliceXY(U2,U3(1,1,k),z(k),nxp1,nyp1,nz2p1,ZZZ, 0.0,spval)
          call sliceXY(V2,V3(1,1,k),z(k),nxp1,nyp1,nz2p1,ZZZ, 0.0,spval)
         end do
         !--- 3D fields ---
         do k=1,nzta
         do j=1,nyp1
         do i=it1(j),it2(j)
            aU(i,j,k) = aU(i,j,k) + U3(i,j,k)
            aV(i,j,k) = aV(i,j,k) + V3(i,j,k)
            aT(i,j,k) = aT(i,j,k) + T3(i,j,k)
         end do
         end do
         end do
         !--- surface fields ---
         do j=1,nyp1
         do i=it1(j),it2(j)
            aTAUx(i,j) = aTAUx(i,j) + TAUx(i,j)
            aTAUy(i,j) = aTAUy(i,j) + TAUy(i,j)
            aQtot(i,j) = aQtot(i,j) + Qtot(i,j)
         end do
         end do
      ELSE IF (cFlag .eq. 'avg' ) THEN
         write(*,9001) 'forming time-avg using ',ns,' samples'
         if (ns .le. 0) STOP 'time_avg'
         fact=1.0/float(ns)
         !--- 3D fields ---
         do k=1,nzta
         do j=1,nyp1
         do i=it1(j),it2(j)
            aU(i,j,k) = aU(i,j,k) * fact
            aV(i,j,k) = aV(i,j,k) * fact
            aT(i,j,k) = aT(i,j,k) * fact
         end do
         end do
         end do
         !--- surface fields ---
         do j=1,nyp1
         do i=it1(j),it2(j)
            aTAUx(i,j) = aTAUx(i,j) * fact
            aTAUy(i,j) = aTAUy(i,j) * fact
            aQtot(i,j) = aQtot(i,j) * fact
         end do
         end do
      ELSE
         write(*,9000) ' ERROR: invalid cFlag=',cFlag
         STOP 'cFlag'
      END IF

      RETURN
 9000 format('(time_avg) ',3a)
 9001 format('(time_avg) ',a,i4,a)
      END

c=======================================================================

      SUBROUTINE nc_histavg(ncFn,info)

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'control.h'      ! global time step & control info
      include 'netcdf.inc'     ! netCDF lib include file
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields
      include 'histavg.h'      ! time-avg data

      !----- arguments -----
      character*(*) ncFN       ! name of netCDF file to create
      integer       info       ! return error code

      !----- local -----
      character*80 str         ! restart file ID string
      integer      strlen      ! string length
      integer      eday        ! elapsed days since y-m-d = 0-1-1, 0sec
      integer      ncID        ! netCDF file ID
      integer      varid_x     ! variable ID: x coordinate
      integer      varid_y     ! variable ID: y coordinate
      integer      varid_z     ! variable ID: z coordinate
      integer      varid_time  ! variable ID: t coordinate
      integer      dimid_x     ! dimension ID: x coordinate
      integer      dimid_y     ! dimension ID: y coordinate
      integer      dimid_z     ! dimension ID: z coordinate
      integer      dimid_time  ! dimension ID: t coordinate
      integer      vdims(4)    ! vector of dimension IDs
      integer      beg  (4)    ! vector of starting indicies
      integer      end  (4)    ! vector of ending   indicies
      real         time(30)    ! t coordinates (elapsed years)
      character*8  dstr        ! real world date string
      character*8  tstr        ! real world time string


c-----------------------------------------------------------------------
c PURPOSE:
c   create & write netCDF history data file
c
c NOTE:
c   This file can readily be modified to add more data fields to
c   the resulting netCDF history data file.  For example, if you
c   un-comment the lines that are commented-out, you will add the
c   z(i,j,k) coordinate data (layer-depth) to the history file.
c-----------------------------------------------------------------------

      write(*,9000) 'creating file   : ',ncFN
      ncID = ncCRE(ncFn,ncCLOB,info)

      !-----------------------------------------------------------------
      ! create global attributes
      !-----------------------------------------------------------------
      str   = case_desc
      strlen=len(str)
      call ncAPTC(ncID,ncGLOBAL,'title'       ,ncCHAR,strlen,str,info)
      
      call   datetime(dstr,tstr)
      write(str,9003) dstr,tstr
      strlen=index(str,'$')-1
      call ncAPTC(ncID,ncGLOBAL,'history'     ,ncCHAR,strlen,str,info)
     
      str   ='NCAR/CSM Pacific Basin Model$'
      strlen=index(str,'$')-1
      call ncAPTC(ncID,ncGLOBAL,'Conventions' ,ncCHAR,strlen,str,info)
    
      str   ='NCAR/CSM Pacific Basin Model$'
      strlen=index(str,'$')-1
      call ncAPTC(ncID,ncGLOBAL,'Source'      ,ncCHAR,strlen,str,info)

      !-----------------------------------------------------------------
      ! define coordinate data
      !-----------------------------------------------------------------

      if (month.ne.1) call ymd2eday(year  ,month-1,15,eday)
      if (month.eq.1) call ymd2eday(year-1,  12   ,15,eday)
      time(1) = float(eday)/365.0 

      call nc_ddim(ncID,'x',nxp1  ,xxx,'longitude', 'degrees_east',
     &             dimid_x,varid_x,info)
      call nc_ddim(ncID,'y',nyp1  ,yyy,'latitude' , 'degrees_north',
     &             dimid_y,varid_y,info)
      call nc_ddim(ncID,'z',nzta  ,z  ,'depth'     ,'m',
     &             dimid_z,varid_z,info)
      !--------------------------
      dimid_time = ncDDEF(ncID,'time',ncUNLIM,info)
      varid_time = ncVDEF(ncID,'time',ncFLOAT,1,dimid_time,info)
      str   ='years'
      strlen=index(str,' ')-1
      call ncAPTC(ncID,varid_time,'units'    ,ncCHAR,strlen,str,info)
      str   ='years elapsed since 0000-Jan-01 00:00:00 GMT $'
      strlen=index(str,'$')-1
      call ncAPTC(ncID,varid_time,'long_name',ncCHAR,strlen,str,info)
      call ncENDF(ncID,info)

      !-----------------------------------------------------------------
      ! put coordinate data
      !-----------------------------------------------------------------

      call nc_pdim(ncID,dimid_x   ,1,nxp1,xxx ,info)
      call nc_pdim(ncID,dimid_y   ,1,nyp1,yyy ,info)
      call nc_pdim(ncID,dimid_z   ,1,nzta,z   ,info)
      call nc_pdim(ncID,dimid_time,1,1   ,time,info)
      call ncREDF (ncID,info)

      !-----------------------------------------------------------------
      ! define 3d data
      !-----------------------------------------------------------------

      vdims(1) = dimid_x
      vdims(2) = dimid_y
      vdims(3) = dimid_z
      vdims(4) = dimid_time
      call nc_dvar(ncID,'U'   ,4,vdims,
     &    'u velocity'                  ,'m/s'    , 'monthly avg',
     &     spval,info)
      call nc_dvar(ncID,'V'   ,4,vdims,
     &    'v velocity'                  ,'m/s'    , 'monthly avg',
     &     spval,info)
      call nc_dvar(ncID,'T'   ,4,vdims,
     &    'temperature'                 ,'celcius', 'monthly avg',
     &     spval,info)

      vdims(3) = dimid_time
      call nc_dvar(ncID,'Qtot',3,vdims,
     &    'heat flux: net'              ,'W/m^2'  , 'monthly avg',
     &     spval,info)
      call nc_dvar(ncID,'taux',3,vdims,
     &    'surface stress: zonal'       ,'N/m^2'  , 'monthly avg',
     &     spval,info)
      call nc_dvar(ncID,'tauy',3,vdims,
     &    'surface stress: meridional'  ,'N/m^2'  , 'monthly avg',
     &     spval,info)
      call ncENDF (ncID,info)

      !-----------------------------------------------------------------
      ! put 3d/4d data
      !-----------------------------------------------------------------

      beg(1) = 1
      beg(2) = 1
      beg(3) = 1
      beg(4) = 1
      end(1) = nxp1
      end(2) = nyp1
      end(3) = nzta
      end(4) = 1
      call nc_pvar(ncID,'U'    ,beg,end,aU   ,info)
      call nc_pvar(ncID,'V'    ,beg,end,aV   ,info)
      call nc_pvar(ncID,'T'    ,beg,end,aT   ,info)
      end(3) = 1
      call nc_pvar(ncID,'Qtot' ,beg,end,aQtot,info)
      call nc_pvar(ncID,'taux' ,beg,end,ataux,info)
      call nc_pvar(ncID,'tauy' ,beg,end,atauy,info)

      !-----------------------------------------------------------------
      ! done, close netCDF file
      !-----------------------------------------------------------------
      call ncCLOS(ncID,info)

      RETURN
 9000 format('(nc_histavg) ',2a)
 9003 format('Data created ',a,1x,a,' $')
      END

c=======================================================================
c                support routines for plane extraction
c=======================================================================

      SUBROUTINE calcZZZ (  H,ZZZ,iMin,iMax,nx,ny,nz,spval)


      INTEGER   nx,ny,nz, iMin(ny),iMax(ny)
      REAL      H(nx,ny,nz), ZZZ(nx,ny,2*nz+1), spval

      LOGICAL   FirstCall
      DATA      FirstCall / .true. /
      SAVE      FirstCall
c-----------------------------------------------------------------------
c  PURPOSE:
c    compute the z-coordinates (ie. the depth) for all grids point
c  NOTE:
c    sum{n=1 to N of H(i,j,n)} = depth at the BOTTOM of layer N
c    ZZZ(i,j,k) increases (not decreases) with depth
c-----------------------------------------------------------------------

      IF ( FirstCall ) THEN !--- set special values & surface depth ---
        DO k=1,2*nz+1
        DO j=1,ny
          DO i=1,iMin(j)-1
            ZZZ(i,j,k) = spval
          END DO
          DO i=iMax(j)+1,nx
            ZZZ(i,j,k) = spval
          END DO
        END DO
        END DO
        DO j=1,ny
        DO i=iMin(j),iMax(j)
           ZZZ(i,j,1) = 0.0
        END DO
        END DO
        FirstCall = .false.
      END IF

      DO k=1,nz
      DO j=1,ny
      DO i=iMin(j),iMax(j)
         ZZZ(i,j,2*k  )=ZZZ(i,j,2*k-1)+H(i,j,k)*0.5
         ZZZ(i,j,2*k+1)=ZZZ(i,j,2*k-1)+H(i,j,k)
      END DO
      END DO
      END DO

      RETURN
      END

c=======================================================================

      SUBROUTINE Expand(Ain,Aout,cFlag,dpv,sigma,iMin,iMax,nx,ny,nz,spv)


      INTEGER        nx,ny,nz, iMin(ny),iMax(ny)
      REAL           Ain(nx,ny,nz), Aout(nx,ny,2*nz+1)
      REAL           dpv,sigma(nz),spv
      CHARACTER*(*)  cFlag

c-----------------------------------------------------------------------
c  PURPOSE:
c    interpolate a (nx,ny,nz) array onto a (nx,ny,2*nz+1) array
c  NOTE:
c    different rules apply for different variables
c-----------------------------------------------------------------------


      !--- set special values ------------------------------------------
      DO k=1,2*nz+1
      DO j=1,ny
        DO i=1,iMin(j)-1
           Aout(i,j,k) = spv
        END DO
        DO i=iMax(j)+1,nx
           Aout(i,j,k) = spv
        END DO
      END DO
      END DO

      IF (cFlag .eq. 'T') THEN !--- temperature ------------------------
        DO j=1,ny
        DO i=iMin(j),iMax(j)
c         --- surface BC ---
          Aout(i,j,1) =  Ain(i,j,1)
c         --- bottom BC ---
          Aout(i,j,2*nz  ) = Ain(i,j,nz)
          Aout(i,j,2*nz+1) = dpv
        END DO
        END DO
        DO k=1,nz-1
          sk0=sigma(k)
          sk1=sigma(k+1)
          DO j=1,ny
          DO i=iMin(j),iMax(j)
c           --- interior region ---
            Aout(i,j,2*k)  = Ain(i,j,k)
            Aout(i,j,2*k+1)=(Ain(i,j,k)*sk0+Ain(i,j,k+1)*sk1)/(sk0+sk1)
          END DO
          END DO
        END DO
      ELSE IF (cFlag .eq. 'W') THEN !--- We: entrainment ---------------
        DO j=1,ny
        DO i=iMin(j),iMax(j)
c         --- surface BC ---
          Aout(i,j,1) = 0.0
          Aout(i,j,2) = Ain(i,j,1)*0.5
          Aout(i,j,3) = Ain(i,j,1)
c         --- bottom BC ---
          Aout(i,j,2*nz  ) = Ain(i,j,nz-1)*0.5
          Aout(i,j,2*nz+1) = 0.0
        END DO
        END DO
        DO k=2,nz-1
        DO j=1,ny
        DO i=iMin(j),iMax(j)
c         --- interior region ---
          Aout(i,j,2*k  ) = (Ain(i,j,k) + Ain(i,j,k-1))*0.5
          Aout(i,j,2*k+1) =  Ain(i,j,k)
        END DO
        END DO
        END DO
      ELSE IF (cFlag .eq. 'R') THEN !--- Richardson number -------------
        DO j=1,ny
        DO i=iMin(j),iMax(j)
c         --- surface BC ---
          Aout(i,j,1) = Ain(i,j,1)*3.0
          Aout(i,j,2) = Ain(i,j,1)*2.0
          Aout(i,j,3) = Ain(i,j,1)
        END DO
        END DO
        DO k=2,nz
        DO j=1,ny
        DO i=iMin(j),iMax(j)
c         --- bottom BC & interior region ---
          Aout(i,j,2*k  ) = (Ain(i,j,k) + Ain(i,j,k-1))*0.5
          Aout(i,j,2*k+1) =  Ain(i,j,k)
        END DO
        END DO
        END DO
      ELSE !--- default interpolation ----------------------------------
        DO j=1,ny
        DO i=iMin(j),iMax(j)
c         --- surface BC ---
          Aout(i,j,1) =  Ain(i,j,1)
          Aout(i,j,2) =  Ain(i,j,1)
          Aout(i,j,3) = (Ain(i,j,1)+Ain(i,j,2))*0.5
c         --- bottom BC ---
          Aout(i,j,2*nz  ) = Ain(i,j,nz)
          Aout(i,j,2*nz+1) = dpv
        END DO
        END DO
        DO k=2,nz-1
        DO j=1,ny
        DO i=iMin(j),iMax(j)
c         --- interior region ---
          Aout(i,j,2*k  ) =  Ain(i,j,k)
          Aout(i,j,2*k+1) = (Ain(i,j,k)+Ain(i,j,k+1))*0.5
        END DO
        END DO
        END DO
      END IF

      RETURN
      END 

c=======================================================================

      SUBROUTINE sliceXY(Ain,Aout,    z0,      nx,ny,nz,ZZZ,dpval,spval)


      INTEGER   nx,ny,nz
      REAL      Ain(nx,ny,nz),Aout(nx,ny),ZZZ(nx,ny,nz)
      REAL      z0,dpval,spval
      REAL      frac
      INTEGER   i,j,k
c-----------------------------------------------------------------------
c  PURPOSE:
c     o extract (by interpolation) a layer of constant z (an XY-plane)
c  ASSUMPTIONS: 
c     o ZZZ(i,j,k-1) .lt. ZZZ(i,j,k)
c     o a column consists of (all spval's) xor (no spval's)
c-----------------------------------------------------------------------

      DO i=1,nx
      DO j=1,ny
        IF       ( Ain(i,j,1) .eq. spval) THEN
          Aout(i,j)=spval
        ELSE IF (z0 .lt. ZZZ(i,j,1) ) THEN
          Aout(i,j)=0.0
        ELSE IF (z0 .gt. ZZZ(i,j,nz)) THEN
          Aout(i,j)=dpval
        ELSE
          k=2
  10      IF ( z0 .gt. ZZZ(i,j,k) ) THEN
            k=k+1
            GO TO 10
          END IF
          frac      = (z0-ZZZ(i,j,k-1))/(ZZZ(i,j,k)-ZZZ(i,j,k-1))
          Aout(i,j) = (1.-frac)*Ain(i,j,k-1) + frac*Ain(i,j,k)
        END IF
      END DO
      END DO

      RETURN
      END

c=======================================================================

      SUBROUTINE UnifGrid(x0,x1,nx,Xunif)

      implicit  none

      integer   nx
      real      x0,x1,Xunif(nx)

      integer   i
      real      dx

c-----------------------------------------------------------------------
c  PURPOSE:
c    make a uniform grid
c-----------------------------------------------------------------------

      dx=(x1-x0)/float(nx-1)
      Xunif( 1)=x0
      do i=2,nx-1
        Xunif(i)=x0+float(i-1)*dx
      end do
      Xunif(nx)=x1

      RETURN
      END

c=======================================================================
