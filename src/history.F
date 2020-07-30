c=======================================================================
c CVS $Id: history.F,v 1.7 1998/08/27 21:29:03 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/history.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE history ()

      implicit none

      !---- global ----
      include 'control.h'      ! global time step & control info

      !----- local -----
      integer      i,j         ! character string index
      character*32 loc_fn      ! local      file name
      character*64 mss_fn      ! mass-store file name
      integer      ier         ! mswrite error code
      logical      file_exists ! true => a file has been created

      DATA         file_exists / .false. /

      SAVE     

c-----------------------------------------------------------------------
c PURPOSE:
c   if hist_now = TRUE, then create history file & save to mass store
c
c NOTE:
c   history file directory, password, retention period, and mswrite 
c   options are the same as for restart files
c-----------------------------------------------------------------------

      IF ( hist_now .eq. 0 ) RETURN

      if(file_exists .and. mss_rmlf.ne.0) call systemcall('rm '//loc_fn)

      !oc_fn='123456789+1234567.nc '
      loc_fn='hyyyy-mm-dd.sssss.nc '
      write(loc_fn( 2: 5),'(i4.4)') year
      write(loc_fn( 7: 8),'(i2.2)') month
      write(loc_fn(10:11),'(i2.2)') day
      write(loc_fn(13:17),'(i5.5)') sec

      write(*,9002) '(history) date, sec, nstep: ',cdate,sec,nstep
      call nc_history(loc_fn,ier)
      file_exists = .true. 

      if (ier.eq.0) then
        i      = index(mss_dir,' ')-1
        j      = index(loc_fn  ,' ')-1
        mss_fn = mss_dir(1:i)//loc_fn(1:j)//' '
        call MSwrite (ier, loc_fn, mss_fn, mss_pass, mss_rtpd, mss_opts)
        i=index(mss_fn,' ')-1
        j=index(loc_fn,' ')-1
        write(*,9004) ier,loc_fn(1:j),mss_fn(1:i)
      else
        write(*,9000) '(history) ERROR creating netcdf history file'
      end if

      RETURN
 9000 format(2a)
 9002 format(a,i8.8,1x,i5,1x,i7)
 9004 format('(history) return code ',i4, ' on  mswrite ',a,' ',a)
      END

c=======================================================================

      SUBROUTINE nc_history(ncFn,info)

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields
      include 'netcdf.inc'     ! netCDF lib include file

      !----- arguments -----
      character*(*) ncFN       ! name of netCDF file to create
      integer       info       ! return error code

      !----- local -----
      integer      i,j,k       ! generic indicies
      character*80 str         ! restart file ID string
      integer      strlen      ! string length
      integer      eday        ! elapsed days since y-m-d = 0-1-1, 0sec
      integer      ncID        ! netCDF file ID
      integer      varid_x     ! variable ID: x coordinate
      integer      varid_y     ! variable ID: y coordinate
      integer      varid_s     ! variable ID: s coordinate
      integer      varid_time  ! variable ID: t coordinate
      integer      dimid_x     ! dimension ID: x coordinate
      integer      dimid_y     ! dimension ID: y coordinate
      integer      dimid_s     ! dimension ID: s coordinate
      integer      dimid_time  ! dimension ID: t coordinate
      integer      vdims(4)    ! vector of dimension IDs
      integer      beg  (4)    ! vector of starting indicies
      integer      end  (4)    ! vector of ending   indicies
      real         sss (nz)    ! s ("sigma") coordinates
      real         time(30)    ! t coordinates (elapsed years)
      real  temp(nxp1,nyp1,nz) ! temporary work array
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

      do k=1,nz
        sss(k)=float(k)
      end do

      call date2eday(cdate,eday)
      time(1) = float(eday)/365.0 + float(sec)/(3600*24*365)

      call nc_ddim(ncID,'x',nxp1  ,xxx,'longitude', 'degrees_east',
     &             dimid_x,varid_x,info)
      call nc_ddim(ncID,'y',nyp1  ,yyy,'latitude' , 'degrees_north',
     &             dimid_y,varid_y,info)
      call nc_ddim(ncID,'s',nz    ,sss,'vertical level' ,'unitless',
     &             dimid_s,varid_s,info)
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
      call nc_pdim(ncID,dimid_s   ,1,nz  ,sss ,info)
      call nc_pdim(ncID,dimid_time,1,1   ,time,info)
      call ncREDF (ncID,info)

      !-----------------------------------------------------------------
      ! define 3d/4d data
      !-----------------------------------------------------------------

      vdims(1) = dimid_x
      vdims(2) = dimid_y
      vdims(3) = dimid_s
      vdims(4) = dimid_time
      call nc_dvar(ncID,'u'   ,4,vdims,'u velocity'          ,'m/s',
     &              'instantaneous',spval,info)
      call nc_dvar(ncID,'v'   ,4,vdims,'v velocity'          ,'m/s',
     &              'instantaneous',spval,info)
      call nc_dvar(ncID,'h'   ,4,vdims,'layer thickness'     ,'m'  ,
     &              'instantaneous',spval,info)
      call nc_dvar(ncID,'T'   ,4,vdims,'temperature'         ,'celcius',
     &              'instantaneous',spval,info)
      vdims(3) = dimid_time
      call nc_dvar(ncID,'zbot',3,vdims,'model depth'         ,'m'    ,
     &              'instantaneous',spval,info)
      call nc_dvar(ncID,'BL'  ,3,vdims,'boundary layer depth','m'    ,
     &              'instantaneous',spval,info)
      call nc_dvar(ncID,'zeta',3,vdims,'surface elevation'   ,'m'    ,
     &              'instantaneous',spval,info)
      call nc_dvar(ncID,'Q'   ,3,vdims,'net heat flux'       ,'w/m^2',
     &              'instantaneous',spval,info)
      call nc_dvar(ncID,'taux',3,vdims,'u surface stress '   ,'N/m'  ,
     &              'instantaneous',spval,info)
      call nc_dvar(ncID,'tauy',3,vdims,'v surface stress '   ,'N/m'  ,
     &              'instantaneous',spval,info)
      call ncENDF (ncID,info)

      !-----------------------------------------------------------------
      ! put 3d/4d data
      !-----------------------------------------------------------------
      
      do k=1,nz
      do j=1,nyp1
      do i=1,nxp1
         temp(i,j,k) = spval
      end do
      end do
      end do

      beg(1) = 1
      beg(2) = 1
      beg(3) = 1
      beg(4) = 1
      end(1) = nxp1
      end(2) = nyp1
      end(3) = nz
      end(4) = 1

      do k=1,nz
      do j=1,nyp1
      do i=it1(j),it2(j)
         temp(i,j,k) = u(i,j,k)
      end do
      end do
      end do
      call nc_pvar(ncID,'u'    ,beg,end,temp,info)

      do k=1,nz
      do j=1,nyp1
      do i=it1(j),it2(j)
         temp(i,j,k) = v(i,j,k)
      end do
      end do
      end do
      call nc_pvar(ncID,'v'    ,beg,end,temp,info)

      do k=1,nz
      do j=1,nyp1
      do i=it1(j),it2(j)
         temp(i,j,k) = h(i,j,k)
      end do
      end do
      end do
      call nc_pvar(ncID,'h'    ,beg,end,temp,info)

      do k=1,nz
      do j=1,nyp1
      do i=it1(j),it2(j)
         temp(i,j,k) = T(i,j,k)
      end do
      end do
      end do
      call nc_pvar(ncID,'T'    ,beg,end,temp,info)

      !------------------------
      !    put 2d variables
      !------------------------
      end(3) = 1

      do j=1,nyp1
      do i=it1(j),it2(j)
         temp(i,j,1) = h(i,j,1)
         do k=2,nz
            temp(i,j,1) = temp(i,j,1) + h(i,j,k)
         end do
      end do
      end do
      call nc_pvar(ncID,'zbot' ,beg,end,temp,info)

      do j=1,nyp1
      do i=it1(j),it2(j)
         temp(i,j,1) = BLdpth(i,j)
      end do
      end do
      call nc_pvar(ncID,'BL'   ,beg,end,temp,info)

      do j=1,nyp1
      do i=it1(j),it2(j)
         temp(i,j,1) = W(i,j,nz)
      end do
      end do
      call nc_pvar(ncID,'zeta' ,beg,end,temp,info)

      do j=1,nyp1
      do i=it1(j),it2(j)
         temp(i,j,1) = Qtot(i,j)
      end do
      end do
      call nc_pvar(ncID,'Q'    ,beg,end,temp,info)

      do j=1,nyp1
      do i=it1(j),it2(j)
         temp(i,j,1) = taux(i,j)
      end do
      end do
      call nc_pvar(ncID,'taux' ,beg,end,temp,info)

      do j=1,nyp1
      do i=it1(j),it2(j)
         temp(i,j,1) = tauy(i,j)
      end do
      end do
      call nc_pvar(ncID,'tauy' ,beg,end,temp,info)

      !-----------------------------------------------------------------
      ! done, close netCDF file
      !-----------------------------------------------------------------
      call ncCLOS(ncID,info)

      RETURN
 9000 format('(nc_history) ',2a)
 9003 format('Data created ',a,1x,a,' $')
      END

c=======================================================================
c     netCDF generic jackets to define & put dimensions & variables
c=======================================================================

      SUBROUTINE nc_ddim(ncID , name , count, data , 
     $                    lname, units, dimid, varid, info)

      implicit none

      !---- global ----
      include 'netcdf.inc'     ! netCDF lib include file

      !--- input ---
      integer ncID                ! netCDF file id
      character*(*) name          ! name of dimension
      integer       count         ! size of dimension
      real          data(*)       ! the data vector
      character*(*) lname         ! long name of dimension
      character*(*) units         ! name of dimension units

      !--- output ---
      integer       dimid         ! dimension id
      integer       varid         ! corresponding variable id
      integer       info          ! error code

c-----------------------------------------------------------------------
c  define a new dimension (coordinate array) in a netCDF file
c-----------------------------------------------------------------------

      !--- define the dimension ---
      dimid = ncDDEF(ncID,name,count,info)
      if (info.ne.0) write(*,*) '(nc_ddim) ncDDEF info=',info

      !--- create the corresponding coordinate variable ---
      varid = ncVDEF(ncID,name,ncFLOAT,1,dimid,info)
      if (info.ne.0) write(*,*) '(nc_ddim) ncVDEF info=',info

      !--- add the mandatory attributes ---
      call ncAPTC(ncID,varid,"long_name",ncCHAR ,len(lname),lname,info)
      call ncAPTC(ncID,varid,"units"    ,ncCHAR ,len(units),units,info)
      call ncAPT (ncID,varid,"minimum"  ,ncFLOAT,  1, data(1)    ,info)
      call ncAPT (ncID,varid,"maximum"  ,ncFLOAT,  1, data(count),info)

      RETURN
      END

c=======================================================================

      SUBROUTINE nc_pdim(ncID,dimid,start,count,data,info)

      implicit none

      !---- global ----
      include 'netcdf.inc'     ! netCDF lib include file

      !--- input ---
      integer       ncID       ! file id
      integer       dimid      ! dimension id
      integer       start      ! starting index
      integer       count      ! number of data elements
      real          data(*)    ! the data vector

      !--- output ---
      integer info             ! error code

c-----------------------------------------------------------------------
c  put dimension (coordinate array) data into a netCDF file
c-----------------------------------------------------------------------

      call ncVPT(ncID,dimid,start,count,data,info)
      if (info.ne.0) write(*,*) '(nc_pdim) ncVPT info=',info

      RETURN
      END

c=======================================================================

      SUBROUTINE nc_dvar(ncID,name,n,dimid,lname,units,trep,spval,info)
      
      implicit none

      !---- global ----
      include 'netcdf.inc'     ! netCDF lib include file

      !--- input ---
      integer       ncID      ! file id
      integer       n         ! number of dimensions
      character*(*) name      ! variable name
      integer       dimid(n)  ! variable dimension ids
      character*(*) lname     ! variable long name
      character*(*) units     ! variable units
      character*(*) trep      ! instantaneous, daily avg, etc
      real          spval     ! missing value, special value

      !--- output ---
      integer info            ! error code

      !--- local ---
      integer varid           ! variable id

c-----------------------------------------------------------------------
c  define a new nD (x,y,...,t) data array in a netCDF file
c-----------------------------------------------------------------------

      !--- define the variable ---
      varid = ncVDEF(ncID,name,ncFLOAT,n,dimid,info)
      if (info.ne.0) write(*,*) '(nc_dvar) ncVDEF info=',info

      !--- add the mandatory attributes ---
      call ncAPTC(ncID,varid,"long_name"  ,ncCHAR,len(lname),lname,info)
      call ncAPTC(ncID,varid,"units"      ,ncCHAR,len(units),units,info)
      call ncAPTC(ncID,varid,"time_rep"   ,ncCHAR,len(trep) ,trep ,info)
      call ncAPT (ncID,varid,"missing_value", ncFLOAT,1     ,spval,info)
      call ncAPT (ncID,varid,"_FillValue"   , ncFLOAT,1     ,spval,info)

      if (info.ne.0) write(*,*) '(nc_dvar) ncAPT info=',info

      RETURN
      END

c=======================================================================

      SUBROUTINE nc_pvar(ncID,name,start,count,data,info)

      implicit none

      !---- global ----
      include 'netcdf.inc'     ! netCDF lib include file

      !--- input ---
      integer       ncID       ! file id
      character*(*) name       ! the variable name
      integer       start(*)   ! starting indicies
      integer       count(*)   ! number of data elements
      real          data (*)   ! the n-d data

      !--- output ---
      integer info

      !--- local ---
      integer varid

c-----------------------------------------------------------------------
c  put n-d (x,y,...,t) data array into a netCDF file
c-----------------------------------------------------------------------

      varid = ncVID(ncID,name,info)
      call ncVPT(ncID,varid,start,count,data,info)
      if (info.ne.0) write(*,*) '(nc_pvar) ncVPT info=',info

      RETURN
      END

c======================================================================
