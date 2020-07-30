c=======================================================================
c CVS $Id: global.F,v 1.1 1998/08/28 20:52:55 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/Coupled/global.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE global_init

      implicit none

      !----- global -----
#include <global.h>          ! OUTPUT - global grids, mask, SST
#include <netcdf.inc>        ! netCDF defs

      !----- local -----
      integer       i,j,k,n  ! generic indicies
      character*40  fn       ! local      file name
      character*240 str      ! generic character string
      integer       ier      ! msread/mswrite error code
      integer       fid      ! nc file ID
      integer       vid      ! nc variable ID
      integer       rcode    ! nc return code

      !----- STR data ------
      integer      nx,ny     ! STR grid dimensions
      parameter   (nx=181,ny=91)
      real       x(nx),y(ny) ! STR x & y coordinates 
      real     SST(nx,ny,12) ! STR SST data on STR grid

      SAVE

c-----------------------------------------------------------------------
c PURPOSE:
c   initialize global domain specifications
c   initialize global SST climatology data
c
c ASSUMPTIONS
c   o netCDF input data is a 12 month STR SST climatology on an 2x2 grid
c     with domain = [0,360]x[-90,90]

c-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! get global domain info
      !-----------------------------------------------------------------
      write(*,9000) 'reading global domain info'

      fn='data.domain.nc'
      write(*,9000) 'open netCDF file ',fn
      rcode = nf_open(fn,NF_NOWRITE,fid)

      !--- print title info from nc file ---
      rcode = nf_inq_attlen  (fid, NF_GLOBAL, 'title', n)
      if (n.le.240) then
        rcode = nf_get_att_text(fid, NF_GLOBAL, 'title'  , str)
      else
        str = '<too long to print>   '
        n = index(str,'   ')
      end if
      write(*,9000) 'o title   = ',str(1:n)
      rcode = nf_inq_attlen  (fid, NF_GLOBAL, 'history', n)
      if (n.le.240) then
        rcode = nf_get_att_text(fid, NF_GLOBAL, 'history', str)
      else
        str = '<too long to print>   '
        n = index(str,'   ')
      end if
      write(*,9000) 'o history = ',str(1:n)
      rcode = nf_inq_attlen  (fid, NF_GLOBAL, 'Source', n)
      if (n.le.240) then
        rcode = nf_get_att_text(fid, NF_GLOBAL, 'Source' , str)
      else
        str = '<too long to print>   '
        n = index(str,'   ')
      end if
      write(*,9000) 'o source  = ',str(1:n)

      write(*,9000) 'read domain data...'

      rcode = nf_inq_varid(fid,'xc' ,vid)
      rcode = nf_get_var_double(fid,vid,xc_g)
      rcode = nf_inq_varid(fid,'xe' ,vid)
      rcode = nf_get_var_double(fid,vid,xe_g)

      rcode = nf_inq_varid(fid,'yc' ,vid)
      rcode = nf_get_var_double(fid,vid,yc_g)
      rcode = nf_inq_varid(fid,'ye' ,vid)
      rcode = nf_get_var_double(fid,vid,ye_g)

      rcode = nf_inq_varid(fid,'mask',vid)
      rcode = nf_get_var_int(fid,vid,mask_g)

      rcode = nf_close(fid)

      write(*,9001)'mx,my = ',mx,my
      write(*,9002)'xc_g: ',(xc_g(i),i=1,4),' ...',(xc_g(i),i=mx-2,mx)
      write(*,9002)'xe_g: ',(xe_g(i),i=1,4),' ...',(xe_g(i),i=mx-2,mx+1)
      write(*,9002)'yc_g: ',(yc_g(j),j=1,4),' ...',(yc_g(j),j=my-2,my)
      write(*,9002)'ye_g: ',(ye_g(j),j=1,4),' ...',(ye_g(j),j=my-2,my+1)

      !---------------------------------------------------------------
      ! read in SST data from an existing local file
      !---------------------------------------------------------------
      write(*,9000) 'reading global SST climatology...'

      fn='data.sst.nc'

      write(*,9000) 'open SST netCDF file ',fn
      rcode = nf_open(fn,NF_NOWRITE,fid)

      !--- print title info from nc file ---
      rcode = nf_inq_attlen  (fid, NF_GLOBAL, 'title', n)
      if (n.le.240) then
        rcode = nf_get_att_text(fid, NF_GLOBAL, 'title'  , str)
      else
        str = '<too long to print>   '
        n = index(str,'   ')
      end if
      write(*,9000) 'o title   = ',str(1:n)
      rcode = nf_inq_attlen  (fid, NF_GLOBAL, 'history', n)
      if (n.le.240) then
        rcode = nf_get_att_text(fid, NF_GLOBAL, 'history', str)
      else
        str = '<too long to print>   '
        n = index(str,'   ')
      end if
      write(*,9000) 'o history = ',str(1:n)
      rcode = nf_inq_attlen  (fid, NF_GLOBAL, 'Source', n)
      if (n.le.240) then
        rcode = nf_get_att_text(fid, NF_GLOBAL, 'Source' , str)
      else
        str = '<too long to print>   '
        n = index(str,'   ')
      end if
      write(*,9000) 'o source  = ',str(1:n)
 
      write(*,9000) 'read domain data...'
      rcode = nf_inq_varid(fid,'xc' ,vid)
      rcode = nf_get_var_double(fid,vid,x )
      rcode = nf_inq_varid(fid,'yc' ,vid)
      rcode = nf_get_var_double(fid,vid,y )

      if (x( 1) .ne.   0.0) stop 1
      if (x(nx) .ne. 360.0) stop 2
      if (y( 1) .ne. -90.0) stop 3
      if (y(ny) .ne.  90.0) stop 4

      write(*,9000) 'read SST data...'
      rcode = nf_inq_varid(fid,'T',vid)
      rcode = nf_get_var_double(fid,vid,SST)
 
      rcode = nf_close(fid)

      !---------------------------------------------------------------
      ! map SST data from STR grid to global model grid
      !---------------------------------------------------------------
      write(*,9000) 'mapping SST climatology onto model grid...'

      do k=1,12
        call Interp2D(x   ,y   ,SST  (1,1,k),nx,ny,
     +                xc_g,yc_g,T12_g(1,1,k),mx,my,1.0e30,0)
      end do

      RETURN
 9000 format('(global_init) ',3a)
 9001 format('(global_init) ',a,2i4)
 9002 format('(global_init) ',2(a,4f9.3))
      END

c=======================================================================

      SUBROUTINE global_embed

      implicit none

      !----- global -----
#include <dims.h>            ! domain dimensions
#include <domain.h>          ! domain specifications
#include <control.h>         ! date/time, input data file names, tmpunit
#include <constant.h>        ! physical constants
#include <data.h>            ! INPUT  - regional ocn fields
#include <global.h>          ! OUTPUT - global   ocn fields

      !----- local -----
      integer      i,j       ! generic indicies
      integer      ii,jj     ! generic indicies into global data
      integer      month1    ! two months to interpolate between
      integer      month2
      real         f1, f2    ! weights for the two months data
      integer      eday1     ! elapsed days since y-m-d=0-1-1 & 0s
      integer      eday2     ! elapsed days since y-m-d=0-1-1 & 0s

      SAVE

c-----------------------------------------------------------------------
c PURPOSE:
c   embed Pacific region state variables into the global state variable 
c   fields to be sent to the Flux Coupler.
c
c NOTE:
c   The Coupler expects one state variable field value (or one area
c   average flux field value) for each active grid cell.  The Coupler 
c   does not recognize the Arakawa A scheme, thus this routine must take
c   data from the Pacific Model's Arakawa A grid and create one value 
c   for each grid cell.
c  
c   State variable locations are specified by cell "center coordinates",
c   grid cell areas are specified by cell "edge coordinates", and
c   active grid cells are specified by non-zero "domain mask" values.
c   Grid cell edge coordinates, center coordinates, and domain mask 
c   values are set within the ocean model and sent to the Coupler.  
c
c   Global ocean state variable data (eg. global SST) is derived from
c   global monthly climatology fields.
c  
c-----------------------------------------------------------------------

      if (info_dbug.ge.3) write(*,9000) 'creating global data fields'

      !----------------------------------------------------------------
      ! set global state variable fields
      !----------------------------------------------------------------
      if (day.lt.16) then
         !--- day=1  & sec=0     => f1=0.5 & f2=0.5
         !--- day=15 & sec=86400 => f1=0.0 & f2=1.0
         month1 = month - 1
         month2 = month
         if (month1.eq.0) month1=12
         f2 = 0.5 + 0.5*(float(day-1) + sec/86400.0)/15.0
         f1 = 1.0 - f2
      else
         !--- day=16 & sec=0                => f1=1.0 & f2=0.0
         !--- last day of month & 86400 sec => f1=0.5 & f2=0.5
         month1 = month
         month2 = month + 1
         if (month2.le.12) then
           call ymd2eday(year,month1,16,eday1)
           call ymd2eday(year,month2, 1,eday2)
         else ! month2 = 13
           call ymd2eday(year  ,12,16,eday1)
           call ymd2eday(year+1, 1, 1,eday2)
           month2=1
         end if
         f1 = 1.0 - 0.5*(float(day-16) + sec/86400.0)/(eday2-eday1)
         f2 = 1.0 - f1
      end if

      DO j=1,my
      DO i=1,mx
           T_g(i,j) = f1*T12_g(i,j,month1)+f2*T12_g(i,j,month2) + Tfrz
           U_g(i,j) = 0.0
           V_g(i,j) = 0.0
        dhdx_g(i,j) = 0.0
        dhdy_g(i,j) = 0.0
           q_g(i,j) = 0.0
      END DO
      END DO

      !-----------------------------------------------------------------
      ! embed (or blend) regional model values into global domain
      !-----------------------------------------------------------------
      DO j=1,nyp1
      DO i=it1(j),it2(j)
	 ii = i + i_os
	 jj = j + j_os
         if ( YYY(j) .gt. 20.0) then
           f1 = max(0.0,(30.0-YYY(j))/10.0)
           f2 = 1.0 - f1
         else if ( YYY(j) .lt. -20.0) then
           f1 = max(0.0,(YYY(j)+30.0)/10.0)
           f2 = 1.0 - f1
         else 
           f1 = 1.0
           f2 = 0.0
         end if
         T_g(ii,jj) = f1*(T(i,j,1)+Tfrz) + f2*T_g(ii,jj)
         U_g(ii,jj) = f1* U(i,j,1)       + f2*U_g(ii,jj)
         V_g(ii,jj) = f1* V(i,j,1)       + f2*V_g(ii,jj)
      END DO
      END DO

      RETURN
 9000 format('(global_embed) ',a)
      END

c=======================================================================

      SUBROUTINE global_extract

      implicit none

      !----- global -----
#include <dims.h>            ! domain dimensions
#include <domain.h>          ! domain specifications
#include <control.h>         ! date/time, input data file names, tmpunit
#include <constant.h>        ! physical constants
#include <data.h>            ! INPUT  - regional ocn fields
#include <global.h>          ! OUTPUT - global   ocn fields

      !----- local -----
      integer      i ,j      ! generic indicies into regional data
      integer      ii,jj     ! generic indicies into global data

c-----------------------------------------------------------------------
c PURPOSE:
c   extract Pacific region forcing data from the global forcing fields
c   received from the Flux Coupler.
c  
c ASSUMPTIONS:
c   The Coupler provides one area average flux field value for each 
c   active grid cell.  The Coupler does not recognize the Arakawa A 
c   scheme.  Thus this routine must take area average fluxes, one value
c   for each grid cell, and map these values as necessary onto the 
c   Pacific Model's Arakawa A grid.
c  
c   Grid cell areas are specified by the cell's "edge coordinates",
c   active grid cells are specified by non-zero "domain mask" values.
c   Both the edge coordinates and the domain mask values are set within
c   the ocean model and sent to the Coupler.  
c  
c ALGORITHM:
c   mapping data from the global (Coupler-compatible) grid to the 
c   model's regional (Arakawa A) grid is done by simply choosing a valid
c   global data point that is is adjacent to the desired Arakawa A grid
c   point.
c-----------------------------------------------------------------------

      if (info_dbug.ge.3) write(*,9000) 'extracting regional data'

      DO j=1,nyp1
      DO i=it1(j),it2(j)
         !-----------------------------------------------
         ! find corresponding regional & global indicies
         !-----------------------------------------------
	 ii = i + i_os
	 jj = j + j_os
	 if (mask_g(ii,jj).eq.0) then
	   ii = i + i_os - 1
	   jj = j + j_os - 1
	 end if
	 if (mask_g(ii,jj).eq.0) then
	   ii = i + i_os
	   jj = j + j_os - 1
	 end if
	 if (mask_g(ii,jj).eq.0) then
	   ii = i + i_os - 1
	   jj = j + j_os
	 end if
	 if (mask_g(ii,jj).eq.0) STOP ! mask=0 => out of domain
         !-----------------------------------------------
         ! copy global data into regional array
         !-----------------------------------------------
         TAUx(i,j) =  taux_g(ii,jj) ! zonal      stress (N/m^2)
         TAUy(i,j) =  tauy_g(ii,jj) ! meridional stress (N/m^2)
         Qsol(i,j) = netsw_g(ii,jj) ! net sw heat
         Qsen(i,j) = lslwm_g(ii,jj) ! (lat + sen + net lw + melt) heat
         Qlat(i,j) = 0.0
         Qlgw(i,j) = 0.0
         Qtot(i,j) = Qsol(i,j) + Qsen(i,j)
         SST (i,j) = T_g(ii,jj) - Tfrz ! used for N+S temperature BC's
         !-----------------------------------------------
      END DO
      END DO

      RETURN
 9000 format('(global_extract) ',a)
      END
c=======================================================================
