c=======================================================================
c CVS $Id: domain.F,v 1.2 1998/08/24 22:45:46 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/domain.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE domain

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !---- local ----
      integer i,j                ! generic loop indicies

c-----------------------------------------------------------------------
c PURPOSE:
c   o verify valid domain description
c   o compute is1,is2,js1,js2 from it1,it2,jt1,jt2
c   o calculate da of each grid cell & whole basin
c   o write coriolis force & grid info to stdout
c
c NOTE:
c   perhaps only jt1 & jt2 should NOT be on the restart data file, 
c   but rather they should be derived here from it1 & it2
c
c ASSUMPTIONS:
c   o west/east coastlines, it1 & it2, are monotonic
c   o X(i) & Y(j) grids are monotonic
c   o EMx(i,j) & EMy(i,j) are di/dx & dj/dy respectively (units are 1/m)
c   o F(j) is the Coriolis force
c-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! verify valid domain description
      !-----------------------------------------------------------------
      write(*,'(a)') '(domain) checking domain info'

      do j=1,ny
      do i=1,nx
        if (EMx(i,j).le.0.0) stop 'EMx'
        if (EMy(i,j).le.0.0) stop 'EMy'
      end do
      end do

      do i=1,nx
        if (XXX(i+1).le.XXX(i)) stop 'dx/di'
      end do
      do j=1,ny
        if (YYY(j+1).le.YYY(j)) stop 'dy/dj'
      end do

      if (it1(1).ne.1    .and. it1(nyp1).ne.1   ) stop 'it1-min'
      if (it2(1).ne.nxp1 .and. it2(nyp1).ne.nxp1) stop 'it2-max'

      !-----------------------------------------------------------------
      ! compute  jt1,jt2 from it1,it2
      !-----------------------------------------------------------------
      ! for now we'll just get it off the restart data file

      !-----------------------------------------------------------------
      ! compute is1,is2 & js1,js2 from it1,it2 & jt1,jt2
      !-----------------------------------------------------------------
      write(*,'(a)') '(domain) computing more domain info'
      is1(1)=it1(1)
      is2(1)=it2(1)
      do j=2,ny
        is1(j)=max0(it1(j-1),it1(j),it1(j+1))
        is2(j)=min0(it2(j-1),it2(j),it2(j+1))
      end do
      is1(nyp1)=it1(nyp1)
      is2(nyp1)=it2(nyp1)

      js1(1)=jt1(1)
      js2(1)=jt2(1)
      do i=2,nx
        js1(i)=max0(jt1(i-1),jt1(i),jt1(i+1))
        js2(i)=min0(jt2(i-1),jt2(i),jt2(i+1))
      end do
      js1(nxp1)=jt1(nxp1)
      js2(nxp1)=jt2(nxp1)

      !-----------------------------------------------------------------
      ! calculate da of each grid cell & whole basin
      !-----------------------------------------------------------------
      call calc_da(nxp1,nyp1,it1,it2,EMx,EMy,area)

      basin = 0.0
      do j=1,nyp1
      do i=it1(j),it2(j)
        basin = basin + area(i,j)
      end do
      end do

      !-----------------------------------------------------------------
      ! write coriolis force & grid info to stdout
      !-----------------------------------------------------------------
      IF (info_dbug.ge.1) THEN
         write(*,9001) 
         write(*,9901) basin*1.0e-6
         write(*,9903) '  F: ',(  F(j),j=1,nyp1)
         write(*,9904) '  X: ',(XXX(i),i=1,nxp1)
         write(*,9904) '  Y: ',(YYY(j),j=1,nyp1)
         write(*,9001) 
      END IF

      RETURN
 9001 format(80('-'))
 9901 format('(domain) total basin area = ',e14.6,' km**2')
 9903 format(a5,10e12.4 / (5x,10e12.4))
 9904 format(a5,10f12.4 / (5x,10f12.4))
      END

c=======================================================================

      SUBROUTINE calc_da(nx,ny,imin,imax,EMx,EMy,da)

      implicit none

      integer  nx,ny                 ! domain resolution
      integer  imin(ny),imax(ny)     ! west/east coastline
      real     EMx(nx,ny),EMy(nx,ny) ! di/dx & dj/dy
      real     da(nx,ny)             ! area associated with a data point

      integer  i,j              ! generic loop indicies
      integer  iStart,iEnd      ! given j, min/max i for interior region

c-----------------------------------------------------------------------
c  PURPOSE:
c    calculate da(i,j), the area associated with each data point
c  INPUT:
c    o nx,ny, EMx,EMy, imin,imax
c  OUTPUT:
c    o da
c  ASSUMPTIONS:
c    o EMx & EMy are on the grid cell corners of an Arakawa A grid
c    o imin(j) & imax(j) are monotonic
c-----------------------------------------------------------------------

      !------ southern boundary ------
      da(imin(1),1) = 0.25/(EMx(imin(1),1)*EMy(imin(1),1))
      da(imax(1),1) = 0.25/(EMx(imax(1),1)*EMy(imax(1),1))
      do i = imin(1)+1,imax(1)-1
        da(i,1) = 0.5/(EMx(i,1)*EMy(i,1))
      end do
 
      !------ northern boundary ------
      da(imin(ny),ny) = 0.25/(EMx(imin(ny),ny)*EMy(imin(ny),ny))
      da(imax(ny),ny) = 0.25/(EMx(imax(ny),ny)*EMy(imax(ny),ny))
      do i = imin(ny)+1,imax(ny)-1
        da(i,ny) = 0.5/(EMx(i,ny)*EMy(i,ny))
      end do

      do j=2,ny-1
        !------ western boundary ------
        if (imin(j) .LT. imin(j-1)) THEN
           da(imin(j  ),j) = 0.25/(EMx(imin(j  ),j)*EMy(imin(j  ),j))
           da(imin(j-1),j) = 0.75/(EMx(imin(j-1),j)*EMy(imin(j-1),j))
           do i = imin(j)+1,imin(j-1)-1
             da(i,j) = 0.5/(EMx(i,j)*EMy(i,j))
           end do
           iStart = imin(j-1)+1
        else if (imin(j) .LT. imin(j+1)) THEN
           da(imin(j  ),j) = 0.25/(EMx(imin(j  ),j)*EMy(imin(j  ),j))
           da(imin(j+1),j) = 0.75/(EMx(imin(j+1),j)*EMy(imin(j+1),j))
           do i = imin(j)+1,imin(j+1)-1
             da(i,j) = 0.5/(EMx(i,j)*EMy(i,j))
           end do
           iStart = imin(j+1)+1
        else
           da(imin(j),j) = 0.5/(EMx(imin(j),j)*EMy(imin(j),j))
           iStart = imin(j)+1
        end if

        !------ eastern boundary ------
        if (imax(j) .GT. imax(j-1)) THEN
           da(imax(j  ),j) = 0.25/(EMx(imax(j  ),j)*EMy(imax(j  ),j))
           da(imax(j-1),j) = 0.75/(EMx(imax(j-1),j)*EMy(imax(j-1),j))
           do i = imax(j-1)+1,imax(j)-1
             da(i,j) = 0.5/(EMx(i,j)*EMy(i,j))
           end do
           iend = imax(j-1)-1
        else if (imax(j) .GT. imax(j+1)) THEN
           da(imax(j  ),j) = 0.25/(EMx(imax(j  ),j)*EMy(imax(j  ),j))
           da(imax(j+1),j) = 0.75/(EMx(imax(j+1),j)*EMy(imax(j+1),j))
           do i = imax(j+1)+1,imax(j)-1
              da(i,j) = 0.5/(EMx(i,j)*EMy(i,j))
           end do
           iend = imax(j+1)-1
        else
           da(imax(j),j) = 0.5/(EMx(imax(j),j)*EMy(imax(j),j))
           iend = imax(j)-1
        end if

        !------ interior region ------
        do i = iStart,iend
          da(i,j) = 1.0/(EMx(i,j)*EMy(i,j))
        end do
      END DO
 
      RETURN
      END
