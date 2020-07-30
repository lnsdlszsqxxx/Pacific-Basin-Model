c=======================================================================
c CVS $Id: domain.h,v 1.1.1.1 1997/04/17 23:00:00 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/domain.h,v $
c CVS $Name: pac3_e $
c=======================================================================

      !-----------------------------------------------------------------
      ! horizontal grid information         
      !-----------------------------------------------------------------
      real     XXX (nxp1     )     ! X coordinates (degrees east)
      real     YYY (     nyp1)     ! Y coordinates (degrees north)
      real     EMx (nxp1,nyp1)     ! dX/di
      real     EMy (nxp1,nyp1)     ! dY/dj
      real     EMxy(nxp1,nyp1)     ! geometric term, d^Y/(dj*di) != 0
      real     area(nxp1,nyp1)     ! area of grid cell = dX/di * dY/dj
      real     basin               ! area of entire domain
      integer  it1(nyp1),it2(nyp1) ! longitude domain for each latitude
      integer  is1(nyp1),is2(nyp1) ! longitude domain for v filtering
      integer  jt1(nxp1),jt2(nxp1) ! latitude domain for each longitude
      integer  js1(nxp1),js2(nxp1) ! latitude domain for u filtering

      COMMON /grid/ 
     +  XXX,YYY,EMx,EMy,EMxy,area,basin, it1,it2,jt1,jt2,is1,is2,js1,js2

      SAVE   /grid/ 

