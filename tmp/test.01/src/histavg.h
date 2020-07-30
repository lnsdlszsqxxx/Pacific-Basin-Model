c=======================================================================
c CVS $Id: histavg.h,v 1.1.1.1 1997/04/17 23:00:00 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/histavg.h,v $
c CVS $Name: pac3_e $
c=======================================================================

      integer    nzta
      parameter  (nzta=11)
      real      z(nzta)             ! fixed-depth coordinate array
      real       z0,z1              ! min/max fixed depth coordinates
      parameter  (z0=0.0,z1=400.0)

      !-----------------------------------------------------------------
      ! time avg data fields on fixed-depth grid
      !-----------------------------------------------------------------
      real  aT   (nxp1,nyp1,nzta)   ! T                          (deg C)
      real  aU   (nxp1,nyp1,nzta)   ! u velocity                 (m/s)
      real  aV   (nxp1,nyp1,nzta)   ! v velocity                 (m/s)
      real  aTAUx(nxp1,nyp1     )   ! momentum flux - zonal      (N/m^2)
      real  aTAUy(nxp1,nyp1     )   ! momentum flux - meridional (N/m^2)
      real  aQtot(nxp1,nyp1     )   ! heat flux - total          (W/m^2)


      COMMON /data_ta/ 
     &  z , aT , aU , aV , aTAUx , aTAUy , aQtot 

      SAVE   /data_ta/ 
