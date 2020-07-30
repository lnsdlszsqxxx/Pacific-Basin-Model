c=======================================================================
c CVS $Id: global.h,v 1.1 1998/08/28 20:52:55 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/Coupled/global.h,v $
c CVS $Name: pac3_e $
c=======================================================================
c  Sign convention:
c    positive value <=> downward flux
c
c  Unit convention:
c    T             ~ Kelvin
c    u,v           ~ m/s
c    pressure      ~ N/m^2 = Pa
c    heat flux     ~ W/m^2
c    momentum flux ~ N/m^2
c    water flux    ~ (kg/s)/m^2
c-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! global grid resolution
      !-----------------------------------------------------------------
      integer    mx    ,my
      parameter (mx=248,my=160)

      !-----------------------------------------------------------------
      ! regional data (i,j) corresponds to global data (i+i_os,j+j_os)
      !-----------------------------------------------------------------
      integer    i_os   ,j_os
      parameter (i_os=51,j_os=35)

      !-----------------------------------------------------------------
      ! global coordinate arrays and model domain mask
      !-----------------------------------------------------------------
      real       xc_g(mx  )   ! x-coordinates of grid centers
      real       yc_g(my  )   ! y-coordinates of grid centers
      real       xe_g(mx+1)   ! x-coordinates of grid edges
      real       ye_g(my+1)   ! y-coordinates of grid edges
      integer  mask_g(mx,my)  ! domain mask

      !-----------------------------------------------------------------
      ! global monthly mean climatology data
      !-----------------------------------------------------------------
      real    T12_g(mx,my,12) ! SST climatology (degrees K)

      !-----------------------------------------------------------------
      ! output: data sent to coupler
      !-----------------------------------------------------------------

      !----- control flags -----
      integer   nibuff
      parameter(nibuff=100)
      integer   ibuffs(nibuff) ! "info-buffer" sent to coupler

      !----- states -----
      real      T_g(mx,my)     ! SST        (degrees K)
      real      U_g(mx,my)     ! u velocity (m/s)
      real      V_g(mx,my)     ! v velocity (m/s)

      !----- fluxes -----
      real   dhdx_g(mx,my)     ! surface slope, zonal
      real   dhdy_g(mx,my)     ! surface slope, meridional
      real      Q_g(mx,my)     ! heat of fusion    (q > 0), or
                               ! melting potential (q < 0)
      !-----------------------------------------------------------------
      ! input: data received from coupler
      !-----------------------------------------------------------------

      !----- control flags -----
      integer   ibuffr(nibuff) ! "info-buffer" veceived from coupler

      !----- fluxes -----
      real   taux_g(mx,my)     ! zonal      momentum flux
      real   tauy_g(mx,my)     ! meridional momentum flux
      real  netsw_g(mx,my)     ! net shortwave                heat flux
      real  lslwm_g(mx,my)     ! lat + sen + net lw + melting heat flux
      real    pem_g(mx,my)     ! P + E + melting             water flux
      real     runoff(100)     ! runoff                      water flux

      !-----------------------------------------------------------------
      ! make data global
      !-----------------------------------------------------------------
      COMMON /global/ 
     +    T12_g,
     +     xc_g,    yc_g,    xe_g,    ye_g,  mask_g,
     +   ibuffs,     T_g,     U_g,     V_g,  dhdx_g,  dhdy_g,    Q_g,
     +   ibuffr,  taux_g,  tauy_g, netsw_g, lslwm_g,   pem_g, runoff

      SAVE   /global/ 
