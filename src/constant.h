c=======================================================================
c CVS $Id: constant.h,v 1.1.1.1 1997/04/17 23:00:00 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/constant.h,v $
c CVS $Name: pac3_e $
c=======================================================================

      real       SPVAL    ! special (undefined) value
      real       PI       ! the famous mathematical constant
      real       DEGtoRAD ! converts degrees to radians ~ rad/deg 
      real       RADtoDEG ! converts radians to degrees ~ deg/rad 
      real       DAYtoSEC ! converts days    to seconds ~ sec/day 
      real       Omega2   ! earth's speed of rotation   ~ rad/sec 
      real       rEarth   ! radius of the earth         ~ m 
      real       gravity  ! gravitational acceleration  ~ m/s^2 
      real       alphaT   ! thermal expansion coeff     ~ 1/[deg C]     
      real       RHOw     ! density of water            ~ kg/m^3 
      real       RHOa     ! density of air              ~ kg/m^3 
      real       Cd       ! momentum drag coefficient   ~ unitless 
      real       Ce       ! latent   drag coefficient   ~ unitless 
      real       Cpw      ! heat capacity of water      ~ J/([deg C]*kg)
      real       Cpa      ! heat capacity of air        ~ J/([deg C]*kg)
      real       Tbot     ! deep ocean temperature      ~ [deg C]
      real       Tfrz     ! freezing temp of water      ~ [deg K]
      real       WNDmin   ! minimum wind speed for evap ~ m/s 

      parameter (SPVAL    = 1.0e30         )
      parameter (PI       = 3.14159265359  )
      parameter (DEGtoRAD = PI/180.0       )
      parameter (RADtoDEG = 180.0/PI       )
      parameter (DAYtoSEC = 86400.0        )
      parameter (Omega2   = 4.0*PI/DAYtoSEC)
      parameter (rEarth   = 6.378e+6       )
      parameter (gravity  = 9.8            )
      parameter (alphaT   =  0.255e-3      )
      parameter (RHOw     =  1.030e+3      )
      parameter (RHOa     =  1.2           )
      parameter (Cd       =  1.420e-3      )
      parameter (Ce       =  1.200e-3      )
      parameter (Cpw      =  4.00e+3       )
      parameter (Cpa      =  1.10e+3       )
      parameter (Tbot     = 10.0           )
      parameter (Tfrz     = 273.15         )
      parameter (WNDmin   =  4.0           )

