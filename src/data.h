c=======================================================================
c CVS $Id: data.h,v 1.1.1.1 1997/04/17 23:00:00 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/data.h,v $
c CVS $Name: pac3_e $
c=======================================================================

      !-----------------------------------------------------------------
      ! surface forcing fields
      !-----------------------------------------------------------------
      real  WNDs(nxp1,nyp1)        ! wind speed
      real  TAUx(nxp1,nyp1)        ! momentum flux - zonal
      real  TAUy(nxp1,nyp1)        ! momentum flux - meridional
      real  Qsol(nxp1,nyp1)        ! heat flux - solar
      real  Qlat(nxp1,nyp1)        ! heat flux - latent
      real  Qlgw(nxp1,nyp1)        ! heat flux - long wave
      real  Qsen(nxp1,nyp1)        ! heat flux - sensible
      real  Qtot(nxp1,nyp1)        ! heat flux - total
      real  F   (     nyp1)        ! Coriolus force
      real  SST (nxp1,nyp1)        ! SST forcing (obs for fixed BC's)

      !-----------------------------------------------------------------
      ! vertical grid & vertical mixing information         
      !-----------------------------------------------------------------
      real    ZZZ(nxp1,nyp1,nz2p1) ! z coordinate (depth)
      real    Hin(nz)              ! avg layer thickness: <H>
      real  sigma(nz)              ! <H(k)> / (sum <H(k)> for k=2,nz)
      real    Tin(nz)              ! T BC vertical profile 
      real  NU   (nz)              ! background mixing coeff (velocity)
      real  KAPPA(nz)              ! background mixing coeff (heat)
      real cNU   (nz)              !    nu*dk/<dH>
      real cKAPPA(nz)              ! kappa*dk/<dH>
      real BLdpth(nxp1,nyp1   )    ! boundary layer depth
      real rNU   (nxp1,nyp1,nz)    ! kpp mixing coeff (veloctity)
      real rKAPPA(nxp1,nyp1,nz)    ! kpp mixing coeff (heat)
      real rGAMMA(nxp1,nyp1,nz)    ! kpp nonlocal transport

      !-----------------------------------------------------------------
      ! prognostic, n-cycle accumulation, & tendancy arrays
      !-----------------------------------------------------------------
      real  U  (nxp1,nyp1,nz)     ! u velocity      (m/s)
      real  V  (nxp1,nyp1,nz)     ! v velocity      (m/s)
      real  T  (nxp1,nyp1,nz)     ! temperature     (degrees C)
      real  H  (nxp1,nyp1,nz)     ! layer thickness (m)
      real HU  (nxp1,nyp1,nz)     ! H*u
      real HV  (nxp1,nyp1,nz)     ! H*v
      real HT  (nxp1,nyp1,nz)     ! H*T
      real FU  (nxp1,nyp1,nz)     ! n-cycle accumulation array for HU
      real FV  (nxp1,nyp1,nz)     ! n-cycle accumulation array for HV
      real FT  (nxp1,nyp1,nz)     ! n-cycle accumulation array for HT
      real FH1 (nxp1,nyp1   )     ! n-cycle accumulation array for H(1)
      real FHnz(nxp1,nyp1   )     ! n-cycle accumulation array for H(nz)
      real dH1 (nxp1,nyp1   )     ! dH( 1)/dt
      real dHnz(nxp1,nyp1   )     ! dH(nz)/dt
      real W   (nxp1,nyp1,nz)     ! vertical velocity
      real Herr(          nz)     ! depth error at last correction

      COMMON /data/ 
     &  WNDs, TAUx, TAUy, Qsol, Qlat, Qlgw, Qsen, Qtot, F, SST,
     &  ZZZ ,Hin ,sigma ,Tin,
     &  NU, KAPPA, cNU ,cKAPPA ,BLdpth ,rNU ,rKAPPA ,rGAMMA,
     &  U,V,T,H  ,HU,HV,HT , FU,FV,FT , FH1,FHnz, dH1,dHnz, W, Herr

      SAVE   /data/ 
