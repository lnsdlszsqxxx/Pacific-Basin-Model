c=======================================================================
c CVS $Id: kpp_coef.F,v 1.3 1999/08/12 17:30:46 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/kpp_coef.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE kpp_coef( U    , V      , T      , H  
     &                   , taux , tauy   , Qsol   , Qtot 
     &                   , imt  , jmt    , km     , kmp1 , it1  , it2 
     &                   , rhow , cpw    , alphaT , Tbot , gravity
     &                   ,  nu  ,  kappa   
     &                   , rNu  , rKappa , rGamma , hbl_2D )
 
c----------------------------------------------------------------------
c PURPOSE
c     This subroutine is the interface between a general ocean model 
c     and the kpp vertical mixing scheme
c
c HISTORY
c   Written by 
c     bill large   6 jun 94
c   Modified by 
c     jan morzel  30 jun 94
c     bill large  11 aug 94
c     bill large     dec 95
c     Nancy Norton   Feb 95 - interface to Gent Pacific model
c     B. Kauffman    Nov 96 - input stress & constants (not wind speed)
c----------------------------------------------------------------------

c input
      real U    (imt,jmt,km) ! zonal velocity                  (m/s)
      real V    (imt,jmt,km) ! meridional velocity             (m/s)
      real H    (imt,jmt,km) ! layer thicknesses               (m)
      real T    (imt,jmt,km) ! temperature                     (deg. C)
      real TAUx (imt,jmt   ) ! momentum flux, zonal            (N/m^2)
      real TAUy (imt,jmt   ) ! momentum flux, meridional       (N/m^2)
      real Qsol (imt,jmt   ) ! solar heat flux                 (w/m^2)
      real Qtot (imt,jmt   ) ! total heat flux                 (w/m^2)
      real nu           (km) ! vertical viscosity   profile
      real kappa        (km) ! vertical diffusivity profile
      real rhow              ! density of water         ~  kg/m^3
      real cpw               ! heat capacity of water   ~  J/([deg C]*kg
      real alphaT            ! thermal expansion coeff  ~  1/[deg C]
      real Tbot              ! deep ocean temperature   ~  deg C
      real gravity           ! gravitational constant   ~  m/s^2
 
      integer it1  (jmt)     ! set of Westernmost gridpoints
      integer it2  (jmt)     ! set of Easternmost gridpoints
 
c output
      real rNu   (imt,jmt,km)! vertical viscosity coefficient  (m^2/s)
      real rKappa(imt,jmt,km)! vertical scalar diffusivity     (m^2/s)
      real rGamma(imt,jmt,km)! nonlocal transport              (s/m^2)
      real hbl_2D(imt,jmt)   ! boundary-layer depth (2D)       (m)
 
c constants
      parameter (mdiff = 3)   ! number of diffusivities for local arrays
      parameter (jwtype_const = 3) ! Jerlov water type
      parameter (epsln = 1.0E-20)  ! a "small" positive number

c local
      real Ritop(imt,km)     ! numerator of bulk Ri number     (m/s)^2
*          Ritop = (-z -(-zref))* delta bouyancy wrt sfc       (m/s^2)
      real dVsq (imt,kmp1)   ! (velocity shear re sfc)^2       (m/s)^2
      real Shsq (imt,kmp1)   ! (local velocity shear)^2        (m/s)^2
      real ustar(imt)        ! surface friction velocity       (m/s)
      real Bo   (imt)        ! surface turbulent buoy. forcing (m^2/s^3)
      real Bosol(imt)        ! radiative buoyancy forcing      (m^2/s^3)
      real dbloc(imt,km)     ! local delta buoyancy            (m/s^2)
      real dbsfc(imt,km)     ! delta buoyancy wrt sfc          (m/s^2)
      real Coriol            ! Coriolis parameter              (s^{-1})
      real hwide(imt,kmp1)   ! layer thicknesses   (row)       (m)
      real hbl(imt)          ! boundary layer depth (1D)       (m)
      real bfsfc(imt)        ! surface buoyancy forcing        (m^2/s^3)
      real ws(imt)           ! momentum velocity scale
      real wm(imt)           ! scalar   velocity scale 
      real caseA(imt)        ! = 1 in case A; =0 in case B
      real stable(imt)       ! = 1 in stable forcing; =0 in unstable
      real dkm1(imt,mdiff)   ! boundary layer difs at kbl-1 level
      real gat1(imt,mdiff)   ! shape function at sigma=1
      real dat1(imt,mdiff)   ! derivative of shape function at sigma=1
      real blmc(imt,km,mdiff)! boundary layer mixing coefficients
      real sigma(imt)        ! normalized depth (d / hbl)
      real Rib(imt,2)        ! bulk Richardson number
      real zgrid(imt,kmp1)   ! vertical grid (<= 0)            (m)
      real visc(imt,0:kmp1)  ! vertical viscosity coefficient  (m^2/s)
      real difs(imt,0:kmp1)  ! vertical scalar diffusivity     (m^2/s)
      real dift(imt,0:kmp1)  ! vertical temperature diffusivity(m^2/s)
      real ghats(imt,km)     ! nonlocal transport              (s/m^2)

      integer jwtype(imt)    ! Jerlov water type               (1 -- 5)
      integer kbl(imt)       ! index of first grid level below hbl

c-----------------------------------------------------------------------

c     Define constants:
 
      alphag   = alphaT*gravity
      alphagQ  = alphag/(Cpw*RHOw)
      Coriol   = 0.0         ! * 
      half_eps = 0.5*epsln
      j1       = 1
      j2       = jmt  
 
*     Note: references to Coriol have been removed from Large & Morzel
c           subroutines
 
 
      do i = 1, imt  
           jwtype(i)     = jwtype_const
      enddo ! i = 1,imt   
      
      do 1000 j = j1, j2
        i1 = it1(j)
        i2 = it2(j)
     
*       Initialize zgrid
 
 
        do k = 1, kmp1
 
         if (k .eq. 1) then
 
           do i = i1, i2
              zgrid (i,1)   = -0.5*H(i,j,1)
           enddo ! i = i1, i2
 
         else if (k .gt. 1  .and. k .lt. kmp1) then
 
           do i = i1, i2
           zgrid (i,k)   = zgrid(i,k-1) -0.5*( H(i,j,k-1) + H(i,j,k) )
           enddo ! i = i1, i2
 
         else if (k .eq. kmp1) then
 
           do i = i1, i2
              zgrid (i,kmp1) = zgrid(i,km) -0.5*H(i,j,km) - half_eps
           enddo ! i = i1, i2
         endif !  k                   
 
        enddo ! k = 1, kmp1
 
*       Initialize rows of local arrays
 
        do k = 1, km
          if (k .le. km-1) then
            do i = i1, i2
               dVsq (i,k) = (U(i,j,k)-U(i,j,1)  )**2 + 
     1                      (V(i,j,k)-V(i,j,1)  )**2
               Shsq (i,k) = (U(i,j,k)-U(i,j,k+1))**2 +
     1                      (V(i,j,k)-V(i,j,k+1))**2  
               hwide(i,k) = H(i,j,k)
               dbloc(i,k) = alphag*(T(i,j,k) - T(i,j,k+1) )
               dbsfc(i,k) = alphag*(T(i,j,1) - T(i,j,k)   )
               Ritop(i,k) = (-zgrid(i,k) + 0.5*zgrid(i,1))*dbsfc(i,k)
            enddo ! i = i1, i2
          else ! k .eq. km
            do i = i1, i2
               dVsq (i,k) = (U(i,j,k)-U(i,j,1)  )**2 + 
     1                      (V(i,j,k)-V(i,j,1)  )**2
               Shsq (i,k) = (U(i,j,k)           )**2 +
     1                      (V(i,j,k)           )**2  
               hwide(i,k) = H(i,j,k)
               dbloc(i,k) = alphag*(T(i,j,k) - Tbot )
               dbsfc(i,k) = alphag*(T(i,j,1) - T(i,j,k)   )
               Ritop(i,k) = (-zgrid(i,k) + 0.5*zgrid(i,1))*dbsfc(i,k)
            enddo ! i = i1, i2
          endif ! (k .le. km-1)
        enddo ! k = 1, km
 
        do i = i1, i2
           hwide(i,kmp1) = epsln
               
           Bo    (i)     = alphagQ*( Qtot(i,j) - Qsol(i,j) )
           Bosol (i)     = alphagQ*  Qsol(i,j)
           ustar (i) = sqrt( sqrt( TAUx(i,j)**2 + TAUy(i,j)**2 )/RHOw )

 
        enddo ! i = i1, i2
 
* zero out the mixing coefficients
 
      do k=0,km
        do i = i1, i2
           visc (i,k) =  0.0
           difs (i,k) =  0.0
           dift (i,k) =  0.0
        enddo ! i = i1, i2
      enddo ! k = 1,km
 
 
c compute RI and IW interior diffusivities everywhere
 
      call ri_iwmix ( km    , kmp1  , imt   , i1    , i2  ,
     1                Shsq  , dbloc , zgrid ,
     1                visc  , difs  , dift  , nu, kappa )
 
c double diffusion not necessary in Gent and Kauffman
c
c fill  kmp1 coefficients for blmix
c
      do i = i1, i2
         visc (i,kmp1) =  visc(i,km)
         difs (i,kmp1) =  difs(i,km)
         dift (i,kmp1) =  dift(i,km)
      enddo ! i
 
c diagnose the new boundary layer depth
 
 
      call  bldepth (km   , kmp1 , imt   ,i1,i2 , zgrid, hwide, 
     1               dVsq , 
     1               dbloc, Ritop, ustar , Bo   , Bosol, Coriol, jwtype,
     1               hbl  , bfsfc, stable, caseA, kbl  ,
     1               Rib  , sigma, wm    , ws   )
 
 
c
c compute boundary layer diffusivities
c
      call blmix   (km   , kmp1 , imt , i1,i2 , mdiff , zgrid, hwide , 
     $              ustar, bfsfc, hbl , stable, caseA, 
     $              visc , difs , dift, kbl   , 
     $              gat1 , dat1 , dkm1, blmc  , ghats, 
     $              sigma, wm   , ws  )
c
c enhance diffusivity at interface kbl - 1
c
      call enhance (km   , kmp1 , imt , i1, i2, mdiff , dkm1  , visc ,
     $              difs , dift , hbl , kbl   , zgrid , caseA,
     $              blmc , ghats )
c
c combine interior and boundary layer coefficients and nonlocal term
c
      do 200 ki= 1,km
         do 190 i = i1, i2
            if(ki.lt.kbl(i)) then
               rNu   (i,j,ki)=blmc(i,ki,1)
               rKappa(i,j,ki)=blmc(i,ki,2)
             rGamma(i,j,ki)=ghats(i,ki)*(Qtot(i,j)-Qsol(i,j))/(Cpw*RHOw)
            else
               rNu   (i,j,ki)=visc(i,ki)
               rKappa(i,j,ki)=difs(i,ki)
               rGamma(i,j,ki)=0.
            endif
 190     continue   
 200  continue

c fill in 2-D boundary-layer depth
 
      do i = i1, i2
         hbl_2D(i,j) = hbl(i)
      enddo ! i = i1, i2
 
 
 
1000  continue            
 
 
      return
      end
       
c=======================================================================

      SUBROUTINE bldepth (
     $          km   , kmp1 , imt   , i1   , i2   ,
     $          zgrid, hwide, dVsq  ,
     $          dbloc, Ritop, ustar , Bo   , Bosol, Coriol, jwtype,
     $          hbl  , bfsfc, stable, caseA, kbl  ,
     $          Rib  , sigma, wm    , ws   )

c-----------------------------------------------------------------------
c     the oceanic planetray boundary layer depth, hbl, is determined as
c     the shallowest depth where the bulk richardson number is
c     equal to the critical value, Ricr.
c
c     bulk richardson numbers are evaluated by computing velocity and
c     buoyancy differences between values at zgrid(kl) < 0 and surface
c     reference values.
c     in this configuration, the reference values are equal to the
c     values in the surface layer.  
c     when using a very fine vertical grid, these values should be 
c     computed as the vertical average of velocity and buoyancy from 
c     the surface down to epsilon*zgrid(kl).
c
c     when the bulk richardson number at k exceeds Ricr, hbl is
c     linearly interpolated between grid levels zgrid(k) and zgrid(k-1).
c
c     The water column and the surface forcing are diagnosed for 
c     stable/ustable forcing conditions, and where hbl is relative 
c     to grid points (caseA), so that conditional branches can be 
c     avoided in later subroutines.
c-----------------------------------------------------------------------

c  model  
      integer km,kmp1      ! number of vertical levels
      integer imt          ! number of horizontal grid points
      real zgrid(imt,kmp1) ! vertical grid (<= 0)              (m)
      real hwide(imt,kmp1) ! layer thicknesses                 (m)
c
c  input
      real dVsq(imt,kmp1)  ! (velocity shear re sfc)^2      (m/s)^2
      real dbloc(imt,km)   ! local delta buoyancy              (m/s^2)
      real Ritop(imt,km)   ! numerator of bulk Richardson Number (m/s)^2
c          Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc(m/s^2)
      real ustar(imt)      ! surface friction velocity         (m/s)
      real Bo(imt)         ! surface turbulent buoyancy forcing(m^2/s^3)
      real Bosol(imt)      ! radiative buoyancy forcing        (m^2/s^3)
      real Coriol          ! Coriolis parameter                (1/s)
      integer jwtype(imt)  ! Jerlov water type                 (1 to 5)
c
c  output
      real hbl(imt)        ! boundary layer depth              (m)
      real bfsfc(imt)      ! Bo+radiation absorbed to d=hbf*hbl(m^2/s^3)
      real stable(imt)     ! =1 in stable forcing; =0 unstable
      real caseA(imt)      ! =1 in case A, =0 in case B 
      integer kbl(imt)     ! index of first grid level below hbl 
c 
c  local
      real Rib(imt,2)      ! Bulk Richardson number
      real sigma(imt)      ! normalized depth (d/hbl)
      real wm(imt),ws(imt) ! turbulent velocity scales         (m/s)
 
*... Modifications for Gent version:
      logical check_for_HEKman_or_hmonob
      data check_for_HEKman_or_hmonob /.false./
 
      save epsln,Ricr,epsilon,cekman,cmonob,cs,cv,vonk,hbf
c
      data epsln           /  1.e-20 /
      data Ricr            /  0.30   /
      data epsilon         /  0.1    /
      data cekman          /  0.7    /
      data cmonob          /  1.0    /
      data cs              / 98.96   /
      data cv              /  1.8    /
      data vonk            /  0.4    /
      data hbf             /  1.0    /

c-----------------------------------------------------------------------
c find bulk Richardson number at every grid level until > Ric
c
c note: the reference depth is -epsilon/2.*zgrid(i,k), but the reference
c       u,v,t,s values are simply the surface layer values,
c       and not the averaged values from 0 to 2*ref.depth,
c       which is necessary for very fine grids(top layer < 2m thickness)
c note: max values when Ricr never satisfied are
c       kbl(i)=km and hbl(i) -zgrid(i,km)
 
      Vtc =  cv * sqrt(0.2/cs/epsilon) / vonk**2 / Ricr
  
c     indices for array Rib(i,k), the bulk Richardson number.
      ka = 1
      ku = 2
 
c     initialize hbl and kbl to bottomed out values
      do 100 i = i1, i2
         Rib(i,ka) = 0.0
         kbl(i)    = km
         hbl(i)    = -zgrid(i,km)
 100  continue

      do 200 kl = 2,km
 
c        compute bfsfc = sw fraction at hbf * zgrid
         call swfrac(imt,i1, i2, hbf,zgrid(1,kl),jwtype,bfsfc)

         do 190 i = i1, i2
 
c           use caseA as temporary array
            caseA(i)  = -zgrid(i,kl)
 
c           compute bfsfc= Bo + radiative contribution down to hbf * hbl
            bfsfc(i)  = Bo(i) 
     $                  + Bosol(i) * (1. - bfsfc(i))
            stable(i) = 0.5 + SIGN( 0.5, bfsfc(i) )
            sigma(i)  = stable(i) * 1. + (1.-stable(i)) * epsilon

  190    continue

c        compute velocity scales at sigma, for hbl= caseA = -zgrid(i,kl)
         call wscale(imt, i1, i2, sigma, caseA, ustar, bfsfc,   wm, ws)

         do 180 i = i1, i2
 
c           compute the turbulent shear contribution to Rib
            bvsq =0.5*
     $            ( dbloc(i,kl-1) / (zgrid(i,kl-1)-zgrid(i,kl  ))+
     $              dbloc(i,kl  ) / (zgrid(i,kl  )-zgrid(i,kl+1)) )
            Vtsq = - zgrid(i,kl) * ws(i) * sqrt(abs(bvsq)) * Vtc
  
c           compute bulk Richardson number at new level, dunder
            Rib(i,ku) = Ritop(i,kl) / (dVsq(i,kl)+Vtsq+epsln)
            if((kbl(i).eq.km).and.(Rib(i,ku).ge.Ricr)) then
 
c              linear interpolate to find hbl where Rib = Ricr
               hbl(i) = -zgrid(i,kl-1) + (zgrid(i,kl-1)-zgrid(i,kl)) *
     $                  (Ricr - Rib(i,ka)) / (Rib(i,ku)-Rib(i,ka))
               kbl(i) = kl
            endif

 180     continue

         ksave = ka
         ka    = ku
         ku    = ksave

 200  continue

c compare hbl to limits
      call swfrac(imt,i1, i2, -1.0,hbl,jwtype,bfsfc)
 
      do 300 i = i1, i2
         bfsfc(i)  = Bo(i) 
     $               + Bosol(i) * (1. - bfsfc(i))
         stable(i) = 0.5 + SIGN( 0.5, bfsfc(i) )
         bfsfc(i)  = bfsfc(i) + stable(i) * epsln !ensures bfsfc never=0
 300  continue

*.... Modification for Gent code:
      if (check_for_HEKman_or_hmonob) then
 
c             check for HEKman or hmonob (11/22/94)
       do 400 i = i1, i2
          if(bfsfc(i).gt.0.0) THEN
          hekman = cekman * ustar(i) / (abs(Coriol)+epsln)
          hmonob = cmonob * ustar(i)*ustar(i)*ustar(i) 
     &            /vonk / (bfsfc(i)+epsln) 
          hlimit = stable(i)     * AMIN1(hekman,hmonob) +
     &            (stable(i)-1.) * zgrid(i,km)
          hbl(i) = AMIN1(hbl(i),hlimit)
          hbl(i) = AMAX1(hbl(i),-zgrid(i,1))
          ENDIF
          kbl(i) = km
 400   continue
c                find new kbl
          do 405 kl=2,km
          do 415 i =i1, i2
          if((kbl(i).eq.km).and.(-zgrid(i,kl).gt.hbl(i))) then
            kbl(i) = kl
          endif
 415      continue
 405      continue

c find stability and buoyancy forcing for final hbl values
       call swfrac(imt,i1, i2, -1.0,hbl,jwtype,bfsfc)
       do 500 i = i1, i2
          bfsfc(i)  = Bo(i) 
     $               + Bosol(i) * (1. - bfsfc(i))
          stable(i) = 0.5 + SIGN( 0.5, bfsfc(i) )
          bfsfc(i)  = bfsfc(i) + stable(i) * epsln 
  500  continue

      endif ! (check_for_HEKman_or_hmonob) then
 
 
c determine caseA and caseB
      do 600 i = i1, i2
         caseA(i)  = 0.5 + 
     $       SIGN( 0.5,-zgrid(i,kbl(i)) -0.5*hwide(i,kbl(i)) -hbl(i))
 600  continue

      return
      end

c=======================================================================

      SUBROUTINE wscale(imt, i1  , i2  , sigma, hbl, ustar, bfsfc,
     $                  wm , ws  )

c-----------------------------------------------------------------------
c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c-----------------------------------------------------------------------

c lookup table
      parameter ( ni = 890,              ! number of values for zehat
     $            nj = 48)              ! number of values for ustar

      real wmt(0:ni+1,0:nj+1)           ! lookup table for wm
      real wst(0:ni+1,0:nj+1)           ! lookup table for ws
      real deltaz                       ! delta zehat in table
      real deltau                       ! delta ustar in table
      real zmin,zmax                    ! zehat limits for table
      real umin,umax                    ! ustar limits for table
      logical firstf
      save wmt,wst,deltaz,deltau,zmin,zmax,umin,umax,firstf
c
      data zmin,zmax  / -4.e-7, 0.0   / ! m3/s3
      data umin,umax  /  0.   , .04   / ! m/s
      data firstf     / .true.        /

c  model
      integer imt          ! number of horizontal grid points

c  input
      real sigma(imt)      ! normalized depth (d/hbl)
      real hbl(imt)        ! boundary layer depth (m)
      real ustar(imt)      ! surface friction velocity         (m/s)
      real bfsfc(imt)      ! total surface buoyancy flux       (m^2/s^3)

c  output
      real wm(imt),ws(imt) ! turbulent velocity scales at sigma

c local
      real zehat           ! = zeta *  ustar**3
      real zeta            ! = stability parameter d/L

      save epsln,c1,am,cm,c2,zetam,as,cs,c3,zetas,vonk

      data epsln           /   1.0e-20/
      data c1              /   5.0   /
      data am,cm,c2,zetam  /   1.257 ,  8.380 , 16.0 , - 0.2  /
      data as,cs,c3,zetas  / -28.86  , 98.96  , 16.0 , - 1.0  /
      data vonk            /   0.40  /

c-----------------------------------------------------------------------
c
c construct the wm and ws lookup tables
c
      if(firstf) then

         deltaz = (zmax-zmin)/(ni+1) 
         deltau = (umax-umin)/(nj+1)

         do 100 i=0,ni+1
            zehat = deltaz*(i) + zmin
            do 90 j=0,nj+1
               usta = deltau*(j) + umin
               zeta = zehat/(usta**3+epsln)

               if(zehat.ge.0.) then
                  wmt(i,j) = vonk*usta/(1.+c1*zeta)
                  wst(i,j) = wmt(i,j)
               else
                  if(zeta.gt.zetam) then
                     wmt(i,j) = vonk* usta * (1.-c2*zeta)**(1./4.)
                  else
                     wmt(i,j) = vonk* (am*usta**3 - cm*zehat)**(1./3.)
                  endif
                  if(zeta.gt.zetas) then
                     wst(i,j) = vonk* usta * (1.-c3*zeta)**(1./2.)
                  else
                     wst(i,j) = vonk* (as*usta**3 - cs*zehat)**(1./3.)
                  endif
               endif   
 90         continue   
 100     continue
         firstf=.false.
      endif       
 
c use lookup table for zehat < zmax  ONLY;  otherwise use stable formulae
      do 200 i = i1, i2
         zehat = vonk * sigma(i) * hbl(i) * bfsfc(i)

         IF (zehat .le. zmax) THEN

         zdiff  = zehat-zmin
         iz = int( zdiff/deltaz )
         iz = min( iz , ni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar(i)-umin
         ju = int( udiff/deltau)
         ju = min( ju , nj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff/deltaz - float(iz)
         ufrac = udiff/deltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm(i) = (1.-ufrac)* wbm          + ufrac*wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws(i) = (1.-ufrac)* wbs          + ufrac*was

         ELSE

         ucube = ustar(i)**3
         wm(i) = vonk * ustar(i) * ucube / (ucube + c1 * zehat)
         ws(i) = wm(i)

         ENDIF   

  200 continue

      return
      end

c=======================================================================
 
      SUBROUTINE ri_iwmix ( km  , kmp1  , imt   , i1 ,  i2 ,
     $                      Shsq, dbloc , zgrid ,
     $                      visc, difs  , dift  , nu, kappa  )

c-----------------------------------------------------------------------
c     compute interior viscosity diffusivity coefficients due to
c     shear instability (dependent on a local richardson number)
c     and due to background internal wave activity.
c-----------------------------------------------------------------------

c  input
      real Shsq(imt,kmp1)    ! (local velocity shear)^2          (m/s)^2
      real dbloc(imt,km)     ! local delta buoyancy              (m/s^2)
      real zgrid(imt,kmp1)   ! vertical grid (<= 0)              (m)
      integer km,kmp1        ! number of vertical levels
      integer imt            ! number of horizontal grid points
      real nu           (km) ! vertical viscosity   profile
      real kappa        (km) ! vertical diffusivity profile
 
c output
      real visc(imt,0:kmp1)  ! vertical viscosity coefficient    (m^2/s)
      real difs(imt,0:kmp1)  ! vertical scalar diffusivity       (m^2/s)
      real dift(imt,0:kmp1)  ! vertical temperature diffusivity  (m^2/s)
 
c local variables
      real Rigg              ! local richardson number
      real fri               ! function of Rig
 
      save epsln,Riinfty,difm0,difs0,difmiw,difsiw

      data  epsln   / 1.e-20  /
      data  Riinfty / 0.8     /  ! local Ri limit for shear instability
      data  difm0   / 0.005   /  ! max visc due to shear instability
      data  difs0   / 0.005   /  ! max diff ..  .. ..    ..
****** Original values:
*     data  difmiw  / 0.0001  /  ! background/internal waves visc(m^2/s)
*     data  difsiw  / 0.00001 /  ! ..         ..       ..    diff(m^2/s)
*... Modifications for Gent version:
      data  difmiw  / 0.0     /  ! background/internal waves visc(m^2/s)
      data  difsiw  / 0.0     /  ! ..         ..       ..    diff(m^2/s)

      integer    num_v_smooth_Ri ! # of times Ri is vertically smoothed
      parameter (num_v_smooth_Ri=1)

c-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! compute interior gradient Ri at all interfaces, except surface
      ! use visc(i,k) as temporary storage for Ri
      !-----------------------------------------------------------------
      do k =  1,km
      do i = i1,i2
         visc(i,k) = dbloc(i,k) * (zgrid(i,k)-zgrid(i,k+1)) / 
     &                            ( Shsq(i,k) + epsln)
      end do
      end do
   
      !-----------------------------------------------------------------
      ! vertically smooth Ri num_v_smooth_Ri times
      !-----------------------------------------------------------------
      do mr = 1,num_v_smooth_Ri
         call z121(km,imt,i1,i2,visc)
      enddo
 
      !-----------------------------------------------------------------
      ! evaluate diffusivities and viscosity mixing 
      ! due to internal waves, and shear and static instability
      !-----------------------------------------------------------------
      do k =  1,km
      do i = i1,i2
         !--- evaluate function of Ri ---
         Rigg  = MAX( visc(i,k)  , 0.0)
         ratio = MIN( Rigg/Riinfty, 1.0 )
         fri   = (1. - ratio*ratio)
         fri   = fri*fri*fri
   
         !--- sum shear + internal wave diffusivities ---
         visc(i,k) = difm0 * fri  + difmiw + nu(k)
         difs(i,k) = difs0 * fri  + difsiw + kappa(k)
         dift(i,k) = difs(i,k) 
      end do
      end do

      return
      end

c=======================================================================

      subroutine z121 (km,imt,i1,i2,V)
 
      !---- input ----
      integer km
      integer imt
      integer i1,i2

      !---- input/output ----
      real V(imt,0:km+1)  ! 2-D array to be smoothed in k direction
 
      !---- constants ----
      parameter (p25 = 0.25    )
      parameter (p5  = 0.5     )

c-----------------------------------------------------------------------
c    Apply 121 smoothing in k to 2-d array V(i,k=1,km)
c     top (0) value is used as a dummy
c     bottom (km+1) value is set to input value from above.
c-----------------------------------------------------------------------
 
      do i=i1,i2
         V(i,0)    =  p25 * V(i,1)
         V(i,km+1) =        V(i,km)
      enddo
 
      do k=1,km
      do i=i1,i2
          tmp      =  V(i,k)
          V(i,k)   =  V(i,0)  + p5 * V(i,k) + p25 * V(i,k+1)
          V(i,0)   =  p25 * tmp
      enddo
      enddo
 
      return
      end

c=======================================================================

      SUBROUTINE blmix 
     $             (km   , kmp1 , imt , i1    , i2   ,
     $              mdiff, zgrid, hwide,  
     $              ustar, bfsfc, hbl , stable, caseA, 
     $              visc , difs , dift, kbl   , 
     $              gat1 , dat1 , dkm1, blmc  , ghats, 
     $              sigma, wm   , ws  )

c-----------------------------------------------------------------------
c mixing coefficients within boundary layer depend on surface
c forcing and the magnitude and gradient of interior mixing below
c the boundary layer ("matching").
c
c CAUTION if mixing bottoms out at hbl = -zgrid(km) THEN
c fictious layer kmp1 is needed with small but finite width (eg. 1.e-10)
c-----------------------------------------------------------------------

c model
      integer km,kmp1        ! number of vertical levels
      integer imt            ! number of horizontal grid points
      integer mdiff          ! number of viscosities + diffusivities
      real zgrid(imt,kmp1)   ! vertical grid (<=0)               (m)
      real hwide(imt,kmp1)   ! layer thicknesses                 (m)
c
c input
      real ustar(imt)        ! surface friction velocity         (m/s)
      real bfsfc(imt)        ! surface buoyancy forcing        (m^2/s^3)
      real hbl(imt)          ! boundary layer depth              (m)
      real stable(imt)       ! = 1 in stable forcing
      real caseA(imt)        ! = 1 in case A
      real visc(imt,0:kmp1)  ! vertical viscosity coefficient    (m^2/s)
      real difs(imt,0:kmp1)  ! vertical scalar diffusivity       (m^2/s)
      real dift(imt,0:kmp1)  ! vertical temperature diffusivity  (m^2/s)
 
      integer kbl(imt)       ! index of first grid level below hbl
c
c output
      real gat1(imt,mdiff)
      real dat1(imt,mdiff)
      real dkm1(imt,mdiff)   ! boundary layer difs at kbl-1 level
      real blmc(imt,km,mdiff)! boundary layer mixing coefficients(m^2/s)
      real ghats(imt,km)     ! nonlocal scalar transport
c
c  local
      real sigma(imt)        ! normalized depth (d / hbl)
      real ws(imt), wm(imt)  ! turbulent velocity scales         (m/s)
 
      save epsln,epsilon,c1,cs,cstar,vonk

      data epsln             /   1.e-20 /
      data epsilon           /   0.1    /
      data c1                /   5.0    /
      data cs                /  98.96   /
      data cstar             /   5.0    /
      data vonk              /   0.4    /
c
      cg = cstar * vonk * (cs * vonk * epsilon)**(1./3.)
       
c compute velocity scales at hbl
      do 100 i = i1, i2
         sigma(i) = stable(i) * 1.0 + (1.-stable(i)) * epsilon
 100  continue
 
      call wscale(imt, i1, i2, sigma, hbl, ustar, bfsfc,   wm, ws)
      do 200 i = i1, i2 
         kn    = ifix(caseA(i)+epsln) *(kbl(i) -1) +
     $            (1-ifix(caseA(i)+epsln)) * kbl(i)
 
c find the interior viscosities and derivatives at hbl(i) 
         delhat = 0.5*hwide(i,kn) - zgrid(i,kn) - hbl(i)
         R      = 1.0 - delhat / hwide(i,kn)
         dvdzup = (visc(i,kn-1) - visc(i,kn)) / hwide(i,kn) 
         dvdzdn = (visc(i,kn)   - visc(i,kn+1)) / hwide(i,kn+1)
         viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $                        R  * (dvdzdn + abs(dvdzdn)) )
 
         dvdzup = (difs(i,kn-1) - difs(i,kn)) / hwide(i,kn) 
         dvdzdn = (difs(i,kn)   - difs(i,kn+1)) / hwide(i,kn+1)
         difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $                        R  * (dvdzdn + abs(dvdzdn)) )
 
         dvdzup = (dift(i,kn-1) - dift(i,kn)) / hwide(i,kn) 
         dvdzdn = (dift(i,kn)   - dift(i,kn+1)) / hwide(i,kn+1)
         diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $                        R  * (dvdzdn + abs(dvdzdn)) )
 
         visch  = visc(i,kn) + viscp * delhat 
         difsh  = difs(i,kn) + difsp * delhat 
         difth  = dift(i,kn) + diftp * delhat 
c
         f1 = stable(i) * c1 * bfsfc(i) / (ustar(i)**4+epsln) 

         gat1(i,1) = visch / hbl(i) / (wm(i)+epsln)
         dat1(i,1) = -viscp / (wm(i)+epsln) + f1 * visch
         dat1(i,1) = min(dat1(i,1),0.) 

         gat1(i,2) = difsh  / hbl(i) / (ws(i)+epsln)
         dat1(i,2) = -difsp / (ws(i)+epsln) + f1 * difsh 
         dat1(i,2) = min(dat1(i,2),0.) 
  
         gat1(i,3) = difth /  hbl(i) / (ws(i)+epsln)
         dat1(i,3) = -diftp / (ws(i)+epsln) + f1 * difth 
         dat1(i,3) = min(dat1(i,3),0.) 

 200  continue
c
      do 300 ki = 1,km       
c
c     compute turbulent velocity scales on the interfaces
c
         do 290 i  = i1, i2
            sig     = (-zgrid(i,ki) + 0.5 * hwide(i,ki)) / hbl(i)
            sigma(i)= stable(i)*sig + (1.-stable(i))*AMIN1(sig,epsilon)
 290     continue
         call wscale(imt, i1, i2, sigma, hbl, ustar, bfsfc,   wm,  ws)
c
c     compute the dimensionless shape functions at the interfaces
c
         do 280 i = i1, i2
            sig = (-zgrid(i,ki) + 0.5 * hwide(i,ki)) / hbl(i)
            a1 = sig - 2.
            a2 = 3.-2.*sig
            a3 = sig - 1.
c
            Gm = a1 + a2 * gat1(i,1) + a3 * dat1(i,1) 
            Gs = a1 + a2 * gat1(i,2) + a3 * dat1(i,2)
            Gt = a1 + a2 * gat1(i,3) + a3 * dat1(i,3)
c
c     compute boundary layer diffusivities at the interfaces
c
            blmc(i,ki,1) = hbl(i) * wm(i) * sig * (1. + sig * Gm)
            blmc(i,ki,2) = hbl(i) * ws(i) * sig * (1. + sig * Gs)
            blmc(i,ki,3) = hbl(i) * ws(i) * sig * (1. + sig * Gt)
c
c     nonlocal transport term = ghats * <ws>o
            ghats(i,ki) = (1.-stable(i)) * cg / (ws(i) * hbl(i) +epsln)

 280     continue

 300  continue
 
c find diffusivities at kbl-1 grid level 
      do 400 i=i1, i2
         sig      =  -zgrid(i,kbl(i)-1)  / hbl(i)
         sigma(i) =stable(i) * sig + (1.-stable(i)) * AMIN1(sig,epsilon)
 400  continue
c
      call wscale(imt, i1, i2, sigma, hbl, ustar, bfsfc,   wm, ws)
c
      do 500 i = i1, i2
         sig = -zgrid(i,kbl(i)-1) / hbl(i)
         a1= sig - 2.
         a2 = 3.-2.*sig
         a3 = sig - 1.
         Gm = a1 + a2 * gat1(i,1) + a3 * dat1(i,1)
         Gs = a1 + a2 * gat1(i,2) + a3 * dat1(i,2)
         Gt = a1 + a2 * gat1(i,3) + a3 * dat1(i,3)
         dkm1(i,1) = hbl(i) * wm(i) * sig * (1. + sig * Gm)
         dkm1(i,2) = hbl(i) * ws(i) * sig * (1. + sig * Gs)
         dkm1(i,3) = hbl(i) * ws(i) * sig * (1. + sig * Gt)
  500 continue

      return
      end

c=======================================================================

      SUBROUTINE enhance (km   , kmp1  , imt   , i1    , i2    ,
     &                    mdiff, dkm1  , visc  ,
     &                    difs , dift  , hbl   , kbl   , zgrid , caseA,
     &                    blmc , ghats )

c-----------------------------------------------------------------------
c
c enhance the diffusivity at the kbl-.5 interface
c
c-----------------------------------------------------------------------

c input
      integer km,kmp1        ! number of vertical levels
      integer imt            ! number of horizontal grid points
      integer mdiff          ! number of viscosities + diffusivities
      integer kbl(imt)       ! grid above hbl
      real hbl(imt)          ! boundary layer depth             (m)
      real dkm1(imt,mdiff)   ! bl diffusivity at kbl-1 grid level
      real zgrid(imt,kmp1)   ! vertical grid (<= 0)             (m)
      real visc(imt,0:kmp1)  ! enhanced viscosity               (m^2/s) 
      real difs(imt,0:kmp1)  ! enhanced thermal diffusivity     (m^2/s)
      real dift(imt,0:kmp1)  ! enhanced scalar  diffusivity     (m^2/s)
      real caseA(imt)        ! = 1 in caseA, = 0 in case B
 
c input/output
      real ghats(imt,km)     ! nonlocal transport               (s/m**2)
c                              modified ghats at kbl(i)-1 interface
c output
      real blmc(imt,km,mdiff)! enhanced bound. layer mixing coeff.
c
c local
      real delta             ! fraction hbl lies beteen zgrid neighbors

c-----------------------------------------------------------------------

      do 100 ki=1,km-1
        do 90 i = i1, i2

          if(ki .eq. (kbl(i) - 1) ) then

            delta = (hbl(i)+zgrid(i,ki)) / (zgrid(i,ki)-zgrid(i,ki+1))

            dkmp5 = caseA(i) * visc(i,ki) + (1.-caseA(i)) * blmc(i,ki,1)
            dstar = (1.-delta)**2 * dkm1(i,1) + delta**2 * dkmp5      
            blmc(i,ki,1) = (1.-delta) * visc(i,ki) + delta * dstar

            dkmp5 = caseA(i) * difs(i,ki) + (1.-caseA(i)) * blmc(i,ki,2)
            dstar = (1.-delta)**2 * dkm1(i,2) + delta**2 * dkmp5    
            blmc(i,ki,2) = (1.-delta) * difs(i,ki) + delta * dstar

            dkmp5 = caseA(i) * dift(i,ki) + (1.-caseA(i)) * blmc(i,ki,3)
            dstar = (1.-delta)**2 * dkm1(i,3) + delta**2 * dkmp5     
            blmc(i,ki,3) = (1.-delta) * dift(i,ki) + delta * dstar
            
            ghats(i,ki) = (1.-caseA(i)) * ghats(i,ki)

          endif

 90     continue
 100  continue

      return
      end

c=======================================================================
 
      SUBROUTINE swfrac( imt, i1, i2, fact, z, jwtype, swdk )

c-----------------------------------------------------------------------
c     compute fraction of solar short-wave flux penetrating to specified
c     depth (times fact) due to exponential decay in  Jerlov water type
c     reference : two band solar absorption model of simpson and 
c     paulson (1977)
c-----------------------------------------------------------------------

      parameter(nwtype=5) ! max number of different water types 
c
c  model
      integer imt         ! number of horizontal grid points

c  input
      real fact           ! scale  factor to apply to depth array
      real z(imt)         ! vertical height ( <0.) for desired sw 
c                           fraction                                 (m)
      integer jwtype(imt) ! index for jerlov water type

c  output
      real swdk(imt)      !  short wave (radiation) fractional decay

c  local
      real  rfac(nwtype),a1(nwtype),a2(nwtype)
      save  rfac,a1,a2
c
c     jerlov water type :  I       IA      IB      II      III
c                jwtype    1       2       3       4       5
      data rfac         /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
      data a1           /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
      data a2           / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /

c-----------------------------------------------------------------------
      do 100 i = i1, i2

         swdk(i) =      rfac(jwtype(i))  * exp(z(i)*fact/a1(jwtype(i)))
     $            + (1.-rfac(jwtype(i))) * exp(z(i)*fact/a2(jwtype(i)))

  100 continue

      return
      end

c=======================================================================
