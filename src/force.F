c=======================================================================
c CVS $Id: force.F,v 1.2 1998/08/24 22:45:46 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/force.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE force

      implicit none

c-----------------------------------------------------------------------
c PURPOSE:
c   compute the following surface forcing fields 
c   o in force_tau:  TAUx, TAUy, WNDs
c   o in force_q  :  Qtot, Qsol, Qlat, Qsen, Qlgw
c   o in force_sst:  SST (used wrt vertical boundary conditions)
c-----------------------------------------------------------------------

      call force_tau
      call force_q
      call force_sst

      RETURN
      END

c=======================================================================

      SUBROUTINE force_tau

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! OUTPUT - forcing data goes here

      !----- local -----
      integer i,j,k              ! generic indicies
      integer month1,month2      ! two month to interpolate between
      real     fact1, fact2      ! weights for the two months data
      integer  eday1, eday2      ! elapsed days since y-m-d=0-1-1 & 0s
      real    wndx(nxp1,nyp1,12) ! climatological wind data
      real    wndy(nxp1,nyp1,12)
      real    WINDu,WINDv        ! instantaneous wind vector

      character*20 loc_fn      ! local      file name
      character*64 mss_fn      ! mass-store file name
      character*16 opts        ! mswrite options
      character*16 pass        ! mswrite password(s)
      integer      ier         ! msread/mswrite error code
      integer      LstIndex    ! last index funct (see intrinsic index)

      logical      first_call  ! flags one-time initializations
      data         first_call / .true. /

      SAVE         

c-----------------------------------------------------------------------
c PURPOSE:
c   compute TAUx,TAUx,WNDs
c   the wind stress & wind speed, based on 12 monthly mean wind fields
c-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! initializations
      !-----------------------------------------------------------------
      IF ( first_call) THEN
        write(*,9000) 'initializing climatological wind data'
        write(*,9000) 'using an annual cycle thru mean monthly winds'
        write(*,9000) 'Some wind stress parameters:'
        write(*,9008) 'o drag coefficient      : Cd     = ', Cd  
        write(*,9008) 'o density of water      : RHOw   = ', RHOw
        write(*,9008) 'o density of air        : RHOa   = ', RHOa
        write(*,9008) 'o minimum wind speed    : WNDmin = ', WNDmin

        mss_fn = frc_wfile
        i      = LstIndex(mss_fn,'/')+1
        j      =    index(mss_fn,' ')
        loc_fn = mss_fn(i:j)//' '
        pass   = ' '
        opts   = ' '
        call MSread(ier, loc_fn , mss_fn, pass, opts)
        write(*,9010) ier,loc_fn(1:j-i+1),mss_fn(1:j)

        open(tmpunit,file=loc_fn,form='UNFORMATTED')
        DO k=1,12
           read(tmpunit) ((WNDx(i,j,k),i=1,nxp1),j=1,nyp1)
           read(tmpunit) ((WNDy(i,j,k),i=1,nxp1),j=1,nyp1)
        END DO
        close(tmpunit)

        DO j=1,nyp1
        DO i=1,nxp1
           TAUx(i,j) = SPVAL
           TAUy(i,j) = SPVAL
           WNDs(i,j) = SPVAL
        END DO
        END DO
      
        first_call = .false.
      END IF

      !----------------------------------------------------------------
      ! time interpolate between two monthly means
      !----------------------------------------------------------------
      if (day.lt.16) then
         !--- day=1  & sec=0     => fact1=0.5 & fact2=0.5
         !--- day=15 & sec=86400 => fact1=0.0 & fact2=1.0
         month1 = month - 1
         month2 = month
         if (month1.eq.0) month1=12
         fact2 = 0.5 + 0.5*(float(day-1) + sec/86400.0)/15.0
         fact1 = 1.0 - fact2
      else
         !--- day=16 & sec=0                => fact1=1.0 & fact2=0.0
         !--- last day of month & 86400 sec => fact1=0.5 & fact2=0.5
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
         fact1 = 1.0 - 0.5*(float(day-16) + sec/86400.0)/(eday2-eday1)
         fact2 = 1.0 - fact1
      end if

      DO j=1,nyp1
      DO i=it1(j),it2(j)
        WINDu     = fact1*WNDx(i,j,month1)+fact2*WNDx(i,j,month2)
        WINDv     = fact1*WNDy(i,j,month1)+fact2*WNDy(i,j,month2)
        WNDs(i,j) = sqrt(WINDu**2 + WINDv**2)
        TAUx(i,j) = Cd*RHOa*WINDu*WNDs(i,j)
        TAUy(i,j) = Cd*RHOa*WINDv*WNDs(i,j)
        WNDs(i,j) = max(WNDmin,WNDs(i,j))
      END DO
      END DO
      
      RETURN
 9000 format('(force_tau) ',a,a,a)
 9008 format('(force_tau) ',a,1pe10.3)
 9010 format('(force_tau) return code ',i2,' on  msread ',a,a)
      END

c=======================================================================

      SUBROUTINE force_q

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! OUTPUT - forcing data goes here

      !----- heat flux local constants -----
      real       albedo      ! sea surface albedo
      real       aCLD        ! solar    cloud coefficient
      real       bCLD        ! longwave cloud coefficient
      real        ALPHA      ! solar altitude (in degrees)
      real       aALPHA      ! solar altitude coefficient
      real       RH          ! relative humidity
      real       L           ! latent heat of vaporization
      real       gamma       ! molecular weight ratio
      real       CeL         ! Ce*L
      real       RaCeCp      ! RHOa*Ce*Cp
      real       EpsSig      ! Epsilon*Sigma
      real       RaGam       ! RHOa*Gamma/1000mb
      real       DAYtoRAD    ! convert days to rad: assume 365 days/yr
      complex    cPI18       ! e^(0,-PI/18.0)
      parameter (albedo   =  0.06        )
      parameter (aCLD     =  0.75        )
      parameter (bCLD     =  0.6         )
      parameter (aALPHA   =  0.002       )
      parameter (RH       =  0.75        )
      parameter (L        =  2.45e+3     )
      parameter (gamma    =  0.62        )
      parameter (CeL      =  Ce*L        )
      parameter (RaCeCp   =  RHOa*Ce*Cpa )
      parameter (EpsSig   =  5.500e-8    )
      parameter (RaGam    =  RHOa*gamma  )
      parameter (DAYtoRAD =  PI/182.5    )
      parameter (cPI18    = (0.98481,-0.17365))

      !----- heat flux local variables -----
      real cld(nxp1,nyp1,12) ! climatological cloud cover
      real         cloud     ! instantaneous cloud cover
      real         julday    ! Julian day
      real         Q0        ! clear sky solar radiation
      real         RaHum     ! RH*Qs
      real         Tatm      ! atmospheric temperature (degrees C)
      real         Tocn      ! ocean       temperature (degrees C)

      !--- misc vars used to derive clear sky solar radiation ---
      real         PHI
      complex      cexPHI
      complex      cexPHI2
      real         cosPHI ,sinPHI
      real         cosPHI2,sinPHI2
      real         RLAT,cosRLAT,sinRLAT
      complex      cexRLAT
      real         A0,A1,A2,B1,B2

      integer      i,j,k         ! generic loop indicies
      integer      month1,month2 ! two months to interpolate between
      real         fact1 , fact2 ! weights for the two months data
      integer      eday1 , eday2 ! elapsed days since 1 Jan 0000 & 0s
      character*20 loc_fn        ! local      file name
      character*64 mss_fn        ! Mass Store file name
      character*16 opts          ! msread options
      character*16 pass          ! msread password(s)
      integer      ier           ! msread error code
      integer      LstIndex      ! last index function
      logical      first_call    ! flags one-time initializations

      data         first_call / .true. /

      SAVE

      !--- local function ---
      real  Qs  ! function: atm dens * Clausius-Clapeyron spec hum
      real  x   ! function argument
      Qs(x) = 6.4e+08*exp(-5105.0/(x+273.0))

c-----------------------------------------------------------------------
c  PURPOSE:
c    compute surface heat flux
c
c  OUTPUT:
c    heat flux components (watts/meter**2)
c    o Qsol ~ solar 
c    o Qlat ~ latent
c    o Qsen ~ sensible
c    o Qlgw ~ net longwave
c    o Qtot ~ total heat flux = Qsol - Qlat - Qsen - Qlgw 
c
c  METHOD:
c    o time interpolate climatology data to get current cloud cover
c
c    o compute clear sky solar using julian day and the solar angle 
c      formula from Weare 1980 p.8-9.  Note that, using w/m**2, radians,
c      and assuming Y lies in [67s,67n] degrees, the formula becomes:
c      A1    = 9.63 + 192.44*cos(Y + 90)
c            = 9.63 + 192.44*sin(Y     )
c      alpha = arcsin[ cos(Y)cos(Ysun) + sin(Y)sin(Ysun)]
c            = arcsin[ cos( Y - Ysun       ) ]
c            = arcsin[ sin( Y - Ysun + 90. ) ]
c            = 90 - abs( Y - Ysun), for Y in [-67,67] degrees
c
c    o compute net heat flux using bulk fomula by E. Brady (NCAR)
c
c  ASSUMPTIONS:
c    o The array WNDs(i,j) is valid for the current time step
c
c-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! initializations
      !-----------------------------------------------------------------
      if ( first_call) then
        mss_fn = frc_cfile
        i      = LstIndex(mss_fn,'/')+1
        j      =    index(mss_fn,' ')
        loc_fn = mss_fn(i:j)//' '
        pass   = ' '
        opts   = ' '
        call MSread(ier, loc_fn , mss_fn, pass, opts)
        write(*,9010) '(force_q)  ',ier,loc_fn(1:j-i+1),mss_fn(1:j)
 9010   format(a,' return code ',i2,' on  msread ',a,a)

        open(tmpunit,file=loc_fn,form='UNFORMATTED')
        DO k=1,12
           read(tmpunit) ((CLD(i,j,k),i=1,nxp1),j=1,nyp1)
        END DO
        close(tmpunit)

        write(*,9000) 'Some heat flux parameters:'
        write(*,9001) '  o sea surface albedo       : albedo = ', albedo
        write(*,9001) '  o cloud coeff (short-wave) : aCLD   = ', aCLD
        write(*,9001) '  o cloud coeff (long -wave) : bCLD   = ', bCLD
        write(*,9002) '  o solar altitude coeff     : aALPHA = ', aALPHA
        write(*,9001) '  o relative humidity        : RH     = ', RH

        DO j=1,nyp1
        DO i=1,nxp1
            Qsol(i,j) = SPVAL
            Qlat(i,j) = SPVAL
            Qsen(i,j) = SPVAL
            Qlgw(i,j) = SPVAL
            Qtot(i,j) = SPVAL
        END DO
        END DO
        first_call = .false.
      end if

      !----------------------------------------------------------------
      ! get cloud cover by time interpolation between monthly means
      !----------------------------------------------------------------
      if (day.lt.16) then
         !--- day=1  & sec=0     => fact1=0.5 & fact2=0.5
         !--- day=15 & sec=86400 => fact1=0.0 & fact2=1.0
         month1 = month - 1
         month2 = month
         if (month1.eq.0) month1=12
         fact2 = 0.5 + 0.5*(float(day-1) + sec/86400.0)/15.0
         fact1 = 1.0 - fact2
      else
         !--- day=16 & sec=0                => fact1=1.0 & fact2=0.0
         !--- last day of month & 86400 sec => fact1=0.5 & fact2=0.5
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
         fact1 = 1.0 - 0.5*(float(day-16) + sec/86400.0)/(eday2-eday1)
         fact2 = 1.0 - fact1
      end if

      !----------------------------------------------------------------
      ! compute heat flux from SST, wind, cloud cover & solar angle
      !----------------------------------------------------------------
      call ymd2eday(0,month,day,eday1)
      julday = eday1 + sec/86400.0

      PHI     = (julday-21.0)*DAYtoRAD
      cexPHI  = cEXP ((0,1)*PHI)
      cosPHI  = REAL ( cexPHI )
      sinPHI  = aIMAG( cexPHI )
      cexPHI2 = cexPHI**2
      cosPHI2 = REAL (cexPHI2)
      sinPHI2 = aIMAG(cexPHI2)

      DO j=1,nyp1
        !--- compute clear sky solar radiation ---
        RLAT    = YYY(j)/RADtoDEG
        cexRLAT = cEXP ((0,1)*RLAT)
        cosRLAT = REAL ( cexRLAT )
        sinRLAT = aIMAG( cexRLAT )
        A0 = -15.82 + 326.87*cosRLAT
        A1 =   9.63 - 192.44*sinRLAT
        B1 = - 3.27 + 108.70*sinRLAT
        A2 = - 0.64 -   7.80*REAL(  cexRLAT**2        )
        B2 = - 0.50 +  14.42*REAL( (cexRLAT**2)*cPI18 )
        Q0 = A0 + A1*cosPHI + B1*sinPHI + A2*cosPHI2 + B2*sinPHI2
        !--- compute solar angle ---
        alpha = 90. - abs( YYY(j) - 23.45*sin((julday-82.0)*DAYtoRAD) )
        DO i=it1(j),it2(j)
           cloud     = fact1*cld(i,j,month1)+fact2*cld(i,j,month2) 
           Tocn      = T(i,j,1)
           Tatm      = (5.0*Tocn + 22.0)/6.0
           RaHum     = RH*Qs(Tatm) 
           Qsol(i,j) = (1.0-albedo)*Q0*(1.0 - aCLD*cloud + aALPHA*alpha)
           Qlat(i,j) =    CeL*WNDs(i,j)*( Qs(Tocn) - RaHum)
           Qsen(i,j) = RaCeCp*WNDs(i,j)*(    Tocn  - Tatm ) 
           Qlgw(i,j) = EpsSig*(Tocn+273.0)**3
     &               * ( (Tocn+273.0) * (0.39-0.05*Sqrt(RaHum/RaGam))
     &               *   (1.0 - bCLD*cloud**2) + 4.0*(Tocn - Tatm)    )
           Qtot(i,j) = Qsol(i,j) - Qlat(i,j) - Qsen(i,j) - Qlgw(i,j) 
        END DO
      END DO

      RETURN
 9000 format('(force_q) ',a)
 9001 format('(force_q) ',a,   f6.3)
 9002 format('(force_q) ',a,1pe10.3)
      END

c=======================================================================

      SUBROUTINE force_sst

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! OUTPUT - forcing data goes here


      !----- local -----
      real  SSTc(nxp1,nyp1,12) ! climatological SST data

      integer i,j,k            ! generic indicies
      integer month1,month2    ! two month to interpolate between
      real     fact1, fact2    ! weights for the two months data
      integer  eday1, eday2    ! elapsed days since y-m-d=0-1-1 & 0s
      character*20 loc_fn      ! local      file name
      character*64 mss_fn      ! mass-store file name
      character*16 opts        ! mswrite options
      character*16 pass        ! mswrite password(s)
      integer      ier         ! msread/mswrite error code
      integer      LstIndex    ! last index funct (see intrinsic index)

      logical      first_call  ! flags one-time initializations
      data         first_call / .true. /

      SAVE

c-----------------------------------------------------------------------
c PURPOSE:
c   compute SST(i,j)
c   the instantaneous climatology SST based on 12 monthly means
c-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! initializations
      !-----------------------------------------------------------------
      if ( first_call) then
        write(*,9000) 'initializing climatological SST data'
        mss_fn = frc_tfile
        i      = LstIndex(mss_fn,'/')+1
        j      =    index(mss_fn,' ')
        loc_fn = mss_fn(i:j)//' '
        pass   = ' '
        opts   = ' '
        call MSread(ier, loc_fn , mss_fn, pass, opts)
        write(*,9010) ier,loc_fn(1:j-i+1),mss_fn(1:j)

        open(tmpunit,file=loc_fn,form='UNFORMATTED')
        DO k=1,12
           read(tmpunit) ((SSTc(i,j,k),i=1,nxp1),j=1,nyp1)
        END DO
        close(tmpunit)

        DO j=1,nyp1
        DO i=1,nxp1
           SST(i,j) = SPVAL
        END DO
        END DO
        first_call = .false.
      end if

      !----------------------------------------------------------------
      ! time interpolate between two monthly means
      !----------------------------------------------------------------
      if (day.lt.16) then
         !--- day=1  & sec=0     => fact1=0.5 & fact2=0.5
         !--- day=15 & sec=86400 => fact1=0.0 & fact2=1.0
         month1 = month - 1
         month2 = month
         if (month1.eq.0) month1=12
         fact2 = 0.5 + 0.5*(float(day-1) + sec/86400.0)/15.0
         fact1 = 1.0 - fact2
      else
         !--- day=16 & sec=0                => fact1=1.0 & fact2=0.0
         !--- last day of month & 86400 sec => fact1=0.5 & fact2=0.5
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
         fact1 = 1.0 - 0.5*(float(day-16) + sec/86400.0)/(eday2-eday1)
         fact2 = 1.0 - fact1
      end if

      DO j=1,nyp1
      DO i=it1(j),it2(j)
        SST(i,j) = fact1*SSTc(i,j,month1)+fact2*SSTc(i,j,month2)
      END DO
      END DO

      RETURN
 9000 format('(force_sst) ',a)
 9010 format('(force_sst) return code ',i2,' on  msread ',a,a)
      END
c=======================================================================
