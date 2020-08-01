c=======================================================================
c CVS $Id: initial.F,v 1.6 1998/08/25 20:47:29 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/initial.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE initial1()

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local -----
      integer i                ! generic index

      !----- input namelist -----
      NAMELIST /inparm/  
     &   rest_type   ,
     &   rest_bfile  ,
     &   rest_pfile  ,
     &   rest_date   ,
     &   rest_freq   ,
     &   rest_n      ,
     &   rest_odate  ,
     &   case_name   ,
     &   case_desc   ,
     &   stop_option ,
     &   stop_n      ,
     &   stop_date   ,
     &   hist_freq   ,
     &   hist_n      ,
     &   hist_odate  ,
     &   hist_tavg   ,
     &   diag_freq   ,
     &   diag_n      ,
     &   diag_odate  ,
     &   frc_wfile   ,
     &   frc_cfile   ,
     &   frc_tfile   ,
     &   mss_dir     ,
     &   mss_opts    ,
     &   mss_pass    ,
     &   mss_rtpd    ,
     &   mss_rmlf    ,
     &   msg_group   ,
     &   msg_method  ,
     &   msg_maxt    ,
     &   info_dbug   ,
     &   info_date   ,
     &   info_time   ,
     &   dt_ncyc     ,
     &   dt_nday     ,
     &   dt_hcor     ,
     &   Tin  , Hin  , 
     &   kappa, nu   ,
     &   mBCH , mBCU , mBCT  ,
     &   nOrdrH,nOrdrU,nOrdrT,
     &   nShapH,nShapU,nShapT,mShapH,mShapU,mShapT

      DATA Tin   / 1.00,  0.98,  0.94,  0.91,  0.86,
     &             0.81,  0.75,  0.69,  0.62,  0.54,
     &             0.46,  0.37,  0.28,  0.18,  0.06/
      DATA Hin   /10.00, 10.00, 15.00, 15.00, 20.00,
     &            20.00, 25.00, 25.00, 30.00, 30.00,
     &            35.00, 35.00, 40.00, 40.00, 50.00/
      DATA kappa /15*1.0e-5/
      DATA nu    /15*1.0e-4/

c-----------------------------------------------------------------------
c  PURPOSE:
c    o read input parameter namelist
c    o read input parameters
c-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! set input parm defaults (see also "data" statement)
      !-----------------------------------------------------------------
      rest_type   = 'continue'
      rest_bfile  = '/USER/invalid/rest_bfile '
      rest_pfile  = 'null'
      rest_date   = 00000101
      rest_freq   = 'monthly'
      rest_n      = 30
      rest_odate  = -1
      case_name   = 'null'
      case_desc   = 'null'
      stop_option = 'newmonth'
      stop_n      = 365
      stop_date   = 00010101
      hist_freq   = 'monthly'
      hist_n      = 30
      hist_odate  =  1
      hist_tavg   =  0
      diag_freq   = 'monthly'
      diag_n      = 90
      diag_odate  = -1
      info_dbug   =  1
      info_date   =  00000101
      info_time   =  0
      frc_wfile   = '/USER/invalid/winddata '
      frc_cfile   = '/USER/invalid/clouddata '
      frc_tfile   = '/USER/invalid/sstdata '
      mss_dir     = '/USER/invalid/mss_dir/ '
      mss_opts    = 'nowait,nomail '
      mss_pass    = ' '
      mss_rtpd    = 365
      mss_rmlf    =  1
      msg_group   = 'unknown'
      msg_method  = 'unknown'
      msg_maxt    = 300
      dt_nCyc     =  4  ! for N-cycle scheme, n = dt_nCyc
      dt_nday     =  8  ! # of intermediate steps/day = dt_nCyc*dt_nday
      dt_hcor     = 64  ! intermediate steps between H corrections
      mBCH        =  0  ! BC's for H,U,T fields
      mBCU        =  3
      mBCT        =  1
      nOrdrH      = 16  ! order of Shapiro filter
      nOrdrU      = 12
      nOrdrT      = 16
      nShapH      =  8  ! frequency of filtering
      nShapU      =  8
      nShapT      =  8
      mShapH      =  4  ! filter behaviour at boundaries
      mShapU      =  7
      mShapT      =  5

      !-----------------------------------------------------------------
      ! read input parms from stdin
      !-----------------------------------------------------------------
      write(*,9000)   'reading parameters from input namelist: ocn.parm'
      read (*,nml=inparm)

      !--- MS dir name needs trailing slash ---
      i = index(mss_dir,' ')-1
      if (mss_dir(i:i).ne.'/') mss_dir(i+1:i+2) = '/ '

      !----- case name --------------------------------------
      if ( case_name (1:4).eq.'null') then
        write(*,9000) 'ERROR: must input case_name'
        stop
      end if
      if ( case_desc (1:4).eq.'null') then
          write(*,9000) 'WARNING: should input case_desc'
          case_desc = ' '
      end if

      !----- restart options --------------------------------
      if ( rest_pfile(1:4).eq.'null') then
        if ( case_name(1:4).eq.'null') then
          rest_pfile = 'rpointer'
        else
          i = index(case_name,' ')-1
          rest_pfile = '$HOME/pac3.'//case_name(1:i)//'.rpointer'
        endif
      else if ( rest_pfile(1:11).eq.'./rpointer ') then
        rest_pfile = 'rpointer'
      endif

      !-----------------------------------------------------------------
      ! print values of all input parms
      !-----------------------------------------------------------------
      write(*,9000)   'namelist values are now:'
      write(*,nml=inparm)
      write(*,9001)

      RETURN
 9000 format('(initial1) ',3a)
 9001 format(79('-'))
      END

c=======================================================================

      SUBROUTINE initial2()

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local -----
      integer i,j,k          ! loop indicies
      real    hsum           ! sum of H(k) for k=2,nz
      logical valid_date     ! function to check validity of coded date

c-----------------------------------------------------------------------
c PURPOSE
c  o reset starting date (depending on rest_type & rest_date)
c  o initialize dt & dNt (based on dt_ncyc & dt_nday)
c  o initialize stretched horizontal grid
c  o initialize various derived fields & coefficients
c
c REFERENCE:
c   Gent & Cane (1989), pg 457, eq 70
c   Gent & Cane (1989), pg 458, eq 73,74
c
c ASSUMPTIONS:
c   input parameter namelist and IC/restart data files have been read in
c-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! select initial date, sec, nstep
      !-----------------------------------------------------------------
      if (rest_type(1:4).eq.'init' .OR. rest_type(1:4).eq.'bran') then
        if ( rest_date .gt. 0 ) cdate = rest_date
      end if
      if ( .not. valid_date(cdate) ) then
         write(*,*) '(initial2) invalid initial date ',cdate
         STOP
      else
         write(*,'(a,i8.8)') '(initial2) using initial date: ',cdate
      end if
      call date2ymd (cdate,year,month,day)
      nstep=0
      sec  =0

      !-----------------------------------------------------------------
      ! initialize dt & dNt based on the number of steps per day 
      !-----------------------------------------------------------------
      write(*,9002) 'higher order steps per day: ',dt_nday
      write(*,9002) 'intermediate steps per day: ',dt_nday*dt_nCyc
      dNt = 86400.0/float(dt_nday)
      dt  = 86400.0/float(dt_nday*dt_ncyc)

      !-----------------------------------------------------------------
      ! initialize stretched horizontal grid
      !-----------------------------------------------------------------
      write(*,9001)
      call domain

      !-----------------------------------------------------------------
      ! initialize various dervived fields & coefficients
      !-----------------------------------------------------------------

      !--- W & zeta (zeta is stored in W(i,j,nz) ) ---
      do k=1,nz
      do j=1,nyp1
      do i=1,nxp1
         W(i,j,k) = SPVAL
      end do
      end do
      end do

      !--- dH(1)/dt=0 (for fixed depth mix layer) ---
      do j=1,nyp1
      do i=1,nxp1
         dH1   (i,j) = 0.0
         BLdpth(i,j) = SPVAL
      end do
      end do

      !--- sigma(k) ---
      Hsum=0.0
      do k=2,nz
         Hsum=Hsum+Hin(k)
      end do
      do k=1,nz
         sigma(k)=Hin(k)/Hsum
      end do

      !--- nu & kappa related coeffs ---
      do k=1,nz-1
        cNU   (k)=   nu(k)/( 0.5*(Hin(k)+Hin(k+1)) )
        cKAPPA(k)=kappa(k)/( 0.5*(Hin(k)+Hin(k+1)) )
      end do
      cNU   (nz)=   nu(nz)/( 0.5*Hin(nz) )
      cKAPPA(nz)=kappa(nz)/( 0.5*Hin(nz) )

      !--- form (hu,hv) & hT from (u,v) & T ---
      call formHUVT

      write(*,9001)
      RETURN
 9001 format(79('-'))
 9002 format('(initial2) ',a,i5,a)
      END

c=======================================================================
