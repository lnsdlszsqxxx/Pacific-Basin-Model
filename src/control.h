c=======================================================================
c CVS $Id: control.h,v 1.5 1998/08/27 16:49:08 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/control.h,v $
c CVS $Name: pac3_e $
c=======================================================================

      !----- unit number reserved for temporary use --------------------
      integer     tmpunit      ! unit # reserved for temporary usage
      parameter  (tmpunit = 90)

      !----- restart file selection ------------------------------------
      character*32 rest_type   ! run type: initial, cont, regen, branch
      character*64 rest_bfile  ! IC data file for branch runs
      character*64 rest_pfile  ! restart pointer file path & name
      integer      rest_date   ! start date for initial or branch runs
      character*32 rest_freq   ! restart file frequency
      integer      rest_n      ! restart file interval    (nday option)
      integer      rest_odate  ! restart file offset date (nday option)

      !----- case ID info ----------------------------------------------
      character*32 case_name   ! case name
      character*64 case_desc   ! case description

      !----- stopping date selection -----------------------------------
      character*16 stop_option ! specify when to stop models
      integer      stop_date   ! date on which to stop the models
      integer      stop_n      ! number of days to advance the models

      !----- forcing data file selection -------------------------------
      character*64 frc_wfile   ! forcing data file: wind speed vectors
      character*64 frc_cfile   ! forcing data file: cloud cover
      character*64 frc_tfile   ! forcing data file: SST

      !----- history file selection ------------------------------------
      character*16 hist_freq   ! specify when to create hist data
      integer      hist_n      ! history file interval    (nday option)
      integer      hist_odate  ! history file offset date (nday option)
      integer      hist_tavg   ! T => make monthly avg data files

      !----- diagnostic output selection -------------------------------
      character*16 diag_freq   ! specify when to create diag data
      integer      diag_n      ! diagnostic interval    (nday option)
      integer      diag_odate  ! diagnostic offset date (nday option)

      !----- mass store mswrite options --------------------------------
      character*64 mss_dir     ! MSS ouput file directory
      character*16 mss_pass    ! MSS file MS password
      character*16 mss_opts    ! MSS file mswrite options
      integer      mss_rtpd    ! MSS file MS retention period
      integer      mss_rmlf    ! remove local file after an mswrite

      !----- info needed by some msg passing libraries -----------------
      character*8  msg_group   ! msg passing group name
      character*16 msg_method  ! msg passing method choice
      integer      msg_maxt    ! msg passing dt for time out

      !----- dbugging & timing info selection --------------------------
      integer      info_dbug   ! dbug level: > 0 => print extra info
      integer      info_date   ! date when info_dbug > 3 output starts
      integer      info_time   ! T => print extra timing info

      !----- rest, stop, hist, diag control flags ----------------------
      integer      stop_now    ! T => stop model now
      integer      rest_now    ! T => create restart data now
      integer      hist_now    ! T => create history data now
      integer      diag_now    ! T => create diagnostic data now
      integer      stop_eod    ! T => stop model at end-of-day
      integer      rest_eod    ! T => create restart data at eod
      integer      hist_eod    ! T => create history data at eod
      integer      diag_eod    ! T => create diagnostic data at eod

      !----- timestep & model date info --------------------------------
      integer      dt_nCyc     ! order of N-cycle time stepping scheme
      integer      dt_nday     ! number of full N-cycle steps per day
      integer      dt_hcor     ! number of steps between H corrections
      real         dt          !   dt = 86400sec/(dt_nday*dt_ncyc)
      real         dNt         ! N*dt = 86400sec/(dt_nday)
      integer      year        ! model year
      integer      month       ! model month
      integer      day         ! model day
      integer      sec         ! model elapsed seconds during model day
      integer      cdate       ! model year+month+date, coded: yyyymmdd
      integer      nstep       ! model step  (nstep = 0 on restart)

      !----- global error code -----------------------------------------
      integer      ierr        ! 0 <=> no errors or warnings

      !----- boundary condition selection ------------------------------
      integer  mBCU  ,mBCH  ,mBCT     ! toggles fixed/free BC's

      !----- horizontal diffusion filter options  ----------------------
      integer  nOrdrU,nOrdrH,nOrdrT   ! Shapiro filter order for U,H,T
      integer  nShapH,nShapU,nShapT   ! filter time step interval
      integer  mShapH,mShapU,mShapT   ! filter boundary treatment

      !----- global data via common block ------------------------------
      COMMON /control1/
     &   rest_type   , rest_bfile , rest_pfile , 
     &   rest_freq   , rest_n     , rest_odate , rest_date  ,
     &   case_name   , case_desc  ,
     &   stop_option , stop_n     , stop_date  , 
     &   frc_wfile   , frc_cfile  , frc_tfile  , 
     &   hist_freq   , hist_n     , hist_odate , hist_tavg  ,
     &   diag_freq   , diag_n     , diag_odate , 
     &   mss_dir     , mss_pass    , mss_opts   , mss_rtpd   , mss_rmlf,
     &   msg_group   , msg_method , msg_maxt   ,
     &   info_dbug   , info_date  , info_time  ,
     &   stop_now    , rest_now   , hist_now   , diag_now   ,
     &   stop_eod    , rest_eod   , hist_eod   , diag_eod   ,
     &   dt_nCyc     , dt_nday    , dt_hcor    , dt  , dNt  ,
     &   year        , month      , day        , sec        , 
     &   cdate       , nstep      , ierr       ,
     &   mBCU        , mBCH       , mBCT       ,
     &   nOrdrU      ,nOrdrH      , nOrdrT     ,
     &   nShapH      ,nShapU      , nShapT     ,
     &   mShapH      ,mShapU      , mShapT      

      SAVE   /control1/
