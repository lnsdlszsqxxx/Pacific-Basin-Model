c=======================================================================
c CVS $Id: control.F,v 1.4 1998/08/24 23:12:51 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/control.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE control()

      implicit none

      !---- global ----
      include 'control.h'      ! global time step & control info

      !----- local -----
      integer    true,  false
      parameter (true=1,false=0)
      integer    eday                   ! elapsed days
      integer    eday0                  ! elapsed days at startup
      integer    eday1                  ! elapsed days at end of run
      integer    rest_oeday             ! restart    offset eday
      integer    hist_oeday             ! history    offset eday
      integer    diag_oeday             ! diagnostic offset eday
      integer    first_call             ! true <=> 1st call to routine

      logical   valid_date              ! true <=> yymmdd is valid date

      data      first_call / true /

      save      

c-----------------------------------------------------------------------
c PURPOSE:
c   Set control flags wrt stopping & restart & history data creation.
c
c Note:
c   o when integer values are used to implement boolean logic, 
c     zero <=> false, non-zero <=> true
c   o "elapsed days" means the number of elapsed days since some
c     reference date (normally 00-Jan-01, 0sec) determined by
c     the calendar routines.
c-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! do one-time calculations
      !-----------------------------------------------------------------
      if ( first_call .eq. true ) then

        !-------------------------------------------------------
        ! compute stop_date
        !-------------------------------------------------------
        call date2eday(cdate,eday0)
        IF (stop_option(1:4).eq.'date') THEN
           !-- stop_date is given --
           call date2eday(stop_date,eday1)
        ELSE IF (stop_option(1:5).eq.'ndays') THEN
           !-- stop after n days --
           eday1 = eday0 + stop_n   
           call eday2date(eday1,stop_date)
        ELSE IF (stop_option(1:8).eq.'newmonth') THEN
           !-- stop on day 1 of next month -- 
           if (month.eq.12) call ymd2eday(year+1,      1,1,eday1)
           if (month.ne.12) call ymd2eday(year  ,month+1,1,eday1)
           call eday2date(eday1,stop_date)
        ELSE IF (stop_option(1:7).eq.'newyear') THEN
           !-- stop on January 1 of next year -- 
           call ymd2eday(year+1,1,1,eday1)
           call eday2date(eday1,stop_date)
        ELSE IF (stop_option(1:9).eq.'newdecade') THEN
           !-- stop on January 1 of next decade --
           call ymd2eday(((year/10)+1)*10,1,1,eday1)
           call eday2date(eday1,stop_date)
        ELSE IF (stop_option(1:7).eq.'coupler') THEN
           !-- stop when flux coupler sez so --
           write(*,9000) 'WARNING: model stops on cue from coupler only'
           call ymd2eday(9999,1,1,eday1)
           call eday2date(eday1,stop_date)
        ELSE 
           write(*,*) '(control) invalid stop_option = ',stop_option
           STOP
        END IF
        if ( .NOT. valid_date(stop_date) ) then
          write(*,*) '(control) ERROR: bad stop_date=',stop_date
          stop
        end if
        write(*,9002) eday1-eday0,stop_date

        !-------------------------------------------------------
        ! compute restart file offset eday from offset date
        !-------------------------------------------------------
        if (rest_freq(1:5).eq.'ndays') then
          if (rest_odate .eq. 1) then
            !--- 1st of month & every n days therafter ---
            rest_oeday=-1 
          else
            !--- every n days after an offset date ---
            if (rest_odate .le. 0) rest_odate = cdate
            if (.not. valid_date(rest_odate)) then
               write(*,*) '(control) ERROR: bad rest_odate=',rest_odate
               stop 'rest_odate'
             end if
             call date2eday(rest_odate,rest_oeday)
          end if
        end if

        !-------------------------------------------------------
        ! compute history file offset eday from offset date
        !-------------------------------------------------------
        if (hist_freq(1:5).eq.'ndays') then
          if (hist_odate .eq. 1) then
            !--- 1st of month & every n days therafter ---
            hist_oeday=-1 
          else
            !--- every n days after an offset date ---
            if (hist_odate .le. 0) hist_odate = cdate
            if (.not. valid_date(hist_odate)) then
               write(*,*) '(control) ERROR: bad hist_odate=',hist_odate
               stop 'hist_odate'
             end if
             call date2eday(hist_odate,hist_oeday)
          end if
        end if

        !-------------------------------------------------------
        ! compute diagnostic offset eday from offset date
        !-------------------------------------------------------
        if (diag_freq(1:5).eq.'ndays') then
          if (diag_odate .eq. 1) then
            !--- 1st of month & every n days therafter ---
            diag_oeday=-1 
          else
            !--- every n days after an offset date ---
            if (diag_odate .le. 0) diag_odate = cdate
            if (.not. valid_date(diag_odate)) then
               write(*,*) '(control) ERROR: bad diag_odate=',diag_odate
               stop 'diag_odate'
             end if
             call date2eday(diag_odate,diag_oeday)
          end if
        end if


        !-------------------------------------------------------
        ! warn user about the dbug data they've requested
        !-------------------------------------------------------
        if      (info_dbug.ge.3) then
           write(*,9004) ' WARNING: info_dbug = ',info_dbug
           write(*,9004) ' * This implies a lot of extra text output'
        end if
        if (info_dbug.eq.4) then
           write(*,9004) ' * This implies hist & diag info every day'
           write(*,9004) ' * starting cdate = ',info_date
        else if (info_dbug.ge.5) then
           write(*,9004) ' * This implies hist & diag info every step'
           write(*,9004) ' * starting cdate = ',info_date
        end if

      end if

      !-----------------------------------------------------------------
      ! set stop,rest,hist,diag flags
      !-----------------------------------------------------------------
      stop_now  = false         
      rest_now  = false         
      hist_now  = false         
      diag_now  = false         

      if ( sec .eq. 0 ) then
        call date2eday(cdate,eday) !  today's eday??

        !-------------------------------------------------------
        ! stop now ??
        !-------------------------------------------------------
        if ( cdate .ge. stop_date)  stop_now = true         

        !-------------------------------------------------------
        ! make restart data now ??
        !-------------------------------------------------------
        if ( nstep .eq. 0 ) then
          ! model has just started, don't make a restart file
        else if ( rest_freq(1:6).eq.'monthl' ) then
          if (                  day.eq. 1) rest_now = true
        else if ( rest_freq(1:6).eq.'biweek' ) then
          if (                  day.eq. 1) rest_now = true
          if (                  day.eq.15) rest_now = true
        else if ( rest_freq(1:6).eq.'weekly' ) then
          if (                  day.eq. 1) rest_now = true
          if (                  day.eq. 8) rest_now = true
          if (                  day.eq.15) rest_now = true
          if (                  day.eq.22) rest_now = true
        else if ( rest_freq(1:6).eq.'quarte' ) then
          if (month.eq. 1 .and. day.eq. 1) rest_now = true
          if (month.eq. 4 .and. day.eq. 1) rest_now = true
          if (month.eq. 7 .and. day.eq. 1) rest_now = true
          if (month.eq.10 .and. day.eq. 1) rest_now = true
        else if ( rest_freq(1:6).eq.'yearly' ) then
          if (month.eq. 1 .and. day.eq. 1) rest_now = true
        else if ( rest_freq(1:5).eq.'ndays'   ) then 
          if (rest_oeday.lt.0) then
            if ( mod(day -1,rest_n   ).eq.0) rest_now =true
          else
            if (mod(eday-rest_oeday,rest_n   ).eq.0) rest_now =true
          end if
        else if ( rest_freq(1:7).eq.'coupler' ) then
          if (first_call .eq. true) write(*,9000)
     &       'WARNING: restart files on cue from coupler only' 
        else
          write(*,*) '(control) ERROR: bad rest_freq=',rest_freq
          STOP
        end if

        !-------------------------------------------------------
        ! make history data now ??
        !-------------------------------------------------------
        if      ( hist_freq(1:6).eq.'monthl' ) then
          if (                  day.eq. 1) hist_now = true
        else if ( hist_freq(1:6).eq.'biweek' ) then
          if (                  day.eq. 1) hist_now = true
          if (                  day.eq.15) hist_now = true
        else if ( hist_freq(1:6).eq.'weekly' ) then
          if (                  day.eq. 1) hist_now = true
          if (                  day.eq. 8) hist_now = true
          if (                  day.eq.15) hist_now = true
          if (                  day.eq.22) hist_now = true
        else if ( hist_freq(1:6).eq.'quarte' ) then
          if (month.eq. 1 .and. day.eq. 1) hist_now = true
          if (month.eq. 4 .and. day.eq. 1) hist_now = true
          if (month.eq. 7 .and. day.eq. 1) hist_now = true
          if (month.eq.10 .and. day.eq. 1) hist_now = true
        else if ( hist_freq(1:6).eq.'yearly' ) then
          if (month.eq. 1 .and. day.eq. 1) hist_now = true
        else if ( hist_freq(1:5).eq.'ndays'   ) then 
          if (hist_oeday.lt.0) then
            if ( mod(day -1,hist_n   ).eq.0) hist_now =true
          else
            if (mod(eday-hist_oeday,hist_n   ).eq.0) hist_now =true
          end if
        else if ( hist_freq(1:7).eq.'coupler' ) then
          if (first_call .eq. true) write(*,9000)
     &       'WARNING: history files on cue from coupler only' 
        else
          write(*,*) '(control) ERROR: bad hist_freq=',hist_freq
          STOP
        end if

        !-------------------------------------------------------
        ! make diagnostic data ??
        !-------------------------------------------------------
        if      ( diag_freq(1:6).eq.'monthl' ) then
          if (                  day.eq. 1) diag_now = true
        else if ( diag_freq(1:6).eq.'biweek' ) then
          if (                  day.eq. 1) diag_now = true
          if (                  day.eq.15) diag_now = true
        else if ( diag_freq(1:6).eq.'weekly' ) then
          if (                  day.eq. 1) diag_now = true
          if (                  day.eq. 8) diag_now = true
          if (                  day.eq.15) diag_now = true
          if (                  day.eq.22) diag_now = true
        else if ( diag_freq(1:6).eq.'quarte' ) then
          if (month.eq. 1 .and. day.eq. 1) diag_now = true
          if (month.eq. 4 .and. day.eq. 1) diag_now = true
          if (month.eq. 7 .and. day.eq. 1) diag_now = true
          if (month.eq.10 .and. day.eq. 1) diag_now = true
        else if ( diag_freq(1:6).eq.'yearly' ) then
          if (month.eq. 1 .and. day.eq. 1) diag_now = true
        else if ( diag_freq(1:5).eq.'ndays'   ) then 
          if (diag_oeday.le.0) then
            if (mod(       day - 1 ,diag_n   ).eq.0) diag_now =true
          else
            if (mod(eday-diag_oeday,diag_n   ).eq.0) diag_now =true
          end if
        else
          write(*,*) '(control) ERROR: bad diag_freq=',diag_freq
          STOP
        end if
      end if

      !-----------------------------------------------------------------
      ! info_dbug requests extra data??
      !-----------------------------------------------------------------
      if ( info_dbug.eq.4 .and. cdate.ge.info_date .and. sec.eq.0) then
        diag_now = true
        hist_now = true
      end if
      if ( info_dbug.ge.5 .and. cdate.ge.info_date ) then
        diag_now = true
        hist_now = true
      end if

      !-----------------------------------------------------------------
      ! print out the stop, rest, & hist flags??
      !-----------------------------------------------------------------
      if ( info_dbug .eq. 1 .and. sec.eq.0) then
        if (stop_now .ne. 0) write(*,9000) 'stop_now = true'
        if (rest_now .ne. 0) write(*,9000) 'rest_now = true'
        if (hist_now .ne. 0) write(*,9000) 'hist_now = true'
        if (diag_now .ne. 0) write(*,9000) 'diag_now = true'
      else if ((info_dbug.eq.2 .and. sec.eq.0) .or. info_dbug.ge.3) then
        if (stop_now .ne. 0) write(*,9000) 'stop_now = true'
        if (stop_now .eq. 0) write(*,9000) 'stop_now = false'
        if (rest_now .ne. 0) write(*,9000) 'rest_now = true'
        if (rest_now .eq. 0) write(*,9000) 'rest_now = false'
        if (hist_now .ne. 0) write(*,9000) 'hist_now = true'
        if (hist_now .eq. 0) write(*,9000) 'hist_now = false'
        if (diag_now .ne. 0) write(*,9000) 'diag_now = true'
        if (diag_now .eq. 0) write(*,9000) 'diag_now = false'
      end if

      first_call = false

      RETURN
 9000 format('(control) ',2a)
 9002 format('(control) model will advance ',i4,' days, '
     &       'stopping on ',i8.8)
 9004 format('(control) ',a,i8.8)
      END

c=======================================================================
