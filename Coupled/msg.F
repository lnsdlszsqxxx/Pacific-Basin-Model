c=======================================================================
c CVS $Id: msg.F,v 1.1 1998/08/28 20:52:56 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/Coupled/msg.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE msg_pass(cOption)

      implicit none

      !----- arguments -----
      character*(*) cOption

      !----- global -----
#include <global.h>             ! data sent/received from coupler
#include <control.h>            ! control info, date/time
#include <msgid.h>              ! has task & msg IDs

      !----- local -----
      integer    nsnd           ! number of 2d fields to send
      integer    nrcv           ! number of 2d fields to receive
      parameter (nsnd  =  6)
      parameter (nrcv  =  5)
      real    buffs(mx,my,nsnd) ! contiguous array for sending data
      real    buffr(mx,my,nrcv) ! contiguous array for recving data

      integer   i,j             ! generic loop indicies

      integer*8   irtc          ! clock tics since last reboot
      integer*8   irtc_w        ! rtc tics when waiting for msg
      integer*8   irtc_r        ! rtc tics when msg received
      integer*8   irtc_u        ! rtc tics when msg unpacked
      integer*8   irtc_p        ! rtc tics when packing msg
      integer*8   irtc_s        ! rtc tics when sending msg 
      
c-----------------------------------------------------------------------
c PURPOSE:
c   implement various message passing tasks using the MPI library
c
c FUCTIONALITY:
c   cOption=
c     'connect'    => connect to the msg passing system
c     'init'       => recv/send initial info
c     'send'       => send model state & output fluxes
c     'recv'       => receive input fluxes & stop/continue flags
c     'disconnect' => disconnect from the msg passing system
c-----------------------------------------------------------------------

      IF (cOption .eq. 'connect') THEN 
        !---------------------------------------------------------------
        ! connect to msg passing system
        !---------------------------------------------------------------
        write(*,9000) '(connect) connecting to cpl'
        call msg_init
        
      ELSE IF (cOption .eq. 'init') THEN 
        !---------------------------------------------------------------
        ! recv/send initial data
        !---------------------------------------------------------------
        write(*,9000) '(initial) recv msg'

        call msg_recv_i(ibuffr,nibuff    ,MSGID_TID_CPL,MSGID_C2O_INIT)

        ierr      = ibuffr( 1) ! error code
       !cdate     = ibuffr( 4) ! current date
       !sec       = ibuffr( 5) ! elapsed sec on current date
        info_time = ibuffr(11) ! T => turn on msg-passing timing
       !info_dbug = ibuffr(12) ! cpl debug level

        !--- send grid & domain mask ---
        write(*,9000) '(initial) send msg'

        do i=1,nibuff
          ibuffs(i)=0     ! unused
        end do
        ibuffs( 1)=ierr   ! error code
        ibuffs( 7)=mx     ! number of grid points in x-direction
        ibuffs( 8)=my     ! number of grid points in y-direction
        ibuffs( 9)=1      ! number of msg-pairs per day
        
        call msg_send_i(ibuffs,nibuff    ,MSGID_TID_CPL,MSGID_O2C_INIT)
        call msg_send_r(  xc_g,mx        ,MSGID_TID_CPL,MSGID_O2C_INIT)
        call msg_send_r(  yc_g,my        ,MSGID_TID_CPL,MSGID_O2C_INIT)
        call msg_send_r(  xe_g,mx+1      ,MSGID_TID_CPL,MSGID_O2C_INIT)
        call msg_send_r(  ye_g,my+1      ,MSGID_TID_CPL,MSGID_O2C_INIT)
        call msg_send_i(mask_g,mx*my     ,MSGID_TID_CPL,MSGID_O2C_INIT)
       
      ELSE IF (cOption .eq. 'send') THEN 
        !---------------------------------------------------------------
        ! send output field msg
        !---------------------------------------------------------------
        if (info_dbug.ge.3) write(*,9000) '(send) send data to cpl'

        !--- set info buffer flags ---
        irtc_p = irtc()   ! packing msg
        do i=1,nibuff
          ibuffs(i)= 0     ! unused
        end do
        ibuffs( 1) = ierr  ! error code
        ibuffs( 3) = 0     ! stop flag: 0 <=> able to continue
        ibuffs( 4) = cdate ! model date, coded: yyyymmdd
        ibuffs( 5) = sec   ! elapsed seconds on model date
        ibuffs( 7) = mx    ! number of grid points in x-direction
        ibuffs( 8) = my    ! number of grid points in y-direction
        ibuffs( 9) = 1     ! number of msg-pairs per day
        ibuffs(31) = 0     ! precip adjustment factor (* 1.0e+6)
        ibuffs(33) = 0     ! number of runoff basins in ocn model

        do j=1,my
        do i=1,mx
           buffs(i,j, 1) =    T_g(i,j)  ! temperature
           buffs(i,j, 2) =    U_g(i,j)  ! velocity, zonal
           buffs(i,j, 3) =    V_g(i,j)  ! velocity, meridional
           buffs(i,j, 4) = dhdx_g(i,j)  ! surface slope, zonal
           buffs(i,j, 5) = dhdy_g(i,j)  ! surface slope, meridional
           buffs(i,j, 6) =    q_g(i,j)  ! heat of fusion xor melt pot
        end do
        end do

        irtc_s = irtc()   ! sending msg
        call msg_send_i(ibuffs,nibuff    ,MSGID_TID_CPL,MSGID_O2C)
        call msg_send_r( buffs,mx*my*nsnd,MSGID_TID_CPL,MSGID_O2C)

        if (info_time.ne.0) then
           write(*,9099) irtc_p,'o->c packing'
           write(*,9099) irtc_s,'o->c sending'
        end if

      ELSE IF (cOption .eq. 'recv') THEN
        !---------------------------------------------------------------
        ! recv input field msg
        !---------------------------------------------------------------
        if (info_dbug.ge.3) write(*,9000) '(recv) recv data to cpl'
        
        irtc_w = irtc()   ! waiting for msg
        call msg_recv_i(ibuffr,nibuff    ,MSGID_TID_CPL,MSGID_C2O)
        irtc_r = irtc()   ! received msg
        call msg_recv_r( buffr,mx*my*nrcv,MSGID_TID_CPL,MSGID_C2O)
        call msg_recv_r(runoff,100       ,MSGID_TID_CPL,MSGID_C2O)
        
        do j=1,my
        do i=1,mx
            taux_g(i,j) = buffr(i,j, 1) ! surface stress, zonal
            tauy_g(i,j) = buffr(i,j, 2) ! surface stress, merid
           netsw_g(i,j) = buffr(i,j, 3) ! net sw rad
           lslwm_g(i,j) = buffr(i,j, 4) ! lat+sen+(net lw)+melt
             pem_g(i,j) = buffr(i,j, 5) ! P+E+melting
        end do
        end do

        irtc_u = irtc()   ! unpacked msg

        ierr      = ibuffr( 1) ! cpl error code (0 => no error)
        stop_eod  = ibuffr( 2) ! T => stop at end-of-day
        stop_now  = ibuffr( 3) ! T => stop now
        !           ibuffr( 4) ! current date
        !           ibuffr( 5) ! elapsed sec on current date
        !           ibuffr(11) ! T => turn on msg-passing timing
        !           ibuffr(12) ! cpl debug level
        rest_eod  = ibuffr(21) ! T => create restart data at end-of-day
        rest_now  = ibuffr(22) ! T => create restart data now
        hist_eod  = ibuffr(23) ! T => create history data at end-of-day
        hist_now  = ibuffr(24) ! T => create history data now
        hist_tavg = ibuffr(25) ! T => create create monthly avg data
        diag_eod  = ibuffr(26) ! T => create diagnotic data at eod
        diag_now  = ibuffr(27) ! T => create diagnotic data now
        !           ibuffr(33) ! number of runoff basins in lnd model

        if (info_time.ne.0 .and. stop_now.eq.0) then
           write(*,9099) irtc_w,'c->o waiting'
           write(*,9099) irtc_r,'c->o received'
           write(*,9099) irtc_u,'c->o unpacked'
        end if

      ELSE IF (cOption .eq. 'disconnect') THEN 
        !---------------------------------------------------------------
        ! disconnect from msg passing
        !---------------------------------------------------------------
        write(*,9000) '(disconnect) disconnect from cpl'
        call msg_finalize

      ELSE 
        !---------------------------------------------------------------
        ! this should never happen
        !---------------------------------------------------------------
        write(*,9000) 'ERROR: invalid cOption: ',cOption
        stop
      END IF

      RETURN
 9000 format('(msg_pass) ',3a)
 9099 format('[mp timing]  irtc = ',i20,' ',a)
      END

c=======================================================================
