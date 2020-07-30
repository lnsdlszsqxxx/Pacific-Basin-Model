c=======================================================================
c CVS $Id: misc.F,v 1.4 1998/08/25 20:46:13 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/misc.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      integer FUNCTION LstIndex(string,substr)

      character*(*)   string,substr
      character       str*240,char*1

c-----------------------------------------------------------------------
c  PURPOSE:
c    find the index of the last occurance of substr within str
c-----------------------------------------------------------------------

      str=' '
      str=string

      if (substr(1:1).eq.'x') char='y'
      if (substr(1:1).ne.'x') char='x'

      LstIndex = index(str,substr)
      i = LstIndex
 100  IF (i .gt. 0) THEN
          LstIndex = i
          str(i:i) = char
          i =  index(str,substr)
          GO TO 100
      END IF

      RETURN
      END

c=======================================================================
!liang
      subroutine datetime(cdate,ctime)

      !--- arguments ---
      character*8 cdate       ! Output:  Current date
      character*8 ctime       ! Output:  Current time

      !--- local ---
      integer mon,day,yr      ! calendar date
      integer hr,min,sec      ! time of day
      integer iarray(3)       ! time of day temporary for SGI
      character*10 temptime    ! Current time

c-----------------------------------------------------------------------
c PURPOSE: returns real world date & time strings
c-----------------------------------------------------------------------

      call date_and_time(DATE=cdate,TIME=temptime) !cdate is yyyymmdd
                                                   !temptime is
                                                   !hhmmss.ms
      ctime(1:6)=temptime(1:6) !hhmmss
      ctime(7:8)=temptime(8:9) !millisecond

      write(*,*) "(datetime) ", cdate,", ",ctime
      return
      end

c=======================================================================
!liang
      SUBROUTINE systemcall(str)

      IMPLICIT none

      !----- arguments ---
      character*(*) str !command to be run

      !----- external routines -----
      integer      ier         ! function return error code
      integer      ishell      ! function to envoke a c-shell command
      external     ishell

c-----------------------------------------------------------------------
c PURPOSE: make a system call (invoke the system command)
c-----------------------------------------------------------------------

      ier=system(str)
      if (ier.ne.0) write(*,*) '(systemcall) ierr=',ier,' for ',str

      RETURN
      END

c=======================================================================
!liang
      integer FUNCTION irtc()

c-----------------------------------------------------------------------
c PURPOSE: a stand-in for Cray/SGI irtc function
c-----------------------------------------------------------------------

      irtc = -999 ! clock ticks since...?

      RETURN
      END

c=======================================================================
!liang

      SUBROUTINE MSwrite (ier,loc_fn,mss_fn,mss_pass,mss_rtpd,mss_opts)

      IMPLICIT none

      integer       ier
      character*(*) loc_fn
      character*(*) mss_fn
      character*(*) mss_pass
      integer       mss_rtpd
      character*(*) mss_opts

c-----------------------------------------------------------------------
c PURPOSE: a stand-in for NCAR's MSS routine
c-----------------------------------------------------------------------

      ier=0
      call systemcall('cp '//loc_fn//' '//mss_fn)

      RETURN
      END

c=======================================================================
!liang

      SUBROUTINE MSread  (ier,loc_fn,mss_fn,mss_pass,mss_opts)

      IMPLICIT none

      integer       ier
      character*(*) loc_fn
      character*(*) mss_fn
      character*(*) mss_pass
      character*(*) mss_opts

c-----------------------------------------------------------------------
c PURPOSE: a stand-in for NCAR's MSS routine
c-----------------------------------------------------------------------

      ier=0
      call systemcall('cp '//mss_fn//' '//loc_fn)

      RETURN
      END

c=======================================================================
