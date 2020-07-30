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

      subroutine datetime(cdate,ctime)

      !--- arguments ---
      character*8 cdate       ! Output:  Current date
      character*8 ctime       ! Output:  Current time

      !--- local ---
      integer mon,day,yr      ! calendar date
      integer hr,min,sec      ! time of day
      integer iarray(3)       ! time of day temporary for SGI

c-----------------------------------------------------------------------
c PURPOSE: returns real world date & time strings
c-----------------------------------------------------------------------

#if   ( defined CRAY )
      write(cdate,'(a8)') date()
      write(ctime,'(a8)') clock()
#endif
#if ( defined SGI )
      call idate(mon,day,yr)
      yr = mod(yr,100)
      call itime(iarray)
      hr  = iarray(1)
      min = iarray(2)
      sec = iarray(3)
      write(cdate,'(i2.2,a1,i2.2,a1,i2.2)') mon,'/',day,'/',yr
      write(ctime,'(i2.2,a1,i2.2,a1,i2.2)') hr ,':',min,':',sec
#endif
#if ( defined SX4 )
      call datim(cdate,ctime)
#endif
#if ( !defined CRAY && !defined SGI && !defined SX4 )
      error wrt cpp: must define CRAY, SGI, or SX4
#endif

      return
      end

c=======================================================================
#if ( defined CRAY )

      SUBROUTINE system(str)

      IMPLICIT none

      !----- arguments ---
      character*(*) str

      !----- external routines -----
      integer      ier         ! function return error code
      integer      ishell      ! function to envoke a c-shell command
      external     ishell

c-----------------------------------------------------------------------
c PURPOSE: make a system call
c-----------------------------------------------------------------------

      ier=ishell(str)
      if (ier.ne.0) write(*,*) '(system) ierr=',ier,' for ',str

      RETURN
      END

#endif
c=======================================================================
#if ( !defined CRAY && !defined SGI )

      integer FUNCTION irtc()

c-----------------------------------------------------------------------
c PURPOSE: a stand-in for Cray/SGI irtc function
c-----------------------------------------------------------------------

      irtc = -999 ! clock ticks since...?

      RETURN
      END

#endif
c=======================================================================
#if (!defined MSS)

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
      call system('cp '//loc_fn//' '//mss_fn)

      RETURN
      END

#endif
c=======================================================================
#if (!defined MSS)

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
      call system('cp '//mss_fn//' '//loc_fn)

      RETURN
      END

#endif
c=======================================================================
