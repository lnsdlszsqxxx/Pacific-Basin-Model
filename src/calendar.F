c=======================================================================
c CVS $Id: calendar.F,v 1.2 1998/08/24 22:52:31 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/calendar.F,v $
c CVS $Name: pac3_e $
c=======================================================================
c  Note: 
c   these calendar routines to conversions between
c   o the integer number of elapsed days, 
c   o the integer triple (year,month,day)
c   o the integer coded calendar date: yyyymmdd
c
c  Assumptions:
c   o there is a year 0
c   o all years have 365 days (no leap years)
c   o elapsed days = 0 <=> start of  1 Jan, year 0
c
c=======================================================================

      SUBROUTINE eday2date(eday,date)

      implicit none

      integer  eday,date

      integer  calendar(12),k,year,month,day
      data     calendar /0,31,59,90,120,151,181,212,243,273,304,334/
      save     calendar

c-----------------------------------------------------------------------
c  PURPOSE:
c    compute the calendar date: year/month/day
c  INPUT:
c    an integer number of elapsed days
c  OUTPUT:
c    coded (yyyymmdd) calendar date
c  NOTE:
c    this calendar has a year zero (but no day or month zero)
c-----------------------------------------------------------------------

      year = eday/365       ! calandar year (note: Fortran truncation)
      day  = mod(eday,365)  ! elapsed days within current year
      DO k=1,12
        IF (day .ge. calendar(k)) month=k   ! calendar month
      END DO
      day = day-calendar(month) + 1         ! calendar day

      date = year*10000 + month*100 + day  ! coded calendar date

      RETURN
      END

c=======================================================================

      SUBROUTINE eday2ymd (eday,year,month,day)

      implicit none

      integer  eday,year,month,day

      integer  calendar(12),k
      data     calendar /0,31,59,90,120,151,181,212,243,273,304,334/
      save     calendar

c-----------------------------------------------------------------------
c  PURPOSE:
c    compute the calendar date: year/month/day
c  INPUT:
c    an integer number of elapsed days
c  OUTPUT:
c    uncoded calendar date, integer year, month, & day
c  NOTE:
c    this calendar has a year zero (but no day or month zero)
c-----------------------------------------------------------------------

      year = eday/365       ! calandar year (note: Fortran truncation)
      day  = mod(eday,365)  ! elapsed days within current year
      DO k=1,12
        IF (day .ge. calendar(k)) month=k   ! calendar month
      END DO
      day = day-calendar(month) + 1         ! calendar day

      RETURN
      END

c=======================================================================

      SUBROUTINE date2ymd (date,year,month,day)

      implicit none

      integer  date,year,month,day
      logical  valid_date

c-----------------------------------------------------------------------
c  PURPOSE:
c    decode the calendar date
c  INPUT:
c    calendar date in (integer) yyyymmdd format
c  OUTPUT:
c    calendar year,month,day
c-----------------------------------------------------------------------

      if (.not. valid_date(date)) stop

      year =int(     date       /10000)
      month=int( mod(date,10000)/  100)
      day  =     mod(date,  100) 

      RETURN
      END
c=======================================================================

      SUBROUTINE date2eday(date,eday)

      implicit none

      integer  date,eday
      logical  valid_date

      integer  calendar(12),year,month,day
      data     calendar /0,31,59,90,120,151,181,212,243,273,304,334/
      save     calendar

c-----------------------------------------------------------------------
c  PURPOSE:
c    derive elapsed days from the calendar date
c  INPUT:
c    calendar date in (integer) yyyymmdd format
c  OUTPUT:
c    elapsed days since yy-mm-dd = 00-01-01, with 0 elapsed seconds
c-----------------------------------------------------------------------

      if (.not. valid_date(date)) stop 

      year =int(     date       /10000)
      month=int( mod(date,10000)/  100)
      day  =     mod(date,  100) 

      eday = year*365 + calendar(month) + (day-1)

      RETURN
      END

c=======================================================================

      SUBROUTINE  ymd2date(year,month,day,date)

      implicit none

      integer  year,month,day,date
      logical  valid_ymd

c-----------------------------------------------------------------------
c  PURPOSE:
c    encode the calendar date
c  INPUT:
c    year, month, & date
c  OUTPUT:
c    coded (yyyymmdd) calendar date
c  NOTE:
c    this calendar has a year zero (but no day or month zero)
c-----------------------------------------------------------------------

      if (.not. valid_ymd(year,month,day)) stop 

      date = year*10000 + month*100 + day  ! coded calendar date

      RETURN
      END

c=======================================================================

      SUBROUTINE  ymd2eday(year,month,day,eday)

      implicit none

      integer  year,month,day,eday
      logical  valid_ymd

      integer  calendar(12)
      data     calendar /0,31,59,90,120,151,181,212,243,273,304,334/
      save     calendar

c-----------------------------------------------------------------------
c  PURPOSE:
c    derive elapsed days from the calendar date
c  INPUT:
c    calendara year, month, & date
c  OUTPUT:
c    elapsed days since yy-mm-dd = 00-01-01, with 0 elapsed seconds
c-----------------------------------------------------------------------

      if (.not. valid_ymd(year,month,day)) stop 

      eday = year*365 + calendar(month) + (day-1)

      RETURN
      END

c=======================================================================

      logical FUNCTION valid_date(date)

      implicit none

      integer  date

      integer  year,month,day,dpm(12)
      data     dpm /31,28,31,  30,31,30,  31,31,30,  31,30,31/
      save     dpm ! J  F  M    A  M  J    J  A  S    O  N  D

c-----------------------------------------------------------------------
c  PURPOSE:
c    determine if a valid year, month & day can be decoded 
c    from the coded calendar date
c  INPUT:
c    calendar date in (integer) yyyymmdd format
c  RETURNS:
c    true or false
c-----------------------------------------------------------------------

      year =int(     date       /10000)
      month=int( mod(date,10000)/  100)
      day  =     mod(date,  100) 

      valid_date = .true.
      if (year  .lt.0) valid_date = .false.
      if (month.lt. 1) valid_date = .false.
      if (month.gt.12) valid_date = .false.
      if (day  .lt. 1) valid_date = .false.
      if (.not. valid_date ) return
      if (day .gt. dpm(month)) valid_date = .false.

      RETURN
      END

c=======================================================================

      logical FUNCTION valid_ymd(year,month,day)

      implicit none

      integer  year,month,day,dpm(12)
      data     dpm /31,28,31,  30,31,30,  31,31,30,  31,30,31/
      save     dpm ! J  F  M    A  M  J    J  A  S    O  N  D

c-----------------------------------------------------------------------
c  PURPOSE:
c    determine if a given year, month & day constitute a valid date
c  INPUT:
c    calendar year, month, and day
c  RETURNS:
c    true or false
c-----------------------------------------------------------------------

      valid_ymd = .true.
      if (year  .lt.0) valid_ymd = .false.
      if (month.lt. 1) valid_ymd = .false.
      if (month.gt.12) valid_ymd = .false.
      if (day  .lt. 1) valid_ymd = .false.
      if (.not. valid_ymd ) return
      if (day .gt. dpm(month)) valid_ymd = .false.

      RETURN
      END

c=======================================================================
