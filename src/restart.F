c=======================================================================
c CVS $Id: restart.F,v 1.5 1998/08/25 20:47:29 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/restart.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE restart (cOption)

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- arguments -----
      character*(*) cOption    ! selects read xor write option

      !----- local -----
      character*80 str         ! restart file ID string
      character*8  dstr        ! real world date string
      character*8  tstr        ! real world time string
      character*20 loc_fn      ! local      file name
      character*64 mss_fn      ! mass-store file name
      character*16 opts        ! mswrite options
      character*16 pass        ! mswrite password(s)
      integer      rtpd        ! mswrite retention period
      integer      ier         ! msread/mswrite error code
      integer      i,j         ! character string indicies
      integer      LstIndex    ! last index funct (see intrinsic index)
      logical      file_exists ! a local restart file exists
      integer      mx,my,mz    ! dimensions of data on restart file

      DATA     file_exists /.false. /

      SAVE     loc_fn,file_exists

c-----------------------------------------------------------------------
c PURPOSE:
c   IF (cOption='read')
c     Acquire ocn restart/IC data (read a data file xor fabricate)
c     Restart data is the model output and daily avg fluxes
c     computed during the previous day of integration.
c   IF (cOption='write')
c     AND  the flag rest_now is non-zero (ie. is "true"), 
c     THEN create restart & pointer files & send them to the mass store
c
c DETAILS 
c    wrt cOption='read':
c    IF (rest_type = 'initial' )
c       IC's are fabricated by this routine
c       initial date is specified by input parameter "rest_date"
c    IF (rest_type = 'continue' )
c       IC's are in the data file pointed to by the pointer file
c       initial date is specified by the IC data file
c    IF (rest_type = 'branch' )
c       IC's are in the file specified by "rest_bfile"
c       if   rest_date > 0, initial date is specified by "rest_date"
c       else initial date is specified by the IC data file
c    IF (rest_type = 'regeneration' )
c       IC's are in the file specified by "rest_bfile"
c       initial date is specified by the IC data file
c
c MACHINE DEPENDENCY:
c    The "mass storage system" (MS) is site specific to NCAR
c-----------------------------------------------------------------------
      
      !-----------------------------------------------------------------
      !  read a restart file
      !-----------------------------------------------------------------
      IF (cOption .eq. 'read') THEN

         !----- get IC file name -----
         IF (rest_type(1:4) .eq. 'cont') THEN
            write(*,9000) 'this is a continuation run...'
            write(*,9000) 'IC data file specified by pointer file'
            if (rest_pfile(1:9).ne.'rpointer ') then
              i = index(rest_pfile,' ')-1
              call system(  'cp '//rest_pfile(1:i)//' rpointer')
              write(*,9000) 'cp '//rest_pfile(1:i)//' rpointer'
            end if
            open (tmpunit,file='rpointer',form='formatted')
            read (tmpunit,'(a)') mss_fn
            close(tmpunit)
         ELSE IF (rest_type(1:4) .eq. 'init') THEN
            write(*,9000) 'this is an initial run...'
            write(*,9000) 'IC data file specified by input parameter'
            mss_fn=rest_bfile
         ELSE IF (rest_type(1:4) .eq. 'bran') THEN
            write(*,9000) 'this is a branch run...'
            write(*,9000) 'IC data file specified by input parameter'
            mss_fn=rest_bfile
         ELSE IF (rest_type(1:4) .eq. 'rege') THEN
            write(*,9000) 'this is a regeneration run...'
            write(*,9000) 'IC data file specified by input parameter'
            mss_fn=rest_bfile
         ELSE
            write(*,*) '(restart) invalid rest_type = ',rest_type
            STOP 'restart'
         END IF

         !----- acquire IC data -----
         write(*,9000) 'IC data file      : ',mss_fn
         i      = LstIndex(mss_fn,'/')+1
         j      =    index(mss_fn,' ')
         loc_fn = mss_fn(i:j)//' '
         pass   = ' '
         opts   = ' '
         call MSread(ier, loc_fn , mss_fn, pass, opts)
         write(*,9004) ier,loc_fn(1:j-i+1),mss_fn(1:j)
         IF (ier.ne.0) write(*,9000) 'WARNING: msread return code'
         file_exists = .true. 
         str = ' '
         open (tmpunit,file=loc_fn,form='unformatted')
         read (tmpunit) str 
         read (tmpunit) cdate
         read (tmpunit) mx,my,mz
         write(*,9000      )           'IC file ID string : ',str
         write(*,'(a,i8.8)') '(restart) IC file date      : ',cdate
         if (mx.ne.nx .or. my.ne.ny .or. mz.ne.nz  ) then
            write(*,*) '(restart) fatal error: incompatible IC file'
            write(*,*) '(restart) input nx,ny,nz = ',mx,my,mz
            write(*,*) '(restart) model nx,ny,nz = ',nx,ny,nz
            stop
         end if
         read (tmpunit) U,V,H,T
         read (tmpunit) XXX,YYY,F
         read (tmpunit) EMx,EMy,EMxy
         read (tmpunit) it1,it2,jt1,jt2
         close(tmpunit)

      !-----------------------------------------------------------------
      !  write a restart file
      !-----------------------------------------------------------------
      ELSE IF (cOption .eq. 'write') THEN

         IF (rest_now .ne. 0) THEN
            !--- remove previous file ---
            IF ( file_exists ) call system('rm '//loc_fn(1:14))
            !--- create a new file name ---
            !oc_fn='123456789+1 '
            loc_fn='ryyyy-mm-dd '
            write(loc_fn( 2: 5),'(i4.4)') year
            write(loc_fn( 7: 8),'(i2.2)') month
            write(loc_fn(10:11),'(i2.2)') day
            !--- create a new file ---
            call datetime(dstr,tstr)
            write(str,9005) dstr,tstr,case_desc
            write(*,'(2a)') '(restart) create restart file : ',loc_fn
            write(*,'(2a)') '(restart) restart file ID str : ',str
            open (tmpunit,file=loc_fn,form='unformatted')
            write(tmpunit) str 
            write(tmpunit) cdate
            write(tmpunit) nx,ny,nz
            write(tmpunit) U,V,H,T
            write(tmpunit) XXX,YYY,F
            write(tmpunit) EMx,EMy,EMxy
            write(tmpunit) it1,it2,jt1,jt2
            close(tmpunit)
            file_exists = .true. 

            !--- write new restart data & pointer files to MS ---
            i      = index(mss_dir,' ')-1
            mss_fn = mss_dir(1:i)//loc_fn//' '
            pass   = mss_pass
            rtpd   = mss_rtpd
            opts   = mss_opts
            call MSwrite ( ier, loc_fn, mss_fn, pass, rtpd, opts )
            i=index(mss_fn,' ')-1
            j=index(loc_fn,' ')-1
            write(*,9002) ier,loc_fn(1:j),mss_fn(1:i)
            IF (ier .eq. 0) THEN
              open (tmpunit,file='rpointer',form='formatted')
              write(tmpunit,'(2a)') mss_fn
              write(tmpunit,'(2a)')'restart pointer for case ',case_name
              write(tmpunit,9005  ) dstr,tstr,case_desc
              close(tmpunit)
              if (rest_pfile(1:9).ne.'rpointer ') then
                i = index(rest_pfile,' ')-1
                call system(  'cp -f rpointer '//rest_pfile(1:i))
                write(*,9000) 'cp -f rpointer '//rest_pfile(1:i)
              end if
            ELSE
              write(*,9000) 'ERROR on mswrite?'
            END IF

         END IF

      !-----------------------------------------------------------------
      !  there aren't any other choices
      !-----------------------------------------------------------------
      ELSE
         write(*,*) '(restart) ERROR: invalid cOption = ',cOption
      END IF

      RETURN
 9000 format('(restart) ',2a)
 9002 format('(restart) return code ',i4, ' on mswrite ',a,' ',a)
 9004 format('(restart) return code ',i4, ' on msread  ',a,' ',a)
 9005 format('(',a,1x,a,') ',a60)
      END

C=======================================================================
