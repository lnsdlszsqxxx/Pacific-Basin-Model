c=======================================================================
c CVS $Id: main.F,v 1.4 1998/08/28 20:46:07 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/main.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      PROGRAM pac3

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local -----
      integer      iStep    ! index for intermediate Lorenz N-cycle step
      integer      eday     ! elapsed days
      real         a,b      ! Lorenz N-cycle coefficients
      integer      i,j,k    ! generic loop indices
      real         sknz     ! sknz=sigma(k)/sigma(nz)
      character*8  dstr     ! real world date string
      character*8  tstr     ! real world time string

c-----------------------------------------------------------------------
c PURPOSE:
c   Impliments a dynamic/thermo-dynamic equitorial ocean simulation 
c   the physics of which are based on a model by M. Cane & P. Gent
c   with significant physics modifications by P. Gent
c   and significant coding/diagnositc modifications by B. Kauffman
c
c HISTORY:
c   pre-SCCS revisions:
c   o EQUD version of 16 Aug 76 by M. CANE
c   o EQUD version of 31 Oct 78 by M. CANE
c   o EQUD version of 05 Dec 85 by P. GENT
c   o EQUD version of 05 Dec 86 by P. GENT
c
c   SCCS revisions circa Jun 90 by B. Kauffman 
c   ocn1 model (aka C&G model), major modifications include:
c   o new separation of code into modules 
c   o new UNICOS compatible file format & I/O routines
c   o new graphics & diagnostic routines
c  
c   SCCS revisions circa Oct 90 by B. Kauffman 
c   ocn2 model (aka G&K model), major modifications include:
c   o non-rectangular domain (Pacific equitorial basin)
c
c   SCCS revisions circa 94/95/96
c   pac3 model (aka ocn3), major modifications include:
c   o new heat flux formulation             (Aug 94 by Brady & Kauffman)
c   o new kpp boundary layer mixing         (Feb 95 by Large & Norton  )
c   o new control logic, calendar, output data (Sep 96 by B. Kauffman)
c   o more uniform use of physical constants   (Nov 96 by B. Kauffman)
c
c   CVS revisions circa Aug 1998 by B. Kauffman
c   o can be coupled to CSM1.2
c
c COPYRIGHT:
c   NCAR/CSM Pacific Basin Model, Copyright (C) 1998
c   University Corporation for Atmospheric Research,
c   All Rights Reserved.
c-----------------------------------------------------------------------

      write(*,9002)
      call datetime(dstr,tstr)
      write(*,9000)'Executing CSM Pacific Basin Model -- ',dstr,' ',tstr
      write(*, *  )
      write(*,9999) ! copyright notice
      write(*,9002)
      write(*,9000) 'CVS ',
     &  '$Id: main.F,v 1.4 1998/08/28 20:46:07 kauff Exp $' 
      write(*,9000) 'CVS ',
     &  '$Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/main.F,v $' 
      write(*,9000) 'CVS $Name: pac3_e $' 
      write(*,9002)

      !-----------------------------------------------------------------
      ! initialization
      !-----------------------------------------------------------------
      call initial1
      call restart('read')
      call initial2

      !--- record initial model date & time ---
      if (sec.le.0 .or. sec.ge.86400 .or. info_dbug.ge.3) then
         write(*,9010) year,month,day,sec,nstep
         call flush(101)
      end if

      !--- initialize control flags ---
      call control

      !-----------------------------------------------------------------
      !              start of MAIN INTEGRATION LOOP
      !-----------------------------------------------------------------
      iStep = 0
      DO WHILE (.true.)

        !--- record wall clock time ---
        if (day.eq.1 .and. sec.eq.0) then
           call datetime(dstr,tstr)
           write(*,9001) 'wall clock time:  ',dstr,'  ',tstr
           call flush(101)
        end if

        !--- Lorenz n-cycle time stepping scheme ---
        a=-iStep/dNT
        b=dNT/(dt_nCyc-iStep)
        DO k=1,nz
        DO j=1,nyp1
        DO i=it1(j),it2(j)
           FU(i,j,k)=FU(i,j,k)*a
           FV(i,j,k)=FV(i,j,k)*a
           FT(i,j,k)=FT(i,j,k)*a
        END DO
        END DO
        END DO
        DO j=1,nyp1
        DO i=it1(j),it2(j)
           FH1 (i,j)=FH1 (i,j)*a
           FHnz(i,j)=FHnz(i,j)*a
        END DO
        END DO

        !--- advance boundary forcing ---
        call force

        !--- calc diagnostic fields ---
        call calcDIV
        call calcW
        call calcZETA

        !--- make restart, diagnostic, history data ---
        if (nstep.gt.0) call restart('write')
        call diagnos
        if (nstep.gt.0 .or. rest_type(1:4).ne.'cont') call history
        call histavg

        !----------------------------------------------------
        ! terminate program OR continue advancing the model?
        !----------------------------------------------------
        IF ( stop_now .ne. 0 ) GO TO 900


        !--- calculate vertical, horizontal, PGF terms ---
        call dVERT
        call dHORIZ
        call dPGF

        !--- advance HU,HV & HT fields ---
        DO k=1,nz
        DO j=1,nyp1
        DO i=it1(j),it2(j)
           FU(i,j,k)=FU(i,j,k)*b
           FV(i,j,k)=FV(i,j,k)*b
           FT(i,j,k)=FT(i,j,k)*b
           HU(i,j,k)=HU(i,j,k)+FU(i,j,k)
           HV(i,j,k)=HV(i,j,k)+FV(i,j,k)
           HT(i,j,k)=HT(i,j,k)+FT(i,j,k)
        END DO
        END DO
        END DO
        !--- advance height fields for k=1 and nz ---
        DO j=1,nyp1
        DO i=it1(j),it2(j)
           FH1 (i,j) = FH1 (i,j)*b
           FHnz(i,j) = FHnz(i,j)*b
           H(i,j, 1) = H(i,j, 1)+FH1 (i,j)
           H(i,j,nz) = H(i,j,nz)+FHnz(i,j)
        END DO
        END DO
        !--- no slip BC's => fix H field in corners ---
        IF (mBCU.EQ.3 .or. mBCU.EQ.2) call Hfix
        !--- filter H field  (iff no lap or bih friction) ---
        if ( nShapH.GT.0 .AND. mod(nStep+1,nShapH).EQ.0) call FilterH
        !--- advance height field for k=2,..,nz-1 ---
        DO k=2,nz-1
           sknz=sigma(k)/sigma(nz)
           DO j=1,nyp1
           DO i=it1(j),it2(j)
              H(i,j,k)=H(i,j,nz)*sknz
           END DO
           END DO
        END DO
        !--- impose boundary conditions ---
        call BCset

        !--- filter HU,HV,HT fields ---
        if ( nShapU.GT.0 .AND. mod(nStep+1,nShapU).EQ.0) call FilterU
        if ( nShapT.GT.0 .AND. mod(nStep+1,nShapT).EQ.0) call FilterT
        !--- form U,V,T from HU,HV,HT ---
        call formUVT
        !--- do vertical mixing (implicit) ---
        IF ( (iStep+1) .EQ. dt_nCyc )  THEN
          call   kpp_coef( U    , V      , T      , H  
     &                   , taux , tauy   , Qsol   , Qtot 
     &                   , nxp1 , nyp1   , nz     , nzp1 , it1  , it2 
     &                   , rhow , cpw    , alphaT , Tbot , gravity
     &                   ,  nu  ,  kappa   
     &                   , rNu  , rKappa , rGamma , bldpth )
           call dRMIX
           call formHUVT
        END IF
        !--- correct H field (non-conservative filter?)---
        if ( dt_hcor.GT.0 .AND. mod(nStep+1,dt_hcor).EQ.0) then
           call correctH
           call formHUVT
        END IF

        !--- advance nCycle index, model step, model time ---
        iStep=mod(iStep+1,dt_nCyc)
        nStep=nStep+1
        if (mod(nStep,dt_ncyc*dt_nday).eq.0 ) then
          call date2eday(cdate,eday)
          eday = eday + 1 
          call eday2date(eday,cdate)
          call eday2ymd (eday,year,month,day)
          sec   = 0
        else
          sec =  dt*mod(nStep,dt_ncyc*dt_nday)
        end if

        !--- record current model date & time ---
        if (sec.le.0 .or. sec.ge.86400 .or. info_dbug.ge.3) then
           write(*,9010) year,month,day,sec,nstep
           call flush(101)
        end if

        !--- update control flags ---
        call control

      END DO
      !-----------------------------------------------------------------
      !              end of MAIN INTEGRATION LOOP
      !-----------------------------------------------------------------
  900 CONTINUE

      write(*,9002)
      call datetime(dstr,tstr)
      write(*,9000)'Terminating Pacific Basin Model -- ',dstr,' ',tstr
      write(*,9002)

      STOP '(main)'
 9000 format(4a)
 9001 format('(main)',4a)
 9002 format(80('='))
 9003 format(80('=')/)
 9004 format(80('-'))
 9005 format(80('-')/)
 9010 format('(main) date=',i4.4,'-',i2.2,'-',i2.2,i6,'s, step=',i6)
 9999 format(
     $ 'NCAR/CSM Pacific Basin Model, Copyright (C) 1998'/,
     $ 'University Corporation for Atmospheric Research, ',
     $ 'All Rights Reserved.')

      END

c=======================================================================
