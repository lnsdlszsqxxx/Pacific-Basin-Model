c=======================================================================
c CVS $Id: diagnos.F,v 1.2 1998/08/24 22:45:46 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/diagnos.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE diagnos

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local variables -----
      real     eTO(0:nz)       ! total     energy
      real     eKE(0:nz)       ! kinetic   energy
      real     ePE(0:nz)       ! potential energy
      real     eHC(0:nz)       ! heat content
      real     eHV(0:nz)       ! heat variance
      real    Umin(0:nz), Umax(0:nz) ! min/max u velocity
      real    Vmin(0:nz), Vmax(0:nz) ! min/max v velocity
      real    Hmin(0:nz), Hmax(0:nz) ! min/max layer depth
      real    Tmin(0:nz), Tmax(0:nz) ! min/max T 
      real    Bmin      , Bmax       ! min/max BL depth
      real    D2(nxp1,nyp1)    ! 2*(total depth at middle of layer k)
      integer i,j,k            ! generic loop indicies
      real    ZZ               ! area*zeta^2
      real    diffT            ! T(i,j,k)-Tbot
      real    volume           ! H*area
      real    Tavg             ! spatial average of T

c-----------------------------------------------------------------------
c  PURPOSE:
c    print diagnostics: energies & hi/low values of U,V,H,T
c-----------------------------------------------------------------------

      if (diag_now.eq.0) RETURN

      write(*,9000) cdate,sec,nstep

      !-----------------------------------------------------------------
      ! this diagnostic is to check bit-for-bit duplication of a run
      !-----------------------------------------------------------------
      Tavg=0.0
      DO j=1,nyp1
      DO i=it1(j),it2(j)
         Tavg=Tavg+T(i,j,1)*area(i,j)
      END DO
      END DO
      Tavg=Tavg/basin
      write(*,*)'(diagnos) bit-check: avg T(i,j,1)=',Tavg

      !-----------------------------------------------------------------
      ! compute global energies
      !-----------------------------------------------------------------
      eKE(0)=0.0
      ePE(0)=0.0
      eHC(0)=0.0
      eHV(0)=0.0
      DO k=1,nz
        IF (k .eq. 1) THEN
          ZZ=0.0
          DO j=1,nyp1
          DO i=it1(j),it2(j)
            ZZ     =ZZ+area(i,j)*W(i,j,nz)**2
            D2(i,j)=H(i,j,1)
          END DO
          END DO
        ELSE
          DO j=1,nyp1
          DO i=it1(j),it2(j)
             D2(i,j)=D2(i,j)+H(i,j,k-1)+H(i,j,k)
          END DO
          END DO
        END IF
        eKE(k)=0.0
        ePE(k)=0.0
        eHC(k)=0.0
        eHV(k)=0.0
        DO j=1,nyp1
        DO i=it1(j),it2(j)
          volume=area(i,j)*H(i,j,k)
          diffT =T(i,j,k)-Tbot
          eKE(k)=eKE(k)+volume*(U(i,j,k)**2 + V(i,j,k)**2 )
          ePE(k)=ePE(k)+volume*diffT*D2(i,j)
          eHC(k)=eHC(k)+volume*diffT
          eHV(k)=eHV(k)+volume*diffT**2
        END DO
        END DO
        eKE(k)=0.5*eKE(k)/basin
        ePE(k)=0.5*ePE(k)*alphaT*gravity/basin
        eHC(k)=    eHC(k)/basin
        eHV(k)=    eHV(k)/basin
        eTO(k)=eKE(k)+ePE(k)
        !-------------------
        eKE(0)=eKE(0)+eKE(k) 
        ePE(0)=ePE(0)+ePE(k)
        eHC(0)=eHC(0)+eHC(k)
        eHV(0)=eHV(0)+eHV(k)
      END DO
      ePE(0)=ePE(0)-0.5*gravity*ZZ/basin
      eTO(0)=eKE(0)+ePE(0)

      !-----------------------------------------------------------------
      ! locate the min/max values
      !-----------------------------------------------------------------
      DO k=1,nz
        Umin(k) =  1.e30
        Umax(k) = -1.e30
        Vmin(k) =  1.e30
        Vmax(k) = -1.e30
        Hmin(k) =  1.e30
        Hmax(k) = -1.e30
        Tmin(k) =  1.e30
        Tmax(k) = -1.e30
        DO j=1,nyp1
        DO i=it1(j),it2(j)
          Umin(k) = min(Umin(k),U(i,j,k))
          Umax(k) = max(Umax(k),U(i,j,k))
          Vmin(k) = min(Vmin(k),V(i,j,k))
          Vmax(k) = max(Vmax(k),V(i,j,k))
          Hmin(k) = min(Hmin(k),H(i,j,k))
          Hmax(k) = max(Hmax(k),H(i,j,k))
          Tmin(k) = min(Tmin(k),T(i,j,k))
          Tmax(k) = max(Tmax(k),T(i,j,k))
        END DO
        END DO
      END DO
      if (nstep.ge.dt_ncyc) then
        Bmin =  1.e30
        Bmax = -1.e30
        DO j=1,nyp1
        DO i=it1(j),it2(j)
            Bmin = min(Bmin,BLdpth(i,j))
            Bmax = max(Bmax,BLdpth(i,j))
        end do
        end do
      else
        Bmin = -1.0
        Bmax = -1.0
      end if

      !-----------------------------------------------------------------
      ! print out the diagnositcs
      !-----------------------------------------------------------------
      write(*,9002) 'total',eTO(0),'KE',eKE(0),'PE',ePE(0),
     &                             'HC',eHC(0),'HV',eHV(0)
      DO k=1,nz
        write(*,9004) k,'KE+PE',eTO(k),'KE',eKE(k),'PE',ePE(k),
     &                     'HC',eHC(k),'HV',eHV(k),'HE',Herr(k)
      END DO
      write(*,9006) 'BL',Bmin,Bmax
      DO k=1,nz
        write(*,9008) k,'U',Umin(k),Umax(k), 'V',Vmin(k),Vmax(k),
     +                  'H',Hmin(k),Hmax(k), 'T',Tmin(k),Tmax(k)
      END DO


      call diageq ! same diagnostics, but confined to very near equator

      RETURN
 9000 format('(diagnos) date (yyyymmdd): ',i8.8,i6,'s, nStep=',i7)
 9002 format('(diagnos) ENERGY    ',2x,7(1x,a, (1x,1pe10.3)))
 9004 format('(diagnos) energy  k=',i2,7(1x,a, (1x,1pe10.3)))
 9006 format('(diagnos) min/max   ',2x,7(   a,2(1x,1pe10.3)))
 9008 format('(diagnos) min/max k=',i2,7(1x,a,2(1x,1pe10.3)))
      END

c=======================================================================

      SUBROUTINE diageq

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local variables -----
      real     eTO(0:nz)       ! total     energy
      real     eKE(0:nz)       ! kinetic   energy
      real     ePE(0:nz)       ! potential energy
      real     eHC(0:nz)       ! heat content
      real     eHV(0:nz)       ! heat variance
      real    Umin(0:nz), Umax(0:nz) ! min/max u velocity
      real    Vmin(0:nz), Vmax(0:nz) ! min/max v velocity
      real    Hmin(0:nz), Hmax(0:nz) ! min/max layer depth
      real    Tmin(0:nz), Tmax(0:nz) ! min/max T 
      real    Bmin      , Bmax       ! min/max BL depth
      real    D2(nxp1,nyp1)    ! 2*(total depth at middle of layer k)
      integer i,j,k            ! generic loop indicies
      integer j1,j2            ! j-indicies for near-equitorial region
      real    ZZ               ! area*zeta^2
      real    diffT            ! T(i,j,k)-Tbot
      real    volume           ! H*area
      real    Tavg             ! spatial average of T
      real    basineq          ! area of equitorial basin

c-----------------------------------------------------------------------
c  PURPOSE:
c    diagnostics done within 5 degress of the equator
c    print diagnostics: energies & hi/low values of U,V,H,T
c-----------------------------------------------------------------------

      if (diag_now.eq.0) RETURN

      j1=nyp1
      do j=nyp1,1,-1
	if (yyy(j) .gt. -5.0) j1=j
      end do
      j2=1
      do j=1,nyp1
	if (yyy(j) .lt.  5.0) j2=j
      end do
      if (j1.gt.j2) write(*,*) '(diageq) WARNING: no equitorial region?'
      if (j1.gt.j2) RETURN

      write(*,9000) cdate,sec,nstep

      !-----------------------------------------------------------------
      ! this diagnostic is to check bit-for-bit duplication of a run
      !-----------------------------------------------------------------
      basineq=0.0
      Tavg=0.0
      DO j=j1,j2
      DO i=it1(j),it2(j)
         Tavg=Tavg+T(i,j,1)*area(i,j)
         basineq=basineq+area(i,j)
      END DO
      END DO
      Tavg=Tavg/basineq
      write(*,*)'(diageq ) bit-check: avg T(i,j,1)=',Tavg

      !-----------------------------------------------------------------
      ! compute global energies
      !-----------------------------------------------------------------
      eKE(0)=0.0
      ePE(0)=0.0
      eHC(0)=0.0
      eHV(0)=0.0
      DO k=1,nz
        IF (k .eq. 1) THEN
          ZZ=0.0
          DO j=j1,j2
          DO i=it1(j),it2(j)
            ZZ     =ZZ+area(i,j)*W(i,j,nz)**2
            D2(i,j)=H(i,j,1)
          END DO
          END DO
        ELSE
          DO j=j1,j2
          DO i=it1(j),it2(j)
             D2(i,j)=D2(i,j)+H(i,j,k-1)+H(i,j,k)
          END DO
          END DO
        END IF
        eKE(k)=0.0
        ePE(k)=0.0
        eHC(k)=0.0
        eHV(k)=0.0
        DO j=j1,j2
        DO i=it1(j),it2(j)
          volume=area(i,j)*H(i,j,k)
          diffT =T(i,j,k)-Tbot
          eKE(k)=eKE(k)+volume*(U(i,j,k)**2 + V(i,j,k)**2 )
          ePE(k)=ePE(k)+volume*diffT*D2(i,j)
          eHC(k)=eHC(k)+volume*diffT
          eHV(k)=eHV(k)+volume*diffT**2
        END DO
        END DO
        eKE(k)=0.5*eKE(k)/basin
        ePE(k)=0.5*ePE(k)*alphaT*gravity/basin
        eHC(k)=    eHC(k)/basin
        eHV(k)=    eHV(k)/basin
        eTO(k)=eKE(k)+ePE(k)
        !-------------------
        eKE(0)=eKE(0)+eKE(k) 
        ePE(0)=ePE(0)+ePE(k)
        eHC(0)=eHC(0)+eHC(k)
        eHV(0)=eHV(0)+eHV(k)
      END DO
      ePE(0)=ePE(0)-0.5*gravity*ZZ/basin
      eTO(0)=eKE(0)+ePE(0)

      !-----------------------------------------------------------------
      ! locate the min/max values
      !-----------------------------------------------------------------
      DO k=1,nz
        Umin(k) =  1.e30
        Umax(k) = -1.e30
        Vmin(k) =  1.e30
        Vmax(k) = -1.e30
        Hmin(k) =  1.e30
        Hmax(k) = -1.e30
        Tmin(k) =  1.e30
        Tmax(k) = -1.e30
        DO j=j1,j2
        DO i=it1(j),it2(j)
          Umin(k) = min(Umin(k),U(i,j,k))
          Umax(k) = max(Umax(k),U(i,j,k))
          Vmin(k) = min(Vmin(k),V(i,j,k))
          Vmax(k) = max(Vmax(k),V(i,j,k))
          Hmin(k) = min(Hmin(k),H(i,j,k))
          Hmax(k) = max(Hmax(k),H(i,j,k))
          Tmin(k) = min(Tmin(k),T(i,j,k))
          Tmax(k) = max(Tmax(k),T(i,j,k))
        END DO
        END DO
      END DO
      if (nstep.ge.dt_ncyc) then
        Bmin =  1.e30
        Bmax = -1.e30
        DO j=j1,j2
        DO i=it1(j),it2(j)
           Bmin = min(Bmin,BLdpth(i,j))
           Bmax = max(Bmax,BLdpth(i,j))
        end do
        end do
      else
        Bmin = -1.0
        Bmax = -1.0
      end if

      !-----------------------------------------------------------------
      ! print out the diagnositcs
      !-----------------------------------------------------------------
      write(*,9002) 'total',eTO(0),'KE',eKE(0),'PE',ePE(0),
     &                             'HC',eHC(0),'HV',eHV(0)
      DO k=1,nz
        write(*,9004) k,'KE+PE',eTO(k),'KE',eKE(k),'PE',ePE(k),
     &                     'HC',eHC(k),'HV',eHV(k)
      END DO
      write(*,9006) 'BL',Bmin,Bmax
      DO k=1,nz
        write(*,9008) k,'U',Umin(k),Umax(k), 'V',Vmin(k),Vmax(k),
     +                  'H',Hmin(k),Hmax(k), 'T',Tmin(k),Tmax(k)
      END DO


      RETURN
 9000 format('(diageq ) date (yyyymmdd): ',i8.8,i6,'s, nStep=',i7)
 9002 format('(diageq ) ENERGY    ',2x,7(1x,a, (1x,1pe10.3)))
 9004 format('(diageq ) energy  k=',i2,7(1x,a, (1x,1pe10.3)))
 9006 format('(diageq ) min/max   ',2x,7(   a,2(1x,1pe10.3)))
 9008 format('(diageq ) min/max k=',i2,7(1x,a,2(1x,1pe10.3)))
      END

c=======================================================================
