c=======================================================================
c CVS $Id: integrate.F,v 1.4 1998/08/26 20:46:40 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/integrate.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE BCset

      implicit none

      !----- global -----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields
      include 'control.h'      ! global time step & control info

      !----- local -----
      integer i,j,k       ! generic loop indicies
      real    temp,tink

c-----------------------------------------------------------------------
c PURPOSE:
c   set boundary conditions on HU,HV,HT, & H
c
c NOTE: 
c   this SUBROUTINE has hard-coded dependencies on
c   a particular shape of domain
c-----------------------------------------------------------------------

      if (nstep.eq.0) write(*,'(a)') 
     &   '(BCset) WARNING: hard-coded domain assumptions'

      !-----------------------------------------------------------------
      ! impose velocity BC's
      !-----------------------------------------------------------------

      !---- always set (u,v) o n = 0 at all boundaries -----
      do k=1,nz
        do j=1,nyp1
          HU(it1(j),j,k)=0.0
          HU(it2(j),j,k)=0.0
        end do
        do i=1,nxp1
          HV(i,jt1(i),k)=0.0
          HV(i,jt2(i),k)=0.0
        end do
      end do
      !---- (u,v)=0 velocity at water boundaries? ----------
      if (mBCU.eq.1 .or. mBCU.eq.3) then
        do k=1,nz
          do i=1,122
            HU(i,nyp1,k)=0.0
          end do
          do i=31,nxp1
            HU(i,   1,k)=0.0
          end do
        end do
      end if
      !---- (u,v)=0 velocity at land boundaries? -----------
      if (mBCU.eq.2 .or. mBCU.eq.3) then
        do k=1,nz
          do i=1,31
            HU(i,js1(i),k)=0.0
          end do
          do i=122,nxp1
            HU(i,js2(i),k)=0.0
          end do
          do j=1,nyp1
            HV(is1(j),j,k)=0.0
            HV(is2(j),j,k)=0.0
          end do
        end do
      end if

      !-----------------------------------------------------------------
      ! impose temperature BC's (on the HT field)
      !-----------------------------------------------------------------
 
      !---- fixed T at water boundaries? -------------------
      if (mBCT.eq.1 .OR. mBCT.eq.3) then
        do k=1,nz
          Tink =   Tin(k)
          temp = (-Tink + 1.0)*Tbot
          do i=31,nxp1
          ! HT(i,   1,k)=(Tin(k)*(SST(i,1)-Tbot) + Tbot)*H(i,1,k)
            HT(i,   1,k)=(Tink  * SST(i,1) + temp      )*H(i,1,k)
          end do
          do i=1,122
            HT(i,nyp1,k)=(Tink*SST(i,nyp1) + temp   )*H(i,nyp1,k)
          end do
        end do
      end if
      !---- fixed T at land boundaries? --------------------
      if (mBCT.eq.2 .OR. mBCT.eq.3) then
        write(*,*) '(BCset) ERROR: no implimentation for mBCT=',mBCT
        STOP 'BCset'
      end if

      !-----------------------------------------------------------------
      ! impose layer depth BC's
      !-----------------------------------------------------------------

      !---- fixed H at water or land boundaries? -----------
      if (mBCH.ne.0) then
        write(*,*) '(BCset) ERROR: no implimentation for mBCH=',mBCH
        STOP 'BCset'
      end if

      RETURN

      ENTRY Hfix
c-----------------------------------------------------------------------
c  PURPOSE: 
c    fix H field in corners, with no-slip BC's, 
c    H would remain constant in corners
c  REFERENCE:
c    Gent & Cane (1989), pg 450, eq 37
c-----------------------------------------------------------------------

      H(it1(1),1, 1) = 0.5*(H(it1(1)+1,1, 1) + H(it1(1),2, 1))
      H(it1(1),1,nz) = 0.5*(H(it1(1)+1,1,nz) + H(it1(1),2,nz))
      H(it2(1),1, 1) = 0.5*(H(it2(1)-1,1, 1) + H(it2(1),2, 1))
      H(it2(1),1,nz) = 0.5*(H(it2(1)-1,1,nz) + H(it2(1),2,nz))
 
      do j=2,ny
        if (it1(j).LT.it1(j-1)) then
          H(it1(j),j, 1) = 0.5*(H(it1(j)+1,j, 1)+H(it1(j),j+1, 1))
          H(it1(j),j,nz) = 0.5*(H(it1(j)+1,j,nz)+H(it1(j),j+1,nz))
        ELSE if (it1(j).LT.it1(j+1)) then
          H(it1(j),j, 1) = 0.5*(H(it1(j)+1,j, 1)+H(it1(j),j-1, 1))
          H(it1(j),j,nz) = 0.5*(H(it1(j)+1,j,nz)+H(it1(j),j-1,nz))
        end if
        if (it2(j).GT.it2(j-1)) then
          H(it2(j),j, 1) = 0.5*(H(it2(j)-1,j, 1)+H(it2(j),j+1, 1))
          H(it2(j),j,nz) = 0.5*(H(it2(j)-1,j,nz)+H(it2(j),j+1,nz))
        ELSE if (it2(j).GT.it2(j+1)) then
          H(it2(j),j, 1) = 0.5*(H(it2(j)-1,j, 1)+H(it2(j),j-1, 1))
          H(it2(j),j,nz) = 0.5*(H(it2(j)-1,j,nz)+H(it2(j),j-1,nz))
        end if
      end do
 
      H(it1(nyp1),nyp1, 1) = 0.5*(H(it1(nyp1)+1,nyp1, 1)
     *                           +H(it1(nyp1)  ,ny  , 1))
      H(it1(nyp1),nyp1,nz) = 0.5*(H(it1(nyp1)+1,nyp1,nz)
     *                           +H(it1(nyp1)  ,ny  ,nz))
      H(it2(nyp1),nyp1, 1) = 0.5*(H(it2(nyp1)-1,nyp1, 1)
     *                           +H(it2(nyp1)  ,ny  , 1))
      H(it2(nyp1),nyp1,nz) = 0.5*(H(it2(nyp1)-1,nyp1,nz)
     *                           +H(it2(nyp1)  ,ny  ,nz))
 
      RETURN
      END

c=======================================================================

      SUBROUTINE calcDIV

      implicit none

      !----- global -----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local -----
      integer i,j,k            ! generic loop indicies
      real    TA(nxp1,nyp1)    ! temporary work array
      real    TB(nxp1,nyp1)    ! temporary work array

c-----------------------------------------------------------------------
c  PURPOSE:
c    calculate divergence for all layers
c  INPUT:
c    HU,HV
c  OUTPUT:
c    W    = div(hU), for k=1,..,nz
c    dHnz = sum{k=1 to nz} of dH(k)/dt
c         = sum{k=1 to nz} of -div(hU)
c-----------------------------------------------------------------------

      do j=1,nyp1
      do i=it1(j),it2(j)
        dHnz(i,j)=0.0
      end do
      end do
 
      do k=1,nz
        call DIFFx(HU(1,1,k),TA) !--- TA = d(HU)/di ---
        call DIFFy(HV(1,1,k),TB) !--- TB = d(HV)/dj ---
        do j=1,nyp1
        do i=it1(j),it2(j)
          W(i,j,k)  = EMx(i,j)*TA(i,j) + EMy(i,j)*TB(i,j)
          dHnz(i,j) = dHnz(i,j) - W(i,j,k)
        end do
        end do
      end do
 
      RETURN
 
      ENTRY calcW
c-----------------------------------------------------------------------
c  PURPOSE:
c    calculates W (entrainment) at base of layers 2,..,nz-1
c  NOTE:
c    W(1) is already set and W(nz) is always 0
c  INPUT:
c    W(1) = entrainment at base of mix layer
c    W(k) = div(hU) for layers k=2,..,nz
c    dHnz = sum{k=2,..,nz} of dH(k)/dt
c  OUTPUT:
c    W(k) at base of layers 2,..,nz-1  (over-writes W)
c  REFERENCE:
c    Gent & Cane (1989), pg 454, eq 56
c-----------------------------------------------------------------------

      do k=2,nz-1
        do j=1,nyp1
        do i=it1(j),it2(j)
          W(i,j,k) = W(i,j,k-1) + W(i,j,k) + sigma(k)*dHnz(i,j)
        end do
        end do
      end do

      RETURN
 

      ENTRY calcZETA
c-----------------------------------------------------------------------
c  PURPOSE:
c    calculate ZETA, the sea surface elevation
c  INPUT:
c    H,T
c  OUTPUT:
c    ZETA  (stored in the lowest layer of the W array)
c  ASSUMPTION:
c    the T array contains temperature, *not* hT
c  REFERENCE:
c    Gent & Cane (1989), pg 454, eq 53
c-----------------------------------------------------------------------

      do j=1,nyp1
      do i=it1(j),it2(j)
         W(i,j,nz)=alphaT*H(i,j,1)*(T(i,j,1)-Tbot)
      end do
      end do

      do k=2,nz
      do j=1,nyp1
      do i=it1(j),it2(j)
         W(i,j,nz)=W(i,j,nz)+alphaT*H(i,j,k)*(T(i,j,k)-Tbot)
      end do
      end do
      end do

      RETURN
      END

c=======================================================================

      real function solarpen(depth,n)

      implicit none

      real    depth   ! depth ( >0 ) for desired short-wave fraction (m)
      integer n       ! index for Jerlov water type

      integer   nwtype
      parameter(nwtype=5)
      real  rfac(nwtype),a1(nwtype),a2(nwtype)

      ! Jerlov water type:   I      IA      IB      II      III
      !           index n    1       2       3       4       5
      !                    -----   -----   -----   -----   -----
      data rfac         /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
      data a1           /  0.35 ,  0.60 ,  1.00 ,  1.50 ,  1.40 /
      data a2           / 23.00 , 20.00 , 17.00 , 14.00 ,  7.90 /

      save  

c-----------------------------------------------------------------------
c  PURPOSE:
c    Compute the fraction of solar short-wave radiation penetrating to
c    the specified depth by considering the exponential short-wave 
c    radiation decay of the specified Jerlov water type.
c
c  REFERENCE: 
c    solar absorption model of Simpson and Paulson (1977)
c-----------------------------------------------------------------------

      solarpen =      rfac(n)  * exp(-depth/a1(n))
     &         + (1.- rfac(n)) * exp(-depth/a2(n))

      RETURN
      END

c=======================================================================

      SUBROUTINE dVERT

      implicit none

      !----- global -----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local -----
      integer    jwtype
      PARAMETER (jwtype=3)

      integer   i,j,k          ! generic loop indicies
      real      FTk,FUk,FVk
      real      QS1,QS2,QS3
      real      Sk0,SK1
      real      temp
      real      SSI

      !----- functions -----
      real     solarpen
      external solarpen

c-----------------------------------------------------------------------
c  PURPOSE:
c    calculate vertical terms & update FU,FV,FT appropriately
c  INPUT:
c    TAUx,TAUy,Q
c    T,u,v
c    dH1  = dH(1)/dt
c    dHnz = sum{k=2,..,nz} of dH(k)/dt
c  OUTPUT:
c    modified FU,FV,FT
c  NOTES:
c    o on output, dHnz = dH(nz)/dt
c  REFERENCE:
c    Gent & Cane (1989), pg 453, eq 48,49
c    Gent & Cane (1989), pg 455, eq 59,64
c    Gent & Cane (1989), pg 458, eq 73,74
c-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! top layer
      !-----------------------------------------------------------------
      Sk0 = Sigma(1)
      Sk1 = Sigma(2)
      SSi = 1.0/(Sk0+Sk1)
      do j=1,nyp1
      do i=it1(j),it2(j)
         !--- surface forcing ---
         FU(i,j,1) = FU(i,j,1) + TAUx(i,j)*(1.0/RHOw)
         FV(i,j,1) = FV(i,j,1) + TAUy(i,j)*(1.0/RHOw)
         Qs3=Qsol(i,j)*solarpen(H(i,j,1)+H(i,j,2),jwtype)
         Qs2=Qsol(i,j)*solarpen(H(i,j,1)         ,jwtype)
         Qs1=Qtot(i,j)-Qs2-Qs3
         FT(i,j,1) = FT(i,j,1) + Qs1*(1.0/(Cpw*RHOw))
         FT(i,j,2) = FT(i,j,2) + Qs2*(1.0/(Cpw*RHOw))
         FT(i,j,3) = FT(i,j,3) + Qs3*(1.0/(Cpw*RHOw))
         !--- H(i,j,1) update ---
         FH1(i,j) = FH1(i,j) + dH1(i,j)
         !--- vertical diffusion terms ---
         FUk =    cNU(1)*(U(i,j,2)-U(i,j,1))
         FVk =    cNU(1)*(V(i,j,2)-V(i,j,1))
         FTk = cKAPPA(1)*(T(i,j,2)-T(i,j,1))
         !--- vertical advection terms ---
         FUk = FUk + W(i,j,1)*0.5*(    U(i,j,1)+    U(i,j,2))
         FVk = FVk + W(i,j,1)*0.5*(    V(i,j,1)+    V(i,j,2))
         FTk = FTk + W(i,j,1)*SSi*(Sk0*T(i,j,1)+Sk1*T(i,j,2))
         !--- advection/diffusion redistribution with next layer ---
         FU(i,j,1) = FU(i,j,1) + FUk
         FV(i,j,1) = FV(i,j,1) + FVk
         FT(i,j,1) = FT(i,j,1) + FTk
         FU(i,j,2) = FU(i,j,2) - FUk
         FV(i,j,2) = FV(i,j,2) - FVk
         FT(i,j,2) = FT(i,j,2) - FTk
         !--- part of horizontal divergence term ---
         temp = 0.5*( 0.0-W(i,j,1)  +  dH1(i,j) )
         FU(i,j,1) = FU(i,j,1) + temp*U(i,j,1)
         FV(i,j,1) = FV(i,j,1) + temp*V(i,j,1)
         FT(i,j,1) = FT(i,j,1) + temp*T(i,j,1)
      end do
      end do

      !-----------------------------------------------------------------
      ! mid layers 
      !-----------------------------------------------------------------
      do k=2,nz-1
         Sk0 = Sigma(k)
         Sk1 = Sigma(k+1)
         SSi = 1.0/(Sk0+Sk1)
         do j=1,nyp1
         do i=it1(j),it2(j)
            !--- part of horizontal divergence term ---
            temp = 0.5*(W(i,j,k-1)-W(i,j,k)  +  dHnz(i,j)*Sk0)
            FU(i,j,k) = FU(i,j,k) + temp*U(i,j,k)
            FV(i,j,k) = FV(i,j,k) + temp*V(i,j,k)
            FT(i,j,k) = FT(i,j,k) + temp*T(i,j,k)
            !--- vertical diffusion terms ---
            FUk =    cNU(k)*(U(i,j,k+1)-U(i,j,k))
            FVk =    cNU(k)*(V(i,j,k+1)-V(i,j,k))
            FTk = cKAPPA(k)*(T(i,j,k+1)-T(i,j,k))
            !--- vertical advection terms ---
            FUk = FUk + W(i,j,k)*0.5*(    U(i,j,k)+    U(i,j,k+1))
            FVk = FVk + W(i,j,k)*0.5*(    V(i,j,k)+    V(i,j,k+1))
            FTk = FTk + W(i,j,k)*SSi*(Sk0*T(i,j,k)+Sk1*T(i,j,k+1))
            !--- advection/diffusion redistribution with next layer ---
            FU(i,j,k  ) = FU(i,j,k  ) + FUk
            FV(i,j,k  ) = FV(i,j,k  ) + FVk
            FT(i,j,k  ) = FT(i,j,k  ) + FTk
            FU(i,j,k+1) = FU(i,j,k+1) - FUk
            FV(i,j,k+1) = FV(i,j,k+1) - FVk
            FT(i,j,k+1) = FT(i,j,k+1) - FTk
         end do
         end do
      end do

      !-----------------------------------------------------------------
      ! bottom layer
      !-----------------------------------------------------------------
      Sk0 = sigma(nz)
      do j=1,nyp1
      do i=it1(j),it2(j)
         !--- H(i,j,nz) update ---
         FHnz(i,j) = FHnz(i,j) +  dHnz(i,j)*Sk0
         !--- part of horizontal divergence term ---
         temp = 0.5*(W(i,j,nz-1)-0.0  +  dHnz(i,j)*Sk0)
         FU(i,j,nz) = FU(i,j,nz) + temp*U(i,j,nz)
         FV(i,j,nz) = FV(i,j,nz) + temp*V(i,j,nz)
         FT(i,j,nz) = FT(i,j,nz) + temp*T(i,j,nz)
         !--- advection/diffusion to deep ocean ---
         FU(i,j,nz) = FU(i,j,nz) +    cNU(nz)*( 0.0-U(i,j,nz))
         FV(i,j,nz) = FV(i,j,nz) +    cNU(nz)*( 0.0-V(i,j,nz))
         FT(i,j,nz) = FT(i,j,nz) + cKAPPA(nz)*(Tbot-T(i,j,nz))
      end do
      end do

      RETURN
      END

c=======================================================================

      SUBROUTINE dHORIZ

      implicit none

      !----- global -----
      include 'dims.h'         ! defines array dimensions
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local -----
      integer i,j,k            ! generic loop indicies
      real    TA(nxp1,nyp1),TB(nxp1,nyp1),TC(nxp1,nyp1)
      real    TD(nxp1,nyp1),TE(nxp1,nyp1),TF(nxp1,nyp1)

c-----------------------------------------------------------------------
c  PURPOSE:
c    calculate horizontal advection & coriolis terms
c  INPUT:
c    u,v,T,hu,hv,f
c  ASSUMPTION:
c    the T array contains temperature, *not* hT
c  REFERENCE:
c    Gent & Cane (1989), pg 455, eq 59,64
c-----------------------------------------------------------------------

      do k=1,nz
         !--- update FU ---
         do j=1,nyp1
         do i=it1(j),it2(j)
           TA(i,j)=HU(i,j,k)*U(i,j,k)
           TB(i,j)=HU(i,j,k)*V(i,j,k)
         end do
         end do
         call DIFFx (U(1,1,k),TC)
         call DIFFx (TA      ,TD)
         call DIFFy (U(1,1,k),TE)
         call DIFFy (TB      ,TF)
         do j=1,nyp1
         do i=it1(j),it2(j)
           FU(i,j,k) = FU(i,j,k) + F(j)*HV(i,j,k) - 0.5*(
     &                 EMx(i,j)*( TD(i,j)+HU(i,j,k)*TC(i,j) )
     &               + EMy(i,j)*( TF(i,j)+HV(i,j,k)*TE(i,j) )    )
           TA(i,j) = HV(i,j,k)*V(i,j,k)
         end do
         end do
         !--- update FV ---
         call DIFFx (V(1,1,k),TC)
         call DIFFx (TB      ,TD)
         call DIFFy (V(1,1,k),TE)
         call DIFFy (TA      ,TF)
         do j=1,nyp1
         do i=it1(j),it2(j)
           FV(i,j,k) = FV(i,j,k) - F(j)*HU(i,j,k) - 0.5*(
     &                 EMx(i,j)*( TD(i,j)+HU(i,j,k)*TC(i,j) )
     &               + EMy(i,j)*( TF(i,j)+HV(i,j,k)*TE(i,j) )    )
           TA(i,j) = HU(i,j,k)*T(i,j,k)
           TB(i,j) = HV(i,j,k)*T(i,j,k)
         end do
         end do
         !--- update FT ---
         call DIFFx (T(1,1,k),TC)
         call DIFFx (TA      ,TD)
         call DIFFy (T(1,1,k),TE)
         call DIFFy (TB      ,TF)
         do j=1,nyp1
         do i=it1(j),it2(j)
           FT(i,j,k) = FT(i,j,k) - 0.5*(
     &                 EMx(i,j)*( TD(i,j)+HU(i,j,k)*TC(i,j) )
     &               + EMy(i,j)*( TF(i,j)+HV(i,j,k)*TE(i,j) )    )
         end do
         end do
      end do

      RETURN
      END

c=======================================================================

      SUBROUTINE dPGF

      implicit none

      !----- global -----
      include 'constant.h'     ! defines physical constants
      include 'dims.h'         ! defines array dimensions
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local -----
      integer i,j,k            ! generic loop indicies
      real    temp
      real    TAL
      real    GAL
      real    TA(nxp1,nyp1),TB(nxp1,nyp1) ! temp work arrays
      real    TC(nxp1,nyp1),TD(nxp1,nyp1) ! temp work arrays
      real    TE(nxp1,nyp1)               ! temp work arrays

c-----------------------------------------------------------------------
c  PURPOSE:
c    calculate pressure gradient force terms
c  INPUT:
c    h,hT,w
c  ASSUMPTION:
c    lowest layer of W array contains zeta
c    T array should contain hT here, *not* T
c  NOTE:
c    hb       = h(k)*alphaT*gravity*(T(k) - Tbot)
c    grad(hb) = alphaT*gravity*grad( hT )
c  REFERENCE:
c    Gent & Cane (1989), pg 454, eq 55
c-----------------------------------------------------------------------

      GAL=gravity*alphaT

      !------------------------------------------
      ! compute h*(g-b)*grad(zeta) term
      !------------------------------------------
      call DIFFx (W(1,1,nz),TA)
      call DIFFy (W(1,1,nz),TB)
      do j=1,nyp1
      do i=it1(j),it2(j)
        TD(i,j) = GAL*EMx(i,j)*TA(i,j)
        TE(i,j) = GAL*EMy(i,j)*TB(i,j)
      end do
      end do
      TAL = Tbot + 1.0/alphaT
      do k=1,nz
      do j=1,nyp1
      do i=it1(j),it2(j)
        temp = HT(i,j,k) - TAL*H(i,j,k)
        FU(i,j,k) = FU(i,j,k) + temp*TD(i,j)
        FV(i,j,k) = FV(i,j,k) + temp*TE(i,j)
      end do
      end do
      end do
      
      !------------------------------------------
      ! set TA=2h, TB=2hT at middle of kth layer
      !------------------------------------------
      do j=1,nyp1
      do i=it1(j),it2(j)
        TA(i,j)=H (i,j,1)
        TB(i,j)=HT(i,j,1)
      end do
      end do
      do k=1,nz
         if (k .gt. 1) then
            do j=1,nyp1
            do i=it1(j),it2(j)
               TA(i,j) = TA(i,j)+H(i,j,k)+H(i,j,k-1)
               TB(i,j) = TB(i,j)+HT(i,j,k)+HT(i,j,k-1)
            end do
            end do
         end if
         !--- calculate x derivatives & update U momentum eq ---
         call DIFFx (TA,TC)
         call DIFFx (TB,TD)
         do j=1,nyp1
         do i=it1(j),it2(j)
           FU(i,j,k)= FU(i,j,k)
     &         -0.5*GAL*EMx(i,j)*(HT(i,j,k)*TC(i,j)-H(i,j,k)*TD(i,j) )
         end do
         end do
         !--- calculate y derivatives & update V momentum eq ---
         call DIFFy (TA, TC)
         call DIFFy (TB, TD)
         do j=1,nyp1
         do i=it1(j),it2(j)
           FV(i,j,k)= FV(i,j,k)
     &         -0.5*GAL*EMy(i,j)*(HT(i,j,k)*TC(i,j)-H(i,j,k)*TD(i,j) )
         end do
         end do
      end do

      RETURN
      END

c=======================================================================

      SUBROUTINE dRMIX

      implicit none

      !----- global -----
      include 'dims.h'         ! defines array dimensions
      include 'constant.h'     ! defines physical constants
      include 'control.h'      ! time stepping info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local -----
      integer i,j,k            ! generic loop indicies
      integer ierru,ierrv,ierrt
      real    FA1(nz),FA2(nz),FA3(nz) 
      real    ANUU(nz,3),BU(nz) 
      real    ANUV(nz,3),BV(nz) 
      real    AKAP(nz,3),BT(nz)
c-----------------------------------------------------------------------
c  PURPOSE:
c    do implicit vertical mixing 
c  INPUT:
c    U,V,T,rNU,rKAPPA,rGAMMA 
c  OUTPUT:
c    U,V,T
c  ASSUMPTIONS:
c    rNU,rKAPPA,rGAMMA are valid for this time step
c  IMPLEMENTATION:
c    diffusion is based on the vertical laplacian 
c    at the *next* time step
c    T(t+dNT) = T(t) + dNT*kappa*lap(T(t+dNT))
c    U(t+dNT) = U(t) + dNT*nu   *lap(U(t+dNT))
c    Note: this requires a matrix inversion, Ax=b 
c          where b=T(t), and x=T(t+dNT)
c-----------------------------------------------------------------------

      do j=1,nyp1
      do i=it1(j),it2(j)
          do k=2,nz
             FA1(k)  = 2.0*dNT/H(i,j,k)
             FA2(k-1)=rNU   (i,j,k-1)/(H(i,j,k)+H(i,j,k-1))
             FA3(k-1)=rKAPPA(i,j,k-1)/(H(i,j,k)+H(i,j,k-1))
             BU (k)  = U(i,j,k)
             BV (k)  = V(i,j,k)
             BT (k)  = T(i,j,k)+H(i,j,k)*(rKAPPA(i,j,k-1)*
     &                 rGAMMA(i,j,k-1)-rKAPPA(i,j,k)*rGAMMA(i,j,k))
          end do 
          FA1(1)  = 2.0*dNT/H(i,j,1) 
          BU (1)  = U(i,j,1) 
          BV (1)  = V(i,j,1) 
          BT (1)  = T(i,j,1)-H(i,j,1)*rKAPPA(i,j,1)*rGAMMA(i,j,1) 
          FA2(nz) = rNU   (i,j,nz)/H(i,j,nz)
          FA3(nz) = rKAPPA(i,j,nz)/H(i,j,nz)
          BT (nz) = BT(nz)+Tbot*FA1(nz)*FA3(nz)
          do k=2,nz 
             ANUu(k,1)=   -FA1(k)* FA2(k-1)
             ANUu(k,2)=1.0+FA1(k)*(FA2(k-1)+FA2(k))
             ANUu(k,3)=   -FA1(k)* FA2(k)
             ANUv(k,1)=ANUu(k,1)
             ANUv(k,2)=ANUu(k,2)
             ANUv(k,3)=ANUu(k,3)
             AKAP(k,1)=   -FA1(k)* FA3(k-1)
             AKAP(k,2)=1.0+FA1(k)*(FA3(k-1)+FA3(k))
             AKAP(k,3)=   -FA1(k)* FA3(k)
          end do
          ANUu( 1,1)=0.0
          ANUu( 1,2)=1.0+FA1(1)*FA2(1)
          ANUu( 1,3)=   -FA1(1)*FA2(1)
          ANUv( 1,1)=0.0
          ANUv( 1,2)=ANUu(1,2)
          ANUv( 1,3)=ANUu(1,3)
          AKAP( 1,1)=0.0
          AKAP( 1,2)=1.0+FA1(1)*FA3(1)
          AKAP( 1,3)=   -FA1(1)*FA3(1)
          ANUu(nz,3)=0.0
          ANUv(nz,3)=0.0
          AKAP(nz,3)=0.0
#if ( defined CRAY )
          call SGTSL(nz,ANUu(1,1),ANUu(1,2),ANUu(1,3),BU,ierrU)
          call SGTSL(nz,ANUv(1,1),ANUv(1,2),ANUv(1,3),BV,ierrV)
          call SGTSL(nz,AKAP(1,1),AKAP(1,2),AKAP(1,3),BT,ierrT)
#else
          call DGTSL(nz,ANUu(1,1),ANUu(1,2),ANUu(1,3),BU,ierrU)
          call DGTSL(nz,ANUv(1,1),ANUv(1,2),ANUv(1,3),BV,ierrV)
          call DGTSL(nz,AKAP(1,1),AKAP(1,2),AKAP(1,3),BT,ierrT)
#endif
          do k=1,nz
             U(i,j,k)=BU(k)
             V(i,j,k)=BV(k)
             T(i,j,k)=BT(k)
          end do
      end do
      end do

      RETURN
      END

c=======================================================================

      SUBROUTINE formHUVT

      implicit none

      !----- global -----
      include 'dims.h'         ! defines array dimensions
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local -----
      integer i,j,k            ! generic loop indicies
c-----------------------------------------------------------------------
c  PURPOSE: 
c    form HU,HV,HT from U,V,T
c-----------------------------------------------------------------------

      do k=1,nz
      do j=1,nyp1
      do i=it1(j),it2(j)
        HU(i,j,k)=U(i,j,k)*H(i,j,k)
        HV(i,j,k)=V(i,j,k)*H(i,j,k)
        HT(i,j,k)=T(i,j,k)*H(i,j,k)
      end do
      end do
      end do

      RETURN
 

      ENTRY formUVT
c-----------------------------------------------------------------------
c  PURPOSE: 
c    form U,V,T from HU,HV,HT 
c-----------------------------------------------------------------------

      do k=1,nz
      do j=1,nyp1
      do i=it1(j),it2(j)
        U(i,j,k)=HU(i,j,k)/H(i,j,k)
        V(i,j,k)=HV(i,j,k)/H(i,j,k)
        T(i,j,k)=HT(i,j,k)/H(i,j,k)
      end do
      end do
      end do

      RETURN
      END

c=======================================================================

      SUBROUTINE DIFFset

      implicit none

      !----- global -----
      include 'dims.h'         ! defines array dimensions
      include 'domain.h'       ! global grid & domain info

      !----- local -----
      integer i,j              ! generic loop indicies

      REAL      G(nxp1,nyp1),dG(nxp1,nyp1)
      REAL      a1,a2,a3
      REAL      b1,b2,b3,b4
      REAL      c1,c2
      PARAMETER ( a1=-14./12., a2=16./12., a3=-2./12.             )
      PARAMETER ( b1=- 6./12., b2=-1./12., b3= 8./12., b4=-1./12. )
      PARAMETER ( c1=  8./12., c2=-1./12.)

c-----------------------------------------------------------------------
c  PURPOSE:
c    initialize & implement the centered difference scheme 
c  REFERENCE:
c    Gent & Cane (1989), pg 452
c  NOTE: 
c    The boundary form has been modified so that 
c    sum(weight(i)*dG(i)) = G(nxp1) - G(1),
c    where weight(1)=weight(nxp1)=1/2, all other weight(i)=1.0 
c-----------------------------------------------------------------------
 
      write(*,*) 
      write(*,*)  'Initializing finite difference coefficients:'
      write(*,*)  '  o interior: 4th order centered'
      write(*,*)  '  o boundary: 3rd order'
      write(*,*) 

      RETURN
 

      ENTRY DIFFx (G,dG)
c-----------------------------------------------------------------------
c  PURPOSE: 
c    calculate dG/di
c-----------------------------------------------------------------------

      do j=1,nyp1
      do i=it1(j)+2,it2(j)-2
        dG(i,j) = c1*( G(i+1,j)-G(i-1,j) ) + c2*( G(i+2,j)-G(i-2,j) )
      end do
      end do
      do j=1,nyp1
        dG(it1(j)  ,j) = + a1*G(it1(j)  ,j) + a2*G(it1(j)+1,j)
     ~                   + a3*G(it1(j)+2,j)
        dG(it1(j)+1,j) = + b1*G(it1(j)  ,j) + b2*G(it1(j)+1,j)
     ~                   + b3*G(it1(j)+2,j) + b4*G(it1(j)+3,j)
        dG(it2(j)  ,j) = - a1*G(it2(j)  ,j) - a2*G(it2(j)-1,j)
     ~                   - a3*G(it2(j)-2,j)
        dG(it2(j)-1,j) = - b1*G(it2(j)  ,j) - b2*G(it2(j)-1,j)
     ~                   - b3*G(it2(j)-2,j) - b4*G(it2(j)-3,j)
      end do
      RETURN
 

      ENTRY DIFFy (G,dG)
c-----------------------------------------------------------------------
c  PURPOSE: 
c    calculate dG/dj
c-----------------------------------------------------------------------

      do i=1,nxp1
      do j=jt1(i)+2,jt2(i)-2
        dG(i,j) = c1*( G(i,j+1)-G(i,j-1) ) + c2*( G(i,j+2)-G(i,j-2) )
      end do
      end do
      do i=1,nxp1
        dG(i,jt1(i)  ) = + a1*G(i,jt1(i)  ) + a2*G(i,jt1(i)+1)
     ~                   + a3*G(i,jt1(i)+2)
        dG(i,jt1(i)+1) = + b1*G(i,jt1(i)  ) + b2*G(i,jt1(i)+1)
     ~                   + b3*G(i,jt1(i)+2) + b4*G(i,jt1(i)+3)
        dG(i,jt2(i)  ) = - a1*G(i,jt2(i)  ) - a2*G(i,jt2(i)-1)
     ~                   - a3*G(i,jt2(i)-2)
        dG(i,jt2(i)-1) = - b1*G(i,jt2(i)  ) - b2*G(i,jt2(i)-1)
     ~                   - b3*G(i,jt2(i)-2) - b4*G(i,jt2(i)-3)
      end do

      RETURN
      END

c=======================================================================

      SUBROUTINE correctH()

      implicit none

      !----- global -----
      include 'dims.h'         ! defines array dimensions
      include 'data.h'         ! global prognostic & diagnostic fields
      include 'domain.h'       ! global grid & domain info

      !----- local -----
      integer i,j,k            ! generic loop indicies
      real    volume           ! integrated volume of a level

c-----------------------------------------------------------------------
c  PURPOSE:
c    o calculate depth errors: eDE(k)
c    o correct depths
c    o check for fatal error: H .le. 0
c  INPUT:
c    o H,Hin,area,basin
c  OUTPUT:
c    o H    (volume corrected)
c    o Herr (for diagnostics only)
c  ASSUMPTIONS:
c    o area  = area of a single grid cell
c    o basin = area of entire basin
c    o Hin   = average layer depth
c-----------------------------------------------------------------------

      do k=1,nz
         volume=0.0
         do j=1,nyp1
         do i=it1(j),it2(j)
            volume=volume+area(i,j)*H(i,j,k)
         end do
         end do
         Herr(k)=volume/basin-Hin(k)
         do j=1,nyp1
         do i=it1(j),it2(j)
            H(i,j,k)=H(i,j,k)-Herr(k)
         end do
         end do
      end do

      do j=1,nyp1
      do i=it1(j),it2(j)
         if (H(i,j,nz) .le. 0.0) then
           write(*,9000) '(correctH) ERROR: H<0',XXX(i),YYY(j)
           stop '(H<0)'
         end if
      end do
      end do

      RETURN
 9000 format(a,' at x=',f9.4,', y=',f9.4)
      END

c=======================================================================
