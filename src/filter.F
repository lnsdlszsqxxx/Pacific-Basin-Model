c=======================================================================
c CVS $Id: filter.F,v 1.2 1998/08/24 22:45:46 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/filter.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE FilterU

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local -----
      integer rFilt,cFilt1,cFilt2,cFilt3,mBCUr,mBCUc
      integer mCode,mFixWater,mFixLand,mRedOrd,mCons    
      logical first_call

      DATA    first_call /.true./

      SAVE    

c-----------------------------------------------------------------------
c PURPOSE:
c   o Initialize the scheme numbers based on input parmaeter mShapU.
c     See routines Filter(), ShapiroR(), ShapiroC() for more details.
c   o Call the filtering routine with appropriate options.
c
c NOTE: 
c   This routine is hard-coded for a particlar boundary shape.
c
c FILTER OPTIONS:
c   rFilt  : type of row filtering    (only 1 region)
c   cFiltX : type of column filtering (for 3 regions)
c   mBCU   : BC's => filtering domain
c   nOrdrU : order of filter (2**n)
c
c Decoding the rFilt & cFilt options (in binary): 
c   nnn1 => max fixed     , nnn0 => filtered
c   nn1n => min fixed     , nn0n => filtered
c   n1nn => reduced order , n0nn => full order
c   1nnn => conservative  , 0nnn => *possibly* non-conservative
c
c   Thus we have:
c   0100 = 4 => cons?, reduced order, min filtered, max filtered
c   0101 = 5 => cons?, reduced order, min filtered, max fixed
c   0110 = 6 => cons?, reduced order, min fixed,    max filtered
c   0111 = 7 => cons?, reduced order, min fixed   , max fixed   
c
c Decoding the mBC options:
c   00 = 0 => land free  , water free
c   01 = 1 => land free  , water fixed  
c   10 = 2 => land fixed , water free
c   11 = 3 => land fixed , water fixed
c    
c   The choice of mBCU (which specifies free/no-slip at land/water 
c   boundaries) dictates the choice of row & column filtering 
c   domains (ie. it1(j) vs. is1(j), etc).  There is hard-coded
c   knowledge of which boundary points represent land and water.
c
c-----------------------------------------------------------------------

      if (first_call) then
        write(*,'(a)') 
     &    '(filterU) WARNING: hard-coded for particular domain shape.'
        first_call = .false.

        !--- decode user requests ---
        mCode     = mShapU
        mFixWater = mod(mCode,2)
        mCode     = mCode/2
        mFixLand  = mod(mCode,2)
        mCode     = mCode/2
        mRedOrd   = mod(mCode,2)
        mCode     = mCode/2
        mCons     = mod(mCode,2)

        !--- set row and column filter option ---
        rFilt =8*mCons + 4*mRedOrd + 2*mFixLand  + mFixLand
        cFilt1=8*mCons + 4*mRedOrd + 2*mFixLand  + mFixWater
        cFilt2=8*mCons + 4*mRedOrd + 2*mFixWater + mFixWater
        cFilt3=8*mCons + 4*mRedOrd + 2*mFixWater + mFixLand

        !--- preserve (u,v) dot n = 0 at boundary ---
        !--- by reducing filter domains           ---
        mBCUr=3
        mBCUc=3
      end if

      call filter(HU,1,nz,rFilt,cFilt1,cFilt2,cFilt3,mBCUr,mBCU ,nOrdrU)
      call filter(HV,1,nz,rFilt,cFilt1,cFilt2,cFilt3,mBCU ,mBCUc,nOrdrU)

      RETURN
      END

c=======================================================================

      SUBROUTINE FilterH

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local -----
      integer rFilt,cFilt1,cFilt2,cFilt3
      integer mCode,mFixWater,mFixLand,mRedOrd,mCons    
      logical first_call
      DATA    first_call /.true./

      SAVE

c-----------------------------------------------------------------------
c PURPOSE:
c   o Initialize the scheme numbers based on input parmaeter mBCT.
c   o Call the filtering routine with appropriate options.
c     See routine FilterU() or Filter() for more details.
c NOTE: 
c   This routine is hard-coded for a particlar boundary shape.
c-----------------------------------------------------------------------

      if (first_call) then
        write(*,'(a)') 
     &    '(filterH) WARNING: hard-coded for particular domain shape.'
        first_call = .false.

        !--- decode user requests ---
        mCode     = mShapH
        mFixWater = mod(mCode,2)
        mCode     = mCode/2
        mFixLand  = mod(mCode,2)
        mCode     = mCode/2
        mRedOrd   = mod(mCode,2)
        mCode     = mCode/2
        mCons     = mod(mCode,2)

        !--- set row and column filter option ---
        rFilt =8*mCons + 4*mRedOrd + 2*mFixLand  + mFixLand
        cFilt1=8*mCons + 4*mRedOrd + 2*mFixLand  + mFixWater
        cFilt2=8*mCons + 4*mRedOrd + 2*mFixWater + mFixWater
        cFilt3=8*mCons + 4*mRedOrd + 2*mFixWater + mFixLand
      end if

ccc   assume  H(i,j,1) is constant
ccc   call filter(H, 1, 1,rFilt,cFilt1,cFilt2,cFilt3,mBCH,mBCH,nOrdrH)
      call filter(H,nz,nz,rFilt,cFilt1,cFilt2,cFilt3,mBCH,mBCH,nOrdrH)
     
      RETURN
      END

c=======================================================================

      SUBROUTINE FilterT

      implicit none

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'control.h'      ! global time step & control info
      include 'domain.h'       ! global grid & domain info
      include 'data.h'         ! global prognostic & diagnostic fields

      !----- local -----
      integer rFilt,cFilt1,cFilt2,cFilt3
      integer mCode,mFixWater,mFixLand,mRedOrd,mCons    
      logical first_call

      DATA    first_call /.true./

      SAVE

c-----------------------------------------------------------------------
c PURPOSE:
c   o Initialize the scheme numbers based on input parmaeter mBCT.
c   o Call the filtering routine with appropriate options.
c   See routine FilterU() or Filter() for more details.
c NOTE: 
c   This routine is hard-coded for a particlar boundary shape.
c   Do you want to filter T or hT?
c-----------------------------------------------------------------------


      if (first_call) then
        write(*,'(a)') 
     &    '(filterT) WARNING: hard-coded for particular domain shape.'
        first_call = .false.

        !--- decode user requests ---
        mCode     = mShapT
        mFixWater = mod(mCode,2)
        mCode     = mCode/2
        mFixLand  = mod(mCode,2)
        mCode     = mCode/2
        mRedOrd   = mod(mCode,2)
        mCode     = mCode/2
        mCons     = mod(mCode,2)

        !--- set row and column filter option ---
        rFilt =8*mCons + 4*mRedOrd + 2*mFixLand  + mFixLand
        cFilt1=8*mCons + 4*mRedOrd + 2*mFixLand  + mFixWater
        cFilt2=8*mCons + 4*mRedOrd + 2*mFixWater + mFixWater
        cFilt3=8*mCons + 4*mRedOrd + 2*mFixWater + mFixLand
      end if

      call filter(HT,1,nz,rFilt,cFilt1,cFilt2,cFilt3,mBCT,mBCT,nOrdrT)

      RETURN
      END

c=======================================================================

      SUBROUTINE Filter(FLD,k0,k1,s1,s2,s3,s4,mBCr,mBCc,nOrdr)

      implicit none 

      !---- global ----
      include 'dims.h'         ! defines array dimensions
      include 'domain.h'       ! global grid & domain info

      !----- arguments -----
      real      FLD(nxp1,nyp1,nz)
      integer   k0,k1
      integer   s1,s2,s3,s4
      integer   mBCr,mBCc,nOrdr

      !----- local -----
      integer   i,j,k
      real      TA(nxp1,nyp1)
      real      TB(nxp1,nyp1)

c-----------------------------------------------------------------------
c PURPOSE:
c   Call ShapiroR & ShapiroC to filter the field FLD
c
c   The row & column domains 
c   (ie. it1(j) vs. is1(j), j=1,nyp1 vs. j=2,ny, etc) 
c   are selected based on input parameter mBCx and 
c   hard-coded knowledge of the boundary of the domain.
c
c Decoding the s1,s2,s3,s4 flags (in binary): 
c   s1        applies to row filtering
c   s2,s3,s4  apply to the 3 column filtering regions
c
c   nnn1 => max fixed     , nnn0 => filtered
c   nn1n => min fixed     , nn0n => filtered
c   n1nn => reduced order , n0nn => full order
c   1nnn => conservative  , 0nnn => *possibly* non-conservative
c 
c   Thus we have:
c   0100 = 4 => cons?, reduced order, min filtered, max filtered
c   0101 = 5 => cons?, reduced order, min filtered, max fixed
c   0110 = 6 => cons?, reduced order, min fixed,    max filtered
c   0111 = 7 => cons?, reduced order, min fixed   , max fixed   
c   (see routine ShapiroR() for more details)
c
c Decoding the mBCr,mBCc flags: 
c   00 = 0 => land free  , water free  
c   01 = 1 => land free  , water fixed  
c   10 = 2 => land fixed , water free  
c   11 = 3 => land fixed , water fixed
c
c-----------------------------------------------------------------------

      !--- filter rows -------------------------------------------------
      IF (mBCr.eq.0) THEN
        DO 100 k=k0,k1
          call ShapiroR(
     +         FLD(1,1,k),TA,TB,nOrdr,s1,nxp1,nyp1,1,nyp1,it1,it2)
           DO 100 j=1,nyp1
           DO 100 i=it1(j),it2(j)
              FLD(i,j,k)=FLD(i,j,k)-TA(i,j)
  100   CONTINUE
      ELSE IF (mBCr.eq.1) THEN
        DO 200 k=k0,k1
          call ShapiroR(
     +         FLD(1,1,k),TA,TB,nOrdr,s1,nxp1,nyp1,2,ny  ,it1,it2)
           DO 200 j=2,ny
           DO 200 i=it1(j),it2(j)
              FLD(i,j,k)=FLD(i,j,k)-TA(i,j)
  200   CONTINUE
      ELSE IF (mBCr.eq.2) THEN
        DO 300 k=k0,k1
          call ShapiroR(
     +         FLD(1,1,k),TA,TB,nOrdr,s1,nxp1,nyp1,1,nyp1,is1,is2)
           DO 300 j=1,nyp1
           DO 300 i=is1(j),is2(j)
              FLD(i,j,k)=FLD(i,j,k)-TA(i,j)
  300   CONTINUE
      ELSE IF (mBCr.eq.3) THEN
        DO 400 k=k0,k1
          call ShapiroR(
     +         FLD(1,1,k),TA,TB,nOrdr,s1,nxp1,nyp1,2,ny  ,is1,is2)
           DO 400 j=2,ny
           DO 400 i=is1(j),is2(j)
              FLD(i,j,k)=FLD(i,j,k)-TA(i,j)
  400   CONTINUE
      ELSE
        write(*,*) 'FILTER: filter() has invalid mBCr = ',mBCr
        STOP 'filter'
      END IF

      !--- filter columns ----------------------------------------------
      IF (mBCc.eq.0 .or. mBCc.eq.1) THEN
        DO 500 k=k0,k1
          call ShapiroC(
     +         FLD(1,1,k),TA,TB,nOrdr,s2,nxp1,nyp1,  1,  30,jt1,jt2)
          call ShapiroC(
     +         FLD(1,1,k),TA,TB,nOrdr,s3,nxp1,nyp1, 31, 122,jt1,jt2)
          call ShapiroC(
     +         FLD(1,1,k),TA,TB,nOrdr,s4,nxp1,nyp1,123,nxp1,jt1,jt2)
          DO 500 i=1,nxp1
          DO 500 j=jt1(i),jt2(i)
             FLD(i,j,k)=FLD(i,j,k)-TA(i,j)
  500   CONTINUE
      ELSE IF (mBCc.eq.2 .or. mBCc.eq.3) THEN 
        DO 600 k=k0,k1
          call ShapiroC(
     +         FLD(1,1,k),TA,TB,nOrdr,s2,nxp1,nyp1,  2,  31,js1,js2)
          call ShapiroC(
     +         FLD(1,1,k),TA,TB,nOrdr,s3,nxp1,nyp1, 32, 121,js1,js2)
          call ShapiroC(
     +         FLD(1,1,k),TA,TB,nOrdr,s4,nxp1,nyp1,122,nx  ,js1,js2)
          DO 600 i=2,nx
          DO 600 j=js1(i),js2(i)
             FLD(i,j,k)=FLD(i,j,k)-TA(i,j)
  600   CONTINUE
      ELSE
        write(*,*) 'FILTER: filter() has invalid mBCc = ',mBCc
        STOP 'filter'
      END IF

      RETURN
      END

c=======================================================================

      SUBROUTINE ShapiroR
     +  (A,FA,Work, nOrder,scheme, nx,ny, jMin,jMax,iMin,iMax)

      implicit none 

      integer   nOrder,scheme,nx,ny,jMin,jMax,iMin(ny),iMax(ny)
      real      A(nx,ny),FA(nx,ny),Work(nx,ny)

      integer   i,j,k,kmax

      real      factor(20)

      DATA      factor/
     + 2.500000E-01, 6.25000000E-02, 1.562500000E-02, 3.90625000000E-03,
     + 9.765625E-04, 2.44140625E-04, 6.103515625E-05, 1.52587890625E-05,
     + 3.814697E-06, 9.53674300E-07, 2.384186000E-07, 5.96046400000E-08,
     + 1.490116E-08, 3.72529000E-09, 9.313226000E-10, 2.32830600000E-10,
     + 5.820766E-11, 1.45519200E-11, 3.637979000E-12, 9.09494700000E-13/

      SAVE factor
 
c-----------------------------------------------------------------------
c   PURPOSE:
c      This routine applies a 1d shapiro filter to rows jMin through 
c      jMax of the array A, and returns the correction array as FA
c
c   INPUT:
c      A      : array to be filtered
c      nOrder : order of the filter = 2*k, for k in {1,2,3,...,20}
c      scheme : flags the type of boundary scheme to be used...
c         nnn1 => max fixed     , nnn0 => filtered
c         nn1n => min fixed     , nn0n => filtered
c         n1nn => reduced order , n0nn => full order
c         1nnn => conservative  , 0nnn => *possibly* non-conservative
c         Thus we have:
c         0100 = 4 => cons?, reduced order, min filtered, max filtered
c         0101 = 5 => cons?, reduced order, min filtered, max fixed
c         0110 = 6 => cons?, reduced order, min fixed,    max filtered
c         0111 = 7 => cons?, reduced order, min fixed   , max fixed   
c      nx,ny  : array dimensions
c      jMin,jMax,iMin(j),iMax(j) : boundaries of region to be filtered
c
c    OUTPUT:
c       FA : correction array (ie. filtered field = A-FA)
c
c    NOTE:
c       o array A is unchanged
c       o array Work is trashed
c       o factor(k) = (.25)**k for k=1,2,3,...,20
c-----------------------------------------------------------------------

      kMax = nOrder/2  ! number of iterations of 1-2-1 filter

      !--- k=1 iteration: filter entire region (including boundaries) --
      DO j=jMin,jMax
         FA(iMin(j),j) = 2.0*(A(iMin(j),j)-A(iMin(j)+1,j))
         FA(iMax(j),j) = 2.0*(A(iMax(j),j)-A(iMax(j)-1,j))
      END DO
      DO j=jMin,jMax
      DO i=iMin(j)+1,iMax(j)-1
        FA(i,j) = 2.0*A(i,j)-A(i-1,j)-A(i+1,j)
      END DO
      END DO

      !--- k=2,..,kMax: higher order filtering on decreasing subsets ---
      DO k=   2,kMax
      DO j=jMin,jMax
        DO i=iMin(j)+k-2,iMax(j)-k+2
          work(i,j) = FA(i,j)
        END DO
        FA(iMin(j)+k-2,j) = factor(k-1)*FA(iMin(j)+k-2,j)
        FA(iMax(j)-k+2,j) = factor(k-1)*FA(iMax(j)-k+2,j)
        DO i=iMin(j)+k-1,iMax(j)-k+1
          FA(i,j) = 2.0*Work(i,j)-Work(i-1,j)-Work(i+1,j)
        END DO
      END DO
      END DO

      !--- constant multiplication on full order subset ----------------
      DO j=jMin,jMax
      DO i=iMin(j)+kMax-1,iMax(j)-kMax+1
         FA(i,j) = factor(kMax)*FA(i,j)
      END DO
      END DO

      !--- set no-change BC's ------------------------------------------
      IF      (scheme .eq. 4 ) THEN
         !--- do nothing ---
      ELSE IF (scheme .eq. 7 ) THEN
         DO j=jMin,jMax
           FA(iMin(j),j) = 0.0
           FA(iMax(j),j) = 0.0
         END DO
      ELSE IF (scheme .eq. 5 ) THEN
         DO j=jMin,jMax
           FA(iMax(j),j) = 0.0
         END DO
      ELSE IF (scheme .eq. 6 ) THEN
         DO j=jMin,jMax
           FA(iMin(j),j) = 0.0
         END DO
      ELSE 
         write(*,*) 'ERROR (ShapiroR): invalid scheme flag =',scheme
         STOP 'ShapiroR'
      END IF


      RETURN
      END

c=======================================================================

      SUBROUTINE ShapiroC
     +  (A,FA,Work, nOrder,scheme, nx,ny, iMin,iMax,jMin,jMax)

      implicit none 

      integer   nOrder,scheme,nx,ny,iMin,iMax,jMin(nx),jMax(nx)
      real      A(nx,ny),FA(nx,ny),Work(nx,ny)

      integer   i,j,k,kmax

      real      factor(20)

      DATA      factor/
     + 2.500000E-01, 6.25000000E-02, 1.562500000E-02, 3.90625000000E-03,
     + 9.765625E-04, 2.44140625E-04, 6.103515625E-05, 1.52587890625E-05,
     + 3.814697E-06, 9.53674300E-07, 2.384186000E-07, 5.96046400000E-08,
     + 1.490116E-08, 3.72529000E-09, 9.313226000E-10, 2.32830600000E-10,
     + 5.820766E-11, 1.45519200E-11, 3.637979000E-12, 9.09494700000E-13/

      SAVE factor
 
c-----------------------------------------------------------------------
c   PURPOSE:
c      This routine applies a 1d shapiro filter to columns iMin through 
c      iMax of the array A, and returns the correction array as FA
c
c   INPUT:
c      A      : array to be filtered
c      nOrder : order of the filter = 2*k, for k in {1,2,3,...,20}
c      scheme : flags the type of boundary scheme to be used...
c         nnn1 => max fixed     , nnn0 => filtered
c         nn1n => min fixed     , nn0n => filtered
c         n1nn => reduced order , n0nn => full order
c         1nnn => conservative  , 0nnn => *possibly* non-conservative
c         Thus we have:
c         0100 = 4 => cons?, reduced order, min filtered, max filtered
c         0101 = 5 => cons?, reduced order, min filtered, max fixed
c         0110 = 6 => cons?, reduced order, min fixed,    max filtered
c         0111 = 7 => cons?, reduced order, min fixed   , max fixed   
c      nx,ny  : array dimensions
c      iMin,iMax,jMin(i),jMax(i) : boundaries of region to be filtered
c
c    OUTPUT:
c       FA : correction array (ie. filtered field = A-FA)
c
c    NOTE:
c       o array A is unchanged
c       o array Work is trashed
c       o factor(k) = (.25)**k for k=1,2,3,...,20
c-----------------------------------------------------------------------


      kMax = nOrder/2  ! number of iterations of 1-2-1 filter

      !--- k=1 iteration: filter entire region (including boundaries) --
      DO i=iMin,iMax
         FA(i,jMin(i)) = 2.0*(A(i,jMin(i))-A(i,jMin(i)+1))
         FA(i,jMax(i)) = 2.0*(A(i,jMax(i))-A(i,jMax(i)-1))
      END DO
      DO i=iMin,iMax
      DO j=jMin(i)+1,jMax(i)-1
        FA(i,j) = 2.0*A(i,j)-A(i,j-1)-A(i,j+1)
      END DO
      END DO

      !--- k=2,..,kMax: higher order filtering on decreasing subsets ---
      DO k=2, kMax
      DO i=iMin,iMax
         DO j=jMin(i)+k-2,jMax(i)-k+2 
           work(i,j)=FA(i,j)
         END DO
         FA(i,jMin(i)+k-2) = factor(k-1)*Work(i,jMin(i)+k-2)
         FA(i,jMax(i)-k+2) = factor(k-1)*Work(i,jMax(i)-k+2)
         DO j=jMin(i)+k-1,jMax(i)-k+1
           FA(i,j) = 2.0*Work(i,j)-Work(i,j-1)-Work(i,j+1)
         END DO
      END DO
      END DO

      !--- constant multiplication on full order subset ----------------
      DO i=iMin,iMax
      DO j=jMin(i)+kMax-1,jMax(i)-kMax+1
         FA(i,j) = factor(kMax)*FA(i,j)
      END DO
      END DO

      !--- set no-change BC's ------------------------------------------
      IF      (scheme .eq. 4) THEN
         !-- do nothing ---
      ELSE IF (scheme .eq. 7) THEN
         DO i=iMin,iMax
           FA(i,jMin(i)) = 0.0
           FA(i,jMax(i)) = 0.0
         END DO
      ELSE IF (scheme .eq. 5) THEN
         DO i=iMin,iMax
           FA(i,jMax(i)) = 0.0
         END DO
      ELSE IF (scheme .eq. 6) THEN
         DO i=iMin,iMax
           FA(i,jMin(i)) = 0.0
         END DO
      ELSE 
         write(*,*) 'ERROR (ShapiroC) invalid scheme flag =',scheme
         STOP 'ShapiroC'
      END IF

      RETURN
      END
c=======================================================================
