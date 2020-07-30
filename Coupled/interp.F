c=======================================================================
c CVS $Id: interp.F,v 1.1 1998/08/28 20:52:56 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/Coupled/interp.F,v $
c CVS $Name: pac3_e $
c=======================================================================

      SUBROUTINE Interp2D(Xin ,Yin ,Ain ,mx,my,
     +                    Xout,Yout,Aout,nx,ny,spv,mFLAG)

      implicit none

      !----- input -----
      integer mx,my              ! input array dimensions
      real    Xin (mx),Yin(my)   ! input coordinate arrays
      real    Ain (mx,my)        ! data on input grid
      real    spv                ! denotes invalid data
      integer nx,ny              ! output array dimensions
      real    Xout(nx),Yout(ny)  ! output coordinate arrays
      integer mFLAG              ! selects extrapolation option

      !----- output -----
      real    Aout(nx,ny)        ! data on output grid

      !----- local -----
      integer i,j,ii,jj       ! loop indicies
      real    a,b,c,d         ! ... as defined in algorithm below
      real    dx,dy           ! ... as defined in algorithm below
      real    alpha,beta      ! ... as defined in algorithm below
      real    Ca,Cb,Cc,Cd     ! bilinear mapping weights 
      real    XL,XM,XH        ! XL=Xin(ii) =< XM=Xout(i) =< XH=Xin(ii+1)
      real    YL,YM,YH        ! YL=Yin(ii) =< YM=Yout(i) =< YH=Yin(ii+1)

c-----------------------------------------------------------------------
c     PURPOSE:
c        map: Ain -> Aout
c        where Ain (i,j) has real world coordinates (Xin (i),Yin (j))
c        and   Aout(i,j) has real world coordinates (Xout(i),Yout(j))
c
c     NOTE:
c        mFLAG = 0 => map only if there are 4 valid values
c                     adjacent to the the output grid point
c        mFLAG = 1 => extrapolate if there is at least one valid value
c                     adjacent to the the output grid point
c        mFLAG = other => assume there are no special values
c
c     ASSUMPTION: 
c        Xin(1) <= Xout(1) < Xout(nx) <= Xin(mx)
c        Yin(1) <= Yout(1) < Yout(ny) <= Yin(my)
c
c     ALGORITHM:
c       for each i & j, 
c         find ii s.t. Xin(ii) < Xout(i) <= Xin(ii) and
c         find jj s.t. Yin(jj) < Yout(j) <= Yin(jj)
c       then 
c         a = F(Xin(ii  ),Yin(jj  ))
c         b = F(Xin(ii+1),Yin(jj  ))
c         c = F(Xin(ii+1),Yin(jj+1))
c         d = F(Xin(ii  ),Yin(jj+1))
c         dx = Xin(ii+1)-Xin(ii)
c         dy = Yin(jj+1)-Yin(jj)
c         alpha = (Xout(i)-Xin(ii))/dx
c         beta  = (Yout(j)-Yin(jj))/dy
c       with this info, interpolate to get the value
c         Aout(i,j) = F(Xout(i),Yout(j))
c
c-----------------------------------------------------------------------
c           ^  d       c
c       b   |    
c       e  dy      f                interpolate f fromm a,b,c,d
c       t   |                       (avoid spval's)
c       a   v  a       b
c              <--dx--->
c               alpha
c-----------------------------------------------------------------------

      IF (mFLAG.eq.0) THEN
        jj=1
        DO j=1,ny
          DO WHILE (Yin(jj+1) .lt. Yout(j))
            jj=jj+1
          END DO
          ii=1
          DO i=1,nx
            DO WHILE (Xin(ii+1) .lt. Xout(i))
              ii=ii+1
            END DO
            a = Ain(ii  ,jj  )
            b = Ain(ii+1,jj  )
            c = Ain(ii+1,jj+1)
            d = Ain(ii  ,jj+1)
            IF (a.eq.spv .or. b.eq.spv .or .c.eq.spv .or. d.eq.spv) THEN
              Aout(i,j) = spv
            ELSE
              alpha = ( Xout(i)-Xin(ii) )/( Xin(ii+1)-Xin(ii) )
              beta  = ( Yout(j)-Yin(jj) )/( Yin(jj+1)-Yin(jj) )
              Ca = (1.-alpha)*(1.-beta)
              Cb = (   alpha)*(1.-beta)
              Cc = (   alpha)*(   beta)
              Cd = (1.-alpha)*(   beta)
              Aout(i,j) = Ca*a + Cb*b + Cc*c + Cd*d
            END IF
          END DO
        END DO
      ELSE IF (mFLAG.eq.1) THEN
        jj=1
        DO j=1,ny
          DO WHILE (Yin(jj+1) .lt. Yout(j))
            jj=jj+1
          END DO
          ii=1
          DO i=1,nx
            DO WHILE (Xin(ii+1) .lt. Xout(i))
              ii=ii+1
            END DO
            a = Ain(ii  ,jj  )
            b = Ain(ii+1,jj  )
            c = Ain(ii+1,jj+1)
            d = Ain(ii  ,jj+1)
            XL = Xin (ii  )
            XM = Xout(i   )
            XH = Xin (ii+1)
            YL = Yin (jj  )
            YM = Yout(j   )
            YH = Yin (jj+1)
            call IntrpBox(a,b,c,d, XL,XM,XH,YL,YM,YH, spv, Ca,Cb,Cc,Cd)
            Aout(i,j) = Ca*a + Cb*b + Cc*c + Cd*d
          END DO
        END DO
      ELSE
        jj=1
        DO j=1,ny
          DO WHILE (Yin(jj+1) .lt. Yout(j))
            jj=jj+1
          END DO
          ii=1
          DO i=1,nx
            DO WHILE (Xin(ii+1) .lt. Xout(i))
              ii=ii+1
            END DO
            a = Ain(ii  ,jj  )
            b = Ain(ii+1,jj  )
            c = Ain(ii+1,jj+1)
            d = Ain(ii  ,jj+1)
              alpha = ( Xout(i)-Xin(ii) )/( Xin(ii+1)-Xin(ii) )
              beta  = ( Yout(j)-Yin(jj) )/( Yin(jj+1)-Yin(jj) )
              Ca = (1.-alpha)*(1.-beta)
              Cb = (   alpha)*(1.-beta)
              Cc = (   alpha)*(   beta)
              Cd = (1.-alpha)*(   beta)
              Aout(i,j) = Ca*a + Cb*b + Cc*c + Cd*d
          END DO
        END DO
      END IF

      RETURN
      END

c=======================================================================

      SUBROUTINE IntrpBox(a,b,c,d, XL,XM,XH, YL,YM,YH, spv, Ca,Cb,Cc,Cd)

      implicit none

      !----- input -----  
      real    a,b,c,d
      real    XL,XM,XH
      real    YL,YM,YH
      real    spv

      !----- output -----  
      real    Ca,Cb,Cc,Cd

      !----- local -----  
      integer ID
      real    dy,dx
      real    alpha ,beta
      real    alpha2,beta2

c-----------------------------------------------------------------------
c           ^  d       c
c       b   |    
c       e  dy      f                interpolate f from a,b,c,d
c       t   |                       (avoiding special values)
c       a   v  a       b
c              <--dx--->
c               alpha
c-----------------------------------------------------------------------
c     PURPOSE:  
c        given XL<=XM<=XH
c              YL<=YM<=YH
c              a=F(XL,YL)
c              b=F(XH,YL)
c              c=F(XH,YH)
c              d=F(XL,YH)
c        find Ca,Cb,Cc,Cd to get an interpolated value for F(XM,YM)
c             F(XM,YM) = Ca*a + Cb*b + Cc*c + Cd*d
c     INPUT:  
c        a,b,c,d,
c        XL,XM,XH
c        YL,YM,YH
c        spv
c     OUTPUT: 
c        Ca,Cb,Cc,Cd
c-----------------------------------------------------------------------

      ID = 0
      IF (a .ne. spv) ID=ID+1
      IF (b .ne. spv) ID=ID+2
      IF (c .ne. spv) ID=ID+4
      IF (d .ne. spv) ID=ID+8

      Ca=0.0
      Cb=0.0
      Cc=0.0
      Cd=0.0

*     IF (ID .gt. 2) THEN
        dx    =  XH-XL
        dy    =  YH-YL
        alpha = (XM-XL)/dx
        beta  = (YM-YL)/dy
*     END IF

      !--- interpolate with four points --------------------------------
      IF (ID .eq. 15) THEN
         Ca = (1.-alpha)*(1.-beta)
         Cb = (   alpha)*(1.-beta)
         Cc = (   alpha)*(   beta)
         Cd = (1.-alpha)*(   beta)
      !--- interpolate with no points ----------------------------------
      ELSE IF (ID .eq. 0) THEN
         Ca=1.0
      !--- interpolate with one point ----------------------------------
      ELSE IF (ID .eq. 1) THEN
         Ca=1.0
      ELSE IF (ID .eq. 2) THEN
         Cb=1.0
      ELSE IF (ID .eq. 4) THEN
         Cc=1.0
      ELSE IF (ID .eq. 8) THEN
         Cd=1.0
      !--- interpolate with two adjacent points ------------------------
      ELSE IF (ID .eq.  3) THEN
         Ca=1.-alpha
         Cb=alpha
      ELSE IF (ID .eq. 12) THEN
         Cd=1.-alpha
         Cc=alpha
      ELSE IF (ID .eq.  6) THEN
         Cb = 1.-beta
         Cc = beta
      ELSE IF (ID .eq.  9) THEN
         Ca = 1.-beta 
         Cd = beta
      !--- interpolate with two opposite points ------------------------
      ELSE IF (ID .eq. 5) THEN
         alpha2 = sqrt((alpha*dx)**2+(beta*dy)**2)/sqrt(dx**2+dy**2)
         Ca = 1.-alpha2
         Cc = alpha2
      ELSE IF (ID .eq. 10) THEN
         alpha2= sqrt(((1.-alpha)*dx)**2+(beta*dy)**2)/sqrt(dx**2+dy**2)
         Cb = 1.-alpha2 
         Cd = alpha2
      !--- interpolate with three points -------------------------------
      ELSE IF (ID .eq. 11) THEN
         Ca = 1. - alpha - beta
         Cb = alpha
         Cd = beta
      ELSE IF (ID .eq. 7) THEN
         Ca = 1. - alpha
         Cb = alpha -beta
         Cc = beta
      ELSE IF (ID .eq. 13) THEN
         Ca = 1. - beta
         Cc = alpha
         Cd = -alpha + beta
      ELSE IF (ID .eq. 14) THEN
         Cc = alpha + beta -1.
         Cb = 1. - beta
         Cd = 1. - alpha
      ELSE
         write(*,*) '(InterpBox) FATAL ERROR'
         STOP 'IntrpBx'
      END IF

      RETURN
      END

c=======================================================================
