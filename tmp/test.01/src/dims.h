c=======================================================================
c CVS $Id: dims.h,v 1.1.1.1 1997/04/17 23:00:00 kauff Exp $
c CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/dims.h,v $
c CVS $Name: pac3_e $
c=======================================================================

      integer      nx    ! number of zonal      grid cells
      integer      ny    ! number of meridional grid cells
      integer      nz    ! number of vertical   grid cells
      integer      nxp1  ! nx+1
      integer      nyp1  ! ny+1
      integer      nzp1  ! nz+1
      integer      nz2p1 ! 2*nz + 1

      parameter  ( nx    =  166   )
      parameter  ( ny    =   80   )
      parameter  ( nz    =   15   )
      parameter  ( nxp1  = nx+1   )
      parameter  ( nyp1  = ny+1   )
      parameter  ( nzp1  = nz+1   )
      parameter  ( nz2p1 = 2*nz+1 )

