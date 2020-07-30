#! /bin/csh -f
#=======================================================================
# CVS $Id: ocn.setup.csh,v 1.2 1998/08/28 22:49:04 kauff Exp $
# CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/Coupled/ocn.setup.csh,v $
# CVS $Name: pac3_e $
#=======================================================================
#  This is a c shell script that builds a binary and creates an
#  input parameter namelist FOR THE COUPLED VERSION
#=======================================================================

echo '================================================================='
echo ' pac3 setup shell script                                         '
echo '================================================================='

echo '-----------------------------------------------------------------'
echo ' (a) Build an executable                                         '
echo '-----------------------------------------------------------------'
#rm -f ocn
if ( -e ocn) then
  echo 'Note: using an existing binary'
  ls -lFt ocn src/Build.log.* | head
else
  #--- create a src code sub-directory ---
  if (! -d src) mkdir src
  cd    src

  #--- document the build ---
  set BLDLOG = Build.log."`date +%y%m%d-%H%M%S`"
  echo 'Note: (re)building a binary'
  echo "See: $BLDLOG"
  echo "Build log "                        >!  $BLDLOG
  echo "Date: `date` "                     >>& $BLDLOG
  echo "Dir : `pwd` "                      >>& $BLDLOG
  echo "User: $LOGNAME "                   >>& $BLDLOG
  echo '---------------------------------' >>& $BLDLOG

  #--- gather source code ---
  if ( ! $?CSMROOT ) set CSMROOT = /fs/cgd/csm
  echo "$CSMROOT/csm_share/share1.3/src   "  >! Filepath
  echo "$CSMROOT/models/ocn/pac3.e        "  >> Filepath
  echo "$CSMROOT/models/ocn/pac3.e/Coupled"  >> Filepath
  foreach SRCDIR (`cat Filepath`)
    echo '---------------------------------' >>& $BLDLOG
    echo "o gathering src code from $SRCDIR" >>& $BLDLOG
    cat $SRCDIR/README                       >>& $BLDLOG
    ls -lF $SRCDIR                           >>& $BLDLOG
    cp -fp $SRCDIR/* .                       >>& $BLDLOG
  end

  #--- create make's include files & invoke make ---
  echo '---------------------------------' >>& $BLDLOG
  Makeprep                                 >>& $BLDLOG
  make EXEC=pac3 ARCH=$ARCH                >>& $BLDLOG || exit 2

  #--- document the source code used ---
  echo '---------------------------------' >>& $BLDLOG
  echo "o contents of `pwd`"               >>& $BLDLOG
  ls -alFt                                 >>& $BLDLOG
  echo '---------------------------------' >>& $BLDLOG
  echo "o revision control info:"          >>& $BLDLOG
  grep 'CVS' *.[hF]                        >>& $BLDLOG

  #--- link binary into ./ directory ---
  cd ..
  rm -f ocn ; ln -s  src/pac3 ocn
endif

echo ' '
echo '-----------------------------------------------------------------'
echo ' (b) document the source code used                               '
echo '-----------------------------------------------------------------'
echo "o contents of /src:"      ; ls -alFt   src        ; echo ' '
echo "o revision control info:" ; grep 'CVS' src/*.[hF] ; echo ' '

echo ' '
echo '-----------------------------------------------------------------'
echo ' (c) Create an input parameter namelist file                     '
echo '-----------------------------------------------------------------'
#et  DATADIR = /insert/data/directory/here
set  DATADIR = $CSMROOT/models/ocn/pac3.d/Coupled
cp -p   $DATADIR/*.nc   .
rm data.sst.nc    ; ln -s  sst.str.nc       data.sst.nc
rm data.domain.nc ; ln -s  domain.pac3g.nc  data.domain.nc

cat >! ocn.parm  << EOF
 &inparm
  case_name   = '$CASE '
  case_desc   = '$CASE $CASESTR '
  rest_type   = '$RUNTYPE'
  rest_pfile  = '$RPTDIR/pac3.$CASE.rpointer '
  rest_bfile  = '/KAUFF/pac3/data/r0000-01-01'
  rest_date   = 00010102
  rest_freq   = 'coupler'
  stop_option = 'coupler'
  hist_freq   = 'monthly'
  hist_tavg   = 1
  diag_freq   = 'monthly'
  mss_dir     =  '$MSSDIR/ocn'
  mss_pass    = ',$MSSPWD'
  mss_rtpd    =   $MSSRPD
 /
EOF

echo "o contents of ocn.parm:" ; cat ocn.parm ; echo ' '
echo "o contents of `pwd`:"    ; ls -alF      ; echo ' '

echo '================================================================='
echo ' end of setup shell script                                       '
echo '================================================================='
