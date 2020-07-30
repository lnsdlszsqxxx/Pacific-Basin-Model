#! /bin/csh -f
#=======================================================================
# CVS $Id: pac3.setup.csh,v 1.4 1998/08/28 22:48:46 kauff Exp $
# CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/pac3.setup.csh,v $
# CVS $Name: pac3_e $
#=======================================================================
#  This is a c shell script that builds a binary and creates an
#  input parameter namelist FOR THE UNCOUPLED VERSION
#=======================================================================

echo '================================================================='
echo ' pac3 setup shell script (pac3.setup.csh)                        '
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
  echo "See: ./src/$BLDLOG"
  echo "Build log "                        >!  $BLDLOG
  echo "Date: `date` "                     >>& $BLDLOG
  echo "Dir : `pwd` "                      >>& $BLDLOG
  echo "User: $LOGNAME "                   >>& $BLDLOG
  echo '---------------------------------' >>& $BLDLOG

  #--- gather source code ---
  if ( ! $?CSMROOT ) set CSMROOT = /Users/liangyu
  echo "$CSMROOT/pac3.e/src        "  >! Filepath #liang, change source dir
  foreach SRCDIR (`cat Filepath`)
    echo '---------------------------------' >>& $BLDLOG
    echo "o gathering src code from $SRCDIR" >>& $BLDLOG
#liang, that's why I moved all source files in ./src
    cat $SRCDIR/README                       >>& $BLDLOG
    ls -lF $SRCDIR                           >>& $BLDLOG
    cp -fp $SRCDIR/* .                       >>& $BLDLOG
  end

  #--- create make's include files & invoke make ---
  echo '---------------------------------' >>& $BLDLOG
  echo 'calling Makeprep-----------------' >>& $BLDLOG
  csh -f Makeprep                          >>& $BLDLOG #liang ,add "csh -f"
  echo 'calling Make---------------------' >>& $BLDLOG
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
cat >! ocn.parm  << EOF
 &inparm
  case_name   = '$CASE '
  case_desc   = '$CASE $CASESTR '
  rest_type   = '$RUNTYPE'
  rest_pfile  = '$RPTDIR/pac3.$CASE.rpointer '
  rest_bfile  = '/KAUFF/pac3/data/r0000-01-01'
  rest_date   = 00000101
  rest_freq   = 'monthly'
  stop_option = 'newmonth'
  hist_freq   = 'monthly'
  hist_tavg   = 1
  diag_freq   = 'monthly'
  mss_dir     =  '$MSSDIR'
  mss_pass    = ',$MSSPWD'
  mss_rtpd    =   $MSSRPD
  frc_wfile   = '/KAUFF/pac3/data/wind1.fsu.mnth'
  frc_cfile   = '/KAUFF/pac3/data/cloud.isccp.mnth'
  frc_tfile   = '/KAUFF/pac3/data/sst.str.mnth'
 /
EOF

echo "o contents of ocn.parm:" ; cat ocn.parm ; echo ' '
echo "o contents of `pwd`:"    ; ls -alF      ; echo ' '

echo '================================================================='
echo ' end of setup shell script                                       '
echo '================================================================='
