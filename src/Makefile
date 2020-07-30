#=======================================================================
# CVS $Id: Makefile,v 1.2 1998/08/21 22:54:03 kauff Exp $
# CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/Makefile,v $
# CVS $Name: pac3_e $
#=======================================================================
# PURPOSE:
#   this is a generic makefile
#
# NOTES:
#   it is assumed that that following include files are available...
#   o Objects        - defines OBJS = list of *.o file
#   o Macros.$(ARCH) - defines CPP,FC,LDR,SLIBS,MLIBS
#   o Depends        - defines *.o file dependencies
#-----------------------------------------------------------------------

EXEC = a.out

include Objects

include Macros.liang

.SUFFIXES:
.SUFFIXES: .F .o

#-----------------------------------------------------------------------

$(EXEC): $(OBJS)
	$(LDR) -o $(EXEC) $(OBJS) $(SLIBS) $(MLIBS)

.F.o:
	$(FC) $*.f

clean:
	rm -f $(OBJS) *.i $(EXEC)

#-----------------------------------------------------------------------

include Depends

#-----------------------------------------------------------------------

