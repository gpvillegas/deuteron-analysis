#
# kin_sig_int.1 uses full 3d interolation
#

FC            = gfortran


FFLAGS        = -ffixed-line-length-none -w -fno-automatic -g -fbounds-check

 
# set HERE to the current directory
HERE := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))

# print info on home directory setting here
$(info Home directory set to $(HERE))


CROSEC        = $(HERE)

CERNLIB       =

DEST	      = 

HDRS	      = 

LIBDIR        = -L $(HERE)/lib/gfortran


LIBS	      =

LINKER	      = gfortran


#suffix rules

%.o: %.f
	$(FC) $(FFLAGS) -c $<


OBJS          = 

PRINT	      = 

PROGRAM	      = heepcheck2025

SRCS          = heepcheck2025.f

all:		$(PROGRAM) 

$(PROGRAM):     $(OBJS)
		$(LINKER) $(LDFLAGS) $(LIBDIR) $(OBJS) $(LIBS) -o $(PROGRAM)  

clean:;		rm $(OBJS) $(PROGRAM)

program:        $(PROGRAM)