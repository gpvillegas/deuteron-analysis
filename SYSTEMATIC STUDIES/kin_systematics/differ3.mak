
#
# differ3
#

FC            =  gfortran

#FFLAGS        = -f -C -g -s 
FFLAGS	       = -g

#suffix rule necessary for compiling with absoft

.f.o :
	$(FC) $(FFLAGS) -c $<

HERE          = ./

CROSEC        = ./

CERNLIB       =

DEST	      = 

HDRS	      = 

# LIBDIR        = -L/data/boeglin.1/lib/gfortran

LDFLAGS       = 

LIBS	      = 

LINKER	      = gfortran



OBJS          = differ3.o phaslib.1.o wlib.o

PRINT	      = 

PROGRAM	      = differ3

SRCS          = differ3.f phaslib.1.f wlib.f

all:		$(PROGRAM) 

$(PROGRAM):     $(OBJS)
		$(LINKER) $(LDFLAGS) $(LIBDIR) $(OBJS) $(LIBS) -o $(PROGRAM)  

clean:;		rm $(OBJS)

program:        $(PROGRAM)










