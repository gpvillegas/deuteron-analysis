#
#  make response function pyrhont module respy
#

#------Makefile for wraping the boris tracker library ----------

# modules
SRC = response_functions_mod.f90  
OBJ = response_functions_mod.o  

F2PY     = f2py

F2PY_F1   = --overwrite-signature -m 
F2PY_F2   = -c --backend meson --fcompiler=gfortran --f90flags="-ffixed-line-length-none -w -fno-automatic "

PROGRAM = respy_fsi

#----------------------------------------------------------
all: $(OBJ)

%.o: %.f90
	$(F2PY) $(F2PY_F1) $(PROGRAM) -h sgn_$(PROGRAM).pyf $<
	# need to modify one line in the signature file to make sure that the bfield_array works correctly, of the code changes this needs to be adjusted
	# cat sgn_fluxpy.pyf | sed s/"real(kind=8), allocatable,dimension(:,:) :: b_array"/"real(kind=8), allocatable,dimension(size(r),4) :: b_array"/ > xx
	# mv xx sgn_$(PROGRAM).pyf
	$(F2PY) $(F2PY_F2) sgn_$(PROGRAM).pyf $(SRC) 

#----------------------------------------------------------
.PHONY : clean

clean:
	rm -f *.so *.pyf $(OBJ)
	rm -rf *.dSYM
