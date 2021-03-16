sources = parameters.f90 constants.f90 #Main f90 sources

objects = parameters.o constants.o #Main object modules

 # Set flag options:
norm_flags = -O3 -mcmodel=medium #Compiler flags
test_flags = -O3 -mcmodel=medium -fbounds-check -fbacktrace  #Compiler flags

 # Set compiler and flag used:
f90=gfortran
flags=$(norm_flags)

#----------------------------------------------------------------------
 # Rules:
all: nearest ranpar

nearest: $(objects) nearest.f90
	skin -pomp $(f90) $(objects) nearest.f90 -o nearest $(flags)

ranpar: $(objects) ranpar.f90
	$(f90) $(objects) ranpar.f90 -o ranpar $(flags)

$(objects): $(sources)
	$(f90) $(sources) -c

.PHONY: clean
clean:
	rm -f *.o *.mod
