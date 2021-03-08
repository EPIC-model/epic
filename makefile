f90=gfortran
flags=-O3

sources = parser.f90 constants.f90 model.f90
objects = parser.o constants.o model.o


all: epic

epic: $(objects) epic.f90
	$(f90) $(flags) $(objects) epic.f90 -o epic

$(objects): $(sources)
	$(f90) $(sources) -c

clean:
	rm -f *.o *.mod
