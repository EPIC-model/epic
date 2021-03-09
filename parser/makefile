sources = types.f90 parser.f90
objects = types.o parser.o


$(objects): $(sources)
	$(FC) $(FFLAGS) $(sources) -c

.PHONY: clean
clean:
	rm -f *.o *.mod
