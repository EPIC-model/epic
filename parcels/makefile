sources = *.f90
objects = *.o


$(objects): $(sources)
	$(FC) $(FFLAGS) $(sources) -c

.PHONY: clean
clean:
	rm -f *.o *.mod
