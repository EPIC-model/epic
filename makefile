export FFLAGS=-O3

sources = types.f90 parser.f90 constants.f90 model.f90
objects = types.o parser.o constants.o model.o


all: epic

SUBDIRS = parcels stepper

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

epic: $(SUBDIRS) $(objects) epic.f90
	$(FC) $(FFLAGS) $(objects) epic.f90 -o epic

$(objects): $(sources)
	$(FC) $(FFLAGS) $(sources) -c

.PHONY: clean
clean:
	$(MAKE) -C $(SUBDIRS) clean
	rm -f *.o *.mod
