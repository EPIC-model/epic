export FFLAGS=-O3

sources = parameters.f90 init.f90 constants.f90 model.f90
objects = parameters.o init.o constants.o model.o

all: epic

SUBDIRS = parser parcels stepper
.PHONY: subdirs $(SUBDIRS)
subdirs: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@

# 9 March 2021
# https://stackoverflow.com/questions/15347543/recursive-clean-in-a-makefile
recursive-clean:
	for s in $(SUBDIRS); do cd "$$s" && make clean && cd ..; done

epic: $(SUBDIRS) $(objects) epic.f90
	$(FC) $(FFLAGS) parser/parser.o parcels/parcels.o stepper/rk4.o $(objects) epic.f90 -o epic

$(objects): $(sources)
	$(FC) $(FFLAGS) -I./parcels -I./parser -I./stepper $(sources) -c

.PHONY: clean
clean:
	rm -f *.o *.mod
	make recursive-clean
