SUBDIRS =		\
	src		\
	unit-tests	\
	mpi-tests	\
 	analyse		\
 	models

if ENABLE_PYTHON
SUBDIRS += python-scripts	\
	   run-script
endif

# 3 November 2021
# https://stackoverflow.com/questions/16657991/glibtoolize-on-macos-tells-me-to-consider-adding-i-m4-to-aclocal-amflags-but
ACLOCAL_AMFLAGS = -I m4

dist_doc_DATA = README.md

dataroot_DATA = 				\
	examples/taylor_green.config		\
	examples/robert.config			\
	examples/straka.config			\
	examples/robert2d-one-gaussian.nml	\
	examples/robert2d-one-uniform.nml	\
	examples/robert2d-two-gaussians.nml	\
	examples/straka2d.nml			\
	examples/taylorgreen2d.nml

if ENABLE_PYTHON
dataroot_DATA += examples/write_moist_setup.py


# make python script an executable script
install-exec-hook:
	chmod 755 $(bindir)/animate-parcels.py
	chmod 755 $(bindir)/plot-parcels.py
	chmod 755 $(bindir)/plot-diagnostics.py
	chmod 755 $(bindir)/plot-timings.py
	chmod 755 $(bindir)/merge-ellipses.py
	chmod 755 $(bindir)/epic-run.py
	chmod 755 $(bindir)/subgrid_2d_generalised.py
	chmod 755 $(bindir)/compare_results.py
	chmod 755 $(bindir)/plot_slices.py
endif
