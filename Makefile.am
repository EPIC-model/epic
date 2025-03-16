SUBDIRS =		\
	src		\
	unit-tests	\
	mpi-tests	\
 	analyse		\
 	models		\
 	python-scripts	\
	run-script

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
	examples/taylorgreen2d.nml		\
	examples/write_moist_setup.py


# make python script an executable script
pybindir=$(bindir)/pybin
install-exec-hook:
	chmod 755 $(pybindir)/animate-parcels.py
	chmod 755 $(pybindir)/plot-parcels.py
	chmod 755 $(pybindir)/plot-diagnostics.py
	chmod 755 $(pybindir)/plot-timings.py
	chmod 755 $(pybindir)/merge-ellipses.py
	chmod 755 $(pybindir)/epic-run.py
	chmod 755 $(pybindir)/subgrid_2d_generalised.py
	chmod 755 $(pybindir)/compare_results.py
	chmod 755 $(pybindir)/plot_slices.py
