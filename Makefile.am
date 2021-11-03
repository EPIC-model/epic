SUBDIRS =		\
	src		\
	unit-tests	\
	python-scripts	\
	models		\
	analyse		\
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
	examples/taylorgreen2d.nml


# make python script an executable script
install-exec-hook:
	chmod 755 $(bindir)/animate-parcels.py
	chmod 755 $(bindir)/plot-parcels.py
	chmod 755 $(bindir)/plot-diagnostics.py
	chmod 755 $(bindir)/plot-timings.py
	chmod 755 $(bindir)/merge-ellipses.py
	chmod 755 $(bindir)/epic-run.py
	chmod 755 $(bindir)/plot-mergers.py
