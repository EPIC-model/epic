SUBDIRS = src unit-tests python-scripts models
dist_doc_DATA = README.md

dataroot_DATA = 				\
	taylor_green.config			\
	robert.config				\
	straka.config				\
	models/robert2d-one-gaussian.nml	\
	models/robert2d-one-uniform.nml		\
	models/robert2d-two-gaussians.nml	\
	models/straka2d.nml			\
	models/taylorgreen2d.nml


# make python script an executable script
install-exec-hook:
	chmod 755 $(bindir)/animate-parcels.py
	chmod 755 $(bindir)/bokeh-parcels.py
	chmod 755 $(bindir)/plot-parcels.py
	chmod 755 $(bindir)/plot-fields.py
	chmod 755 $(bindir)/plot-diagnostics.py
	chmod 755 $(bindir)/plot-timings.py
	chmod 755 $(bindir)/merge-ellipses.py
