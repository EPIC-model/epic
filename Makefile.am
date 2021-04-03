SUBDIRS = src test python-scripts
dist_doc_DATA = README.md

# make python script an executable script
install-exec-hook:
	chmod 755 $(bindir)/animate-ellipses.py
	chmod 755 $(bindir)/plot-ellipses.py
	chmod 755 $(bindir)/merge-ellipses.py
