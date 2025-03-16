SUBDIRS = tools

pybindir=$(bindir)/pybin

if ENABLE_PYTHON_CHECK
animate_parcels_PYTHON = animate-parcels.py
animate_parcelsdir = $(pybindir)

plot_parcels_PYTHON = plot-parcels.py
plot_parcelsdir = $(pybindir)

plot_diagnostics_PYTHON = plot-diagnostics.py
plot_diagnosticsdir = $(pybindir)

plot_timings_PYTHON = plot-timings.py
plot_timingsdir = $(pybindir)

merge_ellipses_PYTHON = merge-ellipses.py
merge_ellipsesdir = $(pybindir)

subgrid_2d_generalised_PYTHON = subgrid_2d_generalised.py
subgrid_2d_generaliseddir = $(pybindir)

compare_results_PYTHON = compare_results.py
compare_resultsdir = $(pybindir)

plot_slices_PYTHON = plot_slices.py
plot_slicesdir = $(pybindir)
else
nobase_pybin_DATA =           \
    animate-parcels.py        \
    plot-parcels.py           \
    plot-diagnostics.py       \
    plot-timings.py           \
    merge-ellipses.py         \
    subgrid_2d_generalised.py \
    compare_results.py        \
    plot_slices.py
endif
