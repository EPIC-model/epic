if ENABLE_UNIT_TESTS

AM_FCFLAGS = -I $(top_builddir)/src/

AM_LDFLAGS =					\
	$(top_builddir)/src/parcels/*.o 	\
	$(top_builddir)/src/taylorgreen.o	\
	$(top_builddir)/src/utils/*.o		\
	$(top_builddir)/src/fields/*.o

unittestsdir = $(prefix)/unit-tests

unittests_PROGRAMS = 				\
	test_ellipse_orientation 		\
	test_free_slip				\
	test_fft_1				\
	test_fft_2				\
	test_trilinear				\
	test_diverge				\
	test_ellipse_split			\
	test_ellipse_bi_merge			\
	test_ellipse_multi_merge_1		\
	test_ellipse_multi_merge_2		\
	test_ellipse_multi_merge_3		\
	test_ellipse_multi_merge_symmetry	\
	test_ellipse_bi_vs_multi_merge		\
	test_tri_inversion

# 7 April 2021
# https://www.gnu.org/software/automake/manual/html_node/Scripts.html
unittests_SCRIPTS = run_unit_tests.sh

test_ellipse_orientation_SOURCES = test_ellipse_orientation.f90 unit_test.f90
test_free_slip_SOURCES = test_free_slip.f90 unit_test.f90
test_fft_1_SOURCES = test_fft_1.f90 unit_test.f90
test_fft_2_SOURCES = test_fft_2.f90 unit_test.f90
test_trilinear_SOURCES = test_trilinear.f90 unit_test.f90
test_diverge_SOURCES = test_diverge.f90 unit_test.f90
test_ellipse_split_SOURCES = test_ellipse_split.f90 unit_test.f90
test_ellipse_bi_merge_SOURCES = test_ellipse_bi_merge.f90 unit_test.f90
test_ellipse_multi_merge_1_SOURCES = test_ellipse_multi_merge_1.f90 unit_test.f90
test_ellipse_multi_merge_2_SOURCES = test_ellipse_multi_merge_2.f90 unit_test.f90
test_ellipse_multi_merge_3_SOURCES = test_ellipse_multi_merge_3.f90 unit_test.f90
test_ellipse_multi_merge_symmetry_SOURCES = test_ellipse_multi_merge_symmetry.f90 unit_test.f90
test_ellipse_bi_vs_multi_merge_SOURCES = test_ellipse_bi_vs_multi_merge.f90 unit_test.f90
test_tri_inversion_SOURCES = test_tri_inversion.f90

endif
