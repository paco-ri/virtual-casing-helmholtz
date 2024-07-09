LIBNAME = libfmm3dbie_matlab
DYNAMICLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a

COM = src/common
COMOBJS = $(COM)/hkrand.o $(COM)/dotcross3d.o \
	$(COM)/dlaran.o \
	$(COM)/rotmat_gmres.o $(COM)/setops.o \
	$(COM)/sort.o $(COM)/sparse_reps.o $(COM)/get_fmm_thresh.o \
	$(COM)/common_Maxwell.o \
	$(COM)/rigidbodies.o $(COM)/polytens.o \
	$(COM)/chebexps.o $(COM)/gmres_routs.o

EM = src/maxwell
EMOBJS = $(EM)/em_mfie_pec.o $(EM)/em_aumfie_pec.o \
	$(EM)/em_nrccie_pec.o $(EM)/em_auCKi_pec.o \
	$(EM)/em_dfie_trans.o $(EM)/em_adpie_pec.o \
	$(EM)/em_sdpie_pec.o $(EM)/em_cfie_rwg_pec.o \
	$(EM)/maxwell_common_evaluators.o \
	$(EM)/incoming_fields.o \
	$(EM)/fix_tri.o $(EM)/analytic_sphere_pw_pec.o

HELM = src/helm_wrappers
HOBJS = $(HELM)/helm_comb_dir.o $(HELM)/helm_rpcomb_neu.o \
	$(HELM)/helm_comb_trans.o $(HELM)/helm_rpcomb_imp.o \
	$(HELM)/helm_s_neu.o $(HELM)/helm_common_evaluators.o

KER = src/kernels
KOBJS = $(KER)/helm_kernels.o $(KER)/lap_kernels.o $(KER)/DPIE_kernels.o \
	$(KER)/yuk_kernels.o $(KER)/stok_kernels.o $(KER)/em_kernels.o

LAP = src/lap_wrappers
LOBJS = $(LAP)/lap_comb_dir.o $(LAP)/lap_s_neu.o

QUAD = src/quadratures
QOBJS = $(QUAD)/far_field_routs.o \
	$(QUAD)/ggq-selfquad-routs.o $(QUAD)/ggq-quads.o \
	$(QUAD)/ggq-selfquad.o \
	$(QUAD)/near_field_routs.o $(QUAD)/near_quad_sub.o

SURF = src/surface_routs
SOBJS = $(SURF)/in_go3.o $(SURF)/surf_routs.o $(SURF)/vtk_routs.o \
	$(SURF)/xtri_parameterizations.o \
	$(SURF)/xtri_plot.o $(SURF)/write_go3.o $(SURF)/in_gidmsh2.o \
	$(SURF)/in_gmsh2.o $(SURF)/patch_basis_routs.o \
	$(SURF)/analytic_geometry_routs.o $(SURF)/analytic_charts.o \
	$(SURF)/xquad_parametrizations.o \

TRIA = src/tria_routs
TOBJS = $(TRIA)/ctriaints_main.o $(TRIA)/koornexps.o \
	$(TRIA)/triaintrouts.o $(TRIA)/dtriaints_main.o \
	$(TRIA)/triasymq.o $(TRIA)/triatreerouts.o $(TRIA)/dtriaintrouts.o

STOK = src/stok_wrappers
STOKOBJS = $(STOK)/stok_comb_vel.o

QUAD2 = src/quad_routs
QOBJS2 = $(QUAD2)/cquadints_main.o \
	$(QUAD2)/cquadintrouts.o $(QUAD2)/dquadints_main.o \
	$(QUAD2)/squarearbq.o $(QUAD2)/quadtreerouts.o $(QUAD2)/dquadintrouts.o

OBJS = $(COMOBJS) $(EMOBJS) $(HOBJS) $(KOBJS) $(LOBJS) $(QOBJS) $(SOBJS) $(TOBJS) $(STOKOBJS) $(QOBJS2)

FC = gfortran
FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy -w
OMPFLAGS = -fopenmp
FFLAGS += $(OMPFLAGS)
FFLAGS += -J .mod/

%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@

FMM_INSTALL = ${HOME}/lib/libfmm3d.a
FMMBIE_INSTALL = ${HOME}/lib/libfmm3dbie.a

install:
	cd lib && rm -rf *
	cd lib-static && rm -rf * && ar -x $(FMMBIE_INSTALL) && cp ../src/magneto-static-routs.o . && cp ../src/surf_routs.o surf_routs2.o && ar rcs libvirtualcasing.a *.o && rm -f *.o
	gfortran -shared -fPIC -O3 -march=native -funroll-loops -std=legacy -w -fopenmp -J .mod/ -Wl,--whole-archive lib-static/libvirtualcasing.a -Wl,--no-whole-archive -o libvirtualcasing.so -lm -lstdc++ -lgomp -lblas -llapack
	mv libvirtualcasing.so lib/

test_gradcurllap:
	gfortran -o test/test_gradcurllap test/test_gradcurllap.f90 lib-static/libvirtualcasing.a -fallow-argument-mismatch -fPIC -O3 -march=native -funroll-loops -std=legacy -w -fopenmp -lm -lstdc++ -lgomp -lblas -llapack
	test/test_gradcurllap

