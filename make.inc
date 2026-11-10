# STRATA-NOVA compiler defaults.
# Command-line overrides are supported, e.g.:
#   make F90=gfortran-14 F90FLAGS='...'
# The default flags reproduce the campaign-ready production build.
export F90 ?= gfortran
export F90FLAGS ?= -O3 -fprotect-parens -march=native -flto=auto -fopenmp
export LDFLAGS ?= $(F90FLAGS)
export LAPACKLIBS ?= -llapack -lblas
export LIBS += $(LAPACKLIBS)
