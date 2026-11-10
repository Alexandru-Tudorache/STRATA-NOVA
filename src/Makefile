.PHONY: all clean snec FORCE

include ../make.inc


SOURCES=observer_projection.F90\
	filesystem_utils.F90\
	mixing_utils.F90\
	startup_log.F90\
	composition_profile_helpers.F90\
	composition_profile_store.F90\
	csm_switch.F90\
	legacy_grid_pattern.F90\
	cadence_control.F90\
	staged_state.F90\
	simulation_manager.F90\
	problem.F90\
	allocate_vars.F90\
	eos_switch.F90\
	eos_content.F90\
	grid.F90\
	read_profile.F90\
	read_profile_compositions_custom.F90\
	output.F90\
	timestep.F90\
	blstep.F90\
	artificial_viscosity.F90\
	opacity.F90\
	bicubic_interpolator.F90\
	hydro.F90\
	hydro_rad.F90\
	outinfo.F90\
	input_parser.F90\
	conservation.F90\
	analysis.F90\
	read_BolCorr.F90\
	opal_opacity.F90\
	nickel.F90\
	simple_saha.F90\
	saha_tables.F90\
	luminosity.F90\
	bomb_profile.F90\
	arrays.F90\
	boxcar.F90\
	shock_capture.F90\
	snec.F90


CLEANSTUFF=rm -f *.o *.mod .build-config*

OBJECTS=$(SOURCES:.F90=.o)

all: ../snec

snec: ../snec

# Rebuild when the selected compiler or flags change, including after debug builds.
.build-config: FORCE
	@set -e; config_file=$$(mktemp "$@.XXXXXX"); \
	printf '%s\n' '$(F90)' '$(F90FLAGS)' '$(LDFLAGS)' '$(DEFS)' '$(MODINC)' '$(EXTRAINCS)' '$(EXTRAOBJECTS)' '$(EXTRALIBS)' '$(LIBS)' > "$$config_file"; \
	if ! cmp -s "$$config_file" $@; then cp "$$config_file" $@; fi; \
	rm -f "$$config_file"

FORCE:

blmod.o $(OBJECTS): .build-config


blmod.o: blmod.F90
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) -c blmod.F90

../snec: blmod.o $(OBJECTS) $(EXTRAOBJECTS) Makefile ../make.inc
	$(F90) $(LDFLAGS) -o $@ blmod.o $(OBJECTS) $(EXTRAOBJECTS) $(EXTRALIBS) $(LIBS)

$(OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@

# Fortran module dependencies. Keep these explicit so parallel builds do not
# race module-file creation.
outinfo.o: blmod.o csm_switch.o
composition_profile_helpers.o: blmod.o
composition_profile_store.o: mixing_utils.o
csm_switch.o: blmod.o filesystem_utils.o
legacy_grid_pattern.o: blmod.o
cadence_control.o: blmod.o
staged_state.o: blmod.o
simulation_manager.o: staged_state.o blmod.o startup_log.o filesystem_utils.o observer_projection.o
problem.o: staged_state.o blmod.o csm_switch.o startup_log.o
allocate_vars.o: blmod.o outinfo.o
eos_switch.o: blmod.o
eos_content.o: blmod.o
grid.o: blmod.o legacy_grid_pattern.o
read_profile.o: blmod.o csm_switch.o legacy_grid_pattern.o
read_profile_compositions_custom.o: blmod.o mixing_utils.o startup_log.o composition_profile_helpers.o composition_profile_store.o
output.o: blmod.o csm_switch.o outinfo.o
timestep.o: blmod.o
blstep.o: blmod.o
artificial_viscosity.o: blmod.o csm_switch.o
opacity.o: blmod.o startup_log.o
bicubic_interpolator.o: blmod.o
hydro.o: blmod.o
hydro_rad.o: blmod.o
input_parser.o: blmod.o startup_log.o filesystem_utils.o
conservation.o: blmod.o
analysis.o: blmod.o
read_BolCorr.o: blmod.o
opal_opacity.o: blmod.o
nickel.o: blmod.o
simple_saha.o: blmod.o
saha_tables.o: blmod.o
luminosity.o: blmod.o
bomb_profile.o: blmod.o
arrays.o: blmod.o
boxcar.o: blmod.o mixing_utils.o
shock_capture.o: blmod.o csm_switch.o
snec.o: blmod.o outinfo.o startup_log.o cadence_control.o csm_switch.o simulation_manager.o


clean:
	$(CLEANSTUFF)
