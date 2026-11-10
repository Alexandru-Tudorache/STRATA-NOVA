.PHONY: all clean portable debug check

all:
	$(MAKE) -C src

# Production build: same numerical flags as the campaign-ready baseline.
# -fprotect-parens is deliberately retained and -ffast-math is deliberately
# avoided because reassociation/fast-math changes floating-point semantics.

portable:
	$(MAKE) clean
	$(MAKE) F90FLAGS="-O3 -fprotect-parens -flto=auto -fopenmp"

debug:
	$(MAKE) clean
	$(MAKE) F90FLAGS="-O0 -g -fcheck=all -fbacktrace -fopenmp"

check:
	bash -n run_campaign.sh
	$(MAKE) -C src

clean:
	$(MAKE) -C src clean
	rm -f snec
