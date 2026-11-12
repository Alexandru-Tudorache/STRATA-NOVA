# STRATA-NOVA

STRATA-NOVA is a one-dimensional Lagrangian grey radiation-hydrodynamics code
for Type Ia supernova ejecta interacting with circumstellar material (CSM).
It is derived from [SuperSNEC](https://github.com/cfremling/SuperSNEC) and
[SNEC 1.01](https://stellarcollapse.org/index.php/SNEC.html).

The code adds Type Ia profile handling, attached `r^-2` and `r^-3` envelopes,
finite constant-density shells, a Badenes inner free-wind prescription, and a
native staged manager for sequential CSM encounters. It retains selected
SuperSNEC radioactive-deposition and nickel-mixing optimizations. The production
CSM workflow uses a fixed mass grid after initialization. Photometry and
observer-time products use the code's grey approximations.

## Quick start

Requirements: GNU Fortran with OpenMP, GNU Make, LAPACK/BLAS, Bash 4+ and standard
Unix utilities. Python 3.10+, NumPy and Matplotlib are optional for plotting.
Run these commands from the repository root:

The repository includes seven paired HESMA-derived DDT inputs: **N1, N3, N5,
N20, N40, N100 and N200**, +W7 hydrodynamic/composition pair.

N100 ARTIS UVOIR and early STELLA UBVR reference curves and a 35-measurement
SN 2017cbv BVR observation extract are included under `Analysis/DATA/`.

Source authors, archive acknowledgements and preparation conventions are
documented in [THIRD_PARTY.md](THIRD_PARTY.md) and
[the comparison-data README](Analysis/DATA/README.md).

```sh
make -j2
./run_campaign.sh --base example/parameters.quickstart --current --name smoke
```

This runs the N20 delayed-detonation profile without CSM for one day on 400
zones. Results are saved under `campaign_runs/<timestamp>/smoke/`. It is a short
installation check. The root `parameters` file is a separate research example:
N100 ejecta with one staged non-constant wind encounter at a physical shell radius, ending at day 16.
The non-constant rate is fixed by the Badenes+2007 time-history model.

```sh
make portable   # build without -march=native
make debug      # bounds/runtime checks and backtrace
make check      # runner syntax and current build
make clean
```

The production build retains `-fprotect-parens` and avoids `-ffast-math`.
Compiler, architecture, flags and thread counts can affect numerical results;
record them alongside the inputs for a reproducible campaign.
The build detects changes to compiler selection and flags. The campaign runner checks
the build before launching, so source edits cannot silently use an old executable.
It uses the default production flags unless you set `F90FLAGS` in its environment;
for example, `F90FLAGS='-O0 -g -fcheck=all -fbacktrace -fopenmp' ./run_campaign.sh ...`
keeps runtime checks enabled throughout a campaign.

To run the solver directly, create the output folder first:

```sh
mkdir -p quickstart_out
SNEC_PARAMETER_FILE=example/parameters.quickstart \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./snec
```

`./snec` reads `parameters` unless `SNEC_PARAMETER_FILE` names another file.
The solver resets its selected output directory when a run starts. The campaign
runner manages fresh case directories for you.

## Choose your model and CSM

Set both `profile_name` and `comp_profile_name` when changing the ejecta model.
Paths are relative to the repository root, or absolute. Use a matching
hydrodynamic/composition pair in the format accepted by the solver.
For bundled DDT models, `--model` selects both files for you:

```sh
./run_campaign.sh --list-models
./run_campaign.sh --model n20 --base example/parameters.quickstart --name n20_smoke
```

For W7, select its matching files explicitly:

```sh
./run_campaign.sh --base example/parameters.quickstart --name w7_control \
  --set profile_name=profiles/moriya_cases/wdprofile_moriya_w7_homologous_20ks.dat \
  --set comp_profile_name=profiles/moriya_cases/wdcomp_moriya_w7_sourceNi_20ks.dat
```

Model names are case-insensitive. Unspecified settings still come from the base
file. Explicit `--set` profile paths and case-file overrides take precedence over
`--model`; custom models continue to use the two profile parameters.

| Control | Meaning |
| --- | --- |
| `simulation_mode` | `single` or `staged_pipeline` |
| `csm` | `none` or `active`; staged runs require `active` |
| `mass_scheme` | `continuous_mass`, `constant_mass_loss`, `non_constant_mass_loss` or `constant_density` |
| `csm_mass_msun` | Direct attached mass for `continuous_mass` and `constant_density` |
| `csm_outer_radius_cm` | Outer radius in single mode |
| `csm_wind_velocity_kms` | Wind/uniform-shell velocity; continuous-mass envelopes inherit ejecta-surface velocity |
| `csm_temperature_k` | Positive fixed temperature; staged runs also accept a non-positive value for equilibrium temperature at attachment |
| `csm_mdot_msun_per_yr` | Steady mass-loss rate, used only by `constant_mass_loss` |

Generated CSM profiles retain their radius–mass geometry on the hydro grid.
The `continuous_mass` scheme uses the Piro-type $r^{-3}$ mass profile, while
`constant_mass_loss` uses a steady $r^{-2}$ wind.

For staged runs, set `staged_n_shells`, `staged_shell_inner_radii_cm`,
`staged_shell_widths_cm`, and `staged_final_day`. Direct-mass schemes also need
`staged_shell_masses_msun`. Each active array contains one value per shell.
The code evolves the actual hydrodynamic exterior to the next configured
physical radius; encounter times are computed by the simulation, not provided
by the user. Each newly attached shell has outer radius equal to the live
attachment radius plus its specified width. For multiple shells, the next
inner radius must lie outside the preceding shell's outer edge.
The top-level `csm_outer_radius_cm` does not set staged geometry.

The non-constant wind uses a fixed Badenes+2007-inspired rise, plateau, and decline
before the supernova. Its local density comes from the delayed launch time.
The `csm_mdot_msun_per_yr` control does not rescale this wind. The guide is at
[the guide](guides+papers/STRATA_NOVA_GUIDE.pdf), with editable
[LaTeX source](guides+papers/STRATA_NOVA_GUIDE.tex).

## Campaign runner

Run `./run_campaign.sh` and enter one `parameter=value` change per prompt. Enter
`?` to view the starting parameter file, then finish with a blank line. Strings
may be unquoted in the runner, and staged arrays keep ordinary comma syntax.
You can change model files, CSM controls and numerical/output settings through
the same input method.

For example, a single constant-density shell can be selected with:

```sh
./run_campaign.sh --name uniform_shell \
  --set simulation_mode=single --set csm=active \
  --set mass_scheme=constant_density --set csm_mass_msun=0.01 \
  --set csm_outer_radius_cm=1.0d15 --set csm_temperature_k=1.0d4
```

Unspecified controls come from the base file. Use `--base FILE` to start from
another complete parameter file. Model selection uses the same overrides:

```sh
./run_campaign.sh --name n20_control \
  --set profile_name=profiles/ddt_cases/wdprofile_ddt_2013_n20.dat \
  --set comp_profile_name=profiles/ddt_cases/wdcomp_ddt_2013_n20.dat \
  --set simulation_mode=single --set csm=none --set tend=8.64d5
```

For several configurations, put one case per line in a text file:

```text
control | simulation_mode=single | csm=none | tend=8.64d5
piro_envelope | simulation_mode=single | csm=active | mass_scheme=continuous_mass | csm_mass_msun=0.01 | csm_outer_radius_cm=1.0d15 | csm_temperature_k=1.0d4
```

```sh
./run_campaign.sh --file example/campaign.example --dry-run --save campaign_runs/check
./run_campaign.sh --file my_cases.txt --jobs 4 --save campaign_runs/my_cases
```

Case names must start with a letter, digit or `_` and contain only letters,
digits, `_`, `.` or `-`. The `|` character separates overrides and cannot be used
inside values. Blank lines and lines starting with `#` are ignored. Repeated
`--set` options apply to all cases; case-file values take precedence.

Each case retains `parameters.used`, `run.log`, `exit_status`, and `output/`.
Staged stitched results are in `output/final/`. The runner writes directly into
each case's output folder. `--save` controls its location; output-directory
parameter overrides are reserved for the runner. Existing case folders are
reused with a "Current folder is being overwritten" notice. A new calculation
replaces the case output, parameters, log and status; a dry run preserves the
previous solver output. Duplicate case names within one campaign are rejected. `campaign_summary.tsv` contains only cases
from the current invocation.

`--dry-run` prepares inputs without building or running and records `dry-run`
as the status. It checks assignment syntax and template keys; the solver checks
values, arrays and model compatibility when executed. Override keys must appear
in the root template or your base file, including commented assignments. The
runner accepts optional solver keys when you include them in a custom base.

`--jobs N` runs independent cases concurrently and keeps solver streams in
`run.log`; `--quiet` does the same for one case. OpenMP, OpenBLAS and MKL default
to one thread per case, with explicit environment settings taking precedence.
A solver failure makes the campaign return nonzero and preserves partial output.
Success also requires `run_summary.dat` and `lum_observed.dat` to be nonempty.

## Outputs and plotting

`output_mode=0` provides full diagnostics, `1` minimal output and `2` compact
science output. Available products depend on the mode and cadence.

| Product | Use |
| --- | --- |
| `lum_observed.dat` | All-frequency grey bolometric luminosity; staged observer-time integration when enabled |
| `magnitudes.dat` | Approximate broadband photometry |
| `rho.xg`, `radius.xg`, `vel.xg`, `temp.xg`, `press.xg` | Time-dependent profiles when enabled |
| `generated_csm/wdprofile_generated.dat` | Attached input profile in the relevant run/stage directory |
| `generated_csm/boundary_diagnostics.txt` | Attachment and mass diagnostics |
| `manifest.json` | Staged encounter metadata in the final directory |

The sole public light-curve file, `lum_observed.dat`, contains the selected
all-frequency grey luminosity.
For staged runs, `staged_equal_arrival_time=1` applies the visible-hemisphere
light-travel integration to luminosity and individual passband fluxes.
With observer integration disabled, or in single-run mode, `lum_observed.dat`
is in the emission frame. The radius is the existing photospheric surface;
this remains a grey, spherical estimate rather than a detailed emergent spectrum.
No peak cap, fitted rescaling or extra temporal smoothing is applied. Approximate colours retain
the 5000 K colour-temperature floor and inherited bolometric-correction limits.
The staged observer integration assumes a visible hemispherical surface at the
photospheric radius. Every attachment starts a fresh solver process; the live
state is passed onward, while timestep/iteration history restarts.

```sh
python3 -m pip install numpy matplotlib
python3 Analysis/plot_staged_smoke.py campaign_runs/my_cases/CASE/output/final \
  --out Analysis/figures/local
```

Replace `CASE` with your case name. The supplied N100 HESMA reference can be
overlaid with `--hesma Analysis/DATA/HESMA_DDT_LCS/ddt_2013_n100_lightcurve.dat`.
It is a band-limited UVOIR reference; `lum_observed.dat` is the selected
all-frequency grey estimate. The early HESMA UBVR curves are separate STELLA
products. `Analysis/DATA/README.md` records the time/magnitude conventions and
the provenance of the SN 2017cbv observations.

The existing guide figures are retained under `Analysis/figures/`. The guide describes
the calculations and inputs used for those panels. The currently
shipped campaign examples are `example/campaign.example` and
`example/parameters.guide`. New results should record the code version and the
complete `parameters.used` file.

## Attribution and distribution

STRATA-NOVA additions are maintained by Alexandru Tudorache. Cite the software
version, Fremling & Hinds (2026) for SuperSNEC, and Morozova et al. (2015) for SNEC.
The [citation metadata](CITATION.cff), [modification notice](NOTICE.md) and
[third-party provenance](THIRD_PARTY.md) give the corresponding details.

The SNEC/SuperSNEC-derived code and original additions retain the inherited
[CC BY-NC-SA 4.0 license](LICENSE). Independently supplied data and components
retain their own conditions. HESMA data require archive acknowledgement and
model-specific citations. For the N1–N1600C suite, cite Seitenzahl et al. (2013)
for hydrodynamic/nucleosynthesis inputs and Sim et al. (2013) for the companion
synthetic light curves and spectra. Their full references and the archive paper
are in `THIRD_PARTY.md` and `CITATION.cff`. Include the acknowledgement:

> This work made use of the Heidelberg Supernova Model Archive (HESMA), https://hesma.h-its.org.

The guide's interaction curves are STRATA-NOVA calculations; separately
identified HESMA points are model references. Credit Seitenzahl et al. (2013)
and the N100 work of Röpke et al. (2012) for those models, Sim et al. (2013) and
Kromer & Sim (2009) for ARTIS products/methods, and Blinnikov et al. (2006) for
STELLA methodology. W7 credits Nomoto, Thielemann & Yokoi (1984) and Thielemann,
Nomoto & Yokoi (1986); Moriya et al. (2023) describes the wind-comparison setup.
For the observations, cite Wee et al. (2018), its CDS/VizieR catalogue DOI and
the VizieR service acknowledgement recorded in `Analysis/DATA/README.md`.

The software licence and the original sources' data terms remain distinct.
Bundling these comparison inputs does not make their original simulations or
measurements STRATA-NOVA work. Attribution and source links accompany them.
