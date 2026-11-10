# Third-party attribution and provenance

STRATA-NOVA is derived from SNEC 1.01 and SuperSNEC. The root `LICENSE`
retains CC BY-NC-SA 4.0 for the derived software and original additions.
Independent scientific datasets and components retain their original source
terms and attribution; inclusion does not relicense them as STRATA-NOVA work.

## Software and inherited resources

| Material | Attribution and distribution basis |
| --- | --- |
| SNEC 1.01 | Viktoriya Morozova, Christian D. Ott and Anthony L. Piro; version 1.01 also credits Brian W. Mulligan. [SNEC webpage](https://stellarcollapse.org/index.php/SNEC.html), Morozova et al. (2015), and the full inherited CC BY-NC-SA 4.0 text in `LICENSE`. |
| SuperSNEC | Christoffer Fremling and K-Ryan Hinds. [Public repository](https://github.com/cfremling/SuperSNEC), CC BY-NC-SA 4.0, and Fremling & Hinds (2026). Its requested SNEC and SuperSNEC citations are retained below and in `CITATION.cff`. |
| Timmes EOS/Saha material, including `src/saha_tables.F90` | Frank Timmes and collaborators. The [SNEC manual](https://stellarcollapse.org/codes/snec_notes-1.00.pdf) records distribution with Frank Timmes's permission and separate component terms. Original source notices remain. Cite Timmes & Arnett (1999) and Timmes & Swesty (2000) as relevant to the inherited material. |
| OPAL interpolation and `tables/codataa`–`codatae` | LLNL OPAL project; Iglesias & Rogers (1996). The SNEC manual records distribution with Carlos Iglesias's permission. Original interpolation/source notices remain. |
| `tables/BolCorr.dat` | Inherited SNEC/SuperSNEC bolometric-correction resource, not an original STRATA-NOVA dataset. |
| BLAS/LAPACK, GNU Fortran/OpenMP, Bash, Python, NumPy and Matplotlib | Installed separately; their binaries are not bundled in the source repository. |

The upstream manual also credits Ferguson et al. (2005) for separate
low-temperature opacity data. Those separate files are not present here; the
shipped `codata*` files identify OPAL tables. The Timmes references credit
inherited EOS material without claiming that every EOS option implements a
complete Helmholtz method.

## Included model and comparison data

| Files | Original work and preparation |
| --- | --- |
| `profiles/ddt_cases/` | Seven paired HESMA-derived DDT inputs: N1, N3, N5, N20, N40, N100 and N200, prepared in SNEC format from the radially averaged ejecta/composition products. Seitenzahl et al. (2013); Röpke et al. (2012) for the earlier N100 presentation; HESMA archive credit. |
| `profiles/moriya_cases/` | Prepared W7 profile/composition pair at 20,000 s, approximately 4 × 10^13 cm exterior radius. Original W7: Nomoto, Thielemann & Yokoi (1984) and Thielemann, Nomoto & Yokoi (1986). Wind-comparison setup: Moriya et al. (2023). |
| `Analysis/DATA/HESMA_DDT_LCS/` | Two unchanged native N100 reference tables: ARTIS UVOIR luminosity and early STELLA UBVR magnitudes. N100/DDT: Seitenzahl et al. (2013), Röpke et al. (2012); synthetic ARTIS products: Sim et al. (2013); ARTIS method: Kromer & Sim (2009); STELLA method: Blinnikov et al. (2006). |
| `Analysis/DATA/OBSERVATIONS/SN2017cbv_Wee2018_BVR_Table5.csv` | 35 BVR measurements with published errors from Wee et al. (2018), Table 5, CDS/VizieR J/ApJ/863/90. Dates converted from JD−2457000 to MJD by adding 56999.5. All measurement values are retained; acknowledgement text is in the data README rather than a CSV data row. |
| `Analysis/figures/` and the STRATA-NOVA guide | Project figures/documentation, with the underlying HESMA/model/observation credits retained. Continuous interaction curves are STRATA-NOVA calculations; separately identified reference points are the original authors' products. |

The distribution contains the model and comparison files listed above.
Source credits accompany the prepared inputs and reference data.

## Archive and catalogue acknowledgements

This work made use of the Heidelberg Supernova Model Archive (HESMA), https://hesma.h-its.org.

This research has made use of the VizieR catalogue access tool, CDS, Strasbourg, France (DOI: 10.26093/cds/vizier).

HESMA requests original-model citations and its archive acknowledgement on the
[current homepage](https://hesma.h-its.org/) and
[legacy archive](https://old-hesma.h-its.org/doku.php?id=start).
[CDS/VizieR](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJ/863/90) supplies
the observational catalogue reference and requests the catalogue/service credit.
Data-specific terms remain those supplied by the original archives and authors;
the software licence is not a separate dataset licence. These source and
acknowledgement records do not assert an additional permission grant or a
restriction beyond the original terms. The HESMA website's
[legal notice](https://hesma.h-its.org/legalnotice) and
[CDS terms](https://cds.unistra.fr/legals/) are linked for their original wording.

## Framework references

- Morozova, V., Piro, A. L., Renzo, M., Ott, C. D., Clausen, D., Couch, S. M.,
  Ellis, J. & Roberts, L. F. (2015), *Light Curves of Core-Collapse Supernovae
  with Substantial Mass Loss Using the New Open-Source SuperNova Explosion
  Code (SNEC)*, ApJ 814, 63,
  [doi:10.1088/0004-637X/814/1/63](https://doi.org/10.1088/0004-637X/814/1/63).
  Also cite the [SNEC code webpage](https://stellarcollapse.org/index.php/SNEC.html)
  as requested by its manual.
- Fremling, C. & Hinds, K-Ryan (2026), *SuperSNEC: Fast and Accurate Light
  Curve Production for Large Hydrodynamic Model Grids Using Adaptive Gridding*,
  [arXiv:2603.05680](https://arxiv.org/abs/2603.05680).

## DDT, HESMA and radiative-transfer references

- Seitenzahl, I. R., Ciaraldi-Schoolmann, F., Röpke, F. K., Fink, M.,
  Hillebrandt, W., Kromer, M., Pakmor, R., Ruiter, A. J., Sim, S. A. &
  Taubenberger, S. (2013), *Three-dimensional delayed-detonation models with
  nucleosynthesis for Type Ia supernovae*, MNRAS 429, 1156–1172,
  [doi:10.1093/mnras/sts402](https://doi.org/10.1093/mnras/sts402).
- Röpke, F. K., et al. (2012), *Constraining Type Ia supernova models:
  SN 2011fe as a test case*, ApJL 750, L19,
  [doi:10.1088/2041-8205/750/1/L19](https://doi.org/10.1088/2041-8205/750/1/L19).
- Kromer, M., Ohlmann, S. T. & Röpke, F. K. (2017), *Simulating the observed
  diversity of Type Ia supernovae – Introducing a model data base*,
  [arXiv:1706.09879](https://arxiv.org/abs/1706.09879).
- Sim, S. A., Seitenzahl, I. R., Kromer, M., Ciaraldi-Schoolmann, F., Röpke,
  F. K., Fink, M., Hillebrandt, W., Pakmor, R., Ruiter, A. J. & Taubenberger,
  S. (2013), *Synthetic light curves and spectra for three-dimensional
  delayed-detonation models of Type Ia supernovae*, MNRAS 436, 333–347,
  [doi:10.1093/mnras/stt1574](https://doi.org/10.1093/mnras/stt1574).
- Kromer, M. & Sim, S. A. (2009), *Time-dependent three-dimensional spectrum
  synthesis for Type Ia supernovae*, MNRAS 398, 1809–1826,
  [doi:10.1111/j.1365-2966.2009.15256.x](https://doi.org/10.1111/j.1365-2966.2009.15256.x).
- Blinnikov, S. I., Röpke, F. K., Sorokina, E. I., Gieseler, M., Reinecke, M.,
  Travaglio, C., Hillebrandt, W. & Stritzinger, M. (2006), *Theoretical light
  curves for deflagration models of type Ia supernova*, A&A 453, 229–240,
  [doi:10.1051/0004-6361:20054594](https://doi.org/10.1051/0004-6361:20054594).

## W7 and wind-comparison references

- Nomoto, K., Thielemann, F.-K. & Yokoi, K. (1984), *Accreting white dwarf
  models of Type I supernovae. III — Carbon deflagration supernovae*, ApJ 286,
  644–658, [doi:10.1086/162639](https://doi.org/10.1086/162639).
- Thielemann, F.-K., Nomoto, K. & Yokoi, K. (1986), *Explosive nucleosynthesis
  in carbon deflagration models of Type I supernovae*, A&A 158, 17–33,
  [ADS record](https://ui.adsabs.harvard.edu/abs/1986A%26A...158...17T/abstract).
- Moriya, T. J., Mazzali, P. A., Ashall, C. & Pian, E. (2023), *Early excess
  emission in Type Ia supernovae from the interaction between supernova ejecta
  and their circumstellar wind*, MNRAS 522, 6035–6042,
  [doi:10.1093/mnras/stad1386](https://doi.org/10.1093/mnras/stad1386).

## Observation and catalogue references

- Wee, J., Chakraborty, N., Wang, J. & Penprase, B. E. (2018), *Optical and
  Infrared Photometry of the Nearby SN 2017cbv*, ApJ 863, 90,
  [doi:10.3847/1538-4357/aacd4e](https://doi.org/10.3847/1538-4357/aacd4e).
- CDS/VizieR **J/ApJ/863/90**,
  [doi:10.26093/cds/vizier.18630090](https://doi.org/10.26093/cds/vizier.18630090).
- Ochsenbein, F., Bauer, P. & Marcout, J. (2000), *The VizieR database of
  astronomical catalogues*, A&AS 143, 23–32,
  [doi:10.1051/aas:2000169](https://doi.org/10.1051/aas:2000169).

The observation CSV contains apparent magnitudes and measurement errors.
Distance/extinction/first-light adjustments are plotting choices and are not
edits to the source measurements. HESMA UVOIR is band-limited and early UBVR
comes from STELLA; these should not be labelled as full-bolometric observations.

## EOS and opacity references

- Timmes, F. X. & Arnett, D. (1999), *The Accuracy, Consistency, and Speed of
  Five Equations of State for Stellar Hydrodynamics*, ApJS 125, 277,
  [doi:10.1086/313271](https://doi.org/10.1086/313271).
- Timmes, F. X. & Swesty, F. D. (2000), *The Accuracy, Consistency, and Speed
  of an Electron–Positron Equation of State Based on Table Interpolation of
  the Helmholtz Free Energy*, ApJS 126, 501,
  [doi:10.1086/313304](https://doi.org/10.1086/313304).
- Iglesias, C. A. & Rogers, F. J. (1996), *Updated OPAL Opacities*, ApJ 464,
  943, [doi:10.1086/177381](https://doi.org/10.1086/177381).

## STRATA-NOVA modifications

STRATA-NOVA retains selected SuperSNEC deposition/nickel-control optimizations,
uses a fixed native mass grid during each production evolution segment, and
adds Type Ia input handling, attached CSM construction, sequential live-state
handoffs, shock handling, grey observable prescriptions and conservation/runtime
diagnostics. `NOTICE.md` records the derivative relationship and modifications.
Original-source credits are retained; no upstream endorsement is implied.
