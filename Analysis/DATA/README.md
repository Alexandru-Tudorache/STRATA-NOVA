# Observational and model-reference data

The repository includes the following comparison tables, with their original
scientific attribution. They are distinct from STRATA-NOVA-generated outputs;
the root software licence does not relicense independent datasets.

## N100 synthetic light curves from HESMA

| File | Content |
| --- | --- |
| `HESMA_DDT_LCS/ddt_2013_n100_lightcurve.dat` | ARTIS UVOIR synthetic luminosity; time in days, luminosity in erg/s |
| `HESMA_DDT_LCS/ddt_2013_n100_lightcurves_early.dat` | Early STELLA U, B, V and R synthetic magnitudes; time in days |

The native numerical reference files are retained unchanged. The ARTIS curve
is band-limited UVOIR and should be labelled that way when compared with the
code's all-frequency grey luminosity. The early STELLA curves are a separate
calculation. B−V is obtained by subtracting simultaneous B and V entries.

Credit the N100/DDT work of Seitenzahl et al. (2013),
[doi:10.1093/mnras/sts402](https://doi.org/10.1093/mnras/sts402), and Röpke et al.
(2012), [doi:10.1088/2041-8205/750/1/L19](https://doi.org/10.1088/2041-8205/750/1/L19);
the synthetic ARTIS light curves of Sim et al. (2013),
[doi:10.1093/mnras/stt1574](https://doi.org/10.1093/mnras/stt1574);
ARTIS methodology from Kromer & Sim (2009),
[doi:10.1111/j.1365-2966.2009.15256.x](https://doi.org/10.1111/j.1365-2966.2009.15256.x);
and STELLA methodology from Blinnikov et al. (2006),
[doi:10.1051/0004-6361:20054594](https://doi.org/10.1051/0004-6361:20054594).
The archive reference is Kromer, Ohlmann & Röpke (2017),
[arXiv:1706.09879](https://arxiv.org/abs/1706.09879).

This work made use of the Heidelberg Supernova Model Archive (HESMA), https://hesma.h-its.org.

## SN 2017cbv observed BVR photometry

`OBSERVATIONS/SN2017cbv_Wee2018_BVR_Table5.csv` contains 35 measurements:
11 B, 12 V and 12 R entries. Its columns are:

| Column | Meaning |
| --- | --- |
| `time` | Observation date, MJD in days |
| `band` | B, V or R filter label |
| `mag` | Apparent magnitude |
| `mag_err` | Published magnitude measurement uncertainty |

The extract comes from Wee et al. (2018), Table 5. Published JD−2457000 dates
were converted to MJD by adding 56999.5 days. No distance, first-light or
reddening corrections are applied in the CSV; those belong in the plotting
configuration. The acknowledgement is recorded here so the CSV remains a
four-column numerical table. All 35 measurement fields are preserved.

- Wee, J., Chakraborty, N., Wang, J. & Penprase, B. E. (2018), *Optical and
  Infrared Photometry of the Nearby SN 2017cbv*, ApJ 863, 90,
  [doi:10.3847/1538-4357/aacd4e](https://doi.org/10.3847/1538-4357/aacd4e).
- CDS/VizieR catalogue **J/ApJ/863/90**,
  [doi:10.26093/cds/vizier.18630090](https://doi.org/10.26093/cds/vizier.18630090).
- Ochsenbein, F., Bauer, P. & Marcout, J. (2000), *The VizieR database of
  astronomical catalogues*, A&AS 143, 23–32,
  [doi:10.1051/aas:2000169](https://doi.org/10.1051/aas:2000169).

This research has made use of the VizieR catalogue access tool, CDS, Strasbourg, France (DOI: 10.26093/cds/vizier).

Original data-source terms remain associated with both the observations and
model references. Full software and model credits are in
[THIRD_PARTY.md](../../THIRD_PARTY.md) and [CITATION.cff](../../CITATION.cff).
