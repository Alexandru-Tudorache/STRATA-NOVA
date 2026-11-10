# STRATA-NOVA attribution and modification notice

STRATA-NOVA is a modified derivative of SNEC 1.01 and SuperSNEC.

## SNEC

SNEC was created by Viktoriya Morozova, Christian D. Ott and Anthony L. Piro;
SNEC 1.01 also includes contributions from Brian W. Mulligan. The SNEC manual
identifies Morozova, Ott and Piro as the SNEC copyright holders and distributes
SNEC under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0
International license.

SNEC reference:

V. Morozova, A. L. Piro, M. Renzo, C. D. Ott, D. Clausen, S. M. Couch,
J. Ellis, and L. F. Roberts, *The Astrophysical Journal* 814, 63 (2015),
DOI 10.1088/0004-637X/814/1/63.

## SuperSNEC

SuperSNEC was developed by Christoffer Fremling and K-Ryan Hinds as a derivative
of SNEC 1.01. Its [public repository](https://github.com/cfremling/SuperSNEC) applies
CC BY-NC-SA 4.0 and documents its
modifications relative to SNEC.

SuperSNEC reference: C. Fremling and K-Ryan Hinds (2026), *SuperSNEC: Fast and
Accurate Light Curve Production for Large Hydrodynamic Model Grids Using
Adaptive Gridding*, [arXiv:2603.05680](https://arxiv.org/abs/2603.05680). The
production STRATA-NOVA CSM workflow retains selected upstream optimizations
and uses a fixed mass grid.

## STRATA-NOVA modifications

STRATA-NOVA original additions and modifications are maintained by Alexandru
Tudorache. Relative to SNEC/SuperSNEC, the project includes substantial changes
for Type Ia ejecta/CSM work, including Type Ia profile ingestion, attached CSM
constructors, Badenes free-wind support, staged sequential encounters, fixed
native mass-grid operation for the production CSM workflow, observer/output
approximations, CSM-aware shock handling, and conservation/runtime diagnostics.

This notice records upstream attribution and STRATA-NOVA modifications. The
full inherited license text is retained in `LICENSE`; component and data
provenance is recorded in `THIRD_PARTY.md`.

## Third-party exceptions

The SNEC/STRATA-NOVA root licence covers the derived software, while certain
inherited components and independent scientific datasets retain their own
attribution and conditions. This includes the Timmes EOS/Saha material, OPAL
opacities, HESMA model products and published observational photometry.
Their respective sources are recorded in `THIRD_PARTY.md`.

## Model-data acknowledgement

This work made use of the Heidelberg Supernova Model Archive (HESMA), https://hesma.h-its.org.

The N1–N200 delayed-detonation hydrodynamic inputs are credited to Seitenzahl
et al. (2013), MNRAS 429, 1156; the companion synthetic radiation-transfer data
are credited to Sim et al. (2013), MNRAS 436, 333. Early UBVR STELLA
reference points also credit Blinnikov et al. (2006), A&A 453, 229, the
archive's STELLA methods reference. See `THIRD_PARTY.md` for the
complete references, archive citation and distinct data conditions.

## Included models and comparison data

The publication set includes the N1, N3, N5, N20, N40, N100 and N200 DDT
profile/composition pairs, the prepared W7 pair, two N100 HESMA light-curve
references and the 35-measurement SN 2017cbv BVR extract. The project README and
`THIRD_PARTY.md` identify the authors, data products and preparation conventions.

W7 is credited to Nomoto, Thielemann & Yokoi (1984), ApJ 286, 644, and
Thielemann, Nomoto & Yokoi (1986), A&A 158, 17. Moriya et al. (2023),
MNRAS 522, 6035, describes the associated wind-comparison setup.

The HESMA synthetic references credit the N100/DDT work of Seitenzahl et al.
(2013) and Röpke et al. (2012), the ARTIS synthetic light curves of Sim et al.
(2013), ARTIS methodology from Kromer & Sim (2009), and STELLA methodology
from Blinnikov et al. (2006).

SN 2017cbv measurements are credited to Wee et al. (2018), ApJ 863, 90,
Table 5, DOI 10.3847/1538-4357/aacd4e, and CDS/VizieR J/ApJ/863/90,
DOI 10.26093/cds/vizier.18630090. The extract retains the measured magnitudes
and reported errors and expresses observation dates as MJD.

This research has made use of the VizieR catalogue access tool, CDS, Strasbourg, France (DOI: 10.26093/cds/vizier).

These datasets retain their original scientific attribution and source terms.
They are separate from the STRATA-NOVA-generated interaction curves, and the
root software licence does not relicense them. Full references are given in
`THIRD_PARTY.md` and `CITATION.cff`.

## Endorsement note

No endorsement by the SNEC, SuperSNEC, HESMA, W7, CDS/VizieR, Timmes EOS,
OPAL or other source contributors is implied.
