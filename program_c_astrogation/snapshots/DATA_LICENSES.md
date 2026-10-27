# Data licences and acknowledgements for the snapshots in this directory

The snapshot files here contain third-party astronomical data. **The data are
not covered by the repository's GPL-3.0 code licence.** Each source keeps its
own terms, given below. This file takes precedence over the licence string
embedded in the snapshot files (see "Erratum").

| File | Gaia DR3 | SIMBAD identifiers | NASA Exoplanet Archive records |
|---|---|---|---|
| `gaia_dr3_5pc.txt` | yes (50 objects) | yes (`simbad_used: true`) | queried, no planet rows present |
| `gaia_dr3_10pc.txt` | yes (315 objects) | yes (`simbad_used: true`) | yes (16 planet rows) |
| `gaia_dr3_50pc.txt` | yes (48,809 objects) | no (`simbad_used: false`) | queried, no planet rows present |

The `.manifest.json` files record the queries and retrieval times
(Gaia retrieval: 2026-07-26, table `gaiadr3.gaia_source`).

### Outputs derived from these snapshots

The following tracked outputs were computed from these snapshots. They contain
Gaia DR3 source identifiers and derived astrometry (distances, positions,
directions), and some contain SIMBAD names or NASA Exoplanet Archive planet
data. They carry the same terms and acknowledgements as the snapshots:

- `../results/verification/program_c0/c0_demonstration_gaia_dr3_10pc.txt`;
- `../results/verification/program_c0_real/` (all files);
- `../results/verification/program_c1/`: `data/*.csv` (for example
  `local_neighborhood.csv` and `apparent_sky_*.csv`, 315 objects with SIMBAD
  names; `interception.csv`, `trajectory_comparison.csv` and
  `worldline_proxima_1g.csv` for Proxima Centauri), the `figures/` (PNG, SVG
  and the `*.plotted.csv` data behind them), `animations/`, `contact_sheets/`
  and `verification/`, all of which were made from the 10 pc snapshot.

Purely synthetic outputs (`data_kind` `SYNTHETIC`, for example
`../results/verification/program_c0/c0_snapshot_50pc.txt`) contain no Gaia data.

## Gaia DR3 (ESA)

- **Licence: CC BY-NC 3.0 IGO.** Quoted from the ESA Gaia licence page,
  https://www.cosmos.esa.int/web/gaia-users/license (retrieved 2026-10-27):
  "Gaia data are distributed under the CC BY-NC 3.0 IGO license. For details and
  guidelines concerning commercial use of the Gaia data, please see the Terms
  and Conditions for the use of data in the ESA space science archives."
- **Terms and credit line.** Quoted from "Terms and Conditions for the use of
  data in the ESA Space Science Archives",
  https://www.cosmos.esa.int/web/esdc/terms-and-conditions (retrieved
  2026-10-27): data "are open and free to use by the User subject to proper
  acknowledgment to be given in accordance with the credit lines provided
  below"; "Prior to any commercial use … a detailed request for
  authorisation/licence shall be made by the User by sending email to
  data.licences@esa.int". The credit line given there for Gaia is:

  > Credit: ESA, Gaia DPAC

- **Attribution for these files** (composed to meet the credit-line, licence-link
  and indication-of-changes requirements; not ESA-supplied wording):

  > Contains data from the European Space Agency (ESA) mission Gaia, Data
  > Release 3, processed by the Gaia Data Processing and Analysis Consortium
  > (DPAC). Credit: ESA, Gaia DPAC. Licence: CC BY-NC 3.0 IGO
  > (https://creativecommons.org/licenses/by-nc/3.0/igo/). Changes: a subset of
  > `gaiadr3.gaia_source` columns was selected and reformatted by this
  > repository's import tool (`../tools/import_gaia.py`), and derived quantities
  > (for example inverse-parallax distances and Cartesian positions) were added.

- **Citation.** The repository cites Gaia Collaboration (2016) and Gaia
  Collaboration (2023). ESA publishes DR3-specific acknowledgement and citation
  instructions at
  https://gea.esac.esa.int/archive/documentation/GDR3/Miscellaneous/sec_acknowl/
  and
  https://gea.esac.esa.int/archive/documentation/GDR3/Miscellaneous/sec_credit_and_citation_instructions/
  (linked from https://www.cosmos.esa.int/web/gaia-users/credits). **Not yet
  checked:** both pages returned HTTP 503 (maintenance) on 2026-10-27, so their
  wording is not reproduced here. Copy it from those pages before any
  publication that uses these data.
- **Consequence.** The non-commercial term applies to these files and to
  adapted versions of them. The GPL-3.0 code licence does not cover them and
  cannot relicense them.

## SIMBAD (CDS, Strasbourg)

Object names and classifications in the 5 pc and 10 pc snapshots were obtained
from SIMBAD. Quoted from the SIMBAD user guide, https://simbad.cds.unistra.fr/guide/
(retrieved 2026-10-27): "If the Simbad database was helpful for your research
work the following acknowledgment would be appreciated":

> This research has made use of the SIMBAD database, operated at CDS,
> Strasbourg, France.

The same page cites: Wenger et al. 2000, "The SIMBAD astronomical database",
A&AS, 143, 9.

The CDS general licence page (https://cds.unistra.fr/) could not be reached on
2026-10-27 (connection timeout), so no further CDS terms are recorded here.

## NASA Exoplanet Archive

Planet records in the 10 pc snapshot were obtained from the NASA Exoplanet
Archive. Quoted from https://exoplanetarchive.ipac.caltech.edu/docs/acknowledge.html
(retrieved 2026-10-27; page "Last updated: 13 August 2025"): "Please include
the following standard acknowledgment in any published material that makes use
of the NASA Exoplanet Archive's services":

> This research has made use of the NASA Exoplanet Archive, which is operated by
> the California Institute of Technology, under contract with the National
> Aeronautics and Space Administration under the Exoplanet Exploration Program.

The page also asks users to cite Christiansen et al. (2025) for the archive and,
"if you use data from a specific literature reference, please acknowledge that
reference directly". The archive service DOI listed there is 10.26133/NEA45.

## Erratum: "CC BY-SA 3.0 IGO" in existing files

Files generated before 2026-10-27 state the Gaia licence as
`ESA/Gaia/DPAC, CC BY-SA 3.0 IGO`. **That string is wrong; the licence is
CC BY-NC 3.0 IGO, as above.** It appears in:

- `gaia_dr3_5pc.txt`, `gaia_dr3_10pc.txt`, `gaia_dr3_50pc.txt` (header
  `licence` line and every provenance row) and their `.manifest.json` files;
- `../results/verification/program_c0/c0_demonstration_gaia_dr3_10pc.txt`;
- `../results/verification/program_c0_real/c0_real_metadata.json`;
- the synthetic records `../results/verification/program_c0/c0_manifest.txt`
  and `c0_snapshot_50pc.txt` ("would be … NOT APPLICABLE here").

These files were **not** edited. Each snapshot carries an FNV-1a checksum over
its own bytes, which the C++ reader verifies, and the 5 pc and 10 pc checksums
are recorded in the program C0-real and C1 verification records and figure
metadata. Editing the string would invalidate that chain of evidence. The
generator (`../tools/import_gaia.py`) and the synthetic fixture
(`../verification/c0_fixtures.hpp`) now use the correct licence, so snapshots
and records regenerated from them will carry it.
