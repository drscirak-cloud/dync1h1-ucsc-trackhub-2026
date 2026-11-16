# Human *DYNC1H1* Dual-Track UCSC Genome Browser Track Hub — Release `v2.0.0` (`GRCh38` / `hg38`)

**Direct Live UCSC Genome Browser Launch Link (`GRCh38 chr14:101,964,573-102,056,443`):**
https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubUrl=https://raw.githubusercontent.com/drscirak-cloud/dync1h1-ucsc-trackhub-2026/main/hub.txt&position=chr14:101964573-102056443

**Track Hub URL (`hubUrl`):**
`https://raw.githubusercontent.com/drscirak-cloud/dync1h1-ucsc-trackhub-2026/main/hub.txt`

## Scientific Revision `v2.0.0` Summary
- **Canonical Reference:** *DYNC1H1* (`ENSG00000197102.14`, `GRCh38 chr14:101,964,573-102,056,443`, `+1` strand), MANE Select `NM_001376.5` / `ENST00000360184.10` (`78` exons, `13,941`-nt CDS, final exon 77/78 junction `c.13812`, `55-nt` NMD boundary `c.13757`), heavy chain `NP_001367.2` / `Q14204` (`4,646` AA). Zero unpublished internal cases and zero cross-gene literature contaminants.
- **Included Binary UCSC Tracks (`GRCh38/`):**
  1. **`DYNC1H1 78 Exons` (`bigBed 9`, RGB):** All `78` exons of `NM_001376.5` (`35` CDS exons divisible by 3 = **`34` internal frame-preserving `3n` exons** [Teal] + **`1` terminal `Exon 78`** escaping `55-nt` NMD downstream of `c.13757` [Red], and **`43` phase-shifting exons** [Slate]).
  2. **`DYNC1H1 Quadrants v2` (`bigBed 9`, RGB):** All `11,892` observed public GRCh38 *DYNC1H1* variants classified by out-of-fold Dual-Track mechanistic class (`Q I` Concordant High, `Q IIa` SpliceAI/PTV/NMD-Driven, `Q IIb` Constrained Coding Missense / Regulatory with Low Pure-3D Score, `Q III` 3D Pocket/Interface High, `Q IV` Concordant Low).
  3. **`Track 1A Pure 3D` (`bigWig`):** Pure 3D homodimer biophysical sensitivity across the polymer-relaxed, clash-free `2 x 4,646 = 9,292`-AA Off-Phi and On-DDA homodimers (`0` chain breaks, `0` non-bonded CA overlaps `< 3.40 A`; **zero evolutionary conservation `S3` and zero supervised `SPDV`**).
  4. **`GNM Hinge v2.0` (`bigWig`):** Validated clash-free `9,292`-node Kirchhoff Gaussian Network Model slow-mode mechanical hinge strain.
  5. **`Track 2 LOCO O/E` (`bigWig`):** Leave-One-Codon-Out (`LOCO`) `31`-codon CpG-calibrated regional missense depletion (`1 - O/E`), strictly excluding the focal codon.
