# *DYNC1H1* Blinded Dual-Track Interactive Web Atlas — Release `v2.0.0`

Eigenstaendige, statisch im Browser lauffaehige `3Dmol.js`-Web-Applikation fuer das *DYNC1H1*-Dual-Track-Projekt (`GRCh38 chr14:101,964,573-102,056,443`, MANE Select `NM_001376.5`, `NP_001367.2` / `Q14204`, `4.646` Aminosaeuren, `78` kodierende Exons: `34` interne rahmenerhaltende `3n`-Exons + `1` terminales `Exon 78` + `43` phasenverschiebende Exons).

## Funktionen (`Release v2.0.0`)

1. **Bruch- und kollisionsfreier `2 x 4.646`-AA-Zweizustands-Homodimer-Viewer (`3Dmol.js`):** Umschaltbar zwischen Zustand A (`autoinhibierter Phi-Partikel`, `9.292` AA) und Zustand B (`aktiver DDA-LIS1-MTBD-Komplex`, `9.292` AA; `0` Kettenbrueche, `0` nicht-gebundene Ueberlappungen `< 3,40 A`), farbkodierbar nach `17` Strukturdomaenen, `Track 1A Pure 3D`-Sensitivitaet (`0` Konservierung), `9.292`-Knoten-Kirchhoff-`GNM`-Scharnierdehnung und `Track 2` Leave-One-Codon-Out (`LOCO`) CpG-Depletion (`1 - O/E`).
2. **Kodon- & Varianten-Inspektor (`v2.0.0`):** Vergleich von `Track 1A (Pure 3D)`, `Track 1B (3D + Konservierung)`, `Track 2-Intrinsic` (`0` fokale Allelfrequenzen) und `M3B Joint` ueber die Fuenf-Klassen-Mechanismusmatrix (`Q I`, `Q IIa` Splicing/PTV/NMD, `Q IIb` sequenzdepletierte Missense, `Q III` 3D-Taschen/Naehte, `Q IV` benigne).
3. **Empirischer `GTEx v10` Homodimer-Stoechiometrie-Rechner:** Simulation von `WT/WT`-, `WT/Mut`- und `Mut/Mut`-Homodimeranteilen unter `55-nt`-NMD (`c.13757`) und empirischer Spenderdispersions-Skalierung (`SD(log2(TPM + 1))`: Skelettmuskel `7,70` pp, peripherer Nerv `4,28` pp, Fibroblasten `4,09` pp, Grosshirnkortex `0,00` pp).
4. **Direktanbindung an den binären UCSC-Track-Hub (`Release v2.0.0`):** Git-Commit `2803d5ed2ac509e983a9f208603b1649a173fcfa` (`https://github.com/drscirak-cloud/dync1h1-ucsc-trackhub-2026/releases/tag/v2.0.0`).
