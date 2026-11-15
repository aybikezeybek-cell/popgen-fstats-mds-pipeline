# Ancient DNA Population Genetics Pipeline: Outgroup f3, f4-Statistics & MDS

**Reproducible Evolutionary Genomics Workflow (`Bash`, `R`, `AdmixTools`, `EIGENSTRAT`, `tidyverse`, `ggrepel`)**

---

### Overview
This repository implements an end-to-end, reproducible population genomics pipeline for quantifying shared genetic drift and testing admixture/gene flow hypotheses across **Ancient Near Eastern Neolithic populations** using the publicly available **Allen Ancient DNA Resource (AADR v66, 1240K panel)**.

The workflow automates subset extraction from large EIGENSTRAT datasets (`convertf`), programmatically generates non-redundant population pairs and quartets for **AdmixTools** (`qp3Pop` and `qpDstat`), converts outgroup $f_3$ statistics into pairwise genetic distance matrices (1 - f3), and produces publication-ready **Multidimensional Scaling (MDS)** and **f4 forest plots** in `R`.

---

### Study Populations & Archaeological Horizons
A total of **244 ancient and modern reference individuals** were extracted from the AADR v66 1240K dataset across three major Near Eastern Neolithic horizons and an African outgroup:

| Population Label | Geographic & Cultural Horizon | Approximate Chronology | Sample Size ($n$) |
| :--- | :--- | :--- | :---: |
| `Turkey_PPN` | Central Anatolia Pre-Pottery Neolithic | ~8500–7000 BCE | 32 |
| `Turkey_Catalhoyuk_EN` | Central Anatolia Early Neolithic (Çatalhöyük) | ~7100–6700 BCE | 17 |
| `Turkey_Catalhoyuk_MN` | Central Anatolia Middle Neolithic (Çatalhöyük) | ~6700–6300 BCE | 46 |
| `Turkey_Catalhoyuk_LN` | Central Anatolia Late Neolithic (Çatalhöyük) | ~6300–5950 BCE | 32 |
| `Turkey_N` | Northwestern Anatolia Pottery Neolithic | ~6500–6000 BCE | 68 |
| `Jordan_PPNB` | Southern Levant Pre-Pottery Neolithic B | ~8500–6700 BCE | 17 |
| `Iran_GanjDareh_N` | Zagros / Western Iran Early Neolithic | ~8000 BCE | 17 |
| `Mbuti` | Central Africa Outgroup (Ituri Forest) | Present-day | 15 |

---

### Theoretical & Analytical Framework

1. **Subset Extraction (`EIGENSTRAT / convertf`):**
   - Filters the target 8 populations from the full 1240K dataset (`v66.p1_1240K.aadr.patch.PUB`) into a lightweight, memory-efficient working subset.
2. **Automated Combination Enumeration (`scripts/01_run_admixtools.sh`):**
   - Programmatically enumerates $\binom{n}{2} = 21$ unique unordered pairs for outgroup $f_3(P_1, P_2; \text{Mbuti})$ and $105$ non-redundant quartets for $f_4(W, X; Y, \text{Mbuti})$ to eliminate symmetric duplicate tests.
3. **Outgroup $f_3$-Statistics (`qp3Pop`):**
   - Measures shared genetic drift between two test populations ($A, B$) relative to an outgroup ($O = \text{Mbuti}$):
     $$f3(A, B; O) = \mathbb{E}\big[(p_O - p_A)(p_O - p_B)\big]$$
   - Higher $f_3$ values indicate greater shared ancestry since divergence from the outgroup (`inbreed: YES` enabled for pseudo-haploid ancient DNA genotypes).
4. **$f4-Statistics (`qpDstat`):**
   - Tests formal treeness and excess allele sharing across population quartets:
     $$f_4(W, X; Y, O) = \mathbb{E}\big[(p_W - p_X)(p_Y - p_O)\big]$$
   - Significant positive values ($Z > 3$) indicate excess allele sharing between $X$ and $Y$ (or $W$ and $O$), whereas $Z < -3$ indicates excess allele sharing between $W$ and $Y$.
5. **Genetic Distance Matrix & Classical MDS (`scripts/02_fstats_mds_analysis.R`):**
   - Transforms pairwise outgroup $f_3$ estimates into a symmetric dissimilarity matrix ($D_{ij} = 1 - f_3(i, j; \text{Mbuti})$) and performs classical Multidimensional Scaling (`cmdscale()`) alongside block-jackknife $\pm 3\text{ SE}$ forest plots in `ggplot2`.

---

### Key Evolutionary Findings

- **Anatolian Neolithic Continuity:** All Anatolian Neolithic groups (`Turkey_PPN`, `Turkey_Catalhoyuk_EN/MN/LN`, and `Turkey_N`) exhibit tightly clustered, high shared genetic drift ($f_3 \approx 0.297 - 0.300$, $Z > 125$) in MDS space, confirming strong regional genetic continuity across the aceramic-to-ceramic transition.
- **Levantine Gene Flow into Pottery Neolithic Anatolia:** The test `f4(Turkey_N, Turkey_PPN; Jordan_PPNB, Mbuti) = +0.000848` ($Z = +6.0$) reveals statistically significant excess allele sharing between Pottery Neolithic Anatolia (`Turkey_N`) and Levantine PPNB (`Jordan_PPNB`) compared to earlier Pre-Pottery Anatolian farmers.
- **Zagros Divergence:** Early Neolithic Zagros (`Iran_GanjDareh_N`) separates distinctly along MDS Dimension 1 ($19.2\%$ variance explained), showing no significant direct gene flow into Anatolian Neolithic populations during the same horizon (`f4(Turkey_N, Turkey_PPN; Iran_GanjDareh_N, Mbuti) = +0.000020`, $Z = +0.2$).

---

### Visualizations

#### Classical MDS of Ancient Near Eastern Populations
<img width="2850" height="1800" alt="mds_anatolia_neolithic png" src="https://github.com/user-attachments/assets/531f9a70-d5aa-43ab-9fc3-a7dae158bc31" />


#### Gene Flow Signals via $f_4$-Statistics (`qpDstat`)
<img width="3150" height="1860" alt="f4_statistics_plot png" src="https://github.com/user-attachments/assets/060a3ed4-b6e3-48f3-a5fc-0b4be34451b7" />


---

### Repository Structure
```text
popgen-fstats-mds-pipeline/
├── README.md
├── LICENSE
├── scripts/
│   ├── 01_run_admixtools.sh       # Subset extraction + automated f3/f4 execution
│   └── 02_fstats_mds_analysis.R   # Distance matrix, cmdscale MDS & f4 forest plots
└── outputs/
    ├── tables/
    │   ├── f3_distance_matrix.csv
    │   └── f4_summary_results.csv
    └── figures/
        ├── mds_anatolia_neolithic.png
        └── f4_statistics_plot.png
