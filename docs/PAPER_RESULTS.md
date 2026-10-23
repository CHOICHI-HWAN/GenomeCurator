# Results for the paper: the first clade Ib genomes from South Korea (draft)

**한국어:** [PAPER_RESULTS.ko.md](PAPER_RESULTS.ko.md) (the fuller text; this is its English mirror)

Scope: **clade Ib only** (the clade IIb context analysis is on hold). The analyses so far, arranged as methods, results, limits and reproducibility information. Every figure comes from this app's runs; reproduce it from the same inputs before it goes into a manuscript. The plan and remaining work are in [PAPER_PLAN.md](PAPER_PLAN.md). In public outputs the private genomes are "Korean genome A, B, C" (the key to the internal ids stays outside the manuscript).

## 1. Data
Clade Ib dataset: 1,142 public genomes (Pathoplexus / NCBI GenBank as curated by the app; input hash `a2104adcfa50`) + 3 private Korean genomes (internal ids KDCA205-207), 1,145 in all. Korean genomes: collected 2026-06-09, 2026-07-03, 2026-07-10; 195,623 / 195,625 / 195,626 bp; N 0 / 50 / 54 (<= 0.03%); sequence statistics "good"; the Nextclade QC says "mediocre" (7-8 frameshifts), which 239 of the 275 other public A.4.1 genomes share (the known Ib deletions), so it is not specific to these genomes. References: KJ642613.1 (mutation calls, clade I build), PX667572 (Squirrel alignment), NC_063383.1 (Nextclade dataset).

## 2. Methods
**Alignment and calls**: Squirrel 1.3.2; SNVs against KJ642613.1 from the alignment (mutation-profile 0.2.1, minimap2 2.31); ITR left out (positions 6,501-189,942 used); genomes with more than 5% N left out; four callers compared per mutation (the run's alignment, snippy + lofreq, Nextclade 3.25.0 CLI moved to the run's coordinates). **Tree and lineages**: Nextstrain mpox clade I build (augur 34.2.1, IQ-TREE 3.1.4); lineages from mpxv-lineages commit 6ad7809 (sh2023); the build's `clades.tsv` has no Ib lineages, so one was made from the definitions (4.4). **Clock**: root-to-tip regression separately for APOBEC3-type (TC>TT, GA>AA) and other changes; TreeTime 0.12.1; BEAST 10.5.0 (BEAGLE), HKY + gamma(4), exponential-growth tree prior, strict clock, default vague priors, plus a log-normal rate prior (mean 3.82e-5, log sd 0.2) in the A.4.1 analysis; temporal signal by BETS and model choice by marginal likelihood of four clock x tree-prior combinations (GSS, PS, SS). Dates: day-precision only. **Distinguishing mutations**: the group's mutations against the rest and against the rest of the same lineage; independent origins counted on the tree (branches with at least two descendant tips). **Adaptive candidates** (exploratory, not a selection test): persistence, independent origins, countries, protein change, frequency, a baseline from synonymous / non-coding mutations, a country-fixed time trend.

## 3. Results
### 3.1 The Korean genomes
All three are clade Ib, lineage A.4.1 (the app, the Nextclade CLI - outbreakLineage sh2023/A.4.1 - and the curated table agree). Each has **103 mutations** against KJ642613.1, 17 of them APOBEC3-type (17%); **the mutation sets of the three are identical** (the sequences differ). 102 of the 103 are carried by the other 275 A.4.1 genomes too (A.4.1 markers), all A.4.1 markers are present. The only mutation distinguishing them is **G173208A** (non-coding, position 173,208, a GA>AA context): in all three, elsewhere only in one German A.2 genome (2025-01), in none of the other 275 A.4.1 genomes; called by all four callers; **two independent origins on the tree** (one settled, three Korean tips; one single German tip) - a recurrence at an APOBEC3-favoured site rather than shared descent. Call it a "recurrent mutation at an APOBEC3-favoured site", not a "novel mutation": public data of other clades were not searched and resequencing was not done.

### 3.2 Position on the tree; distribution of A.4.1
The three hang as a 3-tip clade (one extra mutation) inside a polytomy of the **A.4.1 subclade defined by G156869A (205 genomes: Portugal 97, Germany 77, Ireland 11, Belgium 7, Slovakia 6, USA 4, Korea 3; first sampled 2026-02)**. Equally close genomes (one mutation) come from Germany 2026-02/03, Belgium, Ireland, Slovakia and the USA: **the nearest single genome cannot be named** (by the mutation table 12 genomes tie, Belgium 9 and Germany 3, earliest 2026-04-10). A.4.1 (253 tips) is mostly European (Portugal 105, Germany 97, Ireland 20, Belgium 11, USA 7, Slovakia 6, Korea 3, others 4), first sampled 2025-11-10 (Belgium), arising from A.4 (Uganda 204 and others). Per-lineage distribution:

| lineage | tips | first - last collection | main countries |
|---|---|---|---|
| A | 154 | 2023-10 - 2026-05 | DR Congo 134, Uganda 10 |
| A.1 | 64 | 2024-08 - 2026-06 | DR Congo 57 |
| A.2 | 13 | 2024-11 - 2025 | Uganda, Germany, Canada 4 each |
| A.3 | 3 | 2025-02 - 2026-01 | South Africa 2 |
| A.4 | 239 | 2024-08 - 2026-05 | Uganda 204, Belgium 12 |
| A.4.1 | 253 | 2025-11 - 2026-09 | Portugal 105, Germany 97 |
| A.5 | 25 | 2024-07 - 2025 | Uganda 20 |

Sampling decides the monthly distribution: it shows neither direction nor speed of spread.

### 3.3 The clock ticks in APOBEC3-type changes
Root-to-tip regression (584 dated tips, 2023.8-2026.6):

| class | slope (changes/year/genome, 95%) | R^2 | line meets zero | mean per tip |
|---|---|---|---|---|
| APOBEC3-type | 4.63 (4.47-4.79) | 0.85 | 2023.7 | 6.4 |
| other | -0.78 | 0.07 | - | 1.4 |
| all | 3.85 (3.54-4.17) | 0.50 | 2023.1 | 7.8 |

82% of the changes from root to tip are APOBEC3-type; **the temporal signal comes from them alone**, the slow background rate cannot be measured in this data. By lineage the APOBEC3-type slope is A.4.1 5.7, A.4 6.2, A 5.7, A.1 5.0 (R^2 0.33-0.60). So the tMRCAs of clocks fitted to all changes date **the start of APOBEC3 accumulation in sustained human transmission**, not the virus's deep origin. Intervals are too narrow because tips share branches.

### 3.4 Temporal signal (BETS) and model choice
BETS log Bayes factor (dated - undated; 60 tips, 50 path steps x 250,000 states, seed 1): relaxed + exponential PS 44.9, SS 44.6, GSS 42.4; strict + exponential PS 51.9, SS 50.8, GSS 49.9. All far above 5: **strong temporal signal** (estimators differ by 2-2.5, flagged). Four-model comparison (dated analyses, 60 tips, log BF against the best, GSS):

| model | 250,000 states/step | 500,000 states/step | estimator spread (500,000) |
|---|---|---|---|
| strict-exponential | 0 | 0 | 4.1 |
| strict-constant | -4.5 | -2.5 | 0.75 |
| relaxed-exponential | -6.3 | -6.6 | 1.2 |
| relaxed-constant | -10.3 | -9.5 | 1.4 |

**The strict clock beats the relaxed one clearly** (6.6 or more); strict-exponential leads but only by 2.5 over the second and its estimators differ by 4.1: not decisive. Repeats with other seeds and ESS were not checked.

### 3.5 Rate and tMRCA of Ib (BEAST, strict-exponential, 100 tips)
30 million states, 10,000 samples, about 21.7 hours; most ESS 1,700-4,400.

| quantity | median (95% HPD) | ESS |
|---|---|---|
| clock rate (substitutions/site/year) | 3.8e-5 (3.2-4.4e-5) | 2,340 |
| root (Ib) tMRCA | 2023.08 (2022.66-2023.42) | 4,435 |
| lineage A.5 | 2023.91 (2023.69-2024.08) | 2,078 |
| lineage A.1 | 2024.16 (2023.96-2024.33) | **47 (low)** |
| lineage A.4 | 2024.50 (2024.31-2024.66) | 3,393 |
| **lineage A.4.1** | **2025.24 (2024.93-2025.51)** | 2,725 |

Do not compare the rate directly with the root-to-tip slope (observed changes to the tip vs the rate of a substitution model). **The ESS of A.1 is low**: do not rely on that interval. TreeTime (root 2020.4, rate 1.8e-5, R^2 0.04-0.37) differs strongly; BEAST has more behind it (strict clock preferred, high ESS, agreement with a report of Ib emerging around September 2023 that was seen only as a web search result), and TreeTime's difference may come from mixing APOBEC3-type and other changes (not checked by rerunning TreeTime on APOBEC3-type changes).

### 3.6 Time of introduction into Korea (A.4.1 + Korean genomes, BEAST)
Tips: 3 Korean + 5 near + 30 random A.4.1 + 6 A.4 = 44 (day-precision dates); strict + exponential; 20 million states, 10,000 samples, about 8.3 hours; ESS 2,600 or more in both runs.

| quantity | with rate prior | without |
|---|---|---|
| tMRCA of the Korean clade | 2026.40 (2026.32-2026.44) | 2026.40 (2026.33-2026.44) |
| **split of the Korean lineage (stem age)** | **2026.26 (2026.09-2026.33)** | 2026.27 (2026.11-2026.33) |
| MRCA of Korea + the two nearest Belgian genomes | 2025.98 (2025.79-2026.12) | 2026.00 (2025.82-2026.14) |
| tMRCA of A.4.1 | 2025.18 (2024.91-2025.45) | 2025.21 (2024.93-2025.47) |
| clock rate (/site/year) | 4.45e-5 (3.4-5.6) | 4.93e-5 (3.5-6.6) |

The Korean lineage split from a sampled European lineage between **February and the end of April 2026 (median early April)**; the three Korean genomes share an ancestor around the end of May 2026; the first Korean sample is 2026-06-09. The two settings agree. This is **when the Korean lineage last shared an ancestor with the sampled European lineage, not the time of introduction**: the introduction lies between that and the first Korean sample (June 2026), and intermediate hosts were not sampled.

### 3.7 Adaptive candidates (exploratory) - a negative result
Few of 868 mutations (9) arose independently twice or more on the tree; none that is persistent, protein-changing and in several countries and not APOBEC3-type. Ancestral, lineage-defining and APOBEC3-type mutations are kept apart; 11 of 15 candidates sit on nearly the same genomes, so there are 5 independent pieces of evidence. No mutation is above the neutral baseline at q < 0.05 (simulation: planted five-origin mutations recovered 120/120, false positives 0/12,000). **The result is negative**; the score is not a selection test.

## 4. Details and limits
**4.1 APOBEC3 and the clock**: the signal is in APOBEC3-type changes. Assumed: a constant APOBEC3 rate (slopes 5.0-6.2 by lineage); possible saturation at favoured sites (as G173208A). The partitioned BEAST (rates for APOBEC3-context sites and others separately) is implemented but **has not been run on real data**. **4.2 Sampling bias**: A.4.1 is skewed to Europe, Korea has three genomes; "nearest" and the stem age depend on who sequenced and are not a route. **4.3 Dates**: BEAST used day-precision dates only; Germany 97 and Ireland 20 (nearly half of A.4.1) have month-only dates and were left out; an option to use the middle of the month is implemented, its effect has not been compared with the day-only analysis. **4.4 Ib lineage definitions in Nextstrain**: the build's `clades.tsv` lacks Ib lineages, so the 7 sh2023 lineages with 14 sites were moved to KJ642613.1 (liftover; no site masked); by branch labels **730 of 751 tips (97.2%)** have the lineage of the curated table (Korean genomes: A.4.1); the workflow's `clades_renaming.py` reads unknown names as IIb, so a post-processing step fixes `clade_membership`; each lineage has only two sites, so one missing site moves a node up a level; an upstream proposal draft is `docs/drafts/upstream-proposal-clade-i-lineages.md` (not submitted). **4.5 Cannot be claimed**: route or direction of introduction, that adaptive mutations were found, the virus's deep origin, "first" for the tool (only MpoxRadar and PyR0 were read), "first Ib from Korea" (public data not yet checked directly). **4.6 Not yet done**: BEAST repeats (seeds, sample sizes), model-choice repeats, the APOBEC3 partitioned run and masking of recurrent sites, the month-date comparison, why the three genomes share one mutation set (transmission cluster or duplication in processing: epidemiological information and raw data), detailed sequence quality, resequencing of G173208A, reading the predecessor literature.

## 5. Reproducibility
| item | value |
|---|---|
| input hash | `a2104adcfa50` (all runs: 1,142 public + 3 private) |
| Squirrel alignment | run 4 (squirrel 1.3.2, IQ-TREE 3.1.4) |
| mutation profile | run 5 (mutation-profile 0.2.1, minimap2 2.31, reference KJ642613.1, Squirrel alignment) |
| caller agreement | snippy + lofreq flow run 7, Nextclade 3.25.0 |
| Nextstrain trees | run 1 (default `clades.tsv`), run 21 (Ib lineage `clades.tsv`; augur 34.2.1, Nextstrain CLI 10.5.0) |
| molecular clock | run 11 (TreeTime 0.12.1, Squirrel tree run 4) |
| BETS | relaxed run 14 (GSS + PS/SS), strict run 16 |
| model comparison | 250,000/step run 15, 500,000/step run 18 |
| BEAST (all of Ib) | run 17 (BEAST 10.5.0, 100 tips, 30 million states) |
| BEAST (A.4.1 + Korea) | run 19 (rate prior), run 20 (without); the recorded `n_tips` is the default 150, the real set is the 44 tips of `tip_names` |
| lineage definitions | mpxv-lineages/lineage-designation commit 6ad7809 |
| generated `clades.tsv` | `docs/drafts/clade-i-sh2023-clades.tsv` |
| result files | under `data/nextstrain/runs/Orthopoxvirus_monkeypox/...clock-20261021-075419/beast/` (`summary.json`, `a41-korea-*`, `bets/`) |

Each run's settings are in the app's run history and in each folder's `settings.json`.

## 6. Figures and tables (candidates)
1. Time-scaled Ib tree with A.4.1 and the Korean lineage (BEAST MCC tree built from run 19's `beast.trees`). 2. Root-to-tip regression, APOBEC3-type vs other (redraw the app's card). 3. Country-by-month table of the G156869A subclade (3.2). 4. Table 1: the three genomes (lineage, mutations, APOBEC3 share, quality); Table 2: the distinguishing mutation with caller agreement and tree origins; Table 3: model comparison and rate / tMRCA; Table 4: timing of the Korean clade. Supplement: adaptive-candidate scoring and methods, the caller-agreement table.

## 7. Before submission
- [ ] why the three genomes share one mutation set (epidemiological links, raw data, resequencing) - [ ] read-level confirmation of G173208A - [ ] check "first Ib from Korea" directly in public databases (record the date) - [ ] BEAST repeats (seeds, sample sizes) and the APOBEC3 partitioned sensitivity - [ ] comparison with an analysis using month-only dates - [ ] read the predecessor literature (PMC10256197, VIRUS-MVP, ...) and tidy the citations - [ ] no private id in any public table or supplement
