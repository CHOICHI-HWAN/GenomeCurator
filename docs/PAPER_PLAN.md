# Paper analysis plan: the first clade Ib genomes from South Korea (draft)

**한국어:** [PAPER_PLAN.ko.md](PAPER_PLAN.ko.md)

Aim: report the first clade Ib genomes from Korea (three private genomes, lineage A.4.1, collected 2026-06-09 to 2026-07-10) and show (1) what is characteristic about them, (2) where and when they came from in genomic-epidemiology terms, (3) distinctive and notable mutations **with the evidence**. This document separates what is done, what is not, and what can and cannot be claimed. Figures are from this app's runs (clade Ib dataset, mutation-profile run 5, clock run 11); reproduce them from the same inputs before they go into a manuscript.

## 1. Results so far (done / partial / not done)
| # | Question | Status | Result |
|---|---|---|---|
| A1 | Clade / lineage | done | clade Ib, lineage A.4.1 (all three) |
| A2 | Mutation profile | done | all three have 103 mutations (17 APOBEC3-type, 17%), an identical set |
| A3 | Distinctive mutation | done | all A.4.1 markers present; the only mutation that sets them apart is G173208A (non-coding, APOBEC3 context): called by all four callers, 2 independent origins on the tree, elsewhere in only one German A.2 genome |
| A4 | Nearest genomes | done | on the tree the three Korean genomes (one extra mutation, G173208A) sit inside a polytomy of the A.4.1 subclade defined by G156869A (205 genomes: Portugal 97, Germany 77, Ireland 11, Belgium 7, Slovakia 6, USA 4 + Korea 3; first sampled 2026-02); equally close genomes exist from Germany 2026-02/03, Belgium, Ireland, Slovakia and the USA. The earlier "Belgium 2026-04" was only the date order among ties |
| B1 | Evolutionary rate of Ib | partial | BEAST strict-exponential 3.8e-5 /site/year (3.2 to 4.4e-5); root-to-tip APOBEC3-type 4.6 changes/year/genome (R^2 0.85) |
| B2 | tMRCA of Ib and A.4.1 | partial | root 2023.05 (2022.66-2023.42), A.4.1 2025.23 (2024.93-2025.51); 100 tips, one seed; A.1 has ESS 47 |
| B3 | Model choice | partial | strict-exponential highest (GSS: log BF 4.5 over second; estimators differ 2-5); a rerun with doubled paths is in progress |
| B4 | Temporal signal | done | BETS log BF: relaxed 42-45; strict PS 51.9 / SS 50.8 / GSS 49.9 |
| B5 | APOBEC3 concern | done | the temporal signal comes from APOBEC3-type changes (R^2 0.85); the others (1.4 per tip) show none (R^2 0.07) |
| C1 | Introduction time | partial | rough: the Korean lineage lies inside the G156869A subclade, so it split off after that subclade arose (first sample 2026-02) and before the first Korean sample (2026-06-09): about early 2026 to June 2026. The earlier "late 2025 to April 2026" estimate rested on mistaking Belgium April for the nearest relative and is withdrawn |
| C2 | Novel / adaptive mutations | partial | no strong candidate (few mutations with two or more independent origins; none above the neutral baseline) |
| D1 | Public Korean IIb (166 genomes) | not done | sequences not downloaded |

## 2. What can and cannot be claimed
**Can**: the three genomes are Ib / A.4.1 with an identical profile and all A.4.1 markers; the only distinguishing mutation is G173208A, an independently recurrent APOBEC3-favoured change; the clock of human transmission is carried by APOBEC3-type changes (4.6 per year, R^2 0.85), so the estimated tMRCAs date **the start of APOBEC3 accumulation, i.e. sustained human transmission**; Ib root early 2023 (2022.7-2023.4), A.4.1 early 2025 (2024.9-2025.5); the Korean genomes lie inside the G156869A A.4.1 subclade (mostly Portugal and Germany, first sampled 2026-02) with one extra mutation.
**Cannot**: the route or direction of introduction (three Korean genomes, sampling skewed to Europe, Uganda, Congo, unsampled intermediate hosts); that adaptive mutations were found; that the tool is the first of its kind (only MpoxRadar and PyR0 were read; VIRUS-MVP could not be); the virus's deep origin (the non-APOBEC rate is not measurable here); that these are the first Ib genomes from Korea until public databases are checked directly and the date of that check is written down.

## 3. Analyses to do (priority order)
- **P1 Introduction time, formally.** BEAST on A.4.1 genomes: the three Korean genomes + about 20 nearest genomes (by mutation difference) + a random A.4.1 sample; strict clock, exponential tree prior, rate prior from the APOBEC3 rate. Date and 95% HPD of the MRCA of the Korean clade and of the split from the nearest European genomes ("when the Korean lineage separated from the sampled European lineage", not "when it was introduced"). Sensitivity: number of neighbours, relaxed clock, seeds, alignment with recurrent APOBEC3 sites masked. Optional: discrete-trait phylogeography on country-by-week balanced samples, interpreted weakly.
- **P2 APOBEC3 sensitivity.** (a) is now in the code (`partition: apobec`) and passed a short check only; no real run yet. (a) a partitioned BEAST (APOBEC3-context sites vs others, separate rates) - not in the code yet; (b) mask recurrent APOBEC3 sites; (c) a clock from APOBEC3-type changes only vs all changes. Criterion: do the tMRCA intervals overlap.
- **P3 Stability of rate and tMRCA.** Read the four-model comparison with doubled paths; repeat the chosen model with other seeds (ESS >= 200, longer chains for low-ESS parameters); vary the sample (100 / 150 / 300); one unit system (per site per year and per genome per year); test whether TreeTime's older root (2020.4) comes from the mixed clock by rerunning it on APOBEC3-type changes.
- **P4 Evidence for characteristic / novel mutations.** One table of the 103 mutations (gene, amino-acid change, APOBEC3 or not, lineage marker or not; public tables use "Korean genome A/B/C", not private ids). G173208A: four callers, raw reads if available, context, the German A.2 genome (independent recurrence: two tree origins, absent in the other 275 A.4.1). Compare with the mpxv-lineages A.4.1 definition (all present / extra / missing). To call a mutation novel: search all public data (all clades), confirm with another caller, exclude sequencing error (resequencing).
- **P5 Sequence quality and integrity (before writing).** N fraction, covered length (195.6 kb vs 196.4 kb), ITR handling, Squirrel / Nextclade / CheckV QC in one table; whether the identical profile of the three is a real transmission cluster (epidemiological links) or a duplication in processing; distinguish "absent" from "unknown" at masked sites.
- **P6 Context (optional).** Public Korean IIb (166 genomes, 2022-08 to 2025-07, mostly C.1; about 33 MB) via *Get FASTA* after confirming count and size; for introduction and discussion only.
- **P7 Extras.** Site-wise selection tests outside the app on the HyPhy export; read PMC10256197 and the other Ib genomic-epidemiology reports in full before citing; cite only what has been read.

## 4. Figures and tables (draft)
1. Time-scaled Ib tree with A.4.1 and the Korean lineage (BEAST MCC tree, from P1). 2. Root-to-tip regression, APOBEC3-type vs other (redraw the app's card). 3. Mutations of the Korean lineage and its nearest genomes (position x time view or table). Table 1 the three genomes (lineage, mutations, APOBEC3 share, quality); Table 2 distinguishing mutation(s) with caller agreement and tree origins; Table 3 model comparison and rate / tMRCA sensitivity; Supplement: adaptive-candidate scoring and methods.

## 5. Reproducibility
Report every figure with the app's run ids and input hashes (`analysis_run.members_sha`); keep seeds and `settings.json` of BEAST / BETS; write BEAST version, BEAGLE, chain length, ESS and priors in the methods; record the curated database's intake history and the lineage-designation commit.

## 6. Privacy and release
Private genome ids stay out of public outputs (tables, supplements, public repositories) - use pseudonyms with a key held internally. The app never uploads to a public database; submission is a separate human step (the package in the *Submit / share* tab). The HyPhy export holds public genomes only but contains accessions: check before sharing.

## 7. Decisions needed
1. P1: number of neighbours (about 20 suggested) and whether to add phylogeography. 2. P2: implement the partitioned BEAST, or settle for masking and APOBEC3-only analyses. 3. Include the public Korean IIb genomes. 4. Use the wording "first from Korea" (after the public-data check). 5. Adaptive-candidate analysis in the main text or the supplement (the value of a negative result).

## 8. Expected compute
P1: BEAST with about 50 tips, hours. P2(a): similar after the implementation, at least twice to check priors. P3: seed repeats and sample sizes, a day or more (the 100-tip, 30-million-state chain took 21.7 hours). Run BEAST jobs one after another: several at once saturate the machine.
