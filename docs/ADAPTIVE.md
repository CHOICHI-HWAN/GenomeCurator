# Finding candidate adaptive mutations: design notes

**한국어:** [ADAPTIVE.ko.md](ADAPTIVE.ko.md)

How the *adaptive candidates* of the mutation picture can be made more trustworthy: what is implemented, which approaches the literature uses, the gap between them, and what to do next (in priority order).
The literature was checked at **abstract level only** (no paper was read in full). The citations show that an approach exists; figures and details must be checked against the papers before they are copied.

## 1. What "adaptive" can mean here
The tool **ranks candidates; it does not show adaptation.** Sequences alone cannot tell a mutation that raised fitness from a founder effect (one successful branch), a mutation hotspot or APOBEC3 editing, hitchhiking on a successful branch, or heavier sequencing of one place. So the tool should show several kinds of evidence separately, say which spurious signal each is weak against, and keep the one-line score as a sort aid.

## 2. Approaches in the literature (abstract level)
| Approach | Idea | Source |
|---|---|---|
| APOBEC3 signature | A marker of sustained human-to-human transmission; a host edit, kept apart from adaptation | reports on clades IIb and Ib (e.g. [Canada Ib case](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12207198/)) |
| Convergence | Score by how often the same mutation arose independently | [CoVtRec](https://arxiv.org/pdf/2106.07292v1) |
| Frequency growth | Candidates by how fast a mutation's frequency rises | [Wada et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7409725/) |
| Lineage-fitness regression | Model lineage shares, estimate mutations linked to fitness | [PyR0](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9161372/) |
| Convergence + antibody-binding sites | Rank spike changes | [VOCAL](https://evbc.uni-jena.de/events/vibiom2022/oral-presentation-abstracts/) |
| Distribution of fitness effects | Coding genes of mpox | [PMC10256197](https://pmc.ncbi.nlm.nih.gov/articles/PMC10256197/) (**full text not read**) |
| Site-wise positive selection tests | HyPhy FEL/MEME/FUBAR (codon alignment + tree) | not confirmed in this search; a known tool |

No study testing positive selection of clade Ib genes or recurrent mutations per OPG was found in this search.

## 3. Evidence types and what exists now
| # | Evidence | Now |
|---|---|---|
| E1 | Persistence: seen in >= 3 months and still in the last 6 | implemented |
| E2 | Recurrence: appears in several lineages that do not define it | replaced by E7 when a tree exists; lineage-name approximation only without a tree (on real data it counted ancestral mutations as recurring in 7 lineages) |
| E3 | Geographic spread: >= 3 countries | implemented |
| E4 | Protein change (non-synonymous / stop) | implemented |
| E5 | Frequency: rare early, common now; share holds up; **country-fixed logistic slope** (`trend.py`) | implemented; cannot correct sampling that changes inside a country |
| E6 | Classes: ancestral (nearly all genomes), lineage-defining, APOBEC3-type | implemented, kept apart |
| E7 | Number of independent origins on the tree | **implemented** (`origins.py`; only when the tree's reference is the run's reference) |
| E8 | Neutral baseline: independent-origin distribution of synonymous / non-coding mutations (Poisson / negative binomial), BH q | **implemented** (`origins.neutral_baseline`; the tail is a fitted model, extrapolated) |
| E9 | Linkage: separate mutations that travel together | partly: groups count as one piece of evidence ("N candidates = M independent pieces"); no regularised split |
| E10 | Selection tests: dN/dS, site-wise | exploratory per-gene dN/dS + **HyPhy input export** (`hyphyexport.py`: codon alignment per gene + tree, public genomes only); the tests run outside |

The present score (+2/+1/+1/+1/+1/+2) adds E1-E5 with weights that have no basis. The aim is to replace it with an evidence table strengthened by E7-E10.

## 4. Proposed methods
1. **Independent origins (E7, priority 1) - implemented**: from the dataset's Nextstrain tree and branch mutations, count the branches that carry each mutation; count an event as settled only if it has >= 2 descendant tips (single-tip events are indistinguishable from errors / within-host variants and are flagged separately); lift coordinates to the run reference. Replaces the lineage proxy of E2.
2. **Neutral baseline (E8, priority 2)**: synonymous and non-coding mutations are mostly neutral, so their recurrence is the chance level. Per mutation context (APOBEC3 trinucleotide, transition/transversion), compare the independent-origin count of a non-synonymous mutation with the neutral distribution -> empirical p, BH-corrected q. State "more independent origins than neutral mutations of the same context (q = ...)" instead of "score 6".
3. **Frequency growth corrected for sampling (E5, priority 3)**: binomial/multinomial logistic regression of the carrier share per month/week, stratified by country; reuse `lineagetree.growth`. Report the slope with a 95% interval; an interval spanning 0 means "cannot be called a rise".
4. **Linked mutations (E9, priority 4)**: groups (`co_travel`, Jaccard >= 0.6) count as one piece of evidence; later, PyR0-style regression of lineage fitness on mutation presence with regularisation; with little data report "cannot be separated".
5. **Selection tests (E10, priority 5)**: keep the exploratory per-gene dN/dS; support site-wise tests by exporting input for HyPhy (codon alignment + tree) rather than implementing them.
6. **Function**: gene functional categories (Supplementary Table 1) are labels only, never scored: there is no functional (DMS-type) data for mpox, so scoring by function would be circular.
7. **APOBEC3-type**: a separate class; interpret together with the APOBEC3 share per lineage/period (`outbreak.accumulation`), a measure of the pace of human transmission, not of selection.

## 5. Validation plan (decide before implementing)
1. **Recovery on simulated data** (like the planted signals in `scripts/verify_pipeline.py`): simulate trees, give some mutations a fitness advantage, measure how many are recovered and how many neutral ones come up as false positives.
2. **Negative control**: scoring synonymous / non-coding mutations as if they were candidates should give q < 0.05 about 5% of the time.
3. **Sampling robustness**: subsample by country and period and see whether the candidate list holds.
4. **State the missing positive control**: no mpox mutation has proven selection, so "recovered known adaptive mutations" cannot be tested; say so in the docs and on screen.
5. Re-run 1-3 whenever weights change (add to `./test.sh verify`).

## 6. Wording and display
Keep **"adaptive candidate"** (never shorten to "adaptive mutation"); show per candidate which evidence is present and which is absent (E1-E10 ticks), the one-line score only as a sort aid; count linked groups as one; always show the sampling-bias warning; do not claim "first" in a paper or announcement until the predecessor tools (PyR0, MpoxRadar, VIRUS-MVP) have been read and compared.

## 7. To do (priority order)
1. ~~E7 independent origins (tree-based)~~ - done. On run 5 (clade Ib) few mutations arose independently twice or more and no non-APOBEC mutation was both persistent, protein-changing and recurrent (the lineage proxy had counted ancestral mutations as recurrent in 7 lineages)
2. ~~E8 neutral baseline with q values~~ - done (Poisson / negative-binomial tail); no mutation reaches q < 0.05 on run 5
3. ~~E5 country-stratified growth regression~~ - done (per mutation, shown as "adjusted rise" in the candidate table)
4. ~~Evidence-tick UI; score as a sort aid~~ - done (columns persistence, recurrence, above neutral, spread, protein, rise, adjusted rise)
5. ~~Recovery test on simulated data in `./test.sh verify`~~ - done (`planted_adaptive`: planted five-origin mutations recovered 120/120, null false positives 0/12000)
6. E9 regularised regression for linked mutations - **partly** (groups count as one piece of evidence and the independent number is shown; the regression split is not done: on the current data most candidates form one group, so there is little to separate)
7. ~~E10 HyPhy input export~~ - done (codon alignments are cut by mapping the Squirrel alignment's reference columns onto the run reference's annotation; OPG185 translates to a protein start, MTQLPILL...)
8. Predecessor comparison - **partly done** (section 8: MpoxRadar and PyR0 read in full text; VIRUS-MVP blocked by a verification page and not read; the rest abstract level)

## 8. Comparison with predecessor tools (what was actually read)
| | MpoxRadar | PyR0 | This tool |
|---|---|---|---|
| Scope | mpox, NCBI GenBank daily (1,862 samples as of March 2023) | SARS-CoV-2, 6,466,300 GISAID genomes | mpox, the app's curated datasets |
| Read | full text (PMC10320163) | full text (PMC9161372) | - |
| Does | mutation frequency, linear-regression trend, map, shared/unique mutations of two datasets, three references, filters | hierarchical Bayesian multinomial logistic regression of lineage fitness and amino-acid mutation effects | mutation picture, table, flow diagram, evidence-by-evidence candidate ranking |
| Recurrence / convergence | not described | relies on convergence (a mutation arising in several lineages) | independent origins on the tree and a neutral baseline |
| Selection analysis | not described | fitness regression, not a selection test | exploratory dN/dS, HyPhy input export (tests run outside) |
| Sampling bias | users must account for it | lineage proportions in 14-day windows, regions separated | country-fixed regression, linked groups counted once |
| Limits stated by authors | countries that sequence more show more mutations | linear additive model (no epistasis), not causal, mutation-level estimates less consistent than lineage-level (rho = 0.48 between versions), top RBD mutations did not on average raise infectivity | see below |

VIRUS-MVP (`PMC12501772`) could not be read (a verification page blocked the fetch; it was not bypassed), so its features are unconfirmed; CoVtRec, Wada et al. and VOCAL are abstract level. **Weaker than PyR0**: this tool adds evidence into a simple score on a small mpox dataset, while PyR0 estimates mutation effects with a hierarchical Bayesian model on millions of genomes; PyR0's own limits (not causal, unstable mutation-level estimates, rankings that need not match phenotype) apply here as well. **Possibly different** (unconfirmed elsewhere): mpox clades / lineages / the paper's functional categories together, independent origins counted on the tree against a neutral baseline, evidence shown per column, and a bridge to HyPhy. This is not checked against every tool, so do not write "first".
