# Draft: lineages of clade Ib in nextstrain/mpox (not submitted)

*Status: a draft for the authors of this repository to review. Nothing has been posted to nextstrain/mpox or mpxv-lineages. The numbers come from one local run (a clade Ib dataset of 1,142 public and 3 private genomes, `nextstrain/mpox` at the commit in the run's `versions.json`, `mpxv-lineages/lineage-designation` at 6ad7809).*

## What is missing
1. `phylogenetic/defaults/clade-i/clades.tsv` defines clade I, clade Ia, sh2024 and clade Ib (one site, 48419T) only. The sh2023 lineages that mpxv-lineages has defined for clade Ib (A, A.1 ... A.5) never reach the clade I trees, whereas the hMPXV-1 build has its lineages in `defaults/clades.tsv`.
2. `phylogenetic/scripts/clades_renaming.py` has no case for the lineages of sh2023: a clade name it does not know falls into the final `else` branch and is written as clade `IIb`, outbreak `sh2017`, lineage = that name. With Ib lineages in `clades.tsv`, every node comes out as `clade_membership = IIb`.

## What was done locally
* A `clades.tsv` that keeps the build's own rows and adds, for each sh2023 lineage, a `clade` row (parent = the parent lineage) and the sites that introduced it (`nuc <position> <base>`), moved from NC_063383.1 to KJ642613.1 with a minimap2 genome-to-genome liftover: 7 lineages, 14 sites (+ the 2 sites of sh2023 itself). The file is `docs/drafts/clade-i-sh2023-clades.tsv` of GenomeCurator.
* None of the 14 sites falls into `defaults/clade-i/mask.bed` or the `from_beginning` / `from_end` trimming of the clade I build.
* A post-processing step that reads the branch labels and writes `clade_membership = Ib`, `outbreak = sh2023`, `lineage = <name>` (GenomeCurator, `cladedefs.fix_auspice`), because of point 2.

## Result of the local run (751 tips on the tree)
* By branch labels, **730 of 751 tips (97.2%)** get the lineage that the genome table holds (the lineage Pathoplexus / Nextclade assigned): 18 A.4.1 tips stay at A.4, 2 A.5 and 1 A.1 tips at A. Each lineage has only two defining sites here, so one missing or masked site moves a node up one level.
* Replaying the same rule on SNV calls outside the repeat region gives 91%; the missing part is A.5 and A.4.1, whose definitions rest on sites in the left inverted terminal repeat (3474, 4119), which that table leaves out. The build's alignment keeps those sites, hence the higher figure on the real tree.

## Suggested changes to the repository
1. Add the lineage rows to `phylogenetic/defaults/clade-i/clades.tsv` (generated from mpxv-lineages, with the coordinate move; the generator is small and can be shared).
2. Teach `clades_renaming.py` the sh2023 lineages: for example, a clade name found among the lineages of an Ib outbreak (a list in the config, or a naming convention such as `sh2023/<lineage>` for the rows of `clades.tsv`) is written as clade `Ib`, outbreak `sh2023`, lineage = the lineage.
3. Add a check to CI that every site of `clades.tsv` lies outside the mask and trimming of its build.

## Questions for the maintainers
* Are lineages of clade Ib meant to be defined in this repository, or only in mpxv-lineages and added by Nextclade?
* Should lineage rows use two sites each (as in mpxv-lineages today) or more, to be robust to missing data?
* Is `sh2023/<lineage>` an acceptable naming for the rows, so that the renaming script can tell them from hMPXV-1 lineages?

## Limits to state in a post
One dataset, one commit of each repository, three genomes of our own among the 751 tips; the 97.2% compares with the lineage of the genome table, which is itself an assignment, not a ground truth; the rule was not tested on a second clade I dataset.
