# Reproducible figures

`genomecurator.treefigure.render()` draws a publication-style figure of the summary tree of a BEAST analysis from one function call, and `genomecurator.repro.manifest()` writes
a manifest next to every figure, so that a figure can be re-made from the same inputs and checked.

```python
from genomecurator import beastprep, nextstrain, treefigure
run = nextstrain.get_run(RUN_ID)                       # the clock run of the BEAST analysis
tree = beastprep.annotated_tree(run, "my-analysis")    # TreeAnnotator summary tree with country / lineage / date on the tips
treefigure.render(tree, "figures/tree", focus={"tipA", "tipB"}, labels={"tipA": "alias A"}, run=run,
                  intervals=[{"label": "MRCA", "lo": 2024.30, "hi": 2024.50, "on": "clade"}], title="...", collapse_min=3)
```

- **What is drawn.** The tree on a calendar axis, ordered as a ladder with the clade of interest last; side clades of `collapse_min` or more tips collapsed into triangles that list their countries; a chip with
  the lineage and the collection month for each tip; nodes with posterior >= `support`; the clade of interest shaded; the credible intervals you pass in `intervals` drawn on the clade (MRCA) or on its stem.
- **Private genomes.** A tip that is a private genome (the members of the run starting with `private:`, or the names in `private_ids`) must have an alias in `labels`; without one the figure is refused.
- **Manifest.** `<figure>.manifest.json` holds the SHA-256 and size of the tree and of every other input, the parameters, the program version and git commit (and whether files were modified since),
  the Python, library and BEAST / TreeAnnotator versions, and the SHA-256 of every file written. `repro.verify(manifest)` reports files that differ.
- **Same inputs, same file.** SVG output is byte-identical for the same inputs (no timestamps, fixed hash salt); PNG and PDF carry no creation date.
- **Install.** matplotlib is optional: `pip install genomecurator[figures]`.

What the manifest does not guarantee: a BEAST chain is stochastic. The same seed on the same BEAST and BEAGLE build and thread count should repeat a chain, but another CPU or thread count can change it; compare the
summaries (medians, HPD), and run at least one other seed.
