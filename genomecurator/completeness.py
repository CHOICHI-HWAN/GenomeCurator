"""Genome completeness: Complete / Near complete / Partial.

NCBI only says Complete or Partial. Between them sits a large group of genomes that lack a few hundred bases at the ends (primers, UTRs,
terminal repeats): useful for most analyses, but not "complete". This module adds that middle class.

  Complete       NCBI says Complete (the submitter's claim is trusted, so genomes of a divergent lineage that are legitimately
                 shorter than the reference are not demoted)
  Near complete  not Complete, and the length is at least NEAR_COMPLETE_PCT % of the reference genome (default 90)
  Partial        everything else, including when there is no reference to compare with

The reference of a virus is its longest RefSeq genome among the collected records. The same rule is written in SQL in the analysis
view (`db.ANALYSIS_VIEW`), so Superset and the web screens always agree; keep the two in step."""
from . import config

CLASSES = ["Complete", "Near complete", "Partial"]


def reference_map(con) -> dict[int, tuple[str, int]]:
    """{watch_id: (accession, length)} of the longest RefSeq genome each virus has (not Removed)."""
    out: dict[int, tuple[str, int]] = {}
    for wid, acc, ln in con.execute(
            "SELECT watch_id, accession, genome_length FROM genome WHERE refseq_accessions IS NOT NULL AND refseq_accessions != '' "
            "AND genome_length > 0 AND IFNULL(genome_status,'') != 'Removed' ORDER BY genome_length"):
        out[wid] = (acc, ln)       # ascending length: the last one is the longest
    return out


def classify(status: str | None, length: int | float | None, ref_len: int | None, pct: float | None = None) -> tuple[str, float | None]:
    """-> (class, coverage % of the reference or None)."""
    pct = config.NEAR_COMPLETE_PCT if pct is None else pct
    cov = round(100 * float(length) / ref_len, 1) if length and ref_len else None
    if status == "Removed":
        return "Removed", cov
    if status == "Complete":
        return "Complete", cov
    if cov is not None and cov >= pct:
        return "Near complete", cov
    return "Partial", cov


def describe(ref: tuple[str, int] | None, pct: float | None = None) -> str:
    pct = config.NEAR_COMPLETE_PCT if pct is None else pct
    if not ref:
        return "기준 게놈(RefSeq)이 없어 NCBI의 Complete/Partial만 쓰고, Near complete는 판정하지 않습니다."
    return (f"완성도 기준: 기준 게놈 {ref[0]} ({ref[1]:,} bp). Complete = NCBI가 Complete로 표시 · "
            f"Near complete = Complete가 아니면서 길이가 기준의 {pct:g}% 이상 · Partial = 그 미만")
