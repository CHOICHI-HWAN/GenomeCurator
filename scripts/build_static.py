"""Write the public read-only snapshot to docs/index.html (GitHub Pages: main branch, /docs).

    python scripts/build_static.py [--out docs/index.html] [--repo https://github.com/<user>/<repo>]

The same page is also served live by the local app at http://127.0.0.1:8000/public.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomecurator.static_site import build_html  # noqa: E402

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "docs" / "index.html"))
    ap.add_argument("--repo", default="https://github.com/CHOICHI-HWAN/GenomeCurator")
    a = ap.parse_args()
    try:
        html = build_html(a.repo)
    except ValueError as e:
        sys.exit(str(e))
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    print(f"wrote {out} ({out.stat().st_size / 1e6:.2f} MB)")
