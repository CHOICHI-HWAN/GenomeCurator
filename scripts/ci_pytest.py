"""Run pytest for CI and, when it fails, repeat the failing tests' short summary as a GitHub annotation (the job log itself needs admin rights to download; annotations are public)."""
import subprocess
import sys

r = subprocess.run([sys.executable, "-m", "pytest", "-q", "--tb=short", "-rf"], capture_output=True, text=True)
sys.stdout.write(r.stdout + r.stderr)
if r.returncode:
    tail = (r.stdout + r.stderr)[-1800:].replace("%", "%25").replace("\r", "").replace("\n", "%0A")
    print(f"::error title=pytest failed::{tail}")
sys.exit(r.returncode)
