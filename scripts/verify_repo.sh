#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

echo "[1/3] Python syntax"
while IFS= read -r -d '' f; do
  python -m py_compile "$f"
done < <(find . -type f -name '*.py' -print0)

echo "[2/3] Shell syntax"
while IFS= read -r -d '' f; do
  bash -n "$f"
done < <(find . -type f -name '*.sh' -print0)

echo "[3/3] Public-path hygiene"
if grep -RInE '/scratch/(gautschi|bell|negishi)/|lald@purdue\.edu|/depot/jgreeley/apps' .   --exclude-dir=.git --exclude='PROVENANCE.tsv' --exclude='verify_repo.sh'; then
  echo "Found unsanitized Purdue-specific path/account data." >&2
  exit 1
fi

echo "OK"
