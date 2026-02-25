#!/usr/bin/env bash
set -euo pipefail
mkdir -p build/html
doxygen Doxyfile
echo
echo "Done. Open: docs/build/html/index.html"

