#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"
python3 gdp_to_gtp.py --input-dir input --pattern "*.pdb" --out-dir output --cache gtp_all_rcsb_cache.tsv
