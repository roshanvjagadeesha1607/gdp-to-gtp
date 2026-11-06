#!/usr/bin/env bash
cd "$(dirname "$0")"
python3 gdp_to_gtp.py --input-dir input --pattern "*.pdb" --out-dir output --cache gtp_tubulin_cache.tsv
