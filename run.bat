@echo off
REM  Double-click this file, or run it from a terminal.
REM  Every *.pdb in input\ is scanned; every GDP found is replaced by the
REM  best-matching clash-free GTP.  Results go into output\.
cd /d "%~dp0"
python gdp_to_gtp.py --input-dir input --pattern *.pdb --out-dir output --cache gtp_tubulin_cache.tsv
echo.
echo Finished. Results are in the output folder.
pause
