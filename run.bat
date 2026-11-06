@echo off
REM Convert GDP residues in PDB files placed in input\.
REM Any residual overlaps are reported in output\ CSV files.
cd /d "%~dp0"
python gdp_to_gtp.py --input-dir input --pattern "*.pdb" --out-dir output --cache gtp_all_rcsb_cache.tsv
echo.
echo Finished. Review the output folder and site reports.
pause
