# GDP-to-GTP Structure Preparation

**Automated, structure-guided placement of GTP in GDP-containing PDB structures**

## Purpose

Preparing GTP-containing models of large protein complexes can be challenging when multiple chains contain GDP-bound sites. Manually replacing GDP with GTP at each binding site requires repeated nucleotide placement, structural alignment, and evaluation of steric interactions. In large multichain or repetitive assemblies, performing these steps individually becomes increasingly difficult and may introduce inconsistencies between binding sites.

This tool automates GDP-to-GTP structural replacement by identifying GDP residues in PDB structures and aligning experimentally determined GTP conformations with the existing nucleotide coordinates. Candidate placements are evaluated using RMSD-based ranking and steric clash screening to identify suitable GTP geometries. The approach reduces repetitive manual modeling and provides a systematic, reproducible workflow for generating candidate GTP-containing structures for subsequent structural validation and molecular simulation.

The workflow is not restricted to a particular protein family. It performs **structure-guided nucleotide placement**, not molecular docking or simulation of nucleotide exchange.

## Method

For each GDP residue with the required atoms, the program:

1. Identifies the GDP site using its chain, residue number, insertion code, and segment ID.
2. Superposes experimentally observed GTP conformations onto the **28 atoms shared by GDP and GTP**, considering equivalent phosphate-oxygen labels.
3. Ranks the placements by structural alignment RMSD and checks steric overlaps with the surrounding structure.
4. If needed, explores gamma-phosphate rotations and limited rigid-body adjustments to reduce overlaps.
5. Writes the selected GTP coordinates and reports the reference structure, RMSD, placement method, and remaining overlaps.

All eligible GDP sites are processed in one run. Non-target atomic coordinates are retained; fitted GTP atoms may differ from the original GDP coordinates.

## Requirements

- Python 3.8 or newer
- NumPy: `python3 -m pip install numpy`
- `requests` only when downloading or rebuilding the GTP reference library

Input files must use **PDB format**, with GDP represented as `HETATM` records and standard nucleotide atom names expected by the script.

## Usage

Place GDP-containing `.pdb` files in `input/` and run:

```bash
python3 gdp_to_gtp.py \
  --input-dir input \
  --out-dir output \
  --cache gtp_all_rcsb_cache.tsv
```

Outputs in `output/`:

- `*_GTP.pdb`: structure with successfully processed GDP sites replaced by GTP
- `sites_*.csv`: per-site placement results for each input structure
- `all_sites.csv`: combined per-site results

For further options, run `python3 gdp_to_gtp.py --help`.

## GTP reference library

The recommended `gtp_all_rcsb_cache.tsv` was assembled from experimentally determined RCSB PDB entries containing the ligand `GTP` (search performed 6 November 2026). The search returned 1,922 PDB entries; the resulting cache contains **9,257 complete GTP coordinate sets from 1,894 entries**. The associated `gtp_all_pdb_ids.txt` records the search-result entry IDs. Some entries could not be fetched or did not provide a complete usable GTP conformation.

A geometrically compatible GTP donor is not necessarily the preferred bound conformation for an unrelated protein; inspect each result, including any reported residual overlaps.

To refresh the general cache using the **saved list** of RCSB entries, install `requests` and run:

```bash
python3 -m pip install requests
printf 'HEADER    CACHE BUILD ONLY\nEND\n' > /tmp/gtp_cache_build.pdb
python3 gdp_to_gtp.py \
  --files /tmp/gtp_cache_build.pdb \
  --cache gtp_all_rcsb_cache.tsv \
  --entry-list gtp_all_pdb_ids.txt \
  --refresh
```

This refreshes the stored entry IDs; it does not discover newer entries added after the search. The current `--search-rcsb` option in the script still applies a tubulin keyword filter.

## Limitations

A low RMSD and absence of overlaps under the program's distance-based criteria do not establish a physically stable GTP-binding pose. Some sites may retain overlaps, which are identified in the reports. The tool does not predict the conformational response of the protein to GDP-to-GTP exchange.

Inspect the resulting structures and perform appropriate validation, force-field preparation, and energy minimization before molecular dynamics or quantitative analysis. When reporting results, cite the experimental structures used as GTP references.
