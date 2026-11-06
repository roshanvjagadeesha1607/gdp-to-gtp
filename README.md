# GDP-to-GTP Structure Preparation

**Automated, structure-guided placement of GTP in GDP-containing PDB structures**

## Purpose

Preparing GTP-bound models of large protein complexes can require placing a nucleotide at many GDP-bound sites. Manually rebuilding GTP or docking it separately into each binding site is time-consuming, particularly in multichain or repeating assemblies. This tool automates that structural preparation: it identifies GDP residues in a PDB file and replaces them with candidate GTP conformations, avoiding repeated manual placement at individual sites.

The workflow is applicable to GDP-containing protein structures with one or multiple nucleotide-binding sites; it is not restricted to a particular protein family. It provides **coordinate-based GTP placement**, not a molecular docking calculation or a simulation of nucleotide exchange.

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
  --cache gtp_tubulin_cache.tsv
```

Outputs in `output/`:

- `*_GTP.pdb`: structure with successfully processed GDP sites replaced by GTP
- `sites_*.csv`: per-site placement results for each input structure
- `all_sites.csv`: combined per-site results

For further options, run `python3 gdp_to_gtp.py --help`.

## GTP reference library

The included `gtp_tubulin_cache.tsv` contains GTP conformations collected from experimentally determined structures using a tubulin-related search. **The placement algorithm is general, but this bundled reference library is not comprehensive for every GDP-binding protein.** For other target systems, a relevant GTP reference library may be preferable.

To prepare a custom library, list experimental PDB IDs (one per line) in `gtp_pdb_ids.txt` and run:

```bash
python3 -m pip install requests
python3 gdp_to_gtp.py \
  --input-dir input \
  --out-dir output \
  --cache custom_gtp_cache.tsv \
  --refresh \
  --entry-list gtp_pdb_ids.txt
```

## Limitations

A low RMSD and absence of overlaps under the program's distance-based criteria do not establish a physically stable GTP-binding pose. Some sites may retain overlaps, which are identified in the reports. The tool does not predict the conformational response of the protein to GDP-to-GTP exchange.

Inspect the resulting structures and perform appropriate validation, force-field preparation, and energy minimization before molecular dynamics or quantitative analysis. When reporting results, cite the experimental structures used as GTP references.
