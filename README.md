# `gdp_to_gtp.py` v2 — every GDP replaced by a clash-free experimental GTP

Finds, for each bound GDP, the experimentally determined GTP that superimposes on
it best, and writes the substituted structure. Clash-free is the default.

---

## What is new in v2

| | v1 | v2 |
|---|---|---|
| GDP sites per file | one | **all of them** |
| Identity of a nucleotide | chain + residue number | chain + residue + insertion code + **segID** |
| Clash handling | reported only | **avoided by default** |
| Clash environment | the one protein chain | **the whole structure** — all chains, ions, waters, existing GTP, and GTPs placed earlier in the run |
| Mg²⁺ | counted as a clash | recognised as coordination, not a clash |
| Existing GTP in the file | — | left completely untouched |
| Large structures | serial field overflowed | serials preserved above 99,999; CONECT omitted where serials are not unique |

`--pick lowest-rmsd` reproduces v1 exactly (verified byte-for-byte on your four
single-chain files).

---

## Why segID matters for your files

Your microtubule model uses chain IDs `1 2 3 4`, but each of those appears **13
times**, distinguished only by the segment ID in columns 73–76 (`N`, `O`, `P` …).
So "chain 1 residue 446" is not one nucleotide, it is thirteen. v2 keys on segID
and treats them separately. v1 would have merged all thirteen into one broken
416-atom residue.

---

## How a site is solved

Four stages, stopping at the first that gives zero van der Waals overlap:

1. **Direct.** All 6,875 library GTPs are ranked by symmetry-aware RMSD over the
   28 atoms shared with GDP, then screened in that order against every atom
   within 8 Å. The first clash-free one wins — so it is always the *best-fitting*
   clean donor, never just any clean donor.
2. **γ-phosphate rotated.** If nothing is clean outright, the γ-phosphate is
   rotated on its two freely rotatable bonds (about PB–O3B, and the terminal PO₃
   about O3B–PG). **This leaves all 28 GDP-shared atoms exactly where the fit put
   them** — only the four new atoms move. Offsets are tried in order of
   increasing size, so the change is always the smallest one that works. Only
   candidates whose other 28 atoms are already clean are considered, since
   rotating the γ-phosphate cannot fix a clash elsewhere.
3. **Bounded nudge.** Last resort: a small rigid-body move of the whole ligand,
   capped at `--max-shift` (default 0.40 Å RMSD).
4. **Reported.** If even that fails the site is written anyway and flagged
   `unresolved` with its residual overlap, in the console and in the CSV.

The steric model: N/O···N/O pairs may approach to 2.60 Å (hydrogen bonds),
metal···anything to 1.85 Å (coordination), everything else to the sum of van der
Waals radii minus 0.40 Å.

---

## Measured on your real microtubule

`output.pdb` — 368,368 atoms, **52 GDP** at residue 446, 52 GTP at residue 452,
52 Mg²⁺, 4 chain IDs × 13 segments, 260 TER records.

| | result |
|---|---|
| runtime | **2 min 43 s** |
| clash-free | **52 / 52** — 26 direct, 26 by γ-rotation, 0 needed a nudge |
| worst overlap, re-checked independently against all 368k atoms | **−0.014 Å (none)** |
| non-replaced lines | **367,172 of 367,172 byte-identical** |
| GDP left | 0 |
| pre-existing GTP (52), Mg²⁺ (52), TER (260) | all untouched |
| RMSD to the original GDP | 0.525 – 1.275 Å, mean 0.823 Å |

Stress-tested further on synthetic variants of the same model: 104 GDP sites
across 62 chain IDs, 3 min 22 s, 104/104 clash-free.

Note the RMSD range is wider than v1's because clash-freedom now takes priority
over pure RMSD. `--pick lowest-rmsd` restores the old behaviour.

### What the diagnostic columns tell you about this model

`sites_*.csv` now records, for every site, the RMSD of the *best-fitting* donor
and what stopped it being used:

| chain | donor used | its RMSD | best possible | blocked by |
|---|---|---|---|---|
| 1 | 6QVJ | 0.702 | 0.560 | O1G···Ala97:CB 2.37 Å |
| 2 | 9D5N | **1.275** | 0.478 | **O3G···Glu254:OE2 1.53 Å** |
| 3 | 9E5C | 0.525 | 0.457 | O2G···Asn99:N 2.05 Å |
| 4 | 5KMG | 0.791 | 0.559 | O1G···Thr143:OG1 2.26 Å |

Chains 1, 3 and 4 give up only 0.14–0.23 Å of fit to clear a mild contact.
Chain 2 is different: the blocker is **α-tubulin Glu254 of the neighbouring
subunit at 1.53 Å** — the catalytic glutamate that reaches into the β E-site.
In a GDP-state lattice the α T7 loop has moved in over the vacated
γ-phosphate, so no rigid GTP can be placed there in the canonical orientation
without that loop moving. The script therefore falls back to a donor whose
γ-phosphate points elsewhere, at 1.275 Å RMSD.

That is a property of your model, not of the search. Two reasonable ways to
handle chain 2:

* **Keep the clash-free result** (what you have) and let minimisation settle it.
* **Take the canonical fit instead** — run those sites with `--pick lowest-rmsd`
  and let a restrained minimisation push Glu254 back out, which is what happens
  physically on GDP→GTP. The overlap to relieve is large (1.53 Å), so minimise
  with the protein free in that region.

---

## Options

| flag | default | effect |
|---|---|---|
| `--pick` | `clash-free` | `lowest-rmsd` = v1 behaviour, best fit regardless of sterics |
| `--env-radius` | `8.0` | atoms within this distance of the ligand are used for the clash test |
| `--max-scan` | `0` (all) | how many ranked candidates to screen |
| `--torsion-candidates` | `40` | how many near-miss donors get the γ-rotation search |
| `--max-shift` | `0.40` | cap on the last-resort rigid-body nudge, Å |
| `--refine-steps` | `3000` | effort spent on that nudge |
| `--refine-candidates` | `5` | how many donors get the last-resort nudge |
| `--renumber` | `auto` | renumber serials only if the structure fits the 5-digit field |
| `--conect` | `auto` | write CONECT only when serials are unique |
| `--ligand` | `GDP` | residue name to replace |
| `--refresh` | off | re-download the GTP library from RCSB |
| `--entry-list` | — | CSV/TXT of PDB IDs to search (your `2_ligand_presence_matrix.csv` works directly) |
| `--search-rcsb` | off | with `--refresh`, query RCSB live instead of a fixed list |

---

## Output

- `<name>_GTP.pdb` — the substituted structure
- `sites_<name>.csv` — one row per GDP site: chain, residue, segID, donor PDB ID
  and its chain/residue, RMSD, worst overlap, how it was solved, the γ-torsion
  offsets applied, any rigid shift, the closest contact, **the RMSD of the
  best-fitting donor and what blocked it**, and the donor's resolution and method
- `all_sites.csv` — the same across every input file

---

## Two things to keep in mind

**Native GTP only.** GMPCPP (`G2P`) and GTP-γ-S (`GSP`) are excluded from the
library on purpose. They are geometrically valid γ-phosphate donors if you ever
want them — one line in the cache builder.

**A γ-rotation is chemically free, a nudge is not.** Sites marked
`gamma-rotated` have every GDP-shared atom exactly where the superposition put
it. Sites marked `refined` had the whole ligand moved slightly; the amount is in
the `rigid_shift_A` column. Both are still worth a short restrained
minimisation before production MD, but neither leaves a clash.
