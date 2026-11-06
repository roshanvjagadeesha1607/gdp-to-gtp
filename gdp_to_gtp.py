#!/usr/bin/env python3
"""
gdp_to_gtp.py  (version 2)
==========================================================================
Replace every GDP in a structure with the experimentally observed GTP that
best overlaps it, choosing a donor that leaves NO steric clash.

Version 2 changes
-----------------
* Handles structures of any size and any number of chains.  Every GDP in the
  file is found and replaced, not just one.
* Nucleotides are identified by chain + residue number + insertion code +
  **segment ID**, so 13 copies of "chain 1 residue 446" that differ only by
  segID are treated as 13 separate nucleotides.
* Clash-free is now the DEFAULT.  The whole ranked candidate list is screened
  and the lowest-RMSD donor with no van der Waals overlap is used.  If none is
  clean, the gamma-phosphate is rotated about its two free bonds (which leaves
  all 28 GDP-shared atoms exactly where they were) and, only as a last resort,
  the ligand is nudged by a bounded rigid-body refinement.
* The clash environment is the whole structure: all protein atoms, all other
  ligands, any GTP already present, and the GTP molecules placed at earlier
  sites in the same run.  Mg(2+) coordination is recognised as bonding, not
  clashing.
* Everything that is not a GDP is copied through byte for byte -- protein,
  waters, ions, and any GTP that was already in the file.
* Atom serial numbers are preserved when the file is too large for the 5-digit
  PDB serial field (as yours is), and renumbered cleanly when it is not.

Quick start
-----------
    python gdp_to_gtp.py --input-dir input --out-dir output --cache gtp_tubulin_cache.tsv

Requirements
------------
    python 3.8+, numpy          (always)
    requests                    (only for --refresh / --search-rcsb)
==========================================================================
"""

import argparse
import csv
import glob
import itertools
import os
import sys
import time

import numpy as np

VERSION = "2.0"

# ---------------------------------------------------------------------------
# Chemistry
# ---------------------------------------------------------------------------

SHARED = ["PB", "O1B", "O2B", "O3B", "O3A", "PA", "O1A", "O2A",
          "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
          "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"]

GTP_ORDER = ["PG", "O1G", "O2G", "O3G", "O3B", "PB", "O1B", "O2B", "O3A",
             "PA", "O1A", "O2A", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'",
             "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "O6", "N1",
             "C2", "N2", "N3", "C4"]

GAMMA = ["PG", "O1G", "O2G", "O3G"]      # the four atoms GDP does not have

GTP_BONDS = [("PG", "O1G"), ("PG", "O2G"), ("PG", "O3G"), ("PG", "O3B"),
             ("O3B", "PB"), ("PB", "O1B"), ("PB", "O2B"), ("PB", "O3A"),
             ("O3A", "PA"), ("PA", "O1A"), ("PA", "O2A"), ("PA", "O5'"),
             ("O5'", "C5'"), ("C5'", "C4'"), ("C4'", "O4'"), ("C4'", "C3'"),
             ("O4'", "C1'"), ("C3'", "O3'"), ("C3'", "C2'"), ("C2'", "O2'"),
             ("C2'", "C1'"), ("C1'", "N9"), ("N9", "C8"), ("N9", "C4"),
             ("C8", "N7"), ("N7", "C5"), ("C5", "C6"), ("C5", "C4"),
             ("C6", "O6"), ("C6", "N1"), ("N1", "C2"), ("C2", "N2"),
             ("C2", "N3"), ("N3", "C4")]

# The three terminal beta-phosphate oxygens of GDP are chemically equivalent,
# as are the two terminal alpha oxygens, so their PDB names are arbitrary.
PERMUTATIONS = []
for _pb in itertools.permutations(["O1B", "O2B", "O3B"]):
    for _pa in itertools.permutations(["O1A", "O2A"]):
        _n = list(SHARED)
        _n[1], _n[2], _n[3] = _pb
        _n[6], _n[7] = _pa
        PERMUTATIONS.append(_n)

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "H": 1.20,
       "MG": 1.73, "CA": 2.31, "ZN": 1.39, "NA": 2.27, "K": 2.75, "CL": 1.75,
       "F": 1.47, "FE": 2.00, "MN": 2.00, "SE": 1.90, "BR": 1.85, "I": 1.98}

HBOND_MIN = 2.60      # N/O ... N/O may come this close: it is a hydrogen bond
METAL_MIN = 1.85      # Mg2+/metal ... O is coordination, not a clash
VDW_SLACK = 0.40      # tolerance subtracted from the sum of van der Waals radii
METALS = {"MG", "CA", "ZN", "MN", "FE", "NA", "K", "CO", "NI", "CU", "CD"}
WATERS = {"HOH", "WAT", "DOD", "H2O"}

MODELSERVER = "https://models.rcsb.org/v1/{entry}/atoms?label_comp_id=GTP&encoding=cif"
GRAPHQL = "https://data.rcsb.org/graphql"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"


# ---------------------------------------------------------------------------
# PDB handling
# ---------------------------------------------------------------------------

def element_of(line):
    e = line[76:78].strip().upper()
    if e:
        return e
    nm = line[12:16]
    if nm[0] != " " and nm[:2].strip().upper() in VDW and nm[:2].strip().upper() not in ("CA",):
        return nm[:2].strip().upper()
    s = nm.strip()
    return (s[1] if s[:1].isdigit() else s[:1]).upper()


class Structure:
    """A parsed PDB file.

    Ligand copies are keyed by (resname, chainID, resSeq, iCode, segID) so that
    residues sharing a chain and number but differing in segment -- as in a
    multi-protofilament microtubule model -- stay separate.
    """

    def __init__(self, path, ligand="GDP"):
        self.path = path
        self.name = os.path.basename(path)
        with open(path, "r") as fh:
            self.lines = fh.read().split("\n")
        if self.lines and self.lines[-1] == "":
            self.lines.pop()

        self.coords = []          # xyz of every ATOM/HETATM, in line order
        self.elems = []
        self.line_of = []         # index into self.lines
        self.resname = []
        self.label = []           # (resname, chain, resseq, atomname)
        self.sites = {}           # ligand key -> {atom_name: coord_index}
        self.site_lines = {}      # ligand key -> [line indices]

        for i, line in enumerate(self.lines):
            if len(line) < 54:
                continue
            rec = line[:6]
            if rec not in ("ATOM  ", "HETATM"):
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            k = len(self.coords)
            self.coords.append(xyz)
            self.elems.append(element_of(line))
            self.line_of.append(i)
            rn = line[17:20].strip()
            self.resname.append(rn)
            self.label.append((rn, line[21], line[22:27].strip(), line[12:16].strip()))
            if rec == "HETATM" and rn == ligand:
                key = (rn, line[21], line[22:26], line[26], line[72:76])
                self.sites.setdefault(key, {})[line[12:16].strip()] = k
                self.site_lines.setdefault(key, []).append(i)

        self.coords = np.asarray(self.coords, dtype=float) if self.coords \
            else np.zeros((0, 3))
        self.elems = np.array(self.elems, dtype=object)
        self.n_atoms = len(self.coords)

    def site_coords(self, key):
        return {a: self.coords[i] for a, i in self.sites[key].items()}

    def describe(self, key):
        _, ch, rs, ic, sg = key
        s = "chain %s res %s" % (ch.strip() or "-", rs.strip())
        if ic.strip():
            s += ic.strip()
        if sg.strip():
            s += " seg %s" % sg.strip()
        return s


def rotation_matrix(axis, angle_deg):
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    a = np.radians(angle_deg)
    c, s = np.cos(a), np.sin(a)
    x, y, z = axis
    return np.array([
        [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
        [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
        [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)]])


def rotate_about_bond(points, p1, p2, angle_deg):
    R = rotation_matrix(np.asarray(p2) - np.asarray(p1), angle_deg)
    return (np.asarray(points) - p1) @ R.T + p1


# ---------------------------------------------------------------------------
# Superposition
# ---------------------------------------------------------------------------

def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    V, S, W = np.linalg.svd(H)
    d = np.sign(np.linalg.det(V @ W))
    return V @ np.diag([1.0, 1.0, d]) @ W, pc, qc


def batch_fit(Pstack, ep, Qc, eq):
    """Batched Kabsch. Pstack (n,m,3) already centred; Qc (m,3) centred.

    Returns (rmsd (n,), R (n,3,3), centroids are handled by the caller).
    """
    M = np.einsum("nia,ib->nab", Pstack, Qc)
    U, S, Vh = np.linalg.svd(M)
    d = np.sign(np.linalg.det(M))
    tr = S[:, 0] + S[:, 1] + d * S[:, 2]
    rmsd = np.sqrt(np.maximum(0.0, (ep + eq - 2.0 * tr) / Qc.shape[0]))
    D = np.zeros((len(M), 3, 3))
    D[:, 0, 0] = D[:, 1, 1] = 1.0
    D[:, 2, 2] = d
    R = U @ D @ Vh
    return rmsd, R


# ---------------------------------------------------------------------------
# Steric model
# ---------------------------------------------------------------------------

def pair_limits(lig_elems, env_elems):
    """Minimum acceptable distance for every ligand/environment atom pair."""
    le = np.array([e.upper() for e in lig_elems], dtype=object)
    ee = np.array([e.upper() for e in env_elems], dtype=object)
    lr = np.array([VDW.get(e, 1.70) for e in le])
    er = np.array([VDW.get(e, 1.70) for e in ee])
    lim = lr[:, None] + er[None, :] - VDW_SLACK

    l_pol = np.array([e in ("N", "O") for e in le])
    e_pol = np.array([e in ("N", "O") for e in ee])
    lim = np.where(l_pol[:, None] & e_pol[None, :], HBOND_MIN, lim)

    l_met = np.array([e in METALS for e in le])
    e_met = np.array([e in METALS for e in ee])
    lim = np.where(l_met[:, None] | e_met[None, :], METAL_MIN, lim)
    return lim


def batch_worst(V, env_xyz, env_sq, lim2):
    """Screen a stack of ligand placements V (n,m,3) against the environment.

    Returns the largest squared-overlap per placement; <= 0 means clash-free.
    """
    if env_xyz.shape[0] == 0:
        return np.full(len(V), -9.0)
    d2 = ((V ** 2).sum(-1)[:, :, None]
          - 2.0 * np.einsum("nia,ka->nik", V, env_xyz)
          + env_sq[None, None, :])
    return (lim2[None, :, :] - d2).reshape(len(V), -1).max(1)


def worst_overlap(lig_xyz, env_xyz, lim):
    if env_xyz.shape[0] == 0:
        return -9.0, None
    d = np.linalg.norm(env_xyz[None, :, :] - lig_xyz[:, None, :], axis=2)
    ov = lim - d
    i, j = np.unravel_index(np.argmax(ov), ov.shape)
    return float(ov[i, j]), (int(i), int(j), float(d[i, j]))


# ---------------------------------------------------------------------------
# Candidate library
# ---------------------------------------------------------------------------

def read_cache(path):
    meta, copies = {}, []
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "M":
                meta[f[1]] = dict(res=f[2], method=f[3],
                                  title=f[4] if len(f) > 4 else "")
            elif f[0] == "L":
                at = {}
                for tok in f[4].split(";"):
                    nm, xyz = tok.split(":")
                    at[nm] = [float(v) for v in xyz.split(",")]
                copies.append(dict(entry=f[1], chain=f[2], seq=f[3], atoms=at))
    return meta, copies


def write_cache(path, entries, meta, copies):
    with open(path, "w") as fh:
        fh.write("#GTP_CACHE v1\tsource=RCSB ModelServer\tentries=%d\tcopies=%d\n"
                 % (len(entries), len(copies)))
        fh.write("#M\tpdb_id\tresolution\tmethod\ttitle\n")
        fh.write("#L\tpdb_id\tauth_chain\tauth_seq\tatom:x,y,z;...\n")
        for e in entries:
            m = meta.get(e, {})
            fh.write("M\t%s\t%s\t%s\t%s\n" %
                     (e, m.get("res", ""), m.get("method", ""),
                      m.get("title", "").replace("\t", " ")))
        for c in copies:
            atoms = ";".join("%s:%.3f,%.3f,%.3f" % (n, c["atoms"][n][0],
                                                    c["atoms"][n][1],
                                                    c["atoms"][n][2])
                             for n in c["atoms"])
            fh.write("L\t%s\t%s\t%s\t%s\n" % (c["entry"], c["chain"], c["seq"], atoms))


def parse_modelserver_cif(text):
    lines = text.split("\n")
    cols, i, found = [], 0, False
    while i < len(lines):
        if lines[i].strip() == "loop_":
            j, hdr = i + 1, []
            while j < len(lines) and lines[j].strip().startswith("_"):
                hdr.append(lines[j].strip())
                j += 1
            if hdr and hdr[0].startswith("_atom_site."):
                cols = [h.split(".", 1)[1] for h in hdr]
                i, found = j, True
                break
            i = j
            continue
        i += 1
    if not found:
        return []
    need = ["label_atom_id", "Cartn_x", "Cartn_y", "Cartn_z",
            "auth_asym_id", "auth_seq_id", "label_comp_id"]
    for n in need:
        if n not in cols:
            raise RuntimeError("atom_site loop is missing column %r" % n)
    idx = {n: cols.index(n) for n in need}

    groups = {}
    for line in lines[i:]:
        s = line.strip()
        if s.startswith("#"):
            break
        if not s or s.startswith("_") or s == "loop_":
            continue
        f = s.split()
        if len(f) < len(cols) or f[idx["label_comp_id"]] != "GTP":
            continue
        key = (f[idx["auth_asym_id"]], f[idx["auth_seq_id"]])
        nm = f[idx["label_atom_id"]]
        if nm.startswith('"') and nm.endswith('"'):
            nm = nm[1:-1]
        groups.setdefault(key, {})[nm] = [float(f[idx["Cartn_x"]]),
                                          float(f[idx["Cartn_y"]]),
                                          float(f[idx["Cartn_z"]])]
    return [dict(entry=None, chain=k[0], seq=k[1], atoms=v)
            for k, v in groups.items() if len(v) >= 28]


def download_library(entries, workers=12):
    import requests
    from concurrent.futures import ThreadPoolExecutor
    session = requests.Session()
    copies, failed, t0 = [], [], time.time()

    def one(entry):
        try:
            r = session.get(MODELSERVER.format(entry=entry.lower()), timeout=60)
            r.raise_for_status()
            got = parse_modelserver_cif(r.text)
            for g in got:
                g["entry"] = entry
            return entry, got
        except Exception as exc:                              # noqa: BLE001
            return entry, exc

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for n, (entry, got) in enumerate(pool.map(one, entries), 1):
            if isinstance(got, Exception):
                failed.append(entry)
            else:
                copies.extend(got)
            if n % 50 == 0:
                print("      %4d/%d entries, %d GTP copies, %.0fs"
                      % (n, len(entries), len(copies), time.time() - t0))
    if failed:
        print("      WARNING: %d entries not fetched: %s"
              % (len(failed), ", ".join(failed[:10])))
    return copies


def download_metadata(entries, batch=80):
    import requests
    meta = {}
    for a in range(0, len(entries), batch):
        ids = entries[a:a + batch]
        q = ('{entries(entry_ids:[%s]){rcsb_id rcsb_entry_info'
             '{resolution_combined} exptl{method} struct{title}}}'
             % ",".join('"%s"' % i for i in ids))
        try:
            j = requests.get(GRAPHQL, params={"query": q}, timeout=60).json()
            for e in (j.get("data") or {}).get("entries") or []:
                rc = e.get("rcsb_entry_info") or {}
                res = (rc.get("resolution_combined") or [""])[0]
                meta[e["rcsb_id"]] = dict(
                    res="" if res == "" else str(res),
                    method=((e.get("exptl") or [{}])[0] or {}).get("method", ""),
                    title=(e.get("struct") or {}).get("title", ""))
        except Exception as exc:                              # noqa: BLE001
            print("      metadata batch failed (%s)" % exc)
    return meta


def search_rcsb_tubulin_gtp():
    import requests
    body = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
            "operator": "exact_match", "value": "GTP"}},
        {"type": "terminal", "service": "full_text", "parameters": {"value": "tubulin"}}]},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 10000},
                            "results_content_type": ["experimental"]}}
    r = requests.post(SEARCH, json=body, timeout=120)
    r.raise_for_status()
    return [x["identifier"] for x in r.json().get("result_set", [])]


def entries_from_file(path):
    ids = []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        head = fh.readline()
        fh.seek(0)
        if "," in head:
            for row in csv.DictReader(fh):
                key = next((k for k in row if k and k.lower().strip()
                            in ("pdb_id", "pdbid", "entry_id", "id")), None)
                if not key:
                    continue
                if "GTP" in row and str(row["GTP"]).strip() not in ("1", "true", "True"):
                    continue
                ids.append(row[key].strip().upper())
        else:
            ids = [l.strip().upper() for l in fh
                   if l.strip() and not l.startswith("#")]
    seen, out = set(), []
    for i in ids:
        if i and i not in seen:
            seen.add(i)
            out.append(i)
    return out


# ---------------------------------------------------------------------------
# Placement for one GDP site
# ---------------------------------------------------------------------------

def gamma_scan_orders(step_major=10.0, step_minor=15.0):
    """Torsion offsets ordered by increasing magnitude, so the first solution
    found is the smallest change from the donor's own conformation."""
    major = [0.0]
    a = step_major
    while a <= 180.0:
        major += [a, -a]
        a += step_major
    minor = [0.0]
    b = step_minor
    while b < 120.0:
        minor += [b, -b]
        b += step_minor
    return major, minor


def solve_site(struct, key, copies, cand_stacks, placed, args):
    """Choose and place a GTP at one GDP site.

    Returns a dict describing the result, with 'xyz' the final 32 atom
    coordinates keyed by atom name.
    """
    site = struct.site_coords(key)
    Q = np.array([site[a] for a in SHARED])
    Qc = Q - Q.mean(0)
    qc = Q.mean(0)
    eq = float((Qc ** 2).sum())

    # --- rank every candidate by symmetry-aware RMSD ----------------------
    n = len(copies)
    best_r = np.full(n, np.inf)
    best_R = np.zeros((n, 3, 3))
    best_pc = np.zeros((n, 3))
    for P, ep, pcent in cand_stacks:
        r, R = batch_fit(P, ep, Qc, eq)
        upd = r < best_r
        if upd.any():
            best_r[upd] = r[upd]
            best_R[upd] = R[upd]
            best_pc[upd] = pcent[upd]

    order = np.argsort(best_r)

    # --- environment: every atom near this site, except its own -----------
    # Coarse box first (cheap over hundreds of thousands of atoms), then an
    # exact "within env_radius of any GDP atom" test.
    own = set(struct.sites[key].values())
    centre = Q.mean(0)
    reach = args.env_radius + float(np.abs(Q - centre).max()) + 4.0
    coarse = np.where(np.all(np.abs(struct.coords - centre) < reach, axis=1))[0]
    if len(coarse):
        dd = ((struct.coords[coarse][:, None, :] - Q[None, :, :]) ** 2).sum(2).min(1)
        near = dd < (args.env_radius + 4.0) ** 2
        env_idx = np.array([i for i, k in zip(coarse, near)
                            if k and i not in own], dtype=int)
    else:
        env_idx = np.zeros(0, dtype=int)
    env_xyz = struct.coords[env_idx]
    env_el = list(struct.elems[env_idx])
    # GTP molecules already placed at earlier sites in this run
    for pk, pxyz in placed.items():
        if pk == key:
            continue
        arr = np.array([pxyz[a] for a in GTP_ORDER])
        m = ((arr[:, None, :] - Q[None, :, :]) ** 2).sum(2).min(1) < \
            (args.env_radius + 4.0) ** 2
        if m.any():
            env_xyz = np.vstack([env_xyz, arr[m]]) if len(env_xyz) else arr[m]
            env_el += [("P" if a[0] == "P" else a[0]) for a, keep
                       in zip(GTP_ORDER, m) if keep]
    lig_el = ["P" if a[0] == "P" else a[0] for a in GTP_ORDER]
    lim = pair_limits(lig_el, env_el)
    lim2 = lim ** 2
    env_sq = (env_xyz ** 2).sum(1) if len(env_xyz) else np.zeros(0)

    def place(idx):
        c = copies[idx]
        X = np.array([c["atoms"][a] for a in GTP_ORDER])
        return (X - best_pc[idx]) @ best_R[idx] + qc

    def as_dict(A):
        return {a: A[i] for i, a in enumerate(GTP_ORDER)}

    # What the very best-fitting donor would have been, and what stops it being
    # used.  Purely diagnostic -- it tells you when the structure itself, not
    # the library, is the limiting factor.
    i0 = int(order[0])
    A_best = place(i0)
    w0, _ = worst_overlap(A_best, env_xyz, lim)
    diag = dict(best_rmsd=float(best_r[i0]),
                best_entry=copies[i0]["entry"],
                best_block="" if w0 <= 0 else
                describe_contact(A_best, env_xyz, lim, env_idx, struct))

    # --- 1. screen the ranked list in RMSD order, stopping at the first
    #        candidate that is completely clash-free -------------------------
    scan = order if args.max_scan <= 0 else order[:args.max_scan]
    worst_all = np.full(len(scan), np.inf)     # squared-overlap screen values
    worst_core = np.full(len(scan), np.inf)    # same, ignoring the gamma-phosphate
    CH = 256
    hit = None
    if args.pick != "lowest-rmsd":
        for s0 in range(0, len(scan), CH):
            idx = scan[s0:s0 + CH]
            X = np.stack([np.array([copies[i]["atoms"][a] for a in GTP_ORDER])
                          for i in idx])
            A = np.einsum("nia,nab->nib", X - best_pc[idx][:, None, :],
                          best_R[idx]) + qc
            if len(env_xyz) == 0:
                worst_all[s0:s0 + len(idx)] = -9.0
                hit = 0 if hit is None else hit
                break
            d2 = ((A ** 2).sum(-1)[:, :, None]
                  - 2.0 * np.einsum("nia,ka->nik", A, env_xyz)
                  + env_sq[None, None, :])
            ov = lim2[None, :, :] - d2
            w = ov.reshape(len(idx), -1).max(1)
            worst_all[s0:s0 + len(idx)] = w
            worst_core[s0:s0 + len(idx)] = ov[:, len(GAMMA):, :].reshape(
                len(idx), -1).max(1)
            clean_local = np.where(w <= 0.0)[0]
            if len(clean_local):
                hit = s0 + int(clean_local[0])
                break

    if hit is not None:
        k = int(hit)
        i = int(scan[k])
        A = place(i)
        wtrue, _ = worst_overlap(A, env_xyz, lim)
        return dict(entry=copies[i]["entry"], chain=copies[i]["chain"],
                    seq=copies[i]["seq"], rmsd=float(best_r[i]),
                    worst=float(wtrue), how="direct",
                    tors=(0.0, 0.0), shift=0.0, xyz=as_dict(A),
                    contact=describe_contact(A, env_xyz, lim, env_idx, struct), **diag)

    if args.pick == "lowest-rmsd":
        i = int(order[0])
        A = place(i)
        w, _ = worst_overlap(A, env_xyz, lim)
        return dict(entry=copies[i]["entry"], chain=copies[i]["chain"],
                    seq=copies[i]["seq"], rmsd=float(best_r[i]),
                    worst=float(w), how="lowest-rmsd", tors=(0.0, 0.0),
                    shift=0.0, xyz=as_dict(A),
                    contact=describe_contact(A, env_xyz, lim, env_idx, struct), **diag)

    # --- 2. rotate the gamma-phosphate on its two free bonds --------------
    # This leaves all 28 GDP-shared atoms exactly where the fit put them.
    major, minor = gamma_scan_orders()
    gi = [GTP_ORDER.index(a) for a in GAMMA]
    ti = [GTP_ORDER.index(a) for a in ("O1G", "O2G", "O3G")]
    i_pb, i_o3b, i_pg = (GTP_ORDER.index("PB"), GTP_ORDER.index("O3B"),
                         GTP_ORDER.index("PG"))
    combos = [(t1, t2) for t1 in major for t2 in minor]   # ordered by |t1|,|t2|
    # Rotating the gamma-phosphate can only fix a clash that involves the
    # gamma-phosphate, so only consider candidates whose other 28 atoms are
    # already clean. They stay in RMSD order, so the best fit is tried first.
    fixable = [k for k in range(len(scan)) if worst_core[k] <= 0.0]
    for k in fixable[:args.torsion_candidates]:
        i = int(scan[k])
        A0 = place(i)
        V = np.repeat(A0[None, :, :], len(combos), axis=0)
        for c, (t1, t2) in enumerate(combos):
            g = A0[gi]
            if t1:
                g = rotate_about_bond(g, A0[i_pb], A0[i_o3b], t1)
            V[c, gi] = g
            if t2:
                o3b, pg = V[c, i_o3b], V[c, i_pg]
                V[c, ti] = rotate_about_bond(V[c, ti], o3b, pg, t2)
        w = batch_worst(V, env_xyz, env_sq, lim2)
        good = np.where(w <= 0.0)[0]
        if len(good):
            c = int(good[0])
            B = V[c]
            wtrue, _ = worst_overlap(B, env_xyz, lim)
            return dict(entry=copies[i]["entry"], chain=copies[i]["chain"],
                        seq=copies[i]["seq"], rmsd=float(best_r[i]),
                        worst=float(wtrue), how="gamma-rotated",
                        tors=combos[c], shift=0.0, xyz=as_dict(B),
                        contact=describe_contact(B, env_xyz, lim,
                                                 env_idx, struct), **diag)

    # --- 3. bounded rigid-body nudge of the best candidate ----------------
    # Try the few most promising donors; for each, nudge the whole ligand within
    # the allowed bound, then let the gamma-phosphate rotate on top of that.
    rng = np.random.default_rng(0)
    BATCH = 120
    rounds = max(1, args.refine_steps // BATCH)
    fallback = None
    for kk in np.argsort(worst_all)[:max(1, args.refine_candidates)]:
        i = int(scan[int(kk)])
        A0 = place(i)
        cen = A0.mean(0)
        step_t, step_r = 0.15, 3.0
        cur_w = worst_overlap(A0, env_xyz, lim)[0]
        cur_A, cur_d = A0, 0.0
        for _ in range(rounds):
            V = np.empty((BATCH, len(GTP_ORDER), 3))
            disp = np.empty(BATCH)
            for b in range(BATCH):
                R = rotation_matrix(rng.normal(0, 1, 3), rng.normal(0, step_r))
                V[b] = (cur_A - cen) @ R.T + cen + rng.normal(0, step_t, 3)
                disp[b] = np.sqrt(((V[b] - A0) ** 2).sum(1).mean())
            w = batch_worst(V, env_xyz, env_sq, lim2)
            w[disp > args.max_shift] = np.inf
            b = int(np.argmin(w))
            if w[b] < cur_w:
                cur_w, cur_A, cur_d = float(w[b]), V[b], float(disp[b])
                if cur_w <= 0.0:
                    break
            step_t *= 0.75
            step_r *= 0.75
        used_t = (0.0, 0.0)
        if cur_w > 0.0:                      # gamma rotation on the nudged pose
            V = np.repeat(cur_A[None, :, :], len(combos), axis=0)
            for c, (t1, t2) in enumerate(combos):
                g = cur_A[gi]
                if t1:
                    g = rotate_about_bond(g, cur_A[i_pb], cur_A[i_o3b], t1)
                V[c, gi] = g
                if t2:
                    V[c, ti] = rotate_about_bond(V[c, ti], V[c, i_o3b],
                                                 V[c, i_pg], t2)
            wv = batch_worst(V, env_xyz, env_sq, lim2)
            c = int(np.argmin(wv))
            if wv[c] < cur_w:
                cur_A, cur_w, used_t = V[c], float(wv[c]), combos[c]
        wtrue = worst_overlap(cur_A, env_xyz, lim)[0]
        cand = dict(entry=copies[i]["entry"], chain=copies[i]["chain"],
                    seq=copies[i]["seq"], rmsd=float(best_r[i]),
                    worst=float(wtrue),
                    how="refined" if wtrue <= 0 else "unresolved",
                    tors=used_t, shift=cur_d, xyz=as_dict(cur_A),
                    contact=describe_contact(cur_A, env_xyz, lim, env_idx, struct),
                    **diag)
        if wtrue <= 0.0:
            return cand
        if fallback is None or wtrue < fallback["worst"]:
            fallback = cand
    return fallback


def describe_contact(A, env_xyz, lim, env_idx, struct):
    if env_xyz.shape[0] == 0:
        return ""
    d = np.linalg.norm(env_xyz[None, :, :] - A[:, None, :], axis=2)
    ov = lim - d
    i, j = np.unravel_index(np.argmax(ov), ov.shape)
    if j < len(env_idx):
        rn, ch, rs, an = struct.label[env_idx[j]]
        who = "%s%s:%s" % (rn, rs, an)
    else:
        who = "placed GTP"
    return "%s...%s %.2f A" % (GTP_ORDER[i], who, d[i, j])


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

def write_structure(struct, results, out_path, ligand, renumber, conect, remarks):
    """Emit the structure with every solved GDP site replaced by GTP."""
    replaced_lines = {}
    for key, res in results.items():
        first = struct.site_lines[key][0]
        for ln in struct.site_lines[key]:
            replaced_lines[ln] = None
        replaced_lines[first] = (key, res)

    out = list(remarks)
    serial = 0
    ser_by_key = {}

    def fmt(name, ser, tmpl, xyz):
        el = "P" if name[0] == "P" else name[0]
        a4 = (" " + name).ljust(4) if len(name) < 4 else name
        s = min(ser, 99999) if ser is not None else 99999
        return ("HETATM%5d %s GTP %s%s%s   %8.3f%8.3f%8.3f%s%s      %s%2s  "
                % (s, a4, tmpl[21], tmpl[22:26], tmpl[26],
                   xyz[0], xyz[1], xyz[2], tmpl[54:60], tmpl[60:66],
                   tmpl[72:76], el))

    for i, line in enumerate(struct.lines):
        if i in replaced_lines:
            item = replaced_lines[i]
            if item is None:
                continue
            key, res = item
            tmpl = struct.lines[i]
            base = None
            if not renumber:
                try:
                    base = int(tmpl[6:11])
                except ValueError:
                    base = None
            ser_by_key[key] = []
            for n_i, nm in enumerate(GTP_ORDER):
                if renumber:
                    serial += 1
                    s = serial
                else:
                    s = (base + n_i) if base is not None else 99999
                ser_by_key[key].append(min(s, 99999))
                out.append(fmt(nm, s, tmpl, res["xyz"][nm]))
            continue
        rec = line[:6]
        if rec in ("ATOM  ", "HETATM"):
            if renumber:
                serial += 1
                out.append(line[:6] + "%5d" % serial + line[11:])
            else:
                out.append(line)
        elif rec[:3] == "TER":
            if renumber:
                serial += 1
                out.append("TER   " + "%5d" % serial + line[11:])
            else:
                out.append(line)
        elif rec == "CONECT" or line[:3] == "END":
            continue
        elif not line.strip():
            continue
        else:
            out.append(line)

    if conect:
        for key, sers in ser_by_key.items():
            idx = {nm: sers[i] for i, nm in enumerate(GTP_ORDER)}
            conn = {n: [] for n in GTP_ORDER}
            for a, b in GTP_BONDS:
                conn[a].append(idx[b])
                conn[b].append(idx[a])
            for nm in GTP_ORDER:
                out.append("CONECT" + "%5d" % idx[nm]
                           + "".join("%5d" % c for c in sorted(conn[nm])))
    out.append("END")
    with open(out_path, "w") as fh:
        fh.write("\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Replace every GDP with the best-overlapping, clash-free "
                    "experimental GTP (v%s)." % VERSION,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--input-dir", default=".")
    ap.add_argument("--pattern", default="*.pdb")
    ap.add_argument("--files", nargs="*", default=None)
    ap.add_argument("--out-dir", default="GTP_replaced")
    ap.add_argument("--ligand", default="GDP")
    ap.add_argument("--cache", default="gtp_tubulin_cache.tsv")
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--entry-list", default=None)
    ap.add_argument("--search-rcsb", action="store_true")
    ap.add_argument("--pick", choices=["clash-free", "lowest-rmsd"],
                    default="clash-free",
                    help="'clash-free' (default) uses the lowest-RMSD donor that "
                         "leaves no van der Waals overlap; 'lowest-rmsd' ignores "
                         "sterics entirely")
    ap.add_argument("--env-radius", type=float, default=8.0,
                    help="atoms within this distance of any ligand atom are "
                         "included in the clash test, in Angstrom")
    ap.add_argument("--max-scan", type=int, default=0,
                    help="how many ranked candidates to screen (0 = all)")
    ap.add_argument("--torsion-candidates", type=int, default=40,
                    help="how many near-miss candidates get a gamma-phosphate "
                         "torsion search when none is clean outright")
    ap.add_argument("--refine-steps", type=int, default=3000)
    ap.add_argument("--refine-candidates", type=int, default=5,
                    help="how many donors get the last-resort nudge")
    ap.add_argument("--max-shift", type=float, default=0.40,
                    help="largest RMSD displacement the last-resort rigid-body "
                         "refinement may introduce, in Angstrom")
    ap.add_argument("--renumber", choices=["auto", "yes", "no"], default="auto",
                    help="'auto' renumbers atom serials only when the structure "
                         "fits the 5-digit PDB serial field")
    ap.add_argument("--conect", choices=["auto", "yes", "no"], default="auto",
                    help="'auto' writes CONECT records only when serials are unique")
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()

    print("gdp_to_gtp.py v%s" % VERSION)

    paths = args.files if args.files else sorted(
        glob.glob(os.path.join(args.input_dir, args.pattern)))
    if not paths:
        sys.exit("No input files matched.")

    # ---- library ---------------------------------------------------------
    if os.path.exists(args.cache) and not args.refresh:
        print("Loading GTP library from cache: %s" % args.cache)
        meta, copies = read_cache(args.cache)
    else:
        if args.search_rcsb:
            entries = search_rcsb_tubulin_gtp()
        elif args.entry_list:
            entries = entries_from_file(args.entry_list)
        elif os.path.exists(args.cache):
            entries = sorted(read_cache(args.cache)[0])
        else:
            sys.exit("Nothing to search: give --entry-list, --search-rcsb, or a cache.")
        print("Downloading GTP coordinates for %d entries ..." % len(entries))
        copies = download_library(entries, args.workers)
        meta = download_metadata(entries)
        write_cache(args.cache, entries, meta, copies)
        print("Cache written: %s" % args.cache)

    copies = [c for c in copies if all(a in c["atoms"] for a in GTP_ORDER)]
    if not copies:
        sys.exit("No complete GTP molecules in the library.")
    print("Library: %d GTP molecules from %d PDB entries"
          % (len(copies), len({c["entry"] for c in copies})))

    # pre-centre once per permutation: (centred stack, sum-of-squares, centroid)
    stacks = []
    for names in PERMUTATIONS:
        P = np.array([[c["atoms"][nm] for nm in names] for c in copies])
        cent = P.mean(axis=1)
        Pc = P - cent[:, None, :]
        stacks.append((Pc, (Pc ** 2).sum(axis=(1, 2)), cent))

    os.makedirs(args.out_dir, exist_ok=True)
    grand = []

    for path in paths:
        t0 = time.time()
        print("\n" + "=" * 74)
        print("Reading %s" % path)
        st = Structure(path, args.ligand)
        print("  %d atoms, %d %s site(s)" % (st.n_atoms, len(st.sites), args.ligand))
        n_gtp = sum(1 for rn in st.resname if rn == "GTP")
        if n_gtp:
            print("  %d GTP atoms already present -- these are left untouched"
                  % n_gtp)
        if not st.sites:
            print("  nothing to do (no %s found)" % args.ligand)
            continue

        renumber = (args.renumber == "yes") or (
            args.renumber == "auto" and st.n_atoms + 4 * len(st.sites) <= 99999)
        conect = (args.conect == "yes") or (args.conect == "auto" and renumber)
        if not renumber:
            print("  structure exceeds the 5-digit PDB serial field: original "
                  "atom serials preserved, CONECT omitted")

        results, rows = {}, []
        for n_site, key in enumerate(sorted(st.sites), 1):
            missing = [a for a in SHARED if a not in st.sites[key]]
            if missing:
                print("  [%d/%d] %s -- SKIPPED, missing atoms: %s"
                      % (n_site, len(st.sites), st.describe(key), ",".join(missing)))
                continue
            res = solve_site(st, key, copies, stacks,
                             {k: v["xyz"] for k, v in results.items()}, args)
            results[key] = res
            m = meta.get(res["entry"], {})
            flag = {"direct": "clash-free",
                    "gamma-rotated": "clash-free (gamma rotated %+.0f/%+.0f deg)"
                                     % res["tors"],
                    "refined": "clash-free (nudged %.2f A)" % res["shift"],
                    "lowest-rmsd": "lowest RMSD, sterics ignored",
                    "unresolved": "RESIDUAL OVERLAP %+.2f A" % res["worst"]}[res["how"]]
            print("  [%d/%d] %-28s <- %s %s %s  RMSD %.3f A  %s"
                  % (n_site, len(st.sites), st.describe(key), res["entry"],
                     res["chain"], res["seq"], res["rmsd"], flag))
            if res["rmsd"] > res["best_rmsd"] + 0.15 and res["best_block"]:
                print("           best possible fit was %s at %.3f A, blocked by %s"
                      % (res["best_entry"], res["best_rmsd"], res["best_block"]))
            rows.append([st.name, key[1], key[2].strip(), key[4].strip(),
                         res["entry"], res["chain"], res["seq"],
                         "%.4f" % res["rmsd"], "%.3f" % res["worst"],
                         res["how"], "%.0f" % res["tors"][0],
                         "%.0f" % res["tors"][1], "%.3f" % res["shift"],
                         res["contact"], "%.4f" % res["best_rmsd"],
                         res["best_entry"], res["best_block"],
                         m.get("res", ""), m.get("method", "")])
            grand.append(rows[-1])

        if not results:
            continue

        base = os.path.splitext(st.name)[0]
        for suf in ("_GDP", "_gdp", args.ligand, args.ligand.lower()):
            if base.endswith(suf):
                base = base[:-len(suf)].rstrip("_")
                break
        out_pdb = os.path.join(args.out_dir, base + "_GTP.pdb")
        remarks = [
            "REMARK   1 %d %s LIGAND(S) REPLACED BY GTP BY gdp_to_gtp.py v%s"
            % (len(results), args.ligand, VERSION),
            "REMARK   1 EACH DONOR IS THE LOWEST-RMSD EXPERIMENTAL GTP LEAVING NO VDW OVERLAP",
            "REMARK   1 SUPERPOSITION ON THE 28 ATOMS SHARED BY %s AND GTP" % args.ligand,
            "REMARK   1 PROTEIN, IONS, WATERS AND PRE-EXISTING GTP ARE UNCHANGED",
        ]
        write_structure(st, results, out_pdb, args.ligand, renumber, conect, remarks)
        print("  written: %s   (%.1f s)" % (out_pdb, time.time() - t0))

        with open(os.path.join(args.out_dir, "sites_%s.csv" % base), "w",
                  newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["input_file", "chain", "resseq", "segid", "gtp_pdb_id",
                        "gtp_chain", "gtp_resseq", "rmsd_A", "worst_overlap_A",
                        "method_used", "gamma_tors1_deg", "gamma_tors2_deg",
                        "rigid_shift_A", "closest_contact",
                        "best_possible_rmsd_A", "best_possible_pdb_id",
                        "blocked_by", "resolution_A", "experiment"])
            w.writerows(rows)

    if grand:
        with open(os.path.join(args.out_dir, "all_sites.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["input_file", "chain", "resseq", "segid", "gtp_pdb_id",
                        "gtp_chain", "gtp_resseq", "rmsd_A", "worst_overlap_A",
                        "method_used", "gamma_tors1_deg", "gamma_tors2_deg",
                        "rigid_shift_A", "closest_contact",
                        "best_possible_rmsd_A", "best_possible_pdb_id",
                        "blocked_by", "resolution_A", "experiment"])
            w.writerows(grand)
        bad = [r for r in grand if float(r[8]) > 0]
        print("\n%d site(s) replaced in total. %s"
              % (len(grand),
                 "All clash-free." if not bad
                 else "%d site(s) still have an overlap -- see all_sites.csv" % len(bad)))


if __name__ == "__main__":
    main()
