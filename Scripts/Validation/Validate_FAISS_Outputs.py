#!/usr/bin/env python
"""
Validate that ModelSEEDDatabase changes do not perturb the downstream
FAISS-based reaction similarity index used by Ray16/ModelSEED_FAISS.

Two validation tiers, two operations.

# Tier 1: deterministic input check (default; no extra dependencies)

Hashes the reaction-SMILES string that the FAISS pipeline would build
from each non-obsolete reaction, using the same logic as
0_gen_rxn_fps_rxnfp.py. If the hashes match the committed baseline,
the RXNFP fingerprints (deterministic transformer inference) and the
FAISS index would also match. Runs in ~10 seconds with stdlib only.

  python Validate_FAISS_Outputs.py --mode=smiles
      compares against Biochemistry/Embeddings/rxn_smiles_hashes.tsv;
      exits non-zero if any reaction's hash differs.

  python Validate_FAISS_Outputs.py --mode=smiles --regenerate
      rewrites the baseline file from current data.

# Tier 2: full pipeline check (requires rxnfp conda env + sibling repo)

Runs the actual fingerprint generation step from a sibling clone of
Ray16/ModelSEED_FAISS, hashes each reaction's 256-D fingerprint, and
compares against the committed baseline.

  python Validate_FAISS_Outputs.py --mode=fingerprint
                                   --faiss-repo ../ModelSEED_FAISS
      requires transformers, rxnfp, torch.

  python Validate_FAISS_Outputs.py --mode=fingerprint
                                   --faiss-repo ../ModelSEED_FAISS
                                   --regenerate

# Exit codes

  0   all hashes match the baseline (or baseline regenerated)
  1   one or more hashes differ; per-reaction diff printed to stdout
  2   environment error (missing baseline, missing FAISS repo, missing
      python deps for fingerprint mode)
"""
import argparse
import csv
import hashlib
import math
import os
import sys
from collections import OrderedDict

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
BIOCHEM_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', 'Biochemistry'))
EMBED_DIR    = os.path.join(BIOCHEM_ROOT, 'Embeddings')

SMILES_BASELINE      = os.path.join(EMBED_DIR, 'rxn_smiles_hashes.tsv')
FINGERPRINT_BASELINE = os.path.join(EMBED_DIR, 'fingerprint_hashes.tsv')

sys.path.append(os.path.join(SCRIPT_DIR, '..', '..', 'Libs', 'Python'))


# ---------------------------------------------------------------------------
# Tier 1: rxn_smiles hashing (vendored from Ray16/ModelSEED_FAISS)
# ---------------------------------------------------------------------------

def build_rxn_smiles(rxn_obj, compounds_dict):
    """Return a reaction SMILES string, or None if any compound lacks SMILES.

    Matches Ray16/ModelSEED_FAISS:0_gen_rxn_fps_rxnfp.py exactly so
    that fingerprints downstream are bit-identical.
    """
    rgt_smiles, pdt_smiles = [], []
    for rgt in rxn_obj['stoichiometry']:
        cpd_id = rgt['compound']
        if cpd_id not in compounds_dict:
            return None
        cpd_smiles = compounds_dict[cpd_id]['smiles']
        if cpd_smiles == '':
            return None
        count = math.ceil(abs(rgt['coefficient']))
        if rgt['coefficient'] < 0:
            rgt_smiles.extend([cpd_smiles] * count)
        elif rgt['coefficient'] > 0:
            pdt_smiles.extend([cpd_smiles] * count)
    return '>>'.join(['.'.join(rgt_smiles), '.'.join(pdt_smiles)])


def compute_smiles_hashes():
    """Return [(rxn_id, sha256_of_rxn_smiles), ...] in sorted rxn_id order."""
    from BiochemPy import Compounds, Reactions
    compounds = Compounds().loadCompounds()
    reactions = Reactions().loadReactions()
    rows = []
    for rxn_id in sorted(reactions):
        rxn = reactions[rxn_id]
        if rxn.get('is_obsolete') != 0:
            continue
        smi = build_rxn_smiles(rxn, compounds)
        if smi is None:
            continue
        h = hashlib.sha256(smi.encode('utf-8')).hexdigest()
        rows.append((rxn_id, h))
    return rows


# ---------------------------------------------------------------------------
# Tier 2: fingerprint hashing (delegates to Ray16/ModelSEED_FAISS)
# ---------------------------------------------------------------------------

def compute_fingerprint_hashes(faiss_repo):
    """Run the FAISS pipeline's fingerprint step and hash each result.

    Imports two modules from the sibling repo so we run the exact same
    code path as the production index. Requires the rxnfp conda env.
    Returns [(rxn_id, sha256_of_fingerprint_bytes), ...].
    """
    if not os.path.isdir(faiss_repo):
        die(f'FAISS repo not found at {faiss_repo}', code=2)

    sys.path.insert(0, faiss_repo)
    try:
        from BiochemPy import Compounds, Reactions  # MSDB
    except ImportError:
        die('BiochemPy not importable', code=2)
    try:
        import numpy as np
        from rxnfp.transformer_fingerprints import (
            RXNBERTFingerprintGenerator,
            get_default_model_and_tokenizer,
        )
    except ImportError as e:
        die(f'fingerprint mode requires rxnfp + transformers + torch ({e})', code=2)

    compounds = Compounds().loadCompounds()
    reactions = Reactions().loadReactions()

    rxn_ids   = []
    smi_strs  = []
    for rxn_id in sorted(reactions):
        rxn = reactions[rxn_id]
        if rxn.get('is_obsolete') != 0:
            continue
        smi = build_rxn_smiles(rxn, compounds)
        if smi is None:
            continue
        rxn_ids.append(rxn_id)
        smi_strs.append(smi)

    print(f'  built rxn_smiles for {len(smi_strs)} reactions, generating fingerprints...')
    model, tok = get_default_model_and_tokenizer()
    gen        = RXNBERTFingerprintGenerator(model, tok)
    fps        = []
    BATCH      = 1000
    for i in range(0, len(smi_strs), BATCH):
        fps.extend(gen.convert_batch(smi_strs[i:i+BATCH]))
        print(f'    {len(fps)}/{len(smi_strs)}')
    arr = np.array(fps, dtype=np.float32)

    rows = []
    for i, rid in enumerate(rxn_ids):
        h = hashlib.sha256(arr[i].tobytes()).hexdigest()
        rows.append((rid, h))
    return rows


# ---------------------------------------------------------------------------
# baseline file io + diff
# ---------------------------------------------------------------------------

def load_baseline(path):
    if not os.path.isfile(path):
        return None
    out = {}
    with open(path) as fh:
        # skip leading '#' comment lines so DictReader picks up the
        # real header row.
        lines = [ln for ln in fh if not ln.startswith('#')]
    reader = csv.DictReader(lines, dialect='excel-tab')
    for row in reader:
        out[row['reaction_id']] = row['hash']
    return out


def write_baseline(path, rows, header_extra=None):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        if header_extra:
            for line in header_extra.splitlines():
                fh.write('# ' + line + '\n')
        fh.write('reaction_id\thash\n')
        for rid, h in rows:
            fh.write(f'{rid}\t{h}\n')


def diff_against_baseline(rows, baseline):
    """Return (added, removed, changed) lists."""
    cur = dict(rows)
    added   = sorted(set(cur)      - set(baseline))
    removed = sorted(set(baseline) - set(cur))
    changed = sorted(rid for rid in cur if rid in baseline and cur[rid] != baseline[rid])
    return added, removed, changed


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def die(msg, code=2):
    print(f'ERROR: {msg}', file=sys.stderr)
    sys.exit(code)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mode', choices=['smiles', 'fingerprint'], default='smiles')
    ap.add_argument('--regenerate', action='store_true',
                    help='Update the baseline file with current data')
    ap.add_argument('--faiss-repo', default=os.path.join(SCRIPT_DIR, '..', '..', '..', 'ModelSEED_FAISS'),
                    help='Path to a sibling clone of Ray16/ModelSEED_FAISS (fingerprint mode)')
    args = ap.parse_args()

    print(f'mode: {args.mode}  regenerate: {args.regenerate}')
    if args.mode == 'smiles':
        rows = compute_smiles_hashes()
        baseline_path = SMILES_BASELINE
        header = 'Per-reaction sha256 of the rxn_smiles string built from\nthe ModelSEEDDatabase compound and reaction tables, matching\nthe logic in Ray16/ModelSEED_FAISS:0_gen_rxn_fps_rxnfp.py.'
    else:
        rows = compute_fingerprint_hashes(os.path.abspath(args.faiss_repo))
        baseline_path = FINGERPRINT_BASELINE
        header = 'Per-reaction sha256 of the 256-D RXNBERT fingerprint produced\nby Ray16/ModelSEED_FAISS:0_gen_rxn_fps_rxnfp.py over the current\nModelSEEDDatabase. Regenerate this file from the rxnfp conda env.'

    print(f'  reactions hashed: {len(rows)}')

    if args.regenerate:
        write_baseline(baseline_path, rows, header_extra=header)
        print(f'wrote baseline: {baseline_path}')
        return

    baseline = load_baseline(baseline_path)
    if baseline is None:
        die(f'baseline not found: {baseline_path}\nRun with --regenerate to create it.', code=2)

    added, removed, changed = diff_against_baseline(rows, baseline)
    if not (added or removed or changed):
        print('PASS: all hashes match baseline')
        return

    print(f'FAIL: {len(added)} added, {len(removed)} removed, {len(changed)} changed')
    LIMIT = 30
    if added:
        print('\nAdded reactions (first {}):'.format(min(LIMIT, len(added))))
        for r in added[:LIMIT]:
            print(f'  + {r}')
    if removed:
        print('\nRemoved reactions (first {}):'.format(min(LIMIT, len(removed))))
        for r in removed[:LIMIT]:
            print(f'  - {r}')
    if changed:
        print('\nHash-changed reactions (first {}):'.format(min(LIMIT, len(changed))))
        cur = dict(rows)
        for r in changed[:LIMIT]:
            print(f'  ~ {r}  baseline={baseline[r][:12]}…  now={cur[r][:12]}…')
    sys.exit(1)


if __name__ == '__main__':
    main()
