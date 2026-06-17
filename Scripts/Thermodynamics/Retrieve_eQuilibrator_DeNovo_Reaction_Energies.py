#!/usr/bin/env python
"""Tier-2 eQuilibrator coverage expansion: de-novo reaction energies.

The legacy eQuilibrator pipeline only *looks up* energies for compounds already
present in eQuilibrator's precomputed cache (compounds that overlap MetaNetX by
InChIKey). Compounds with a structure in ModelSEED but absent from the cache get
no energy, so any reaction that uses one is left thermodynamically undefined.

This script estimates energies DE-NOVO from structure for those uncached
compounds, using equilibrator-assets' group decomposition (component
contribution), then scores every reaction whose reagents now all resolve
(cached MNX compound OR de-novo compound). Component-contribution uncertainties
cancel across a reaction, so a reaction can get a finite, usable dG'o even when
some reagents' individual formation energies do not.

Output: Biochemistry/Thermodynamics/eQuilibrator/MetaNetX_DeNovo_Reaction_Energies.tbl
  one row per scored reaction:  rxn <tab> dG'o(kcal/mol) <tab> uncertainty(kcal/mol) <tab> ln_RI
  (reactions whose uncertainty is non-finite / above --max-error are written
   with "Unable to retrieve energy" so re-runs skip them).
These feed Update_Reaction_eQuilibrator_DeNovo_Energies.py, which writes them
ADDITIVELY under thermodynamics['eQuilibrator-deNovo'] (a distinct, lower-
confidence source — never conflated with the cache-derived 'eQuilibrator').

ChemAxon note: proper pH-7 microspecies need a licensed ChemAxon cxcalc. When it
is absent we run in "bypass" mode (decompose the user structure, single default
protonation). Transformed energies are then approximate for compounds whose
ModelSEED structure is not already the pH-7 major species. For full fidelity
WITHOUT ChemAxon, pass the database's own MolGpKa/Marvin pKas as specified_pkas
(see --pka-source; the per-compound pkas already live in the compound records).
"""
import os, sys, csv, math, time, shutil, argparse
import warnings; warnings.filterwarnings('ignore')
import logging; logging.disable(logging.WARNING)
import pandas as pd

THERMO = os.path.dirname(os.path.abspath(__file__)) + "/../../Biochemistry/Thermodynamics/"
STRUCT = os.path.dirname(os.path.abspath(__file__)) + "/../../Biochemistry/Structures/"
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../Libs/Python/")

SENTINEL = 10000000
TRIVIAL = {"cpd00067", "cpd00001"}  # H+, H2O — eQuilibrator references these itself


def install_no_chemaxon_shims():
    """Allow equilibrator-assets to create compounds without a ChemAxon license.

    v0.6.0 calls cxcalc unconditionally (get_compound_mappings) and only honours
    bypass_chemaxon afterwards, so without these two shims compound creation
    crashes outright. They make the ChemAxon path yield 'no protonation', which
    forces the documented 'bypass' decomposer (user structure, default
    protonation). They are no-ops if cxcalc IS installed (then don't call this).
    """
    import equilibrator_assets.chemaxon as cax
    cax.get_dissociation_constants = lambda molecules, error_log, num_acidic=20, num_basic=20, mid_ph=7.0: (
        pd.DataFrame({"id": list(molecules["id"]), "major_ms": [None] * len(molecules)}), [])
    import equilibrator_assets.generate_compound as gc
    _orig = gc._populate_compound_information
    gc._populate_compound_information = lambda row: (None if row.method == "chemaxon" else _orig(row))


def load_cached_seed_map():
    """seed cpd id -> MNX accession, for compounds with a VALID cache eQ energy."""
    valid_mnx = set()
    with open(THERMO + "eQuilibrator/MetaNetX_Compound_Energies.tbl") as fh:
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if len(a) >= 3 and "energy" not in a[1] and a[1] != "nan":
                valid_mnx.add(a[0])
    ik2mnx = {}
    with open(STRUCT + "MetaNetX/Structures_in_ModelSEED_and_eQuilibrator.txt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] in valid_mnx:
                ik2mnx.setdefault(p[1], p[0])
    return valid_mnx, ik2mnx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=os.path.expanduser("~/.cache/equilibrator/compounds.sqlite"),
                    help="base eQuilibrator compound cache (read-only source)")
    ap.add_argument("--augmented", default=os.path.expanduser("~/.cache/equilibrator/denovo_compounds.sqlite"),
                    help="writable augmented cache (base copy + de-novo compounds)")
    ap.add_argument("--max-error", type=float, default=100.0, help="kcal/mol; drop reactions above this uncertainty")
    ap.add_argument("--limit", type=int, default=0, help="cap number of candidate reactions (0 = all) for validation")
    ap.add_argument("--out", default=THERMO + "eQuilibrator/MetaNetX_DeNovo_Reaction_Energies.tbl")
    args = ap.parse_args()

    from BiochemPy import Compounds, Reactions
    C = Compounds(); R = Reactions()
    first = lambda d: list(d.keys())[0]
    ik = C.loadStructures(["InChIKey"], ["ModelSEED"])
    inchi = C.loadStructures(["InChI"], ["ModelSEED"])
    cpd_ik = {c: first(ik[c]["InChIKey"]) for c in ik if "InChIKey" in ik[c]}
    cpd_in = {c: first(inchi[c]["InChI"]) for c in inchi if "InChI" in inchi[c]}

    valid_mnx, ik2mnx = load_cached_seed_map()
    seed_cached = {c: ik2mnx[cpd_ik[c]] for c in cpd_ik if cpd_ik.get(c) in ik2mnx}
    print(f"cached seed compounds (valid eQ): {len(seed_cached)}")

    reactions = R.loadReactions()

    # Candidate reactions: not EMPTY, every non-(H+/H2O) reagent has a structure,
    # and at least one reagent is uncached (so the cache pipeline left it incomplete)
    # with an InChI we can decompose.
    candidates = []
    for rid in sorted(reactions.keys()):
        r = reactions[rid]
        if r.get("status") == "EMPTY":
            continue
        core = [s["compound"] for s in (r.get("stoichiometry") or []) if s["compound"] not in TRIVIAL]
        if not core or not all(c in cpd_ik for c in core):
            continue
        uncached = [c for c in core if c not in seed_cached]
        if not uncached or not all(c in cpd_in for c in uncached):
            continue
        candidates.append(rid)
    if args.limit:
        candidates = candidates[: args.limit]
    print(f"candidate reactions de-novo could complete: {len(candidates)}")

    # Resume: skip reactions already in the output file.
    done = set()
    if os.path.exists(args.out):
        with open(args.out) as fh:
            for line in fh:
                done.add(line.split("\t")[0])
    candidates = [r for r in candidates if r not in done]
    print(f"to score this run: {len(candidates)} ({len(done)} already done)")
    if not candidates:
        return

    # Build / load the augmented cache (base copy + de-novo compounds).
    if not os.path.exists(args.augmented):
        print(f"copying base cache -> {args.augmented} ...")
        shutil.copy(args.cache, args.augmented)

    if shutil.which("cxcalc") is None:
        print("ChemAxon cxcalc unavailable -> installing no-ChemAxon bypass shims "
              "(de-novo energies use default protonation; see header note).")
        install_no_chemaxon_shims()

    from equilibrator_assets.local_compound_cache import LocalCompoundCache
    from equilibrator_api import ComponentContribution, Q_
    from equilibrator_api.phased_reaction import PhasedReaction
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from Estimate_Reaction_Reversibility import reversibility_from_energy  # noqa: F401 (re-used downstream)

    lc = LocalCompoundCache(); lc.load_cache(args.augmented); lc.read_only = False

    # Decompose+add every uncached reagent used by the candidates (batched).
    need = sorted({c for rid in candidates for c in
                   [s["compound"] for s in reactions[rid]["stoichiometry"]
                    if s["compound"] not in TRIVIAL and s["compound"] not in seed_cached]})
    print(f"de-novo compounds to decompose: {len(need)}")
    denovo = {}
    BATCH = 250
    t0 = time.time()
    for i in range(0, len(need), BATCH):
        chunk = need[i:i + BATCH]
        res = lc.get_compounds([cpd_in[c] for c in chunk], mol_format="inchi",
                               bypass_chemaxon=True, save_empty_compounds=True)
        for c, rr in zip(chunk, res):
            if getattr(rr, "method", "") != "empty" and rr.compound is not None:
                denovo[c] = rr.compound
        print(f"  decomposed {min(i+BATCH,len(need))}/{len(need)} "
              f"({len(denovo)} ok) [{time.time()-t0:.0f}s]")

    cc = ComponentContribution(ccache=lc.ccache)
    cc.p_h = Q_(7.0); cc.ionic_strength = Q_("0.25M"); cc.temperature = Q_("298.15K")

    cached_obj = {}
    def compound_for(seed):
        if seed in denovo:
            return denovo[seed]
        if seed in seed_cached:
            if seed not in cached_obj:
                try:
                    cached_obj[seed] = lc.ccache.get_compound(seed_cached[seed])
                except Exception:
                    cached_obj[seed] = None
            return cached_obj[seed]
        return None

    out = open(args.out, "a")
    scored = infinite = unresolved = errored = 0
    t0 = time.time()
    for n, rid in enumerate(candidates, 1):
        sparse = {}
        ok = True
        for s in reactions[rid]["stoichiometry"]:
            comp = compound_for(s["compound"])
            if comp is None:
                ok = False; break
            sparse[comp] = sparse.get(comp, 0.0) + float(s["coefficient"])
        if not ok:
            unresolved += 1; continue
        sparse = {k: v for k, v in sparse.items() if abs(v) > 1e-9}
        if not sparse:
            unresolved += 1; continue
        try:
            rxn = PhasedReaction(sparse)
            m = cc.standard_dg_prime(rxn)
            dg = m.value.to("kcal/mol").magnitude
            err = m.error.to("kcal/mol").magnitude
            if not math.isfinite(err) or err > args.max_error:
                infinite += 1
                out.write(f"{rid}\tUnable to retrieve energy\n")
            else:
                try:
                    ln_ri = cc.ln_reversibility_index(rxn)
                    ln_ri = str(ln_ri.magnitude if hasattr(ln_ri, "magnitude") else ln_ri)
                except Exception:
                    ln_ri = "nan"
                out.write(f"{rid}\t{dg}\t{err}\t{ln_ri}\n")
                scored += 1
        except Exception:
            errored += 1
            out.write(f"{rid}\tUnable to retrieve energy\n")
        if n % 250 == 0:
            out.flush()
            print(f"  scored {n}/{len(candidates)} (usable={scored}) [{time.time()-t0:.0f}s]")
    out.close()
    print(f"\nDONE: usable={scored} | infinite/over-error={infinite} | unresolved={unresolved} | errored={errored}")
    print(f"output: {args.out}")


if __name__ == "__main__":
    main()
